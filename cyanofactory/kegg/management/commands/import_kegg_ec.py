"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from django.core.management.base import BaseCommand
import re
import os
from xml.etree.ElementTree import ElementTree
import kegg.models as models


def extract_ecs(text):
    """Extracts EC numbers out of a string and returns a list with all numbers"""
    return re.findall(r"[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+", text)


def uniqify(seq, idfun=None):
    """Order preserving list uniqifier.
    Source: http://www.peterbe.com/plog/uniqifiers-benchmark
    """
    if idfun is None:
        def idfun(x):
            return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        if marker in seen:
            continue
        seen[marker] = 1
        result.append(item)
    return result


def is_overview(pathway):
    return pathway in ["map01100", "map01110", "map01120"]


class Command(BaseCommand):
    help = 'Takes the data from fetch_kegg and loads it into the database'

    def handle(self, *args, **options):
        # [:-5] removes extension
        html_files = map(lambda y: y[:-5], filter(lambda y: y.endswith(".html"), os.listdir("../kegg/fetch")))
        
        html_files_count = len(html_files)
        
        all_ecs = []
        ec_numbers = {}
        map_names = {}
        
        for i, filename in enumerate(html_files):
            print "Parsing {} ({}/{})".format(filename, i + 1, html_files_count)
        
            with open("../kegg/fetch/{}.html".format(filename)) as infile:
                tree = ElementTree()
                tree.parse(infile)
        
                areas = tree.findall("area")
                
                for area in areas:
                    title = area.get("title")
                    coords = area.get("coords")
                    shape = area.get("shape")
        
                    if title is None or\
                        coords is None or\
                        shape is None:
                        print "skipping"
                        continue
                    
                    if filename in title and ":" in title:
                        map_names[filename] = title[title.index(":") + 2:]
                    
                    ecs = extract_ecs(title)
                    
                    if len(ecs) == 0:
                        continue
                    
                    ec_list = ec_numbers.get(filename, [])
                    ec_list += ecs
                    all_ecs += ecs
                    ec_numbers[filename] = ec_list
        
        print("Importing into Database")
        all_ecs = sorted(uniqify(all_ecs))
        ec_items = [models.EcNumber(name=x) for x in all_ecs]
        if models.EcNumber.objects.filter(name=ec_items[0].name).exists():
            print "Already imported, terminating..."
            return
        else:
            models.EcNumber.objects.bulk_create(ec_items)

        for i, filename in enumerate(sorted(html_files)):            
            print "Importing {} ({}/{})".format(filename, i + 1, html_files_count)

            numbers = ec_numbers.get(filename)
            if not numbers is None:
                title = map_names.get(filename, "")
                map_, _ = models.Map.objects.get_or_create(name=filename)
                map_.title = title
                if is_overview(filename):
                    map_.overview = True
                map_.save()

                for ec in models.EcNumber.objects.filter(name__in=numbers):
                    map_.ec_numbers.add(ec)
