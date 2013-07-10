"""
Takes the data from fetch_kegg and converts it to svg files
(All image map components are detected and replaced by svg elements and color
coding is added)
"""

import re
import os
import urllib2
import shutil
from PIL import Image
from xml.etree.ElementTree import ElementTree, Element, SubElement
from StringIO import StringIO

def extract_ecs(text):
    """Extracts EC numbers out of a string and returns a list with all numbers"""
    return re.findall(r"[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+", text)

def uniqify(seq, idfun=None):
    """Order preserving list uniqifier.
    Source: http://www.peterbe.com/plog/uniqifiers-benchmark
    """
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        if marker in seen: continue
        seen[marker] = 1
        result.append(item)
    return result

html_files = filter(lambda x: x.endswith(".html"), os.listdir("fetch"))
html_files_count = len(html_files)

for i, file in enumerate(html_files):
    print "Parsing {} ({}/{})".format(file, i + 1, html_files_count)
    
    file = file[:-5] # Remove extension
    
    with open("fetch/{}.png".format(file), "rb") as image:
        ####import pdb; pdb.set_trace()
        im = Image.open(image)
        idata = list(im.getdata())
        iwidth = im.size[0]

    with open("fetch/{}.html".format(file)) as infile:
        tree = ElementTree()
        tree.parse(infile)
        root = tree.getroot()

        root.tag = "svg"
        root.set("xmlns", "http://www.w3.org/2000/svg")
        root.set("xmlns:xlink", "http://www.w3.org/1999/xlink")
        root.set("width", "100%")
        root.set("height", "100%")

        script = Element("script")
        script.set("xlink:href", "js/SVGPan.js")
        root.append(script)
        
        graphics = Element("g")
        graphics.set("id", "viewport")
        root.append(graphics)
        
        image = Element("image")
        image.set("x", "0")
        image.set("y", "0")
        image.set("width", str(iwidth))
        image.set("height", str(im.size[1]))
        shutil.copyfile("fetch/{}.png".format(file), "svg/img/{}.png".format(file))
        image.set("xlink:href", "img/{}.png".format(file))
        graphics.append(image)

        areas = tree.findall("area")
        
        # Draw polys before the rest because they are line segments and sometimes
        # overlap other elements
        pending = []
        
        for area in areas:
            shape = area.get("shape")
            coords = area.get("coords")
            url = area.get("href")
            title = area.get("title")

            if shape is None:
                continue
            
            if coords is None:
                continue
            
            if url is None:
                continue
            
            if shape == "poly":
                shape = "polygon"
            
            url = "http://www.kegg.jp" + url
        
            coords = coords.split(",")
        
            area.attrib.pop("shape")
            area.attrib.pop("coords")
            area.attrib.pop("href")
            area.attrib.pop("title")
        
            area.tag = shape
            
            elem = Element("a")
            
            # Pathways are green
            # Something with EC numbers blue
            # Everything else red 
            color_component = "rgb(255,0,0)"
            
            # Repoint pathway links to downloaded versions
            if "show_pathway?" in url:
                elem.set("xlink:href", "{}.html".format(url[url.index("?") + 1:]))
                color_component = "rgb(0,255,0)"
            else:
                elem.set("xlink:href", url)

            elem.set("xlink:title", title)
            root.remove(area)
            elem.append(area)
            
            ## Debug EC value text
            ##text_elem = Element("text")
            ##text_elem.set("x", coords[0])
            ##text_elem.set("y", str(int(coords[1]) - 5))
            ##text_elem.set("font-size", "9")
            ##text_elem.set("font-family", "sans-serif")
            ##text_elem.text = ", ".join(uniqify(extract_ecs(title)))
            ##graphics.append(text_elem)
            
            if len(uniqify(extract_ecs(title))) > 0:
                if color_component == "rgb(0,255,0)":
                    print "Pathway link with EC_Number!"
                color_component = "rgb(0,0,255)"
            
            
            
            if shape == "circle":
                pending.append(elem)

                area.set("cx", coords[0])
                area.set("cy", coords[1])
                area.set("r", str(int(coords[2]) + 1))        
                #print coords[0] + " " + coords[1]
                #print idata[int(coords[1]) * iwidth + int(coords[0])]
                ##area.set("fill", "rgb(" + ",".join(str(x) for x in idata[int(coords[1]) * iwidth + int(coords[0])][:3]) + ")")
            elif shape == "rect":
                pending.append(elem)

                area.set("x", coords[0])
                area.set("y", coords[1])
                area.set("width", str(int(coords[2]) - int(coords[0])))
                area.set("height", str(int(coords[3]) - int(coords[1])))
            elif shape == "polygon":
                graphics.append(elem)
                   
                points = zip(*2*[iter(coords)])

                area.set("points", " ".join([",".join(x) for x in points]))
            
            area.set("fill-opacity", "0.0")
            area.set("style", "stroke-width:1;stroke:{}".format(color_component))
        
        for elem in pending:
            graphics.append(elem)
        
        out = StringIO()
        out.write('<html><head><script type="text/javascript" src="js/jquery-1.8.0.min.js"></script><script type="text/javascript" src="js/jquery-svgpan.js"></script></head><body>')
        tree.write(out)
        out.write('<script type="text/javascript">$("svg").svgPan("viewport");</script></body></html>')

        with open("svg/{}.html".format(file), "w") as kegg:
            kegg.write(out.getvalue())
