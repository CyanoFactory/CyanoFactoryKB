"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from django.conf import settings
from django.core.urlresolvers import reverse
import base64


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


def extract_ecs(text):
    """Extracts EC numbers out of a string and returns a list with all numbers"""
    import re
    return re.findall(r"[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+", text)


def get_reaction_map(map_id, enzymes, metabolites, export):
    import StringIO

    # Workaround broken PIL installation
    # see: https://code.djangoproject.com/ticket/6054
    try:
        from PIL import Image
    except ImportError:
        import Image

    from xml.etree.ElementTree import ElementTree, Element

    with open("{}/kegg/img/{}.png".format(settings.STATICFILES_DIRS[0], map_id), "rb") as image:
        if export:
            encoded_string = base64.b64encode(image.read())
            image.seek(0)

        im = Image.open(image)
        iwidth = im.size[0]
        iheight = im.size[1]

    with open("{}/kegg/maps/{}.html".format(settings.ROOT_DIR, map_id)) as infile:
        tree = ElementTree()
        tree.parse(infile)
        root = tree.getroot()

        root.tag = "g"
        root.set("id", "viewport")

        image = Element("image")
        image.set("x", "0")
        image.set("y", "0")
        image.set("width", str(iwidth))
        image.set("height", str(iheight))

        if export:
            href = "data:image/png;base64," + encoded_string
        else:
            href = "{}kegg/img/{}.png".format(settings.STATIC_URL, map_id)
        image.set("xlink:href", href)
        root.append(image)

        areas = tree.findall("area")

        # Draw polys before the rest because they are line segments and sometimes
        # overlap other elements
        pending = []

        enzymes_found = []
        metabolites_found = []

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

            # Pathways are blue
            # Something with EC numbers green
            # Everything else red
            fill_opacity = "0.0"

            # Found objects in the organism show a background color

            elem.set("xlink:title", title)
            root.remove(area)
            elem.append(area)

            ecs = uniqify(extract_ecs(title))

            # Repoint pathway links to our versions
            if "show_pathway?" in url:
                pathway_name = url[url.index("?") + 1:]
                color_component = fill_color = "blue"

                ##try:
                ##    pw_obj = Pathway.objects.for_species(species).for_wid(pathway_name)
                elem.set("xlink:href", reverse("kegg.views.map_view", kwargs={"map_id": pathway_name}))

                for metabolite in metabolites:
                    ltitle = title.lower()
                    if metabolite[0].lower() in ltitle:
                        fill_opacity = "0.3"
                        metabolites_found.append(metabolite)
                        if metabolite[1] is not None:
                            color_component = fill_color = metabolite[1]
            else:
                elem.set("xlink:href", url)
                elem.set("target", "_blank")
                # Handle enzymes
                if len(ecs) > 0:
                    color_component = fill_color = "green"
                    for ec in ecs:
                        for enzyme in enzymes:
                            if enzyme[0] == ec:
                                fill_opacity = "0.3"
                                enzymes_found.append(enzyme)
                                if enzyme[1] is not None:
                                    color_component = fill_color = enzyme[1]
                else:
                    # Everything else: metabolite
                    color_component = fill_color = "red"
                    for metabolite in metabolites:
                        ltitle = title.lower()
                        if metabolite[0].lower() in ltitle:
                            fill_opacity = "0.3"
                            metabolites_found.append(metabolite)
                            if metabolite[1] is not None:
                                color_component = fill_color = metabolite[1]

            if shape == "circle":
                pending.append(elem)

                area.set("cx", coords[0])
                area.set("cy", coords[1])
                area.set("r", str(int(coords[2]) + 1))
            elif shape == "rect":
                pending.append(elem)

                area.set("x", coords[0])
                area.set("y", coords[1])
                area.set("width", str(int(coords[2]) - int(coords[0])))
                area.set("height", str(int(coords[3]) - int(coords[1])))
            elif shape == "polygon":
                root.append(elem)

                points = zip(*2*[iter(coords)])

                area.set("points", " ".join([",".join(x) for x in points]))

            area.set("style", "stroke-width:1;stroke:{};fill-opacity:{};fill:{}".format(color_component, fill_opacity, fill_color))

        for elem in pending:
            root.append(elem)


        ##template = loader.get_template("cyano/pathway/sidebar.html")

        out = StringIO.StringIO()
        ##out.write(template.render(Context()))
        out.write("""<svg id="kegg" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="100%" height="100%" pointer-events="visible">""")
        tree.write(out)
        out.write("</svg>")

        # overwrite enzymes and metabolites content
        # Reference preserving list clear
        # (via http://stackoverflow.com/questions/850795)
        del enzymes[:]
        del metabolites[:]
        enzymes += uniqify(enzymes_found, lambda x: x[0])
        metabolites += uniqify(metabolites_found, lambda x: x[0])

        return out.getvalue()


def request_extract(request):
    import re

    items = []
    if "items" in request.POST:
        get_items = request.POST["items"].replace("%23", "#")
        for item in re.split("\s+", get_items):
            if len(item) == 0:
                continue

            if "#" in item:
                item = item.split("#", 2)
                if len(item[0]) == 0:
                    continue
                item.append(item[1])
                try:
                    h = (100 - int(item[2])) * 2.4
                    item[2] = "hsl(" + str(h) + ", 100%, 50%)"
                except ValueError:
                    pass
                items.append([item[0], item[2], item[1]])
            else:
                items.append([item, None, None])

    metabolite_items = []
    enzyme_items = []

    for item in items:
        maybe_enzyme = item[0].split(".")
        if len(maybe_enzyme) == 4:
            try:
                map(lambda x: int(x), maybe_enzyme)
                enzyme_items.append(item)
                continue
            except ValueError:
                # not a valid EC number, maybe a metabolite
                pass

        metabolite_items.append(item)

    return items, enzyme_items, metabolite_items


def items_to_quoted_string(items):
    itemstr = ""

    for name, color in items:
        if color is not None:
            itemstr += name + "%23" + color
        else:
            itemstr += name
        itemstr += "%20"

    return itemstr
