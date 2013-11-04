"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from urllib import quote
from django.conf import settings
from django.core.urlresolvers import reverse


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


def get_reaction_map(map_id, enzymes, items):
    import StringIO

    # Wordaround broken PIL installation
    # see: https://code.djangoproject.com/ticket/6054
    try:
        from PIL import Image
    except ImportError:
        import Image

    from xml.etree.ElementTree import ElementTree, Element

    with open("{}/kegg/img/{}.png".format(settings.STATICFILES_DIRS[0], map_id), "rb") as image:
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
        image.set("xlink:href", "{}kegg/img/{}.png".format(settings.STATIC_URL, map_id))
        root.append(image)

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

            # Pathways are blue
            # Something with EC numbers green
            # Everything else red
            color_component = "red"
            fill_opacity = "0.0"
            fill_color = "green" # not displayed with 0.0

            # Found objects in the organism show a background color

            # Repoint pathway links to our versions
            if "show_pathway?" in url:
                pathway_name = url[url.index("?") + 1:]
                color_component = "blue"

                ##try:
                ##    pw_obj = Pathway.objects.for_species(species).for_wid(pathway_name)
                elem.set("xlink:href", reverse("kegg.views.map_view", kwargs={"map_id": pathway_name}) + "?items=" + quote(items))
                ##    fill_opacity = "0.2"
                ##    fill_color = "blue"
            else:
                elem.set("xlink:href", url)
                elem.set("target", "_blank")

            elem.set("xlink:title", title)
            root.remove(area)
            elem.append(area)

            ecs = uniqify(extract_ecs(title))

            if len(ecs) > 0:
                color_component = "green"
                for ec in ecs:
                    for enzyme in enzymes:
                        if enzyme[0] == ec:
                            fill_opacity = "0.2"
                            if enzyme[1] is not None:
                                color_component = enzyme[1]
                                fill_color = enzyme[1]

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

        tree.write(out)

        return out.getvalue()


def request_extract_ecs(request):
    import re
    ##if not "items" in request.POST:
    ##    items = [["1.1.1.1", None],
    ##             ["2.2.2.2", None],
    ##             ["4.1.2.20", "green"],
    ##             ["1.1.1.20", None],
    ##             ["1.2.1.3", None],
    ##             ["ascorbate", "red"]]
    ##else:
    items = []
    if "items" in request.GET:
        get_items = request.GET["items"].replace("%23", "#")
        for item in re.split("\s+", get_items):
            if len(item) == 0:
                continue

            if "#" in item:
                item = item.split("#", 2)
                items.append([item[0], item[1]])
            else:
                items.append([item, None])

    ##metabolite_items = []
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

        ##metabolite_items.append(item)

    return enzyme_items
