from PIL import Image
from xml.etree.ElementTree import ElementTree, Element, SubElement
from StringIO import StringIO
import re

im = Image.open("map01100_0.3531465.png")
idata = list(im.getdata())
iwidth = im.size[0]

print_map = 0

map = StringIO()

with open("show_pathway@map01100", "r") as kegg:
    for line in kegg:
        if line[:4] == "<map":
            print_map = 1
        
        if print_map:
            # Fixup to valid XML
            line = re.sub(r"coords=([^\"\s]+)", r'coords="\1"', line.strip())
            line = re.sub(r"shape=([^\"\s]+)", r'shape="\1"', line.strip())
            map.write(line.strip() + "\n")
        
        if line[:6] == "</map>":
            print_map = 0

map.seek(0)
tree = ElementTree()
tree.parse(map)
root = tree.getroot()

root.tag = "svg"
root.set("xmlns", "http://www.w3.org/2000/svg")
root.set("xmlns:xlink", "http://www.w3.org/1999/xlink")
root.set("width", "100%")
root.set("height", "100%")

areas = tree.findall("area")

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
    elem.set("xlink:href", url)
    elem.set("xlink:title", title)
    root.remove(area)
    elem.append(area)
    root.append(elem)
    
    if shape == "circle":
        area.set("cx", coords[0])
        area.set("cy", coords[1])
        area.set("r", coords[2])
        #print coords[0] + " " + coords[1]
        #print idata[int(coords[1]) * iwidth + int(coords[0])]
        area.set("fill", "rgb(" + ",".join(str(x) for x in idata[int(coords[1]) * iwidth + int(coords[0])][:3]) + ")")
    elif shape == "rect":
        area.set("x", coords[0])
        area.set("y", coords[1])
        area.set("width", str(int(coords[2]) - int(coords[0])))
        area.set("height", str(int(coords[3]) - int(coords[1])))
    elif shape == "polygon":
        a = zip(*2*[iter(coords)])
        area.set("points", " ".join([",".join(x) for x in zip(*2*[iter(coords)])]))
        area.set("fill", "rgb(" + ",".join(str(x) for x in idata[int(a[0][1]) * iwidth + int(a[0][0])][:3]) + ")")

tree.write('kegg_out.html')

    