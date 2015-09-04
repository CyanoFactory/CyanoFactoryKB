"""
Converts KEGG pathway elements to a SVG graphic.

This was more of an experiment and is not suitable for all pathways
because they lack information about some components on the image map
"""
from PIL import Image
from xml.etree.ElementTree import ElementTree, Element, SubElement
from io import StringIO
import re

def luminance(param):
    return 0.2126*int(param[0]) + 0.7152*int(param[1]) + 0.0722*int(param[2])

im = Image.open("map01100.png")
idata = list(im.getdata())
iwidth = im.size[0]

print_map = 0

img_map = StringIO()

with open("map01100.html", "r") as kegg:
    for line in kegg:
        if line.startswith("<map"):
            print_map = 1
        
        if print_map:
            # fix up to valid XML
            line = re.sub(r"coords=([^\"\s]+)", r'coords="\1"', line.strip())
            line = re.sub(r"shape=([^\"\s]+)", r'shape="\1"', line)
            # escape non-escaped ampersands
            line = re.sub(r"&(?!(?:apos|quot|[gl]t|amp);|#)", r'&amp;', line)
            img_map.write(line + "\n")
        
        if line.endswith("</map>"):
            print_map = 0

img_map.seek(0)
tree = ElementTree()
tree.parse(img_map)
root = tree.getroot()

root.tag = "svg"
root.set("xmlns", "http://www.w3.org/2000/svg")
root.set("xmlns:xlink", "http://www.w3.org/1999/xlink")
root.set("width", "100%")
root.set("height", "100%")

script = Element("script")
script.set("xlink:href", "SVGPan.js")
root.append(script)

graphics = Element("g")
graphics.set("id", "viewport")
root.append(graphics)

image = Element("image")
image.set("x", "0")
image.set("y", "0")
image.set("width", str(iwidth))
image.set("height", str(im.size[1]))
image.set("xlink:href", "map01100.png")
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
    elem.set("xlink:href", url)
    elem.set("xlink:title", title)
    root.remove(area)
    elem.append(area)
    
    if shape == "circle":
        pending.append(elem)
        
        area.set("cx", coords[0])
        area.set("cy", coords[1])
        area.set("r", str(int(coords[2]) + 1))        
        #print coords[0] + " " + coords[1]
        #print idata[int(coords[1]) * iwidth + int(coords[0])]
        area.set("fill", "rgb(" + ",".join(str(x) for x in idata[int(coords[1]) * iwidth + int(coords[0])][:3]) + ")")
    elif shape == "rect":
        continue # FIXME: Bad font rendering
        pending.append(elem)
        
        area.tag = "text"
        area.set("x", coords[0])
        area.set("y", coords[1])
        area.set("font-size", "6")
        area.set("font-family", "sans-serif")
        if ":" in title:
            title = title.split(":", 2)[1].strip()
        #area.text = title
        width = int(coords[2]) - int(coords[0])
        height = int(coords[3]) - int(coords[1])
        
        text = title.split(" ")
        split_text = []
        
        cur_width = 0
        
        first = True
        #if title == "Glycosphingolipid biosynthesis - lacto and neolacto series":
        #    import pdb; pdb.set_trace()
        print title

        while len(text) > 0:
            cur_width = len(text[0]) * 6
            if cur_width > width:
                elem = Element("tspan")
                if len(split_text) == 0:
                    elem.text = text[0]
                else:
                    split_text = [text[0]]
                elem.set("x", coords[0])
                if first:
                    first = False
                else:
                    elem.set("dy", "1.2em")
                area.append(elem)
            else:
                split_text.append(text[0])
            
            text = text[1:]
            print(" ".join(text))
        
        elem = Element("tspan")
        elem.text = " ".join(split_text)
        elem.set("x", coords[0])
        if not first:
            elem.set("dy", "1.2em")
        area.append(elem)

        #area.set("width", str(int(coords[2]) - int(coords[0])))
        #area.set("height", str(int(coords[3]) - int(coords[1])))
    elif shape == "polygon":
        graphics.append(elem)
        
        points = zip(*2*[iter(coords)])
        
        if points[0] != points[-1]:
            points.append(points[0])
            
        a = 0
        
        # A = 1/2 Sum(x_i * y_i+1 - x_i+1 * y_i)
        for i, p in enumerate(points):
            if i == len(points) - 1:
                break
            x, y = p
            x2 = int(points[i + 1][0])
            y2 = int(points[i + 1][1])
            x = int(x)
            y = int(y)
            a += (x * y2 - x2 * y)
        
        a /= 2
        
        if a == 0:
            a = 1
        
        # Centroid X
        # cx = 1/6A Sum(x_i + x_i+1)*(x_i * y_i+1 - x_i+1 * y_i)
        cx = 0
        
        for i, p in enumerate(points):
            if i == len(points) - 1:
                break
            x, y = p
            x2 = int(points[i + 1][0])
            y2 = int(points[i + 1][1])
            x = int(x)
            y = int(y)
            
            cx += (x + x2) * (x * y2 - x2 * y)
        
        cx /= 6 * a
        
        # Centroid Y
        # cy = 1/6A Sum(y_i + y_i+1)*(x_i * y_i+1 - x_i+1 * y_i)
        cy = 0
        
        for i, p in enumerate(points):
            if i == len(points) - 1:
                break
            x, y = p
            x2 = int(points[i + 1][0])
            y2 = int(points[i + 1][1])
            x = int(x)
            y = int(y)
            
            cy += (y + y2) * (x * y2 - x2 * y)
        
        cy /= 6 * a
        
        area.set("points", " ".join([",".join(x) for x in points]))
        
        #print points
        #print str(cx) + " " + str(cy)
        
        #if title == "2.3.2.-, 2.3.2.-, R05207":
        #    import pdb; pdb.set_trace()

        try:
            # Find pixel with highest luminance in 8-Neighbourhood (Moore)
            # 123
            # 456
            # 789
            
            colors = []
            colors.append(idata[(cy - 1) * iwidth + (cx - 1)]) # 1
            colors.append(idata[(cy - 1) * iwidth + cx]) # 2
            colors.append(idata[(cy - 1) * iwidth + (cx + 1)]) # 3
            colors.append(idata[cy * iwidth + (cx - 1)]) # 4
            colors.append(idata[cy * iwidth + cx]) # 5
            colors.append(idata[cy * iwidth + (cx + 1)]) # 6
            colors.append(idata[(cy + 1) * iwidth + (cx - 1)]) # 7
            colors.append(idata[(cy + 1) * iwidth + cx]) # 8
            colors.append(idata[(cy + 1) * iwidth + (cx + 1)]) # 9

            lum = map(lambda x: luminance(x), colors)
            best_color = colors[lum.index(min(lum))]
            
            area.set("fill", "rgb(" + ",".join(str(x) for x in best_color[:3]) + ")")
                    
        except IndexError:
            print(str(cx) + " " + str(cy))

for elem in pending:
    graphics.append(elem)

out = StringIO()
out.write('<html><head><script type="text/javascript" src="jquery-1.8.0.min.js"></script><script type="text/javascript" src="jquery-svgpan.js"></script></head><body>')
tree.write(out)
out.write('<script type="text/javascript">$("svg").svgPan("viewport");</script></body></html>')

with open("kegg_out.html", "w") as kegg:
    kegg.write(out.getvalue())
