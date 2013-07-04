import re
import urllib2
response = urllib2.urlopen('http://www.kegg.jp/kegg-bin/download_htext?htext=br08901.keg&format=htext&filedir=')
html = response.read()

m = re.findall(r'(?<=\s)([0-9]{5})(?=\s)', html)
mlen = len(m)

for i, match in enumerate(m):
    response = urllib2.urlopen('http://www.kegg.jp/kegg/pathway/map/map' + match + ".png")
    
    print "Downloading Item {} of {}".format(i + 1, mlen)
    image = "map" + match + ".png"
    site = "map" + match + ".html"
    
    print "Downloading " + image
    with open("fetch/" + image, "wb") as imgout:
        imgout.write(response.read())
    
    response = urllib2.urlopen('http://www.kegg.jp/kegg-bin/show_pathway?map' + match)
    
    print_map = 0
    print "Downloading " + site
    with open("fetch/" + site, "w") as kegg:
        for line in response:
            if line[:4] == "<map":
                print_map = 1
            
            if print_map:
                # Fixup to valid XML
                line = re.sub(r"coords=([^\"\s]+)", r'coords="\1"', line.strip())
                line = re.sub(r"shape=([^\"\s]+)", r'shape="\1"', line.strip())
                kegg.write(line.strip() + "\n")
            
            if line[:6] == "</map>":
                print_map = 0
