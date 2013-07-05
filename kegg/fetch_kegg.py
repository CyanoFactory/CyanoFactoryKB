import re
import urllib2
response = urllib2.urlopen('http://www.kegg.jp/kegg-bin/download_htext?htext=br08901.keg&format=htext&filedir=')
html = response.read()

# Find all map ids
m = re.findall(r'(?<=\s)([0-9]{5})(?=\s)', html)
mlen = len(m)

for i, match in enumerate(m):
    print "Downloading {} ({}/{})".format(match, i + 1, mlen)
    image = "map{}.png".format(match)
    site = "map{}.html".format(match)
    
    response = urllib2.urlopen('http://www.kegg.jp/kegg-bin/show_pathway?map' + match)
    
    print_map = 0
    with open("fetch/" + site, "w") as kegg:
        for line in response:
            if line.startswith("<img"):
                last_image = line.strip()
            
            if line.startswith("<map"):
                download_image = last_image
                print_map = 1
            
            if print_map:
                # fix up to valid XML
                line = re.sub(r"coords=([^\"\s]+)", r'coords="\1"', line.strip())
                line = re.sub(r"shape=([^\"\s]+)", r'shape="\1"', line)
                # escape non-escaped ampersands
                line = re.sub(r"&(?!(?:apos|quot|[gl]t|amp);|#)", r'&amp;', line)
                kegg.write(line + "\n")
            
            if line.endswith("</map>"):
                print_map = 0

    # The image appears before the <map>-tag
    # 99% of the images are available under
    # http://www.genome.jp/kegg/pathway/map/ but we also want
    # to handle 1% of the corner cases correctly
    with open("fetch/" + image, "wb") as imgout:
        download_image = re.search(r'(?<=src\=")(.*?)(?=")', download_image).group(1)
        imgdl = urllib2.urlopen('http://www.kegg.jp' + download_image)
        imgout.write(imgdl.read())
