import subprocess, os

# CIRCOS binary
CIRCOS = "circos"

CIRCOS_TEMPLATE = """
karyotype = circos/karyotype.human.hg19.txt

<ideogram>
<spacing>
default = 0.005r
</spacing>
radius    = 0.8r
thickness = 50p
fill      = yes
show_label       = yes
label_radius     = dims(ideogram,radius) + 0.02r
label_size       = 30
label_parallel   = yes
show_bands = yes
fill_bands = yes
</ideogram>

<image>
<<include image.conf>>
#dir = %(imagedir)s
#file = %(imagefn)s
</image>
<<include colors_fonts_patterns.conf>>
<<include housekeeping.conf>>

<links>
<link>
file           = %(linkfn)s
radius         = 0.78r
bezier_radius  = 0r
color          = black
thickness      = 3

</link>
</links>

<plots>
<plot>
type  = histogram
file  = %(histogramfn)s
r1    = 0.99r
r0    = 0.80r
max   = 1
min   = 0
extend_bin = no
fill_color = white
background_color = white
</plot>

<plot>
type  = text
file  = %(textfn)s
r1    = 1.3r
r0    = 1.05r
label_font = normal
label_size = 18p
rpadding = 5p
</plot>
</plots>

"""

def conv_chr(x):
    if "|" in x:
        x = x.replace("|", "_")
    if x.startswith("chr"):
        x = x.lstrip("chr")
    if not x.startswith("hs"):
        x = "hs%s" % (x)
    return x

class Link:
    def __init__(self, chr1, pos1, chr2, pos2, color = "black", thickness = 1, label = "", z = None):
        self.chr1 = chr1
        self.pos1 = pos1
        self.chr2 = chr2
        self.pos2 = pos2
        self.color = color
        self.thickness = thickness
        self.label = label
        self.z = z
    def __str__(self):
        attr = ["color=%s" % (self.color), 
                "thickness=%d" % (self.thickness)]
        if self.z != None:
            attr.append("z=%d" % (self.z))
        return "%s %d %d %s %d %d %s" % (conv_chr(self.chr1), self.pos1, self.pos1, conv_chr(self.chr2), self.pos2, self.pos2, ",".join(attr))

class Region:
    def __init__(self, chrom, pos1, pos2, value, color = "black", z = 1):
        self.chrom = chrom
        self.pos1 = pos1
        self.pos2 = pos2
        self.color = color
        self.value = value
        self.z = z
    def __str__(self):
        attr = ["color=%s,z=%d" % (self.color, self.z)]
        return "%s %d %d %f %s" % (conv_chr(self.chrom), self.pos1, self.pos2, self.value, ",".join(attr))

class Annotation:
    def __init__(self, chrom, pos1, pos2, label = "", color = "black", size = 1):
        self.chrom = chrom
        self.pos1 = pos1
        self.pos2 = pos2
        self.color = color
        self.size = size
        self.label = label
    def __str__(self):
        attr = ["color=%s" % (self.color), 
                "label_size=%s" % (self.size)]
        return "%s %d %d %s %s" % (conv_chr(self.chrom), self.pos1, self.pos2, self.label, ",".join(attr)) 

def __write_links(o, links):
#    for chr1, pos1, chr2, pos2, color, thickness, label in links:
#        o.write("%s %d %d %s %d %d color=%s,thickness=%dp\n" % (conv_chr(chr1), pos1, pos1, conv_chr(chr2), pos2, pos2, color, thickness))
#        o.write("%s %d %d %s %d %d thickness=%dp,url=%s\n" % (conv_chr(chr1), pos1, pos1, conv_chr(chr2), pos2, pos2, thickness, color))
    for l in links:
        o.write("%s\n" % (str(l)))

def __write_histogram(o, histogram):
#    for chr1, pos1, pos2, value, color in histogram:
#        o.write("%s %d %d %s color=%s\n" % (conv_chr(chr1), pos1, pos2, value, color))
    for r in histogram:
        o.write("%s\n" % (str(r)))

def __write_text(o, text):
#    for chrom, pos1, pos2, label, color, size in text:
#        o.write("%s %d %s %s color=%s,label_size=%s\n" % (conv_chr(chrom), pos1, pos2, label, color, size))
    for t in text:
        o.write("%s\n" % (str(t)))

def plot(ofn, links, histogram, text):
    o = open("circos.conf", "w")
    params = {"linkfn" : "links.txt", "histogramfn" : "histogram.txt", "textfn" : "text.txt", "imagefn" : ofn, "imagedir" : "."}
    o.write(CIRCOS_TEMPLATE % params)
    o.close()

    o = open("links.txt", "w")
    __write_links(o, links)
    o.close()

    o = open("histogram.txt", "w")
    __write_histogram(o, histogram)
    o.close()

    o = open("text.txt", "w")
    __write_text(o, text)
    o.close()

    subprocess.call("%s -conf circos.conf -outputdir %s -outputfile %s" % (CIRCOS, os.path.dirname(ofn), os.path.basename(ofn)), shell = True)

if __name__ == "__main__":
    links = [("hs1", 100000, "hs2", 200000), ("hs4", 100000, "hs5", 200000), ("hs6", 100000, "hsX", 200000)]
    histogram = []
    text = []
    plot("test.png", links, histogram, text)
