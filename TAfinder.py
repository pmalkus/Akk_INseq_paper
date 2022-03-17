## This script finds and reports the position of all TA motifs in a genome
## input:  python TAfinder.py argv[1] > output_file.txt
## argv[1] is the fasta sequence of the organism downloaded form NCBI
## The output file has 2 columns:  1) misc_feature; 2) genome position of the TA site on the Forward (+) strand

import sys

fasta={}
fasta[1]=[]

for line in open(sys.argv[1]):
    split=line.split('\t')
    if split[0][0]!=">":
        read=split[0].rstrip()
        fasta[1].append(read)

genome="".join(fasta[1])

TAs={}
TAposition=-1
length=len(genome)
newposition=0

for x in range(0,length):
    position=genome[newposition:].find("TA")
    TAposition = position + TAposition + 1
    TAs[int(TAposition)+1]=[]
    newposition = TAposition + 1

TA_keys=TAs.keys()
TA_keys.sort()

for i in TA_keys:
    print "misc_feature\t%d" %(i)
    
    




    


    
        
