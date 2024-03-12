#!/usr/bin/python3
#program to extract sequences from genome
import os
import subprocess
from datetime import datetime
import argparse, sys

parser = argparse.ArgumentParser(description="A program to extract sequences from genome")
parser.add_argument("-i", "--input", help="File with gene list from edgeR. [Format: mstrg id+tab+gene name]", action="store")
parser.add_argument("-g", "--gtf", help="Stringtie merge gtf file", action="store")
parser.add_argument("-r", "--ref", help="Path of reference files (fasta, fai and chromosome_sizes files). Alternatively, use directory \"path/to/FHSpipe/DB/species/seq_extract\" after running FHSpipe.", action="store")
parser.add_argument("-o", "--output", help="Prefix of output files [Default: out]", action="store", default="out")
args = parser.parse_args()

if args.input == None or args.gtf == None:
    print("\nArgument error. Please check your command.\nUse \"-h\" or \"--help\" for help\n")
    sys.exit(1)

if not os.path.exists(args.ref+"/chicken.fasta") or not os.path.exists(args.ref+"/chicken.chromosome_sizes"):
    print("\n Reference files not found at Default location. Please use \"-r\" or \"--ref\" option\n")
    sys.exit(1)
    
time = datetime.now().strftime("%A %-d %B, %Y %I:%M:%S %p")
print("Sequence extraction started at:", time+"\n")

ref = args.ref
L = []

print("Preparing Bed File")
gene = open(args.input).read().splitlines() #edgeR gene list
gtf = open(args.gtf).read().splitlines() #merge gtf file
bed = open(args.output+".bed", "w+")
for gind in range(len(gene)):
    g = gene[gind]
    w = g.split("\t")[1] #genename
    a = g.split("\t")[0] #mstrg id
    S = [] #start
    E = [] #end
    chrm = ""
    for b in gtf:
        if w.startswith("MSTRG.") and b.startswith("N") and b.split("\t")[2] == "transcript":
            name = a #name col. in bed for fasta header
            mid = b.split("\"")[1] #MSTRG id in GFF
            if a == mid:
                s = int(b.split("\t")[3])
                e = int(b.split("\t")[4])
                S.append(s)
                E.append(e)
                if not chrm.startswith("N"):
                    chrm += b.split("\t")[0]
        if not w.startswith("MSTRG.") and b.startswith("N") and b.split("\t")[2] == "transcript" and "gene_name" in b:
            name = w
            mid = b.split("\"")[1]
            ge = b.split("\"")[5]
            if a == mid and w == ge:
                s = int(b.split("\t")[3])
                e = int(b.split("\t")[4])
                S.append(s)
                E.append(e)
                if not chrm.startswith("N"):
                    chrm += b.split("\t")[0]
    print(chrm+"\t"+str(min(S))+"\t"+str(max(E))+"\t"+name, file = bed)
    if w == "-":
        L.append(a+"\t"+a)
    else:
        L.append(w+"\t"+a)
bed.close()
print("Extracting the Sequences")
sbed = open(args.output+".slop.bed", "w+")
slop = "bedtools slop -i "+args.output+".bed -g "+ref+"/chicken.chromosome_sizes -b 50"
print("\tslop -->\t"+slop)
slp = slop.split()
s1 = subprocess.run(slp, stdout = sbed)
if not s1.returncode == 0:
    print("Error running Bedtools slop.\nCheck manually with command\n"+slop)
    sys.exit(1)
getfa = "bedtools getfasta -fi "+ref+"/chicken.fasta -bed "+args.output+".slop.bed -name+ -fo "+args.output+".fasta"
print("\tgetfa -->\t"+getfa)
gfa = getfa.split()
s2 = subprocess.run(gfa)
if not s2.returncode == 0:
    print("Error extracting sequences.\nCheck manually with command\n"+getfa)
    sys.exit(1)
print("Preparing the List File")
lst = open(args.output+".list", "w+")
slp = open(args.output+".slop.bed").read().splitlines()
for lno in range(len(L)):
    cord = slp[lno].split("\t")
    L[lno] += "\t"+cord[0]+":"+cord[1]+"-"+cord[2]
for line in L:
    print(line, file = lst)
lst.close()

end = datetime.now().strftime("%A %-d %B, %Y %I:%M:%S %p")
print("Sequence extraction finished at:", end)

print("Output is written to "+args.output+".fasta")
