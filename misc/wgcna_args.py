#!/usr/bin/python3
#code to output FPKMs of all genes in all samples for wgcna
'''
File formats:
input = list of genes (gene names / mstrg ids) or lncRNAs (TCONS ids) [Gene/MSTRG.##/TCONS.#####]
sam = list of paths of stringtie 2 gtf files [/path/to/stringtie/sample/new/sample_new.gtf] same as pde input file column 2
'''

from datetime import datetime
import argparse, sys
from tqdm import tqdm

parser = argparse.ArgumentParser(description="A program to output FPKMs of all genes in all samples for wgcna")
parser.add_argument("-g", "--genin", help="File with list of genes (gene names).", action="store")
parser.add_argument("-l", "--lncin", help="File with list of lncRNA (lnc TCONS).", action="store")
parser.add_argument("-s", "--gsam", help="File with list of paths of stringtie2 (-e) gtf files of DEGs (must end with_new.gtf)[Get from pde.txt col 1]", action="store")
parser.add_argument("-S", "--lncsam", help="File with list of paths of stringtie2 (-e) gtf files of DElncRNAs (must end with_new.gtf)[Get from pde.txt col 1]", action="store")
parser.add_argument("-o", "--output", help="Output file name prefix [Default: out]", action="store", default="out")
args = parser.parse_args()

if args.genin == None or args.lncin == None or args.gsam == None or args.lncsam == None:
	print("\nError. Input files not provided.\nPlease check your command.\nUse \"-h\" or \"--help\" for help\n")
	sys.exit(1)

time = datetime.now().strftime("%A %-d %B, %Y %I:%M:%S %p")
print("FPKM extraction started at:", time+"\n")

out = open(args.output+"_wgcna.csv", "w+")
ann = open(args.output+"_wgcna_ann.csv", "w+")
sam = open(args.gsam).read().splitlines()
pre = "gene_id"
print("gene_id,gene_name", file = ann)
for x in range(len(sam)):
	pre += ","+sam[x].split("/")[-1].split("_new")[0]
print(pre, end = "", file = out)

def wgcna(infile, sample, des):
	inpt = open(infile).read().splitlines()
	sam = open(sample).read().splitlines()
	for gind in tqdm(range(len(inpt)), desc = des):
		gene = inpt[gind]
		gid = gind+1
		print("\n"+des+"_"+str(gid), end = "", file = out)
		print(des+"_"+str(gid)+","+gene, file = ann)
		for samind in range(len(sam)):
			sample = sam[samind]
			gtf = open(sample).read().splitlines()
			A = []
			for line in gtf:
				if gene in line and "FPKM" in line:
					if gene == line.split(';')[-5].split("\"")[-2] or gene == line.split(';')[-6].split("\"")[-2]:
							A.append(line.split(";")[-3].split("\"")[-2])
			print(","+str(max(A)), end = "", file = out)
	print("")

wgcna(args.genin, args.gsam, "Gene")
wgcna(args.lncin, args.lncsam, "lncRNA")

out.close()
ann.close()

end = datetime.now().strftime("%A %-d %B, %Y %I:%M:%S %p")
print("FPKM extraction finished at:", end)

print("Output is written to", args.output+"_wgcna.csv")
