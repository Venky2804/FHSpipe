#!/usr/bin/python3
#program to extract 5p sequences of degs from genome
import argparse, sys
import subprocess
import os

##########     PYTHON VERSION CHECK     ##########
if sys.version_info < (3,5,0):
    print("Required python version newer than 3.5.0")
    sys.exit(1)

##########     PARSER ARGUMENTS     ##########
parser = argparse.ArgumentParser(description="A program to extract 5p upstream sequences of DEGs from genome")
parser.add_argument("-i", "--input", help="BED file generated in DEG extraction", action="store")
parser.add_argument("-r", "--ref", help="Path of reference files for extraction (chr sizes, ref fasta and fai). [Default: db/gallus_gallus/seq_extract]", action="store", default="/home/haitoshailesh/Niabdata/venkat/db/gallus_gallus/seq_extract")
parser.add_argument("-o", "--output", help="Prefix of output files [Default: out]", action="store", default="out")
args = parser.parse_args()

##########     CHECK FOR INPUT, REFERENCE AND PROGRAM     ##########
if args.input == None:
	print("\nInput BED file not provided.\nUse \"-h\" or \"--help\" for help\n")
	sys.exit(1)

if not os.path.exists(args.ref+"/chicken.fasta") or not os.path.exists(args.ref+"/chicken.chromosome_sizes"):
	print("\n Reference files not found at Default location. Please use \"-r\" or \"--ref\" option\n")
	sys.exit(1)

rc = subprocess.run(["which", "bedtools"], stdout=subprocess.DEVNULL)
if not rc.returncode == 0:
	print("Unable to find bedtools. Please check the path")
	sys.exit(1)
        
##########     PROGRAM STARTS HERE     ##########
from datetime import datetime

time = datetime.now().strftime("%A %-d %B, %Y %I:%M:%S %p")
print("Sequence extraction started at:", time+"\n")

############STEPS FOR MAKING BED FILE############################################################################################
slop = "bedtools slop -i "+args.input+" -g "+args.ref+"/chicken.chromosome_sizes -l 5000 -r 0"
slop1 = slop.split()
sbed = open("slop.txt", "w+")
s1 = subprocess.run(slop1, stdout = sbed)
if not s1.returncode == 0:
	print("Error running bedtools slop.\nThere must be a problem with input files.\nCheck manually with\n"+slop)
	sys.exit(1)
############CHANGING BED FILE COORDINATES########################################################################################
fp = open("slop.txt").read().splitlines()
tf = open(args.output+"_5p.bed", "w+")
inp = open(args.input).read().splitlines()
print("Generating BED file\n")
for i in inp:
	for j in fp:
		if i.split("\t")[0] == j.split("\t")[0] and i.split("\t")[2] == j.split("\t")[2] and i.split("\t")[3] == j.split("\t")[3]:
			chrm = j.split("\t")[0] #chromosome id from step 1 txt
			strt = j.split("\t")[1] #start from step 1 txt (-5000 of input deg.bed start)
			end = i.split("\t")[1] #start of input deg.bed
			name = i.split("\t")[3] #name of DEG in input deg.bed
			tf.write(chrm+"\t"+strt+"\t"+end+"\t"+name+"\n")
			break
os.remove("slop.txt")
tf.close()
############EXTRACTION OF SEQUENCES##############################################################################################
print("Extracting sequences\n")
getfa = "bedtools getfasta -fi "+args.ref+"/chicken.fasta -name+ -bed "+args.output+"_5p.bed -fo "+args.output+".fasta"
getfa1 = getfa.split()
s2 = subprocess.run(getfa1)
if not s2.returncode == 0:
	print("Error running bedtools getfasta.\nThere must be a problem with input files.\nCheck manually with\n"+getfa)
	sys.exit(1)

end = datetime.now().strftime("%A %-d %B, %Y %I:%M:%S %p")
print("Sequence extraction finished at:", end+"\n")
