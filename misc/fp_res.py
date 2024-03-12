#!/usr/bin/python3
import os, json
import argparse, sys

##########     PYTHON VERSION CHECK     ##########
if sys.version_info < (3,5,0):
    print("Required python version newer than 3.5.0")

##########     PARSER ARGUMENTS     ##########
parser = argparse.ArgumentParser(description="A program to extract results of FASTp as a table")
parser.add_argument("-p", "--path", help="Path of fastp json files", action="store")
parser.add_argument("-o", "--out", help="Prefix of output files [Default: out]", action="store", default="out")
args = parser.parse_args()

##########     CHECK FOR INPUT, DATABASE AND PROGRAMS     ##########
if args.path == None:
    print("\nFastp json files path not provided. Click -h or --help for help\n")
    sys.exit(1)

##########     PROGRAM STARTS HERE     ##########
fpres = open(args.out+".txt", "w+")
#print("Sample name\tReads passed filter\tLow quality reads\tReads with too many N\tToo short reads\tToo long reads\tQ30 bases\tGC content", file = fpres)
print("Sample name\tReads passed filter\tLow quality reads\tReads with too many N\tToo short reads\tToo long reads\tQ30 base content\tGC content", file = fpres)

for x in os.listdir(args.path):
    if x.endswith(".json"):
        x = x.split(".")[0]
        print(x)
        jfile = open(args.path+"/"+x+".json")
        jdict = json.load(jfile)
        bfr = jdict["summary"]["before_filtering"]["total_reads"]
        aft = jdict["summary"]["after_filtering"]["total_reads"]
        tot = str(round(aft/bfr*100, 4))
        lowqual = jdict["filtering_result"]["low_quality_reads"]
        lowqual = str(round(lowqual/bfr*100, 4))
        tooN = jdict["filtering_result"]["too_many_N_reads"]
        tooN = str(round(tooN/bfr*100, 4))
        short = jdict["filtering_result"]["too_short_reads"]
        short = str(round(short/bfr*100, 4))
        toolong = jdict["filtering_result"]["too_long_reads"]
        toolong = str(round(toolong/bfr*100, 4))
        q30a = jdict["summary"]["after_filtering"]["q30_rate"]
        q30a = str(round(q30a*100, 4))
        gc = jdict["summary"]["after_filtering"]["gc_content"]
        gc = str(round(gc*100, 4))
        print(x+"\t"+tot+"\t"+lowqual+"\t"+tooN+"\t"+short+"\t"+toolong+"\t"+q30a+"\t"+gc, file = fpres)

fpres.close()

print("DONE")
