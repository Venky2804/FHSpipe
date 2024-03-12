#!/usr/bin/python3
import os
import argparse, sys

##########     PYTHON VERSION CHECK     ##########
if sys.version_info < (3,5,0):
    print("Required python version newer than 3.5.0")

##########     PARSER ARGUMENTS     ##########
parser = argparse.ArgumentParser(description="A program to extract maping percentages from results of HISAT2 as a table")
required = parser.add_argument_group("required arguments")
required.add_argument("-p", "--path", help="Path of hisat summary files", action="store")
parser.add_argument("-o", "--out", help="Prefix of output files [Default: out]", action="store", default="out")
groups_order = {'positional arguments': 0, 'required arguments': 1, 'optional arguments': 2 }
parser._action_groups.sort(key=lambda g: groups_order[g.title])
args = parser.parse_args()

##########     CHECK FOR INPUT, DATABASE AND PROGRAMS     ##########
if args.path == None:
    print("\nHisat summary files path not provided. Click -h or --help for help\n")
    sys.exit(1)

##########     PROGRAM STARTS HERE     ##########
hsres = open(args.out+".txt", "w+")
print("Sample name\tMapping percentage", file = hsres)

for x in os.listdir(args.path):
    if x.endswith("_align_summary"):
        x = x.split("_align_summary")[0]
        aln = open(args.path+"/"+x+"_align_summary").read().splitlines()
        for lins in aln:
            if "overall" in lins:
                print(x+"\t"+lins.split()[0], file = hsres)

hsres.close()
print("DONE")
