#!/usr/bin/python3
import os
import argparse, sys

##########     PYTHON VERSION CHECK     ##########
if sys.version_info < (3,5,0):
    print("Required python version newer than 3.5.0")

##########     PARSER ARGUMENTS     ##########
parser = argparse.ArgumentParser(description="A program to extract transcript abundance from results of Stringtie as a table")
required = parser.add_argument_group("required arguments")
required.add_argument("-p1", "--path1", help="Path of stringtie 1 results folders", action="store")
required.add_argument("-p2", "--path2", help="Path of stringtie 2 results folders", action="store")
parser.add_argument("-o", "--out", help="Prefix of output files [Default: out]", action="store", default="out")
groups_order = {'positional arguments': 0, 'required arguments': 1, 'optional arguments': 2 }
parser._action_groups.sort(key=lambda g: groups_order[g.title])
args = parser.parse_args()

##########     CHECK FOR INPUT, DATABASE AND PROGRAMS     ##########
if args.path1 == None:
    print("\nStringtie 1 results path not provided. Click -h or --help for help\n")
    sys.exit(1)

if args.path2 == None:
    print("\nStringtie 2 results path not provided. Click -h or --help for help\n")
    sys.exit(1)

##########     PROGRAM STARTS HERE     ##########
stres = open(args.out+".txt", "w+")
print("Sample name\tNumber of transcripts before merging\tNumber of transcripts after merging", file = stres)

ST = []

for x in os.listdir(args.path1):
    abun = open(args.path1+"/"+x+"/"+x+".abund").read().splitlines()
    cnt = -1
    for lin in abun:
        cnt += 1
    ST.append(x+"\t"+str(cnt))

for x in os.listdir(args.path2):
    abun_n = open(args.path2+"/"+x+"/"+x+"_new.abund").read().splitlines()
    cn = -1
    for lins in abun_n:
        cn += 1
    for lines in ST:
        if lines.split("\t")[0] == x:
            print(lines+"\t"+str(cn), file = stres)


stres.close()
print("DONE")
