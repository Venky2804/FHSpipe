#!/usr/bin/python3
#To run hisat2
import subprocess, sys

def run(seq, hisat2, t, hs2db, x, hsres):
    if seq == "S":
            hs = hisat2+" -q -p "+t+" --dta -x "+hs2db+" -U ../qc/"+x+".fastq -S ../hs2/"+x+".sam --summary-file ../hs2/"+x+"_align_summary"
    elif seq == "P":
        hs = hisat2+" -q -p "+t+" --dta -x "+hs2db+" -1 ../qc/"+x+"_R1.fastq -2 ../qc/"+x+"_R2.fastq -S ../hs2/"+x+".sam --summary-file ../hs2/"+x+"_align_summary"
    hs1 = hs.split()
    s2 = subprocess.run(hs1, stderr=subprocess.DEVNULL)
    if not s2.returncode == 0:
        print("Error running Hisat2 on sample", x, "Check manually with command\n"+hs)
        sys.exit(1)
    #####     RESULTS     #####
    aln = open("../hs2/"+x+"_align_summary").read().splitlines()
    for lins in aln:
        if "overall" in lins:
            print(x+"\t"+lins.split()[0], file = hsres)
