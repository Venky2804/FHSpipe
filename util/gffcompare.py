#!/usr/bin/python3
#To run gffcompare

import subprocess, sys

def run(gffcompare, fld, rgff, rfa):
    print("Anootating the transcripts with classcodes", end="", flush=True)
    gcom = gffcompare+" -T -i "+fld+"/stm_gtf.txt -r "+rgff+" -s "+rfa
    gcom1 = gcom.split()
    s6 = subprocess.run(gcom1, stderr=subprocess.DEVNULL)
    if not s6.returncode == 0:
        print("Error running GFFcompare. Check manually with command\n"+gcom)
        sys.exit(1)
    print(" -- DONE", flush=True)
    log = "Classcode annotation DONE"
    return(log)