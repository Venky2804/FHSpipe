#!/usr/bin/python3
#To run prepDE.py

import subprocess

def run(prepDE, i):
    pdec = "python "+prepDE+" -i "+i
    pdes = pdec.split()
    s8 = subprocess.run(pdes)
    if not s8.returncode == 0:
        print("Unable to run prepDE.py. Check manually with command\n"+pdec)
    print(" -- DONE", flush=True)
    log = " -- DONE"
    return(log)