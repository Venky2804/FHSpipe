#!/usr/bin/python3
#To download PFam database
import os, subprocess,sys

def run(pfamdir):
    l="https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Pfam_LE.tar.gz"
    if not os.path.exists(pfamdir):
        os.mkdir(pfamdir)
    w = "wget "+l+" -P "+pfamdir+" -nv --show-progress"
    w = w.split()
    a = subprocess.run(w)
    if not a.returncode == 0:
        sys.exit(1)
    os.chdir(pfamdir)
    z = "tar -xvzf "+pfamdir+"/Pfam_LE.tar.gz -C "+pfamdir
    z = z.split()
    subprocess.run(z, stdout=subprocess.DEVNULL)
    print("Downloaded PFam database")
