#!/usr/bin/python3
#To run sam2bam

import subprocess, sys

def run(samtools, x, t):
    srt = samtools+" sort ../hs2/"+x+".sam -@ "+t+" -O BAM -o ../stringtie/"+x+"/"+x+".bam"
    srt1 = srt.split()
    s3 = subprocess.run(srt1, stderr=subprocess.DEVNULL)
    if not s3.returncode == 0:
        print("Error Sorting SAM file of sample", x, "Check manually with command\n"+srt)
        sys.exit(1)