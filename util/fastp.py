#!/usr/bin/python3
#To run fastp
import subprocess, sys, json

def run(seq, fastp, x, fpres):
    if seq == "S":
        fp = fastp+" -h ../qc/"+x+".html -j ../qc/"+x+".json -a auto -w 16 -e 25 -l 30 -p -i "+x+".fastq -o ../qc/"+x+".fastq"
    elif seq == "P":
        fp = fastp+" --detect_adapter_for_pe --correction -h ../qc/"+x+".html -j ../qc/"+x+".json -a auto -w 16 -e 25 -l 30 -p -i "+x+"_R1.fastq -I "+x+"_R2.fastq -o ../qc/"+x+"_R1.fastq -O ../qc/"+x+"_R2.fastq"
    fp1 = fp.split()
    s1 = subprocess.run(fp1, stderr=subprocess.DEVNULL)
    if not s1.returncode == 0:
        print("Error running Fastp on sample", x, "Check manually with command\n"+fp)
        sys.exit(1)
    #####     RESULTS     #####
    jfile = open("../qc/"+x+".json")
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
    q30 = jdict["summary"]["after_filtering"]["q30_rate"]
    q30 = str(round(q30*100, 4))
    gc = jdict["summary"]["after_filtering"]["gc_content"]
    gc = str(round(gc*100, 4))
    print(x+"\t"+tot+"\t"+lowqual+"\t"+tooN+"\t"+short+"\t"+toolong+"\t"+q30+"\t"+gc, file = fpres)
    jfile.close()
