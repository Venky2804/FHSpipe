#!/usr/bin/python3
#To run fastp, hisat2, stringtie, merge on the given fastq file
import os,shutil
import subprocess
import argparse, sys
import json, yaml

##########     PYTHON VERSION CHECK     ##########
if sys.version_info < (3,5,0):
    sys.exit("Required python version newer than 3.5.0")

##########     FILE PATHS     ##########
realpath = os.path.realpath(__file__)
basename = os.path.basename(__file__)
srcdir = realpath.replace("/"+basename, "")

##########     PARSER ARGUMENTS     ##########
parser = argparse.ArgumentParser(description="A program to run FHS pipeline i.e. FASTp, HISAT2, STRINGTIE, MERGE, GFFCOMPARE, STRINGTIE 2, PREPDE, extract potential lncRNAs (IUX)")
required = parser.add_argument_group("required arguments")
required.add_argument("-p", "--path", help="Path of sample fastq files", action="store")
required.add_argument("-s", "--spcs", help="Reference species name. Ex: Gallus_gallus,Bos_taurus", action="store")
required.add_argument("-e", "--seq", help="Sequencing type - S (single end) or P (paired end)", action="store")
parser.add_argument("-r", "--ref", help="Path of reference genome fasta file if species is not available", action="store")
parser.add_argument("-g", "--gff", help="Path of reference genome gtf/gff file if species is not available", action="store")
parser.add_argument("--len", help="Length filter cut off [Default: 200]", type=int, action="store", default="200")
parser.add_argument("--orfl", help="ORF Length filter cut off [Default: 300]", type=int, action="store", default="300")
parser.add_argument("-o", "--out", help="Prefix of output files [Default: out]", action="store", default="out")
parser.add_argument("-t", "--threads", type=int, help="Number of threads to use in Hisat2, Samtools, Stringtie. [Default: 1]", action="store", default="1")
parser.add_argument("--lnc", help="Use to only extract lncRNAs and skip file processing for differential expression", action="store_true")
groups_order = {'positional arguments': 0, 'required arguments': 1, 'optional arguments': 2 }
parser._action_groups.sort(key=lambda g: groups_order[g.title])
args = parser.parse_args()

##########     CHECK FOR INPUT, DATABASE     ##########
if args.path == None:
    print("\nSample file path not provided. Check below for help\n")
    print(parser.print_help())
    sys.exit(1)

A=[]
for files in sorted(os.listdir(args.path)):
    if files.endswith(".fastq"):
        A.append(files)
if len(A) > 1:
    print(str(len(A))+" FASTQ sample files were given as input")
else:
    print("Sample files at "+args.path+" are not in FASTQ format. Please check.")
    sys.exit(1)

if args.spcs == None:
    print("\nSpecies name not provided. Please use \"-h\" or \"--help\" for help\n")
    sys.exit(1)
else:
    spcs = srcdir+"/DB/"+args.spcs

if args.seq == None:
    print("Sequencing type not provided. Please use \"-h\" or \"--help\" for help\n")
    sys.exit(1)
elif not args.seq == "S" and not args.seq == "P":
    print("Error recognising sequencing type. Please use \"S\" for single end and \"P\" for paired end.")
    print(args.seq)
    sys.exit(1)

pfam = srcdir+"/DB/pfamdb/Pfam"
pfamdir = srcdir+"/DB/pfamdb"
if not os.path.exists(pfamdir):
    print("\nPFAM database not found in the tool directory. Please check the path", pfamdir)
    sys.exit(1)

##########     TOOLS FROM CONFIG AND TEST     ##########
with open(srcdir+"/config.yaml", "r") as c:
    try:
        conf = yaml.safe_load(c)
    except:
        print("Error loading config\n", yaml.YAMLError)
        sys.exit(1)

fastp = conf["fastp"]
hisat2 = conf["hisat2"]
samtools = conf["samtools"]
stringtie = conf["stringtie"]
gffcompare = conf["gffcompare"]
prepDE = conf["prepDE"]
bedtools = conf["bedtools"]
OrfPredictor = conf["OrfPredictor"]
CPC2 = conf["CPC2"]
rpsblast = conf["rpsblast"]

tools = (fastp, hisat2, samtools, stringtie, gffcompare, prepDE, bedtools, OrfPredictor, CPC2, rpsblast)
A = (fastp, hisat2, samtools, gffcompare, bedtools)
B = (stringtie, rpsblast)
C = (prepDE, CPC2)
D = (OrfPredictor)

keys = list(conf.keys())
values = list(conf.values())
errors = []
E = []
for tool in tools:
    if not os.path.exists(tool):
        errors.append(keys[values.index(tool)])
    if tool in A:
        cmd = tool+" --help"
        cmd = cmd.split()
        x = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if x.returncode > 0:
            E.append(keys[values.index(tool)])
    elif tool in B:
        cmd = tool+" -h"
        cmd = cmd.split()
        x = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if x.returncode > 0:
            E.append(keys[values.index(tool)])
    elif tool in C:
        cmd = "python "+tool+" -h"
        cmd = cmd.split()
        x = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if x.returncode > 0:
            E.append(keys[values.index(tool)])
    elif tool in D:
        cmd = "perl "+tool+" -h"
        cmd = cmd.split()
        x = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if x.returncode > 0:
            E.append(keys[values.index(tool)])

if len(errors) > 0:
  print("Unable to find the program(s) ->", ", ".join(errors), "\nPlease check the path in config.yaml file at\n"+srcdir+"/config.yaml")
  sys.exit(1)
else:
  print("All the tools are found at the provided paths")

if len(E) > 0:
  print("But the program(s) ->", ",".join(E), "are not working properly. Please check and reinstall")
  sys.exit(1)
else:
  print("And are working properly")

##########     PROGRAM STARTS HERE     ##########
os.chdir(args.path)
os.chdir("..")
fld = os.getcwd()

log = open(fld+"/fhs_glm_pipeline.log", "a+")
from datetime import datetime
time = datetime.now().strftime("%A %-d %B, %Y %I:%M:%S %p")
print("FHS Pipeline started at:", time+"\n")
print("FHS Pipeline started at:", time+"\n", file = log, flush=True)

t = str(args.threads)
folder = args.path #path of all samples
print(len([name for name in os.listdir(folder) if name.endswith(".fastq")]), "FASTQ sample files were given as input", file=log)
done = open(fld+"/pipeline_done.txt", "a+")
fin = open(fld+"/pipeline_done.txt").read().splitlines()
d = []
for line in fin:
    d.append(line)

##########     BUILDING DATABASE     ##########
BUILD = 0
if not os.path.exists(spcs):
    print("\nSpecies not found in Database")
    print("\nSpecies not found in Database", file = log, flush=True)
    if args.ref == None:
        print("Provide reference genome file to build Database. Please use \"-r\" or \"--ref\" option\n")
        sys.exit(1)
    else:
        print("Building reference database using provided genome file")
        print("Building reference database using provided genome file", file = log, flush=True)
        BUILD += 1
else:
    print("\nSpecies found in the Database. Using reference files from the Database.")
    print("\nSpecies found in the Database. Using reference files from the Database.", file = log, flush=True)
    rfa = spcs+"/"+args.spcs+".fasta"
    rgff = spcs+"/"+args.spcs+".gff"
    hs2db = spcs+"/hs2/"+args.spcs
    chrsize = spcs+"/seq_extract/"+args.spcs+".fasta.chromosome_sizes"

if BUILD == 1:
    if os.path.exists(spcs):
        if os.listdir(spcs):
            shutil.rmtree(spcs)
            os.mkdir(spcs)
    else:
        os.mkdir(spcs)
    shutil.copy(args.ref, spcs+"/"+args.spcs+".fasta")
    shutil.copy(args.gff, spcs+"/"+args.spcs+".gff")
    rfa = spcs+"/"+args.spcs+".fasta"
    rgff = spcs+"/"+args.spcs+".gff"
    #hisat2 build
    if not os.path.exists(spcs+"/hs2"):
        os.mkdir(spcs+"/hs2")
    bld = hisat2+"-build -q -p "+t+" "+rfa+" "+spcs+"/hs2/"+args.spcs
    bld2 = bld.split()
    b1 = subprocess.run(bld2)
    if not b1.returncode == 0:
        print("Error building Hisat2 index. Open new terminal and run command\n"+bld)
        input("\n\nPress any key to continue after running command succefully")
        #hisat2 inspect
        insp = hisat2+"-inspect -s "+spcs+"/hs2/"+args.spcs
        insp = insp.split()
        ins1 = subprocess.run(insp)
        if not ins1.returncode == 0:
            print("Error!!! Index building failed. Exiting")
            shutil.rmtree(spcs)
            sys.exit(1)
    #seq extract files chromosome_sizes
    if not os.path.exists(spcs+"/seq_extract"):
        os.mkdir(spcs+"/seq_extract")
    idx = samtools+" faidx "+spcs+"/"+args.spcs+".fasta"
    idx1 = idx.split()
    idxc = subprocess.run(idx1)
    if not idxc.returncode == 0:
        print("Error fasta index. Open new terminal and run command\n"+idx)
        input("\n\nPress any key to continue after running command succefully")
        if not os.path.isfile(spcs+"/"+args.spcs+".fasta.fai"):
            print("Error!!! Index building failed. Exiting")
            shutil.rmtree(spcs)
            sys.exit(1)
    shutil.copy(spcs+"/"+args.spcs+".fasta", spcs+"/seq_extract")
    shutil.copy(spcs+"/"+args.spcs+".fasta.fai", spcs+"/seq_extract")
    chrsze = open(spcs+"/seq_extract/"+args.spcs+".fasta.chromosome_sizes", "w+")
    cut = "cut -f 1,2 "+spcs+"/"+args.spcs+".fasta.fai"
    cut = cut.split()
    cutcmd = subprocess.run(cut, stdout=chrsze)
    if not cutcmd.returncode == 0:
        print("Error fasta index. Open new terminal and run command\n"+cut)
        sys.exit(1)
    chrsze.close()
    chrsize = spcs+"/seq_extract/"+args.spcs+".fasta.chromosome_sizes"
    hs2db = spcs+"/hs2/"+args.spcs
    print("Database built succesfully. Use the species name \""+args.spcs+"\" with -s option while analysing next set of sample datasets of the same organism.\n")
    print("Database built succesfully. Use the species name \""+args.spcs+"\" with -s option while analysing next set of sample datasets of the same organism.\n", file = log, flush=True)

##########     FILES AND DIRECTORY GENERATION     ##########
stm = open(fld+"/stm_gtf.txt", "a+")
fold = ["qc", "hs2", "stringtie", "fhs_results", "gffcom", "stringtie_2"]
for r in fold:
    if not os.path.exists(fld+"/"+r):
        os.mkdir(fld+"/"+r)

fpres = open("fhs_results/fastp_res.txt", "a+")
if "fpres" not in d :
    print("fpres", file = done, flush=True)
    print("Sample name\tReads passed filter\tLow quality reads\tReads with too many N\tToo short reads\tToo long reads\tQ30 bases\tGC content", file = fpres)
hsres = open("fhs_results/hisat_res.txt", "a+")
if "hsres" not in d :
    print("hsres", file = done, flush=True)
    print("Sample name\tMapping percentage", file = hsres)
st1 = open("fhs_results/st1_temp.txt", "a+")
os.chdir(folder)
folder = os.getcwd() #to get full path

for f in sorted(os.listdir(folder)):
    if args.seq == "S":
        x = f.split(".")[0]
    elif args.seq == "P":
        x = f.split("_")[0]
    if x not in d and f.endswith(".fastq"):
        print("Working with", str(x), end = "\t", flush=True)
        print("Working with", str(x), end = "\t", file = log, flush=True)
################            F A S T P            ####################
        if args.seq == "S":
            fp = fastp+" -h ../qc/"+x+".html -j ../qc/"+x+".json -a auto -w 16 -e 25 -l 30 -p -i "+x+".fastq -o ../qc/"+x+".fastq"
        elif args.seq == "P":
            fp = fastp+" --detect_adapter_for_pe --correction -h ../qc/"+x+".html -j ../qc/"+x+".json -a auto -w 16 -e 25 -l 30 -p -i "+x+"_1.fastq -I "+x+"_2.fastq -o ../qc/"+x+"_1.fastq -O ../qc/"+x+"_2.fastq"
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
        ###########################
        print("Fastp--", end = "", flush=True)
        print("Fastp--", end = "", file = log, flush=True)
        jfile.close()
################           H I S A T 2           ####################
        if args.seq == "S":
            hs = hisat2+" -q -p "+t+" --dta -x "+hs2db+" -U ../qc/"+x+".fastq -S ../hs2/"+x+".sam --summary-file ../hs2/"+x+"_align_summary"
        elif args.seq == "P":
            hs = hisat2+" -q -p "+t+" --dta -x "+hs2db+" -1 ../qc/"+x+"_1.fastq -2 ../qc/"+x+"_2.fastq -S ../hs2/"+x+".sam --summary-file ../hs2/"+x+"_align_summary"
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
        ###########################
        print("Hisat--", end = "", flush=True)
        print("Hisat--", end = "", file = log, flush=True)
################          S A M 2 B A M          ####################
        if not os.path.exists("../stringtie/"+x):
            os.mkdir("../stringtie/"+x)
        srt = samtools+" sort ../hs2/"+x+".sam -@ "+t+" -O BAM -o ../stringtie/"+x+"/"+x+".bam"
        srt1 = srt.split()
        s3 = subprocess.run(srt1, stderr=subprocess.DEVNULL)
        if not s3.returncode == 0:
            print("Error Sorting SAM file of sample", x, "Check manually with command\n"+srt)
            sys.exit(1)
        os.remove("../hs2/"+x+".sam")
        print("SAM2BAM--", end = "", flush=True)
        print("SAM2BAM--", end = "", file = log, flush=True)
################       S T R I N G T I E 1       ####################
        os.chdir("../stringtie")
        st1a = stringtie+" "+x+"/"+x+".bam -G "+rgff+" -o "+x+"/"+x+".gtf -p "+t+" -A "+x+"/"+x+".abund"
        st1b = st1a.split()
        s4 = subprocess.run(st1b, stderr=subprocess.DEVNULL)
        if not s4.returncode == 0:
            print("Error running Stringtie 1 on sample", x, "Check manually with command\n"+st1a)
            sys.exit(1)
        #####     RESULTS     #####
        abun = open(x+"/"+x+".abund").read().splitlines()
        cnt = -1
        for lin in abun:
            cnt += 1
        print(x+"\t"+str(cnt), file=st1, flush=True)
        ###########################
        print(x, file = done, flush=True)
        d.append(x)
        os.chdir(folder)
        print("Stringtie--DONE", flush=True)
        print("Stringtie--DONE", file = log, flush=True)
        stm.write(fld+"/stringtie/"+x+"/"+x+".gtf\n")
    elif x in d and f.endswith(".fastq") and not f.endswith("_2.fastq"):
        print(x, "analysis already done")

print("GTF list file for st merge -->\tstm_gtf.txt")
fpres.close()
hsres.close()
stm.close()
st1.close()
os.chdir(fld)
################   S T R I N G T I E - M E R G E   ##################
if not args.lnc:
    if "st-merge" not in d:
        print("Merging samples with Stringtie --merge")
        stmrg = stringtie+" --merge -G "+rgff+" -o merge.gtf stm_gtf.txt"
        stmrg1 = stmrg.split()
        s5 = subprocess.run(stmrg1)
        if not s5.returncode == 0:
            print("Error running Stringtie merge. Check manually with command\n"+stmrg)
            sys.exit(1)
        print("st-merge", file = done, flush=True)
        print("Stringtie merge DONE", flush=True)
        print("Stringtie merge DONE", file = log, flush=True)
    else:
        print("Stringtie merge DONE already", flush=True)
################       G F F C O M P A R E       ####################
os.chdir("gffcom")
if "gffcom" not in d:
    print("Classifying samples with GFFCompare")
    gcom = gffcompare+" -T -i "+fld+"/stm_gtf.txt -r "+rgff+" -s "+rfa
    gcom1 = gcom.split()
    s6 = subprocess.run(gcom1, stderr=subprocess.DEVNULL)
    if not s6.returncode == 0:
        print("Error running GFFcompare. Check manually with command\n"+gcom)
        sys.exit(1)
    print("gffcom", file = done, flush=True)
    print("GFFCompare DONE", flush=True)
    print("GFFCompare DONE", file = log, flush=True)
else:
    print("GFFCompare DONE already", flush=True)
os.chdir(fld)
################   S T R I N G T I E - 2   ##################
if not args.lnc:
    pdef = open("pde_gtf.txt", "a+")
    stres = open("fhs_results/stringtie_res.txt", "a+")
    if "stres" not in d :
        print("stres", file = done, flush=True)
        print("Sample name\tNumber of transcripts before merging\tNumber of transcripts after merging", file = stres)

    ST = open("fhs_results/st1_temp.txt").read().splitlines()

    for x in sorted(os.listdir(fld+"/stringtie")):
        if x+"_st2" not in d:
            print("Reassembling", str(x), end = "\t", flush=True)
            print("Reassembling", str(x), end = "\t", file=log, flush=True)
            if not os.path.exists("stringtie_2/"+x):
                os.mkdir("stringtie_2/"+x)
            st2a = stringtie+" stringtie/"+x+"/"+x+".bam -G merge.gtf -o stringtie_2/"+x+"/"+x+"_new.gtf -p "+t+" -e -A stringtie_2/"+x+"/"+x+"_new.abund"
            st2b = st2a.split()
            s7 = subprocess.run(st2b)
            if not s7.returncode == 0:
                print("Error running Stringtie 2 on sample", x, "Check manually with command\n"+st2a)
                sys.exit(1)
            #####     R E S U L T S     #####
            abun_n = open("stringtie_2/"+x+"/"+x+"_new.abund").read().splitlines()
            cn = -1
            for lins in abun_n:
                cn += 1
            for lines in ST:
                if lines.split("\t")[0] == x:
                    print(lines+"\t"+str(cn), file = stres)
            #################################
            print(x+"\t"+fld+"/stringtie_2/"+x+"/"+x+"_new.gtf", file = pdef)
            print(x+"_st2", file = done, flush=True)
            print("DONE", flush=True)
            print("DONE", file = log, flush=True)
        else:
            print(x, "Stringtie 2 DONE already", flush=True)

    stres.close()
    pdef.close()
os.remove("fhs_results/st1_temp.txt")
################       P R E P D E . P Y       ####################
if not args.lnc:
    if "pde" not in d:
        print("Generating read count matrix for DEGs", end= "\t", flush=True)
        pdec = "python "+prepDE+" -i pde_gtf.txt"
        pdes = pdec.split()
        s8 = subprocess.run(pdes)
        if not s8.returncode == 0:
            print("Unable to run prepDE.py. Check manually with command\n"+pdec)
        print("pde", file = done, flush=True)
        print("PrepDE (genes) DONE", file = log, flush=True)
        print("DONE")
    else:
        print("PrepDE already DONE")
###############     I U X - E X T R A C T I O N     ###############
if not os.path.exists("gffcom/extract_iux"):
    os.mkdir("gffcom/extract_iux")
os.chdir("gffcom/extract_iux")
if "IUX extraction" not in d:
    print("Extracting I,U,X sequences")
    print("Extracting I,U,X sequences", file = log, flush=True)
    cc = [" class_code \"i\"", " class_code \"u\"", " class_code \"x\""]
    gtf = open("../gffcmp.combined.gtf").read().splitlines() #gffcom merge gtf
    xt = []
    bed = open(args.out+".bed", "w+")
    for g in gtf:
        t1 = g.split(";")
        for T in t1:
            if T.startswith(" class_code"):
                tt = T
        if tt in cc:
            w = g.split("\t")[0]+'\t'+g.split("\t")[3]+'\t'+g.split("\t")[4]+'\t'+g.split(";")[0].split("\t")[-1].split("\"")[-2]
            if w not in xt:
                xt.append(w)
                print(w, file = bed)
    print("Number of coordinates in bed file are:\t", len(xt), file = log, flush = True)
    print("Number of coordinates in bed file are:\t", len(xt))
    bed.close()
    slop = bedtools+" slop -i "+args.out+".bed -g "+chrsize+" -b 50"
    slop1 = slop.split()
    sbed = open("slop.bed", "w+")
    xs1 = subprocess.run(slop1, stdout = sbed)
    if not xs1.returncode == 0:
        print("Error running bedtools slop.\nThere must be a problem with input files.\nCheck manually with\n"+slop)
        sys.exit(1)
    sbed.close()
    getfa = bedtools+" getfasta -fi "+rfa+" -name+ -bed slop.bed -fo "+args.out+".fasta"
    getfa1 = getfa.split()
    xs2 = subprocess.run(getfa1)
    if not xs2.returncode == 0:
        print("Error running bedtools getfasta.\nThere must be a problem with input files.\nCheck manually with\n"+getfa)
        sys.exit(1)
    print("IUX extraction", file = done, flush=True)
else:
    print("IUX sequences already extracted")
###############     L N C - F I L T E R - S T E P S     ###############
if "len-filter" not in d:
    print("Length filter - "+str(args.len), file = log, flush=True)
    print("Length filter - "+str(args.len))
    inseq = open(args.out+".fasta").read().splitlines() #Extracted IUX fasta file
    total = 0
    for line in inseq:
        if line.startswith(">"):
            total +=1
    print("Total IUX sequences are", total)
    print("Total IUX sequences are", total, file = log)
    ##########     LENGTH FILTER     ##########
    print("Filtering sequences with Length less than", args.len)
    delen = open("delete_lenfilter", "w+") #file with headers of filtered sequences
    lf = open(args.out+"_lenfil.fasta", "w+") #len filtered fasta
    count = 0
    for lin in range(len(inseq)):
        if inseq[lin].startswith(">"):
            if len(inseq[lin+1]) > args.len:
                print(inseq[lin], file = lf)
                print(inseq[lin+1], file = lf)
                count += 1
            else:
                print(inseq[lin], file = delen)
    delen.close()
    lf.close()
    print("Sequences retained after length filtering are", count)
    print("Sequences retained after length filtering are", count, file = log)
    print("len-filter", file = done, flush=True)
else:
    print("Length filteration already done")
    ##########     ORF PREDICTOR     ##########
if "orf-filter" not in d:
    print("ORF filter - "+str(args.orfl), file = log, flush=True)
    print("ORF filter - "+str(args.orfl))
    print("Predicting ORFs of the sequences using OrfPredictor")
    open("empty", "w+") #empty file for OrfPred
    cmd1a = "perl "+OrfPredictor+" "+args.out+"_lenfil.fasta empty 0 both 1e-3 "+args.out+"_orfout.fasta"
    cmd1 = cmd1a.split()
    orf = subprocess.run(cmd1)
    if not orf.returncode == 0:
        os.remove("empty")
        print("Error running OrfPredictor.\nThere must be a problem with input files.\nCheck manually with\n"+cmd1a)
        print("Error running OrfPredictor.\nThere must be a problem with input files.\nCheck manually with\n"+cmd1a, file = log)
        sys.exit(1)
    os.remove("empty")
    ##########     ORF DNA FILTER     ##########
    print("Filtering sequences with ORFs less than", args.orfl, "length")
    lenfil = open(args.out+"_lenfil.fasta").read().splitlines() #Reads Len fil file
    orf6 =  open("ORF6frame.txt").readlines()
    dna = [] #DNA header index list
    pro = [] #ORF header index list
    discard = [] #ORF headers to be discarded
    delete = [] #Contigs to be discarded
    for line in range(len(lenfil)):
        if lenfil[line].startswith(">"):
            dna.append(line)
    for index in range(len(orf6)):
        if orf6[index].startswith(">"):
            pro.append(index)
    pro.append(len(orf6))
    for c in range(len(pro)-1):
        for dd in range(pro[c]+1,pro[c+1]):
            for e in orf6[dd].split("*"):
                if len(e) > args.orfl:
                    discard.append(orf6[pro[c]].strip())
                    break
    for f in discard:
        for g in f.split("\t"):
            if g.startswith(">"):
                if g not in delete:
                    delete.append(g)
    print("Sequences retained after ORF filtering:", (len(dna)-len(delete)))
    print("Sequences retained after ORF filtering:", (len(dna)-len(delete)), file = log)
    orffil = open(args.out+"_orffil.fasta", "w+") #orf filtered dna fasta
    for w in dna:
        if lenfil[w] not in delete:
            print(lenfil[w], file=orffil)
            print(lenfil[w+1], file=orffil)
    dele2 = open("delete_orffilter", "w+")
    for x in delete:
        print(x, file = dele2)
    orffil.close()
    dele2.close()
    ##########     ORF PROTEIN FILTER     ##########
    orfout = open(args.out+"_orfout.fasta").read().splitlines()
    dele2 = open("delete_orffilter").read().splitlines()
    orfptn = open(args.out+"_ptn_orffil.fasta", "w+") #orf filtered ptn fasta
    dis = []
    for line in dele2:
        dis.append(line)
    for line in range(len(orfout)):
        new = orfout[line].split("\t")
        for x in new:
            if x.startswith(">") and x not in dis:
                print(orfout[line], file = orfptn)
                print(orfout[line+1], file = orfptn)
    orfptn.close()
    print("orf-filter", file = done, flush=True)
else:
    print("ORF filteration already done")
    ##########     RPSBLAST (PFAM)     ##########
if "RPS-filter" not in d:
    print("Running RPSBlast against Pfam database")
    cmd2a = rpsblast+" -db "+pfam+" -query "+args.out+"_ptn_orffil.fasta -out "+args.out+"_pfam_rps -evalue 0.001 -outfmt 7"
    cmd2 = cmd2a.split()
    rps = subprocess.run(cmd2)
    if not rps.returncode == 0:
        print("Error running RPSBLAST (PFAM).\nThere must be a problem with input files.\nCheck manually with\n"+cmd2a)
        print("Error running RPSBLAST (PFAM).\nThere must be a problem with input files.\nCheck manually with\n"+cmd2a, file = log)
        sys.exit(1)
    rpsfil = open(args.out+"_rpsfil.fasta", "w+") #RPS (PFAM) filtered DNA fasta
    rpsout = open(args.out+"_pfam_rps").read().splitlines() #RPS output file
    norf = open("noOrf.txt").read().splitlines() #No orf file from OrfPredictor.pl
    orffil = open(args.out+"_orffil.fasta").read().splitlines() #ORF filtered DNA sequences
    RPS = []
    for index in range(len(rpsout)):
        if rpsout[index].startswith("# Query"):
            if rpsout[index+2].startswith("# 0 hits"):
                new = rpsout[index].split(": ")
                for x in new:
                    if not x.startswith("#"):
                        RPS.append(">"+x)
    print("Number of sequences with no hits against Pfam database (RPS) are:", len(RPS))
    print("Number of sequences with no hits against Pfam database (RPS) are:", len(RPS), file = log)
    for nof in norf:
        if not nof.startswith("Sequence"):
            RPS.append(">"+nof.split("\t")[0])
    for shind in range(len(orffil)):
        if orffil[shind].startswith(">"):
            if orffil[shind] in RPS:
                print(orffil[shind], file=rpsfil)
                print(orffil[shind+1], file=rpsfil)
    rpsfil.close()
    print("RPS-filter", file = done, flush=True)
else:
    print("RPS filteration already done")
    ##########     CPC2     ##########
if "CPC-filter" not in d:
    print("Running Coding Potential Calculator 2")
    cmd3a = "python3 "+CPC2+" -i "+args.out+"_rpsfil.fasta -o "+args.out+"_cpc"
    cmd3 = cmd3a.split()
    cpc = subprocess.run(cmd3, stderr = subprocess.DEVNULL)
    if not cpc.returncode == 0:
        print("Error running CPC2.\nThere must be a problem with input files.\nCheck manually with\n"+cmd3a)
        sys.exit(1)
    cpcout = open(args.out+"_cpc.txt").read().splitlines()
    rpsfil = open(args.out+"_rpsfil.fasta").read().splitlines()
    final = open(args.out+"_final_lnc.fasta", "w+")
    F = []
    for lines in cpcout:
        if lines.split("\t")[7] == "noncoding":
            F.append(">"+lines.split("\t")[0])
    print("Total number of final lncRNA sequences after filtering:", len(F))
    print("Total number of final lncRNA sequences after filtering:", len(F),file = log)
    for hind in range(len(rpsfil)):
        if rpsfil[hind].startswith(">"):
            if rpsfil[hind] in F:
                print(rpsfil[hind], file=final)
                print(rpsfil[hind+1], file=final)
    final.close()
    print("CPC-filter", file = done, flush=True)
else:
    print("CPC filteration already done")
if not os.path.exists("delnc"):
    os.mkdir("delnc")
###############     IUX gffcom merge gtf     ###############
if not args.lnc:
    if "IUX-merge" not in d:
        print("Preparing IUX merge gtf file", end="\t", flush=True)
        gff = open(fld+"/gffcom/gffcmp.combined.gtf").read().splitlines() #gffcom merge gtf
        flnc = open(args.out+"_final_lnc.fasta").read().splitlines()
        head = [] #final lnc header
        fout = open("delnc/iux_merge.gtf", "w+")
        for line in flnc:
            if line.startswith(">"):
                head.append(line)
        for lin in gff:
            ggg = lin.split("\t")[8].split("\"")[1]
            for line in head:
                tcons = line.split(">")[1].split("::")[0]
                if tcons == ggg:
                    print(lin, file = fout)
                    break
        fout.close()
        print("DONE", flush=True)
        print("IUX-merge", file = done, flush=True)
    else:
        print("IUX merge gtf prepared already")
    os.chdir("delnc")
###############     L N C - S T R I N G T I E 2     ###############
    pdef_lnc = open("pde_lnc.txt", "a+")
    lnc_stres = open("stringtie_res.txt", "a+")
    if "lnc_stres" not in d :
        print("lnc_stres", file = done, flush=True)
        print("Sample name\tNumber of transcripts", file = lnc_stres)
    if not os.path.exists("stringtie"):
        os.mkdir("stringtie")
    for x in sorted(os.listdir(fld+"/stringtie")):
        if x+"_lnc_st" not in d:
            print("Assembling", str(x), "for lnc", end = "\t", flush=True)
            print("Assembling", str(x), "for lnc", end = "\t", file=log, flush=True)
            if not os.path.exists("stringtie/"+x):
                os.mkdir("stringtie/"+x)
            lsta = stringtie+" "+fld+"/stringtie/"+x+"/"+x+".bam -G iux_merge.gtf -o stringtie/"+x+"/"+x+"_new.gtf -p "+t+" -e -A stringtie/"+x+"/"+x+"_new.abund"
            lstb = lsta.split()
            lst = subprocess.run(lstb)
            if not lst.returncode == 0:
                print("Error running Stringtie (lnc) on sample", x, "Check manually with command\n"+lsta)
                sys.exit(1)
            #####     R E S U L T S     #####
            abun_n = open("stringtie/"+x+"/"+x+"_new.abund").read().splitlines()
            cn = -1
            for lins in abun_n:
                cn += 1
            print(x+"\t"+str(cn), file = lnc_stres)
            #################################
            print(x+"\t"+fld+"/gffcom/extract_iux/delnc/stringtie/"+x+"/"+x+"_new.gtf", file = pdef_lnc)
            print(x+"_lnc_st", file = done, flush=True)
            print("DONE", flush=True)
            print("DONE", file = log, flush=True)
        else:
            print(x, "Stringtie (lnc) DONE already", flush=True)
    lnc_stres.close()
    pdef_lnc.close()
################       L N C - P R E P D E . P Y       ####################
    if "pde_lnc" not in d:
        print("Generating read count matrix for DElncs", end="\t", flush=True)
        lpdec = "python "+prepDE+" -i pde_lnc.txt"
        lpdes = lpdec.split()
        s8 = subprocess.run(lpdes)
        if not s8.returncode == 0:
            print("Unable to run prepDE.py. Check manually with command\n"+lpdec)
            sys.exit(1)
        print("pde_lnc", file = done, flush=True)
        print("PrepDE (lncRNAs) DONE", file = log, flush=True)
        print("DONE")
    else:
        print("PrepDE already DONE")

done.close()
end = datetime.now().strftime("%A %-d %B, %Y %I:%M:%S %p")

print("\nFHS Pipeline finished at:", end, file = log, flush=True)
print(
'''
=================================================
|| Run edgeR (GLM) for Differential expression ||
=================================================
'''
, file = log, flush=True)
print("The file \""+fld+"/gffcom/extract_iux/"+args.out+"_final_lnc.fasta\" contains final potential lncRNAs", file = log, flush=True)
print("Use", fld+"/gene_count_matrix.csv for DEGs", file = log, flush=True)
print("Use", fld+"/gffcom/extract_iux/delnc/transcript_count_matrix.csv for DElncRNAs", file = log, flush=True)
log.close()
print("\nFHS Pipeline finished at:", end)
print(
'''
=================================================
|| Run edgeR (GLM) for Differential expression ||
=================================================
'''
)
print("The file \""+fld+"/gffcom/extract_iux/"+args.out+"_final_lnc.fasta\" contains final potential lncRNAs")
print("Use", fld+"/gene_count_matrix.csv for DEGs")
print("Use", fld+"/gffcom/extract_iux/delnc/transcript_count_matrix.csv for DElncRNAs")