#!/usr/bin/python3
#Main pipeline file
import os,shutil, argparse, sys
from util import db_build as build, test_tools, fastp as fp, hisat2 as hs2, sam2bam as s2b, stringtie as st, gffcompare as gffcom, pde, lnc_extract as lncxt, iux_merge as iuxmerge

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
required.add_argument("-p", "--path", help="Path of sample fastq files", action="store", required=True)
required.add_argument("-s", "--spcs", help="Reference species name. Ex: Gallus_gallus,Bos_taurus. Check DB folder at FHSpipe path", action="store", required=True)
required.add_argument("-e", "--seq", help="Sequencing type - S (single end) or P (paired end)", action="store", required=True)
optional = parser.add_argument_group("Optional arguments")
optional.add_argument("-r", "--ref", help="Path of reference genome fasta file if species is not available", action="store")
optional.add_argument("-g", "--gff", help="Path of reference genome gtf/gff file if species is not available", action="store")
optional.add_argument("--len", help="Length filter cut off (in number of nucleotides) [Default: 200]", type=int, action="store", default="200")
optional.add_argument("--orfl", help="ORF Length filter cut off (in number of amino acids) [Default: 100]", type=int, action="store", default="100")
optional.add_argument("-o", "--out", help="Prefix of output files [Default: out]", action="store", default="out")
optional.add_argument("-t", "--threads", type=int, help="Number of threads [Default: 1]", action="store", default="1")
optional.add_argument("--lnc", help="Use to only extract lncRNAs and skip file processing for differential expression", action="store_true")
args = parser.parse_args()

##########     CHECK FOR INPUT, DATABASE     ##########
check = test_tools.check(args.path, parser, args.spcs, srcdir, args.seq)
spath = check[0]
pfam = check[1]
##########     TOOLS FROM CONFIG AND TEST     ##########
config = test_tools.test(srcdir+"/config.yaml")
fastp = config[0]
hisat2 = config[1]
samtools = config[2]
stringtie = config[3]
gffcompare = config[4]
prepDE = config[5]
bedtools = config[6]
OrfPredictor = config[7]
CPC2 = config[8]
rpsblast = config[9]
conf = config[10]

##########     PROGRAM STARTS HERE     ##########
if not args.ref == None:
    args.ref = os.path.realpath(args.ref)
if not args.gff == None:
    args.gff = os.path.realpath(args.gff)
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
numsamp = len([name for name in os.listdir(folder) if name.endswith(".fastq")]) 
if numsamp > 0 :
    print(numsamp, "FASTQ sample files were given as input", file=log)
else:
    print("Input sample folder empty or files not in fastq format. Unzip the files if in zip format")
    sys.exit(0)
done = open(fld+"/pipeline_done.txt", "a+")
fin = open(fld+"/pipeline_done.txt").read().splitlines()
d = []
for line in fin:
    d.append(line)

##########     BUILDING DATABASE     ##########
rfa = spath+"/"+args.spcs+".fasta"
rgff = spath+"/"+args.spcs+".gff"
s0 = build.check(spath, log, args.ref)
BUILD = s0
if BUILD == 0:
    hs2db = spath+"/hs2/"+args.spcs
    chrsize = spath+"/seq_extract/"+args.spcs+".fasta.chromosome_sizes"

elif BUILD == 1:
    if os.path.exists(spath):
        if os.listdir(spath):
            shutil.rmtree(spath)
            os.mkdir(spath)
    else:
        os.mkdir(spath)
    shutil.copy(args.ref, spath+"/"+args.spcs+".fasta")
    shutil.copy(args.gff, spath+"/"+args.spcs+".gff")
    s1 = build.build(spath, hisat2, rfa, t, args.spcs, samtools)
    hs2db = s1[0]
    chrsize = s1[1]
    print(s1[2], file=log, flush=True)

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
st1temp = open("fhs_results/st1_temp.txt", "a+")
os.chdir(folder)
folder = os.getcwd() #to get full path

for f in sorted(os.listdir(folder)):
    if args.seq == "S":
        x = f.split(".")[0]
    elif args.seq == "P":
        x = f.split("_R")[0]
    if x not in d and f.endswith(".fastq"):
        print("Assembling", str(x), end = "", flush=True)
        print("Assembling", str(x), end = "", file = log, flush=True)
################            F A S T P            ####################
        s2 = fp.run(args.seq, fastp, x, t, fpres)
################           H I S A T 2           ####################
        s3 = hs2.run(args.seq, hisat2, t, hs2db, x, hsres)
################          S A M 2 B A M          ####################
        if not os.path.exists("../stringtie/"+x):
            os.mkdir("../stringtie/"+x)
        s4 = s2b.run(samtools, x, t)
        os.remove("../hs2/"+x+".sam")
################       S T R I N G T I E 1       ####################
        os.chdir("../stringtie")
        s5 = st.st1(stringtie, x, rgff, t)
        print(x, file = done, flush=True)
        print(s5[0], file=st1temp, flush=True)
        print(s5[1], file=log, flush=True)
        os.chdir(folder)
        stm.write(fld+"/stringtie/"+x+"/"+x+".gtf\n")
        d.append(x)
    elif x in d and f.endswith(".fastq") and not f.endswith("_R2.fastq"):
        print(x, "analysis already done")

fpres.close()
hsres.close()
stm.close()
st1temp.close()
os.chdir(fld)
################   S T R I N G T I E - M E R G E   ##################
if not args.lnc:
    if "st-merge" not in d:
        s6 = st.merge(stringtie, rgff)
        print("st-merge", file = done, flush=True)
        print(s6, file=log, flush=True)
    else:
        print("Stringtie merge DONE already", flush=True)
################       G F F C O M P A R E       ####################
os.chdir("gffcom")
if "gffcom" not in d:
    s7 = gffcom.run(gffcompare, fld, rgff, rfa)
    print("gffcom", file = done, flush=True)
    print(s7, file=log, flush=True)
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
            abun_n = "stringtie_2/"+x+"/"+x+"_new.abund"
            s8 = st.st2(x, stringtie, t, ST, abun_n, stres)
            print(x+"\t"+fld+"/stringtie_2/"+x+"/"+x+"_new.gtf", file = pdef)
            print(x+"_st2", file = done, flush=True)
            print(s8, file=log, flush=True)
        else:
            print(x, "Stringtie 2 DONE already", flush=True)
    stres.close()
    pdef.close()
    os.remove("fhs_results/st1_temp.txt")
################       P R E P D E . P Y       ####################
if not args.lnc:
    if "pde_DEGs" not in d:
        print("Generating read count matrix for DEGs", end= "", flush=True)
        print("Generating read count matrix for DEGs", end= "", file =log, flush=True)
        i = "pde_gtf.txt"
        s9 = pde.run(prepDE, i)
        print("pde_DEGs", file = done, flush=True)
        print(s9, file=log, flush=True)
    else:
        print("PrepDE already DONE for DEGs")
###############     L N C - E X T R A C T I O N     ###############
if not os.path.exists("gffcom/extract_iux"):
    os.mkdir("gffcom/extract_iux")
os.chdir("gffcom/extract_iux")
if "LNC extraction" not in d:
    gtf = open("../gffcmp.combined.gtf").read().splitlines()
    s10 = lncxt.run(d, log, args.out, gtf, bedtools, chrsize, rfa, done, args.len, args.orfl, OrfPredictor, rpsblast, pfam, CPC2)
    if s10 == 1:
        print("LNC extraction", file = done, flush=True)
    else:
        print("Some error during lncRNA extraction. Rerun the pipeline using same command")
        sys.exit()
else:
    print("LNCRNA sequences already extracted")

if not os.path.exists("delnc"):
    os.mkdir("delnc")
###############     IUX gffcom merge gtf     ###############
if not args.lnc:
    if "IUX-merge" not in d:
        flnc = open(args.out+"_final_lnc.fasta").read().splitlines()
        fout = open("delnc/iux_merge.gtf", "w+")
        s11 = iuxmerge.run(fld, flnc, fout)
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
            print("Assembling", str(x), "for lnc", end = "\t", file=log, flush=True)
            s12 = st.lncst(x, stringtie, fld, t, lnc_stres, pdef_lnc)
            print(x+"_lnc_st", file = done, flush=True)
            print(s12, file = log, flush=True)
        else:
            print(x, "Stringtie (lnc) DONE already", flush=True)
    lnc_stres.close()
    pdef_lnc.close()
################       L N C - P R E P D E . P Y       ####################
    if "pde_lnc" not in d:
        print("Generating read count matrix for DElncs", end="\t", flush=True)
        print("Generating read count matrix for DElncs", end="\t", file=log, flush=True)
        i = "pde_lnc.txt"
        s13 = pde.run(prepDE, i)
        print("pde_lnc", file = done, flush=True)
        print(s13, file=log, flush=True)
    else:
        print("PrepDE already DONE for lncRNAs")

done.close()
end = datetime.now().strftime("%A %-d %B, %Y %I:%M:%S %p")
print("\nFHS Pipeline finished at:", end, file = log, flush=True)
print(
'''
=================================================
|| Run edgeR (GLM) for Differential expression ||
|| Tutorial at https://rpubs.com/bman/79395    ||
=================================================
'''
, file = log, flush=True)
print("Check the file \""+srcdir+"/Output_guide.txt\" for details of output of the pipeline", file = log, flush=True)
log.close()
print("\nFHS Pipeline finished at:", end)
print(
'''
=================================================
|| Run edgeR (GLM) for Differential expression ||
|| Tutorial at https://rpubs.com/bman/79395    ||
=================================================
'''
)
print("Check the file \""+srcdir+"/Output_guide.txt\" for details of output of the pipeline")
