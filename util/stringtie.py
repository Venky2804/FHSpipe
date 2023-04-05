#!/usr/bin/python3
#To run stringtie1

import subprocess, sys, os

def st1(stringtie, x, rgff, t):
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
    st1temp = x+"\t"+str(cnt)
    print(" -- DONE", flush=True)
    log = " -- DONE"
    return(st1temp, log)

def st2(x, stringtie, t, ST, abn, stres):
    print("Reassembling", str(x), end = "", flush=True)
    if not os.path.exists("stringtie_2/"+x):
        os.mkdir("stringtie_2/"+x)
    st2a = stringtie+" stringtie/"+x+"/"+x+".bam -G merge.gtf -o stringtie_2/"+x+"/"+x+"_new.gtf -p "+t+" -e -A stringtie_2/"+x+"/"+x+"_new.abund"
    st2b = st2a.split()
    s7 = subprocess.run(st2b)
    if not s7.returncode == 0:
        print("Error running Stringtie 2 on sample", x, "Check manually with command\n"+st2a)
        sys.exit(1)
    #####     R E S U L T S     #####
    abun_n = open(abn).read().splitlines()
    cn = -1
    for lins in abun_n:
        cn += 1
    for lines in ST:
        if lines.split("\t")[0] == x:
            print(lines+"\t"+str(cn), file = stres)
    #################################
    print(" -- DONE", flush=True)
    log = "Reassembling "+str(x)+"\t"+"DONE"
    return(log)

def lncst(x, stringtie, fld, t, lnc_stres, pdef_lnc):
    print("Assembling", str(x), "for lnc", end = "", flush=True)
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
    print(" -- DONE", flush=True)
    log = "DONE"
    return(log)

def merge(stringtie, rgff):
    print("Merging the assembled samples", end = "", flush = True)
    stmrg = stringtie+" --merge -G "+rgff+" -o merge.gtf stm_gtf.txt"
    stmrg1 = stmrg.split()
    s5 = subprocess.run(stmrg1)
    if not s5.returncode == 0:
        print("Error running Stringtie merge. Check manually with command\n"+stmrg)
        sys.exit(1)
    print(" -- DONE", flush=True)
    log = "Stringtie merge DONE"
    return(log)