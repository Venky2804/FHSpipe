#!/usr/bin/python3
#To run lncRNA extraction steps

import os, subprocess, sys

def run(d, log, out, gtf, bedtools, chrsize, rfa, done, lent, orfl, OrfPredictor, rpsblast, pfam, CPC2):
    if "IUX extraction" not in d:
        print("Extracting I,U,X sequences")
        print("Extracting I,U,X sequences", file = log, flush=True)
        cc = [" class_code \"i\"", " class_code \"u\"", " class_code \"x\""]
        xt = []
        bed = open(out+".bed", "w+")
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
        bed.close()
        slop = bedtools+" slop -i "+out+".bed -g "+chrsize+" -b 50"
        slop1 = slop.split()
        sbed = open("slop.bed", "w+")
        xs1 = subprocess.run(slop1, stdout = sbed)
        if not xs1.returncode == 0:
            print("Error running bedtools slop.\nThere must be a problem with input files.\nCheck manually with\n"+slop)
            sys.exit(1)
        sbed.close()
        getfa = bedtools+" getfasta -fi "+rfa+" -name+ -bed slop.bed -fo "+out+".fasta"
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
        inseq = open(out+".fasta").read().splitlines() #Extracted IUX fasta file
        total = 0
        for line in inseq:
            if line.startswith(">"):
                total +=1
        print("Total IUX sequences are", total)
        print("Total IUX sequences are", total, file = log)
        ##########     LENGTH FILTER     ##########
        print("Filtering sequences with Length less than", lent, end="\t")
        delen = open("delete_lenfilter", "w+") #file with headers of filtered sequences
        lf = open(out+"_lenfil.fasta", "w+") #len filtered fasta
        count = 0
        for lin in range(len(inseq)):
            if inseq[lin].startswith(">"):
                if len(inseq[lin+1]) > lent:
                    print(inseq[lin], file = lf)
                    print(inseq[lin+1], file = lf)
                    count += 1
                else:
                    print(inseq[lin], file = delen)
        delen.close()
        lf.close()
        print("-\tRetained", count)
        print("Sequences retained after length filtering are", count, file = log)
        print("len-filter", file = done, flush=True)
    else:
        print("Length filteration already done")
        ##########     ORF PREDICTOR     ##########
    if "orf-filter" not in d:
        print("Predicting ORFs of the sequences using OrfPredictor")
        open("empty", "w+") #empty file for OrfPred
        cmd1a = "perl "+OrfPredictor+" "+out+"_lenfil.fasta empty 0 both 1e-3 "+out+"_orfout.fasta"
        cmd1 = cmd1a.split()
        orf = subprocess.run(cmd1)
        if not orf.returncode == 0:
            os.remove("empty")
            print("Error running OrfPredictor.\nThere must be a problem with input files.\nCheck manually with\n"+cmd1a)
            print("Error running OrfPredictor.\nThere must be a problem with input files.\nCheck manually with\n"+cmd1a, file = log)
            sys.exit(1)
        os.remove("empty")
        ##########     ORF DNA FILTER     ##########
        print("Filtering sequences with ORFs less than", orfl, "nucleotides", end="\t")
        orfl = int(orfl/3)
        lenfil = open(out+"_lenfil.fasta").read().splitlines() #Reads Len fil file
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
                    if len(e) > orfl:
                        discard.append(orf6[pro[c]].strip())
                        break
        for f in discard:
            for g in f.split("\t"):
                if g.startswith(">"):
                    if g not in delete:
                        delete.append(g)
        print("-\tRetained", (len(dna)-len(delete)))
        print("Sequences retained after ORF filtering:", (len(dna)-len(delete)), file = log)
        orffil = open(out+"_orffil.fasta", "w+") #orf filtered dna fasta
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
        orfout = open(out+"_orfout.fasta").read().splitlines()
        dele2 = open("delete_orffilter").read().splitlines()
        orfptn = open(out+"_ptn_orffil.fasta", "w+") #orf filtered ptn fasta
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
        print("Running RPSBlast against Pfam database", end="\t")
        cmd2a = rpsblast+" -db "+pfam+" -query "+out+"_ptn_orffil.fasta -out "+out+"_pfam_rps -evalue 0.001 -outfmt 7"
        cmd2 = cmd2a.split()
        rps = subprocess.run(cmd2)
        if not rps.returncode == 0:
            print("Error running RPSBLAST (PFAM).\nThere must be a problem with input files.\nCheck manually with\n"+cmd2a)
            print("Error running RPSBLAST (PFAM).\nThere must be a problem with input files.\nCheck manually with\n"+cmd2a, file = log)
            sys.exit(1)
        rpsfil = open(out+"_rpsfil.fasta", "w+") #RPS (PFAM) filtered DNA fasta
        rpsout = open(out+"_pfam_rps").read().splitlines() #RPS output file
        norf = open("noOrf.txt").read().splitlines() #No orf file from OrfPredictor.pl
        orffil = open(out+"_orffil.fasta").read().splitlines() #ORF filtered DNA sequences
        RPS = []
        for index in range(len(rpsout)):
            if rpsout[index].startswith("# Query"):
                if rpsout[index+2].startswith("# 0 hits"):
                    new = rpsout[index].split(": ")
                    for x in new:
                        if not x.startswith("#"):
                            RPS.append(">"+x)
        print("-\tNo hits =", len(RPS))
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
        print("Running Coding Potential Calculator 2", end="\t")
        cmd3a = "python3 "+CPC2+" -i "+out+"_rpsfil.fasta -o "+out+"_cpc"
        cmd3 = cmd3a.split()
        cpc = subprocess.run(cmd3, stderr = subprocess.DEVNULL)
        if not cpc.returncode == 0:
            print("Error running CPC2.\nThere must be a problem with input files.\nCheck manually with\n"+cmd3a)
            sys.exit(1)
        cpcout = open(out+"_cpc.txt").read().splitlines()
        rpsfil = open(out+"_rpsfil.fasta").read().splitlines()
        final = open(out+"_final_lnc.fasta", "w+")
        F = []
        for lines in cpcout:
            if lines.split("\t")[7] == "noncoding":
                F.append(">"+lines.split("\t")[0])
        print("DONE\nTotal number of final lncRNA sequences after filtering:", len(F), flush=True)
        print("Total number of final lncRNA sequences after filtering:", len(F), file = log, flush = True)
        for hind in range(len(rpsfil)):
            if rpsfil[hind].startswith(">"):
                if rpsfil[hind] in F:
                    print(rpsfil[hind], file=final)
                    print(rpsfil[hind+1], file=final)
        final.close()
        print("CPC-filter", file = done, flush=True)
    else:
        print("CPC filteration already done")
