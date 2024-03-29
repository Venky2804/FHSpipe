#!/usr/bin/python3
#To extract lncRNA iux merge

def run(fld, flnc, fout):
    print("Preparing IUX merge gtf file", end="\t", flush=True)
    gff = open(fld+"/gffcom/gffcmp.combined.gtf").read().splitlines() #gffcom merge gtf
    head = [] #final lnc header
    for line in flnc:
        if line.startswith(">"):
            head.append(line.split(">")[1].split("::")[0])
    for lin in gff:
        ggg = lin.split("\t")[8].split("\"")[1]
        if ggg in head:
            print(lin, file = fout)