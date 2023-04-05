#!/usr/bin/python3
#HS2 database build script
import os, subprocess, shutil, sys

def check(spath, log, ref):
    BUILD = 0
    if not os.path.exists(spath):
        print("Species not found in Database")
        print("Species not found in Database", file = log, flush=True)
        if ref == None:
            print("Provide reference genome file to build Database. Please use \"-r\" or \"--ref\" option\n")
            sys.exit(1)
        else:
            print("Building reference database using provided genome file")
            print("Building reference database using provided genome file", file = log, flush=True)
            BUILD += 1
    else:
        print("\nSpecies found in the Database. Using reference files from the Database.")
        print("\nSpecies found in the Database. Using reference files from the Database.", file = log, flush=True)
    return(BUILD)


def build(spath, hisat2, rfa, t, spcs, samtools):
    #hisat2 build
    if not os.path.exists(spath+"/hs2"):
        os.mkdir(spath+"/hs2")
    bld = hisat2+"-build -q -p "+t+" "+rfa+" "+spath+"/hs2/"+spcs
    bld2 = bld.split()
    b1 = subprocess.run(bld2)
    if not b1.returncode == 0:
        print("Error building Hisat2 index. Open new terminal and run command\n"+bld)
        input("\n\nPress any key to continue after running command succefully")
    #hisat2 inspect
    insp = hisat2+"-inspect -s "+spath+"/hs2/"+spcs
    insp = insp.split()
    ins1 = subprocess.run(insp, stdout=subprocess.DEVNULL)
    if not ins1.returncode == 0:
        print("Error!!! Index building failed. Exiting")
        shutil.rmtree(spath)
        sys.exit(1)
    #seq extract files chromosome_sizes
    if not os.path.exists(spath+"/seq_extract"):
        os.mkdir(spath+"/seq_extract")
    idx = samtools+" faidx "+spath+"/"+spcs+".fasta"
    idx1 = idx.split()
    idxc = subprocess.run(idx1)
    if not idxc.returncode == 0:
        print("Error fasta index. Open new terminal and run command\n"+idx)
        input("\n\nPress any key to continue after running command succefully")
        if not os.path.isfile(spath+"/"+spcs+".fasta.fai"):
            print("Error!!! Index building failed. Exiting")
            shutil.rmtree(spath)
            sys.exit(1)
    shutil.copy(spath+"/"+spcs+".fasta", spath+"/seq_extract")
    shutil.copy(spath+"/"+spcs+".fasta.fai", spath+"/seq_extract")
    chrsze = open(spath+"/seq_extract/"+spcs+".fasta.chromosome_sizes", "w+")
    cut = "cut -f 1,2 "+spath+"/"+spcs+".fasta.fai"
    cut = cut.split()
    cutcmd = subprocess.run(cut, stdout=chrsze)
    if not cutcmd.returncode == 0:
        print("Error fasta index. Open new terminal and run command\n"+cut)
        sys.exit(1)
    chrsze.close()
    chrsize = spath+"/seq_extract/"+spcs+".fasta.chromosome_sizes"
    hs2db = spath+"/hs2/"+spcs
    print("Database built succesfully. Use the species name \""+spcs+"\" with -s option while analysing next set of sample datasets of the same organism.\n", flush=True)
    log = "Database built succesfully. Use the species name \""+spcs+"\" with -s option while analysing next set of sample datasets of the same organism.\n"
    return(hs2db, chrsize, log)