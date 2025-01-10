#!/usr/bin/python3
#To test tools given in config file
import os, yaml, sys, subprocess

##########     CHECK FOR INPUT, DATABASE     ##########
def check(path, parser, spcs, srcdir, seq):
    if path == None:
        print("\nSample file path not provided. Check below for help\n")
        print(parser.print_help())
        sys.exit(1)

    A=[]
    for files in sorted(os.listdir(path)):
        if files.endswith(".fastq"):
            A.append(files)
    if len(A) > 1:
        print(str(len(A))+" FASTQ sample files were given as input")
    else:
        print("Sample files at "+path+" are not in FASTQ format. Please check.")
        sys.exit(1)

    if spcs == None:
        print("\nSpecies name not provided. Please use \"-h\" or \"--help\" for help\n")
        sys.exit(1)
    else:
        spath = srcdir+"/DB/"+spcs

    if seq == None:
        print("Sequencing type not provided. Please use \"-h\" or \"--help\" for help\n")
        sys.exit(1)
    elif not seq == "S" and not seq == "P":
        print("Error recognising sequencing type. Please use \"S\" for single end and \"P\" for paired end.")
        print(seq)
        sys.exit(1)

    pfam = srcdir+"/DB/pfamdb/Pfam"
    pfamdir = srcdir+"/DB/pfamdb"
    if os.path.exists(pfamdir):
        if not os.listdir(pfamdir):
          print("\nPFAM database not found in the DB directory. Downloading PFam database")
          import util.download_pfam as dl
          dl.run(pfamdir)
    else:
      print("\nPFAM database not found in the DB directory. Downloading PFam database")
      import util.download_pfam as dl
      dl.run(pfamdir)
    return(spath, pfam)

##########     TOOLS FROM CONFIG     ##########
def test(config):
  with open(config, "r") as c:
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
      continue
    if tool in A:
      cmd = tool+" --help"
      cmd = cmd.split()
      x = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
      if x.returncode > 0:
        E.append(keys[values.index(tool)])
    if tool in B:
      cmd = tool+" -h"
      cmd = cmd.split()
      x = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
      if x.returncode > 0:
        E.append(keys[values.index(tool)])
    if tool in C:
      cmd = "python "+tool+" -h"
      cmd = cmd.split()
      x = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
      if x.returncode > 0:
        E.append(keys[values.index(tool)])
    if tool in D:
      cmd = "perl "+tool+" -h"
      cmd = cmd.split()
      x = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
      if x.returncode > 0:
        E.append(keys[values.index(tool)])

  if len(errors) > 0:
    print("Unable to find the program(s) ->", ", ".join(errors), "\nPlease check the path in config.yaml file at\n", config)
    sys.exit(1)

  if len(E) > 0:
    print("The program(s) ->", ",".join(E), "are not working properly. Please check and reinstall")
    sys.exit(1)
  return(fastp, hisat2, samtools, stringtie, gffcompare, prepDE, bedtools, OrfPredictor, CPC2, rpsblast, conf)
