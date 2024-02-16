# FHSpipe
Transcriptome analysis pipeline using Fastp, Hisat2 and Stringtie. Generates read count files for differential expression analysis using popular R programs like edgeR or DEseq2

![Figure representing the FHSpipe flowchart](https://github.com/Venky2804/FHSpipe/assets/110457541/af63d05b-a8ae-49cf-85c8-1b9ee63fbd24 "Figure representing the FHSpipe flowchart")
### _<div align="center"> <ins> Figure representing the FHSpipe flowchart </ins> </div>_

## R E Q U I R E M E N T S
Requires python version newer than 3.5.0 

### Python modules:
os,shutil, subprocess, argparse, sys, json, yaml 

### TOOLS NEEDED:
- Fastp (https://github.com/OpenGene/fastp)  
- Hisat2 (http://daehwankimlab.github.io/hisat2/)  
- Samtools (http://www.htslib.org/download/)  
- Stringtie (https://ccb.jhu.edu/software/stringtie)  
- PrepDE.py script (https://ccb.jhu.edu/software/stringtie/dl/prepDE.py - Place along with Stringtie)  
- Gffcompare (https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)  
- Bedtools (https://github.com/arq5x/bedtools2/releases)  
- OrfPredictor (http://proteomics.ysu.edu/publication/tools_download/OrfPredictor.zip)  
- CPC2 (http://cpc2.gao-lab.org/download.php)  
- RPSBlast (NCBI Blast+ suite) (http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)  

### FILES REQUIRED:
- Reference genome fasta file (.fasta/.fa)  
- Reference genome annotation file (.gff/.gtf)  
- All Sample fastq files in a single directory (no extra files). Single end sequenced files should be named as Sample.fastq. Pair end sequenced files should be names as Sample_R1.fastq and Sample_R2.fastq to avoid errors. Unzip the sample files if in zip format (.gz)  



## U S A G E
 
```bash
fhs_glm_pipeline.py [-h] [-p PATH] [-s SPCS] [-r REF] [-g GFF] [--len LEN] [--orfl ORFL] [-o OUT] [-t THREADS] [--lnc]  
```
> A program to run FHS pipeline i.e. FASTp, HISAT2, STRINGTIE, MERGE, GFFCOMPARE, STRINGTIE 2, PREPDE, extract potential lncRNAs (IUX)  
>
> ### Required arguments:
> -p PATH,&emsp;&emsp;--path PATH&emsp;&emsp;&emsp;&ensp;-->&ensp;Path of sample fastq files  
> -s SPCS,&emsp;&emsp;--spcs SPCS&emsp;&emsp;&emsp;&emsp;-->&ensp;Reference species name  
>  
> ### Optional arguments:  
> -h,&emsp;&emsp;&emsp;&emsp;&ensp;&nbsp;--help&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;-->&ensp;Show this help message and exit  
> -r REF,&emsp;&emsp;&emsp;--ref REF&emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;-->&ensp;Path of reference genome fasta file if species is not available  
> -g GFF,&emsp;&emsp;&ensp;&nbsp;--gff GFF&emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;-->&ensp;Path of reference genome gtf/gff file if species is not available  
> &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;--len LEN&emsp;&emsp;&emsp;&emsp;&emsp;-->&ensp;Length filter cut off (In number of nucleotides) [Default: 200]  
> &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;--orfl ORFL&emsp;&emsp;&emsp;&emsp;&nbsp;-->&ensp;ORF Length filter cut off (In number of amino acids) [Default: 100]  
> -o OUT,&emsp;&emsp;&ensp;--out OUT&emsp;&emsp;&emsp;&emsp;&ensp;&nbsp;-->&ensp;Prefix of output files [Default: out]  
> -t THREADS,&nbsp;--threads THREADS&ensp;-->&ensp;Number of threads to use in Hisat2, Samtools, Stringtie. [Default: 1]  
> &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;--lnc&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;&nbsp;-->&ensp;Use to only extract lncRNAs and skip file processing for differential expression  



## OUTPUT FILES

Check the file "Output_guide.txt" in pipeline directory for details of output of the pipeline  
