# FHSpipe
Transcriptome analysis pipeline using Fastp, Hisat2 and Stringtie

## R E Q U I R E M E N T S
Required python version newer than 3.5.0 

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
- All Sample fastq files in a single directory (no extra files)  



## U S A G E
 
```bash
fhs_glm_pipeline.py [-h] [-p PATH] [-s SPCS] [-r REF] [-g GFF] [--len LEN] [--orfl ORFL] [-o OUT] [-t THREADS] [--lnc]  
```
> A program to run FHS pipeline i.e. FASTp, HISAT2, STRINGTIE, MERGE, GFFCOMPARE, STRINGTIE 2, PREPDE, extract potential lncRNAs (IUX)  
>
> ### Required arguments:
> -p PATH, --path PATH --> Path of sample fastq files  
> -s SPCS, --spcs SPCS --> Reference species name  
>  
> ### Optional arguments:  
> -h, --help --> show this help message and exit  
> -r REF, --ref REF --> Path of reference genome fasta file if species is not available  
> -g GFF, --gff GFF --> Path of reference genome gtf/gff file if species is not available  
> --len LEN --> Length filter cut off [Default: 200]  
> --orfl ORFL --> ORF Length filter cut off [Default: 300]  
> -o OUT, --out OUT --> Prefix of output files [Default: out]  
> -t THREADS, --threads THREADS --> Number of threads to use in Hisat2, Samtools, Stringtie. [Default: 1]  
> --lnc --> Use to only extract lncRNAs and skip file processing for differential expression  



## OUTPUT FILES

Check the file "Output_guide.txt" in pipeline directory for details of output of the pipeline  
