#! /usr/bin/env python 


###############################################################
# metabolisHMM - A tool for exploring and visualizing the distribution and evolutionary histories of metabolic markers
# create-genome-phylogeny: to create a phylogeny based on ribosomal proteins
# Written by Elizabeth McDaniel emcdaniel@wisc.edu
# November 2018
# This program is free software under the GNU General Public License version 3.0
###############################################################

import glob, argparse, subprocess, os, sys, tempfile, re
from Bio import BiopythonExperimentalWarning
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SearchIO, SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import UnknownSeq, Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, Counter

# Arguments
parser = argparse.ArgumentParser(description = "Create ribosomal phylogenies using specific ribosomal markers for archaea and/or bacteria")
parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
optional = parser.add_argument_group("optional arguments")
required.add_argument('--input', metavar='INPUT', help='Directory where genomes to be screened are held')
required.add_argument('--output', metavar='OUTPUT', help='Directory to store results and intermediate files')
required.add_argument('--domain', metavar='DOMAIN', help='archaea, bacteria, all')
required.add_argument('--phylogeny', metavar='PHY', help='fastree, raxml')
optional.add_argument('--threads',metavar='THREADS',help='Optional: number of threads for calculating a tree using RAxML. This is not taken into account using Fastree')
optional.add_argument('--loci', metavar='LOCI', default='12', help='Output genomes with less than x number of loci. By default prints genomes that have less than 12 ribosomal loci markers.')

# if no arguments given, print help message
if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

# version to print
def version():
    versionFile = open('VERSION')
    return versionFile.read()
VERSION = version()

args = parser.parse_args()
GENOMEDIR = args.input
GENOMEFILES = args.input + "/*.faa"
DOMAIN = args.domain
PHYTOOL = args.phylogeny
THREADS = args.threads
OUTPUT = args.output
out_intm = OUTPUT + "/out"
out_results = OUTPUT + "/results"

# Setup
genomes = glob.glob(GENOMEFILES)
os.makedirs(out_intm)
os.makedirs(out_results)
FNULL = open(os.devnull, 'w')

# different ribosomal markers for archaea/bacteria/all
bacteria_list = ['rpL14','rpL15','rpL16','rpL18','rpL22','rpL24','rpL2','rpL3','rpL4','rpL5','rpL6','rpS10','rpS17','rpS19','rpS3','rpS8']
archaea_list = ['rpL14','rpL15','rpL18','rpL22','rpL24','rpL2','rpL3','rpL4','rpL5','rpL6','rpS17','rpS19','rpS3','rpS8']
all_list = ['rpL14','rpL15','rpL18','rpL22','rpL24','rpL2', 'rpL3','rpL4','rpL5','rpL6','rpS17','rpS19','rpS3','rpS8']

if DOMAIN == 'archaea':    
    prot_list=archaea_list
elif DOMAIN == 'bacteria':
    prot_list=bacteria_list
elif DOMAIN == 'all': 
    prot_list=all_list

# Beginning message
print('')
print('#############################################')
print('metabolisHMM v' + VERSION)

# setup hmmsearch run depending on HMM list
print("Running ribosomal protein HMM searches...")
for genome in genomes:
    name=os.path.basename(genome).replace(".faa", "").strip().splitlines()[0]
    dir=name
    os.makedirs(OUTPUT+"/out/"+dir)
    for prot in prot_list:
        marker ="ribosomal_markers/"+prot+"_bact.hmm"
        outname= OUTPUT + "/out/"+dir+"/"+name + "-" + prot + ".out"
        cmd = ["hmmsearch", "--tblout="+outname, marker, genome]
        subprocess.call(cmd, stdout=FNULL)

# Parse HMM outputs
print("Parsing results...")
if DOMAIN == 'archaea':    
    prot_list=archaea_list
elif DOMAIN == 'bacteria':
    prot_list=bacteria_list
elif DOMAIN == 'all': 
    prot_list=all_list
result_dirs = os.walk(OUTPUT +"/out/")
for prot in prot_list:
    for path, dirs, files in result_dirs: 
        for file in files:
            genome=file.split("-")[0]
            marker=file.replace(".out", "").split("-")[1]
            result=OUTPUT + "/out/"+genome+"/"+file
            outfasta=OUTPUT + "/results/"+marker+".faa"
            genome_file = GENOMEDIR+genome+".faa"
            with open(outfasta, "a") as outf:
                with open(genome_file, "r") as input_fasta:
                    with open(result, "r") as input:
                        for qresult in SearchIO.parse(input, "hmmer3-tab"):
                            hits=qresult.hits
                            num_hits=len(hits)
                            if num_hits >0:
                                for i in range(0,1):
                                    hit_id=hits[i].id
                                for record in SeqIO.parse(input_fasta, "fasta"):
                                    if record.id in hit_id:
                                        outf.write(">"+genome+"\n"+str(record.seq)+"\n")
                       
# Make alignment file for each marker
print("Aligning ribosomal hits...")
prots = OUTPUT + "/results/*.faa"
fastas = glob.glob(prots)
for fasta in fastas:
    outname = os.path.basename(fasta).replace(".faa", "").strip().splitlines()[0]
    output= OUTPUT + "/results/"+outname+".aln"
    mafft_cmd = "mafft --quiet "
    mafft_cmd += fasta+" > "+output
    os.system(mafft_cmd)

# Concatenate alignments
print("Concatenating alignments...")
# referred to the biopython cookbook for concatenating multiple sequence alignments 
# adds question marks for missing loci for a given genome
prot_alignments = OUTPUT + "/results/*.aln"
infiles = glob.glob(prot_alignments)
concatout = OUTPUT + "/out/"+DOMAIN+"-concatenated-ribosomal-alignment.fasta"
alignments = [AlignIO.read(open(f, "r"), "fasta") for f in infiles]
all_labels = set(seq.id for aln in alignments for seq in aln)
tmp = defaultdict(list)
for aln in alignments:
    length = aln.get_alignment_length()
    these_labels = set(rec.id for rec in aln)
    missing = all_labels - these_labels
    for label in missing:
        new_seq = UnknownSeq(length) # prints ? marks for missing
        tmp[label].append(str(new_seq))
    for rec in aln:
        tmp[rec.id].append(str(rec.seq))
msa = MultipleSeqAlignment(SeqRecord(Seq(''.join(v)), id=k)
            for (k,v) in tmp.items())
AlignIO.write(msa,concatout,"fasta")

# check number of loci for a genome and print if less than a certain number to alert the user
counts = dict()
for file in infiles:
    with open(file) as f:
        for line in f:
            if line.startswith(">"):
                id = line.strip('\n').strip('>')
                counts[id] = counts.get(id, 0) + 1

# add argument for this so the user can know if they want a different value, but set a default value
x = int(args.loci)
for (k,v) in counts.items():
    if v < x:
        print("\t" + 'Genome '+ k + ' has fewer than ' + str(x) + ' hits!')

# get rid of "unknown description" descriptor in alignment file that biopython adds and I haven't figured out how to remove
reformatout=OUTPUT + "/results/"+DOMAIN+"-concatenated-ribosomal-alignment-reformatted.fasta"
filter_cmd = "sed 's/<unknown description>//g' "+concatout+" > "+reformatout
os.system(filter_cmd)

# Create tree
if PHYTOOL == 'fastree':
    print("Calculating tree using FastTree...")
    fileIn=OUTPUT + "/results/"+DOMAIN+"-concatenated-ribosomal-alignment-reformatted.fasta"
    outname = OUTPUT + "/results/"+DOMAIN+"-fastTree-ribosomal-tree.tre"
    fastCmd = "FastTree -quiet -nopr "+fileIn+" > "+outname
    os.system(fastCmd)
elif PHYTOOL == "raxml":
    print("Calculating tree with RaxML... be patient...")
    outname= DOMAIN+"-raxml-ribo"
    fileIn=OUTPUT + "/results/"+DOMAIN+"concatenated-ribosomal-alignment-reformatted.fasta"
    raxCmd = "raxmlHPC-PTHREADS -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s "+fileIn+" -T "+THREADS+" -n "+outname
    os.system(raxCmd)

# end message
print("Done! Find your results in "+ OUTPUT + "/results/")
print('#############################################')