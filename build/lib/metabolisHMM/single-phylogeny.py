#! /usr/bin/env python3 

###############################################################
# metabolisHMM - A tool for exploring and visualizing the distribution and evolutionary histories of metabolic markers
# single-marker-phylogeny : to create a phylogeny of a single HMM marker
# Written by Elizabeth McDaniel emcdaniel@wisc.edu
# November 2018
# This program is free software under the GNU General Public License version 3.0
###############################################################

import os, sys, glob, subprocess, argparse
import pandas as pd
from distutils.spawn import find_executable
from subprocess import Popen, DEVNULL
from Bio import BiopythonExperimentalWarning
from collections import defaultdict, Counter
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SearchIO, SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import UnknownSeq, Seq
from Bio.SeqRecord import SeqRecord

# Arguments and Directory setup
parser = argparse.ArgumentParser(description = "Create phylogeny of single marker")
parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
optional = parser.add_argument_group("optional arguments")
ribosomal = parser.add_argument_group("corresponding ribosomal tree arguments")
metadata = parser.add_argument_group("metadata output files for ITOL")
required.add_argument('--input', metavar='INPUT', help='Directory where genomes to be screened are held')
required.add_argument('--output', metavar='OUTPUT', help='Directory to store results and intermediate files')
required.add_argument('--marker', metavar='MARKER', help="Location of single marker to run analysis on")
optional.add_argument('--list', default='hit_list.txt', help="Output list of hits locus tags")
required.add_argument('--phylogeny', metavar='PHY', help="fastree or raxml, choose one")
optional.add_argument("--threads",metavar='THREADS',help="number of threads for tree making")
optional.add_argument("--kofam", metavar='KOFAM', help="Use KEGG HMMs from the KofamKOALA set. Options = ON or OFF")
optional.add_argument("--ko_list", metavar='KOLIST', help="Point to location of the KofamKoala ko_list file if using the KofamKOALA KEGG HMMs")
ribosomal.add_argument("--ribo_tree", metavar='RIBO', default='OFF', help="Make corresponding ribosomal phylogeny of genomes containing hits of the provided single marker. " )
ribosomal.add_argument("--domain", metavar='DOMAIN', help="If constructing corresponding ribosomal tree, select the domain for which your hits belong to. Options: bacteria, archaea")
ribosomal.add_argument('--loci', metavar='LOCI', default='12', help='Output genomes with less than x number of loci. By default prints genomes that have less than 12 ribosomal loci markers.')
metadata.add_argument('--metadata', metavar='METADATA',help='Option for outputting ITOL formatted metadata files. ON or OFF')
metadata.add_argument('--names', metavar='NAMES', help="Provided .csv formatted metadata file of filenames and corresponding taxonomical or group names")
metadata.add_argument('--itol_file', metavar='ITOL', default="itol_metadata.txt", help="Output iTOL formatted metadata file for changing leaf labels to taxonomical or group names")

# check for required installed dependencies
# prodigal
if find_executable('prodigal') is not None:
    pass
else:
    print('You do not have prodigal installed in your path. Please fix this.')
    sys.exit()
# hmmsearch
if find_executable('hmmsearch') is not None:
    pass
else:
    print('You do not have the hmmsearch executable from HMMER in your path. Please fix this.')
    sys.exit()
# mafft
if find_executable('mafft') is not None:
    pass
else: 
    print('You do not have MAFFT installed in your path. Please fix this.')
    sys.exit()
# fasttree, later check if phytool is raxml and check for correct raxml argument
if find_executable('FastTree') is not None:
    pass
else: 
    print('You do not have FastTree installed in your path. Please fix this.')
    sys.exit()


# if no arguments given, print help message
if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

# version to print
VERSION = '2.21'

# Beginning message
print('')
print('#############################################')
print('metabolisHMM v' + VERSION)

args = parser.parse_args()
GENOMEDIR = args.input
GENOMEFILES = GENOMEDIR + "/**"
OUTPUT = args.output
out_intm = OUTPUT + "/out"
out_results = OUTPUT + "/results"
out_genomes = OUTPUT + "/genomes"
MARKER = args.marker
PHYTOOL = args.phylogeny
THREADS = args.threads
RIBO = args.ribo_tree
DOMAIN = args.domain
KOFAM = args.kofam
KO_LIST = args.ko_list
METADATA = args.metadata
NAMES = args.names
ITOL_FILE = args.itol_file

if PHYTOOL == 'raxml':
    if find_executable('raxmlHPC-PTHREADS') is not None:
        pass
    else:
        print('You have opted to use raxml for building trees, but do not have the correct executable, raxmlHPC-PTHREADS installed in your path. Please fix this.')
        sys.exit()

# check if directory exists
if os.path.isdir(OUTPUT) == True:
    print("Directory "+ OUTPUT +" already exists! Please create different directory or remove the existing one.")
    sys.exit()

# if ribosomal option on, check that ribosomal markers directory exists
if RIBO == 'ON':
    if os.path.isdir("ribosomal_markers/") == False:
        print("     The ribosomal markers directory could not be found."+"\n"+"     Please either download the markers from https://github.com/elizabethmcd/metabolisHMM/releases/download/v2.0/metabolisHMM_v2.0_markers.tgz and decompress the tarball, or move the directory to where you are running the workflow from.")
        sys.exit()

# if kofam argument turned on, check if the corresponding ko_list metadata is provided for getting TC scores from
if KOFAM == 'ON':
    if os.path.exists("ko_list") == False:
        print("     You have opted to use a marker from the KofamKOALA distribution, but have not supplied the path to the ko_list.")
        sys.exit()
elif KOFAM == 'OFF':
    if os.path.exists(MARKER) == False: 
        print("     You have not provided an HMM marker to search with.")
        sys.exit()
    if os.path.exists(MARKER) == True: 
        with open(MARKER) as hmmfile:
            if 'TC' not in hmmfile.read():
                print("     You have not provided a threshold cutoff score for your HMM profile. Please add the line TC score-range to your HMM profile.")
                sys.exit()       
# check contents of ko_list to make sure it's the correct file
if KOFAM == 'ON':
    with open(KO_LIST) as ko:
        first_line = ko.readlines()[0]
        elem = 'knum'
        if elem not in first_line:
            print("     This does not look like the ko_list file from the KofamKOALA distribution.")
            sys.exit()

# make directories
os.makedirs(out_intm)
os.makedirs(out_results)
os.makedirs(out_genomes)
genomes=glob.glob(GENOMEFILES)
marker=MARKER
# turns off printing to stdout
DEVNULL = open(os.devnull, 'wb')
FNULL = open(os.devnull, 'w')
prot=os.path.basename(marker).replace(".hmm", "").strip().splitlines()[0]
dir=prot
os.makedirs(OUTPUT + "/out/"+dir)

# optional write out hit list
OUT_LIST_PATH = OUTPUT + "/results/" + args.list
OUT_LIST = open(OUT_LIST_PATH, "w")
OUT_LIST.write ("genome\tlocus_tag\n")

# if metadata option on, check the names file and itol path provided
if METADATA == 'ON':
    if os.path.exists(NAMES) == True:
    # header for itol output
        ITOL_PATH = OUTPUT + "/results/" + ITOL_FILE
        OUT_ITOL = open(ITOL_PATH, "w")
        OUT_ITOL.write("LABELS\nSEPARATOR TAB\nDATA\n")
        with open(NAMES, "r") as f:
            for record in f:
                filename, labelname = record.rstrip().split(',')
                OUT_ITOL.write('%s\t%s\n' % (filename, labelname))
elif METADATA == 'OFF':
    pass


# if .fna predict CDS and reformat header names because prodigal makes them stupid
# if .faa reformat the headers just in case contains weirdness
# if the user didn't provide the right files tell them
n = 0
print("Reformatting fasta files...")
for genome in genomes:
    if genome.endswith('.fna'):
        name = os.path.basename(genome).replace(".fna", "").strip().splitlines()[0]
        out_prot = OUTPUT + "/genomes/" + name + ".faa"
        out_gbk = OUTPUT + "/genomes/" + name + ".gbk"
        out_reformatted = OUTPUT + "/genomes/" + name + ".reformatted.faa"
        prodigal_cmd = "prodigal -q -i "+genome+" -a "+out_prot +" -o "+out_gbk
        os.system(prodigal_cmd)
        for seq_record in SeqIO.parse(out_prot, "fasta"):
            n = n + 1
            a = str(n).zfill(5)
            with open(out_reformatted, "a") as outre:
                outre.write(">" + name + "_" + str(a) + "\n")
                outre.write(str(seq_record.seq).replace("*","") + "\n")
    elif genome.endswith('.faa'):
        name = os.path.basename(genome).replace(".faa", "").strip().splitlines()[0]
        out_reformatted = OUTPUT + "/genomes/" + name + ".reformatted.faa"
        for seq_record in SeqIO.parse(genome, "fasta"):
            n = n + 1
            a = str(n).zfill(5)
            with open(out_reformatted, "a") as outre:
                outre.write(">" + name + "_" + str(a) + "\n")
                outre.write(str(seq_record.seq).replace("*","") + "\n")
    else:
        print("These do not look like fasta files that end in .fna or .faa. Please check your genome files.")
        sys.exit()
reformatted_path = OUTPUT + "/genomes/" + "*.reformatted.faa"
reformatted_genomes = glob.glob(reformatted_path)


# Run HMM for a single marker
print("Searching for " + prot + " marker in genome set...")
for genome in reformatted_genomes: 
    name=os.path.basename(genome).replace(".reformatted.faa", "").strip().splitlines()[0]
    outname= OUTPUT + "/out/"+dir+"/"+name + ".out"
    if KOFAM == 'ON':
        kofam_name = os.path.basename(marker).replace(".hmm", "").strip().splitlines()[0]
        scores = pd.read_csv(KO_LIST, delimiter="\t")
        scoredic = scores.set_index('knum')['threshold'].to_dict()
        tc = scoredic.get(kofam_name)
        cmd = ["hmmsearch","-T"+tc,"--tblout="+outname,marker,genome]
        print(cmd)
        subprocess.call(cmd, stdout=FNULL)
    else:
        cmd = ["hmmsearch","--cut_tc","--tblout="+outname, marker, genome]
        subprocess.call(cmd, stdout=FNULL)

# Parse HMM file 
print("Parsing all results...")
result_dir = os.walk(OUTPUT + "/out/"+dir)
for path, dirs, files in result_dir:
    for file in files:
        genome = file.replace(".out", "").strip().splitlines()[0]
        result = OUTPUT + "/out/"+dir+"/"+file
        output= OUTPUT + "/results/"+dir+".faa"
        genome_file=OUTPUT+"/genomes/"+genome+".reformatted.faa"
        with open(output, "a") as outf:
            with open(genome_file, "r") as input_fasta:
                with open(result, "r") as input:
                    for qresult in SearchIO.parse(input, "hmmer3-tab"):
                        hits = qresult.hits
                        num_hits = len(hits)
                        if num_hits>0:
                            for i in range(0,1):
                                hit_id=hits[i].id
                            for record in SeqIO.parse(input_fasta, "fasta"):
                                if record.id in hit_id:
                                    outf.write(">"+genome+"\n"+str(record.seq)+"\n")
                                    OUT_LIST.write('%s\t%s\n' % (genome, record.id))
OUT_LIST.close()
outf.close()

# Align hits 
print("Aligning hits...")
prots = OUTPUT + "/results/*.faa"
fastas = glob.glob(prots)
for fasta in fastas:
    outname = os.path.basename(fasta).replace(".faa", "").strip().splitlines()[0]
    output= OUTPUT + "/results/"+outname+".aln"
    mafft_cmd = "mafft --quiet "
    mafft_cmd += fasta+" > "+output
    os.system(mafft_cmd)

# Make tree 
if PHYTOOL == 'fastree':
    print("Calculating tree using FastTree...")
    marker_name = os.path.basename(marker).replace(".hmm", "").strip().splitlines()[0]
    logfile = OUTPUT + "/results/"+marker_name+".fastree.logfile"
    alignment_file = OUTPUT + "/results/"+marker_name+".aln"
    output_tree= OUTPUT + "/results/"+marker_name+".tre"
    tree_cmd = ["FastTree","-quiet","-log",logfile,"-out",output_tree,alignment_file]
    subprocess.Popen(tree_cmd, stdout=DEVNULL, stderr=DEVNULL)
elif PHYTOOL == "raxml":
    print("Calculating tree with RaxML... be patient...")
    marker_name = os.path.basename(marker).replace(".hmm", "").strip().splitlines()[0]
    outname= marker_name+"raxml"
    fileIn= OUTPUT + "/results/"+marker_name+".aln"
    raxCmd = "raxmlHPC-PTHREADS -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s "+fileIn+" -T "+THREADS+" -n "+outname
    os.system(raxCmd)


# Corresponding ribosomal tree of hits and options
# lists
bacteria_list = ['rpL14_bact','rpL15_bact','rpL16_bact','rpL18_bact','rpL22_bact','rpL24_bact','rpL2_bact','rpL3_bact','rpL4_bact','rpL5_bact','rpL6_bact','rpS10_bact','rpS17_bact','rpS19_bact','rpS3_bact','rpS8_bact']
archaea_list = ['rpL14_arch','rpL15_arch','rpL16_arch','rpL18_arch','rpL22_arch','rpL24_arch','rpL2_arch','rpL3_arch','rpL4_arch','rpL5_arch','rpL6_arch', 'rpS10_arch','rpS17_arch','rpS19_arch','rpS3_arch','rpS8_arch']
if DOMAIN == 'archaea':    
    prot_list=archaea_list
elif DOMAIN == 'bacteria':
    prot_list=bacteria_list
# runs phylogeny on specific hits of single marker
if RIBO == 'ON':
    print("Making corresponding ribosomal phylogeny of single marker hits...")
    fasta = OUTPUT + "/results/"+dir+".faa"
    hits = []
    for line in open(fasta):
        header = line.strip()
        if header.startswith(">"):
            hits.append(header[1:])
    for hit in hits: 
        genome = OUTPUT + "/genomes/"+hit+".reformatted.faa"
        os.makedirs(OUTPUT+"/ribo_out/"+hit)
        for prot in prot_list:
            marker ="curated_markers/ribosomal_markers/"+prot+".hmm"
            outname= OUTPUT + "/ribo_out/"+hit+"/"+hit + "-" + prot + ".out"
            cmd = ["hmmsearch", "--tblout="+outname, marker, genome]
            subprocess.call(cmd, stdout=FNULL)
    result_dirs = os.walk(OUTPUT + "/ribo_out/")
    for prot in prot_list:
        for path, dirs, files in result_dirs: 
            for file in files:
                genome=file.split("-")[0]
                marker=file.replace(".out", "").split("-")[1]
                result=OUTPUT + "/ribo_out/"+genome+"/"+file
                outfasta=OUTPUT + "/ribo_out/"+marker+".faa"
                genome_file = OUTPUT + "/genomes/"+genome+".reformatted.faa"
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
    ribos = OUTPUT + "/ribo_out/*.faa"
    fastas = glob.glob(ribos)
    for fasta in fastas:
        outname = os.path.basename(fasta).replace(".faa", "").strip().splitlines()[0]
        output= OUTPUT + "/ribo_out/"+outname+".aln"
        mafft_cmd = "mafft --quiet "
        mafft_cmd += fasta+" > "+output
        os.system(mafft_cmd)
    prot_alignments = OUTPUT + "/ribo_out/*.aln"
    infiles = glob.glob(prot_alignments)
    target = os.path.basename(args.marker).replace(".hmm", "").strip().splitlines()[0]
    concatout = OUTPUT + "/ribo_out/"+DOMAIN+"-"+target+"-concatenated-ribosomal-alignment.fasta"
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
    counts = dict()
    for file in infiles:
        with open(file) as f:
            for line in f:
                if line.startswith(">"):
                    id = line.strip('\n').strip('>')
                    counts[id] = counts.get(id, 0) + 1
    x = int(args.loci)
    under_cutoff = OUTPUT + "/ribo_out/genomes-under-cutoff.txt"
    for (k,v) in counts.items():
        if v < x:
            with open(under_cutoff, 'a') as uc:
                uc.write("Genome "+ k + " has fewer than " + str(x) + ' hits!' + "\n")
    print('Find list of genomes with less than ' + str(x) + ' ribosomal markers in ' + under_cutoff)
    reformatout=OUTPUT + "/ribo_out/"+DOMAIN+"-"+target+"-concatenated-ribosomal-alignment-reformatted.fasta"
    filter_cmd = "sed 's/<unknown description>//g' "+concatout+" > "+reformatout
    os.system(filter_cmd)
    if PHYTOOL == 'fastree':
        fileIn=reformatout
        outname = OUTPUT + "/results/"+DOMAIN+"-"+target+"-fastTree-corresponding-ribosomal-tree.tre"
        logfile = OUTPUT + "/ribo_out/"+"ribosomal-tree"+".fastree.logfile"
        fastCmd = ["FastTree","-quiet","-log",logfile,"-out",outname,fileIn]
        subprocess.Popen(fastCmd, stdout=DEVNULL, stderr=DEVNULL)
    elif PHYTOOL == "raxml":
        outname= DOMAIN+"-raxml-ribo"
        fileIn=reformatout
        raxCmd = "raxmlHPC-PTHREADS -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s "+fileIn+" -T "+THREADS+" -n "+outname
        os.system(raxCmd)
    print("Done! Find your results in "+ OUTPUT + "/results/")
    print('#############################################')
else: # end message
    print("Done! Find your results in "+ OUTPUT + "/results/")
    print('#############################################')

