#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  4 22:14:05 2020

@author: minsookim
"""

import os
import argparse
import pandas as pd
import numpy as np
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq


def variant_calling(datadir,libraryid,reffile,genome,minmapq,minbq,minstrand,workingdir,vepdir,vepcache,resultsdir):
    try:
        os.makedirs(f"{resultsdir}/MuTect2_results")
    except OSError:
        pass
    try:
        os.makedirs(f"{resultsdir}/MTvariant_results")
    except OSError:
        pass

    # Running MTvariantpipeline without matched normal
    print("Running MTvariantpipeline..")
    f1 = open(f"{resultsdir}/MTvp_out.txt", "w")
    subprocess.call("python3 " + workingdir + "/MTvariantpipeline.py -d " + datadir + "/ -v " + workingdir + "/TEMPMAFfiles/ -o " + 
        resultsdir + "/MTvariant_results/ -b " + libraryid + ".bam -g " + genome + " -q " + str(minmapq) + " -Q " + str(minbq) + 
        " -s " + str(minstrand) + " -w " + workingdir + "/ -vd " + vepdir + " -vc " + vepcache, shell=True, stdout=f1, stderr=subprocess.STDOUT)
    f1.close()

    # MuTect2 mitochondrial mode
    print("Running MuTect2..")
    f2 = open(f"{resultsdir}/MuTect2_out.txt", "w")
    subprocess.call("gatk --java-options -Xmx4g Mutect2 -R " + reffile + " --mitochondria-mode true -L MT -mbq " + str(minbq) + 
        " --minimum-mapping-quality " + str(minmapq) + " -I " + datadir + "/" + libraryid + ".bam -tumor " + libraryid.replace("-","_") + 
        " -O " + resultsdir + "/MuTect2_results/" + libraryid + ".bam.vcf.gz", shell=True, stdout=f2, stderr=subprocess.STDOUT)
    f2.close()

    # Left align MuTect2 results
    subprocess.call("bcftools norm -m - -f " + reffile + " " + resultsdir + "/MuTect2_results/" + libraryid + ".bam.vcf.gz" + 
        " -o " + resultsdir + "/MuTect2_results/" + libraryid + ".bam.vcf", shell=True)
    
    # Convert the MuTect2 result from vcf to maf file
    f6 = open(f"{resultsdir}/vcf2maf_out.txt", "w")
    subprocess.call("perl " + workingdir + "/vcf2maf/vcf2maf.pl --vep-data " + vepcache + " --vep-path " + vepdir + " --input-vcf " + 
        resultsdir + "/MuTect2_results/" + libraryid + ".bam.vcf" + " --output-maf " + resultsdir + "/MuTect2_results/" + libraryid + 
        ".bam.maf" + " --ncbi-build " + genome + ' --ref-fasta ' + reffile, shell=True, stdout=f6, stderr=subprocess.STDOUT)
    f6.close()


def variant_processing(libraryid, resultsdir):
    """
    Run MTvariantpipeline and MuTect2 on the filtered cells
    MTvariantpipeline: A simple variant calling and annotation pipeline for mitochondrial DNA variants.
    """
    print("Starting variant processing...")

    # Overlap between MuTect and MTvariantpipeline
    # Read in MTvariantpipeline result
    MTvarfile = pd.read_csv(resultsdir + "/MTvariant_results/" + libraryid + ".bam.maf", sep = "\t", comment='#', low_memory=False)

    # Read in MuTect result
    mutectfile = pd.read_csv(resultsdir + "/MuTect2_results/" + libraryid + ".bam.maf", sep = "\t", header=1, low_memory=False)
    saveasthis = resultsdir + "/" + libraryid + ".fillout"

    # Filter out variants falling in the repeat regions of 302-315, 513-525, and 3105-3109 (black listed regions)
    # Make sure End_Position is also not in the region
    rmregions = list(range(301,314)) + list(range(513,524)) + list(range(3105,3109))
    if len(mutectfile['Start_Position'][mutectfile['Start_Position'].isin(rmregions)]) > 0:
        mutectfile = mutectfile[~mutectfile['Start_Position'].isin(rmregions)]
    if len(MTvarfile['Start_Position'][MTvarfile['Start_Position'].isin(rmregions)]) > 0:
        rmthese = MTvarfile['Start_Position'].isin(rmregions)
        MTvarfile = MTvarfile[~rmthese]
    mutectfile.index = range(len(mutectfile.index))
    MTvarfile.index = range(len(MTvarfile.index))

    # Possible solution for filtering blacklited regions
    # rmregions = list(range(301,314)) + list(range(513,524)) + list(range(3105,3109))
    # mutectfile = mutectfile.loc[~mutectfile['Start_Position'].isin(rmregions)]
    # MTvarfile = MTvarfile.loc[~MTvarfile['Start_Position'].isin(rmregions)]
    # mutectfile.index = range(len(mutectfile.index))
    # MTvarfile.index = range(len(MTvarfile.index))
    
    # Output the overlap as final maf file
    combinedfile = pd.merge(mutectfile, MTvarfile, how='inner', on=['Chromosome','Start_Position','Reference_Allele',
        'Tumor_Seq_Allele2','Variant_Classification','Hugo_Symbol','EXON'])
    
    # Fix INDELs in the same position i.e. A:11866:AC and A:11866:ACC
    
    # Final annotation
    final_result = combinedfile.loc[:,['Tumor_Sample_Barcode_y','Matched_Norm_Sample_Barcode_y','Chromosome',
        'Start_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','Hugo_Symbol','EXON',
        'n_depth_y','t_depth_y','t_ref_count_y','t_alt_count_y']]
    final_result.columns = ['Sample','NormalUsed','Chrom','Start','Ref','Alt','VariantClass','Gene','Exon',
        'N_TotalDepth','T_TotalDepth','T_RefCount','T_AltCount']
    
    # output the fillout results
    final_result.to_csv(saveasthis,sep = '\t',na_rep='NA',index=False)


def runhaplogrep(datadir,libraryid,reffile, workingdir, resultsdir):
    """
    Run haplogrep to obtain the haplogroup information from the bam file
    """
    print("Running haplogrep..")
    
    # Filter the bam file for unmapped reads and mapping quality less than 1
    subprocess.call("samtools view -bF 4 -q 1 " + datadir + "/" + libraryid + ".bam > " + resultsdir + "/" + libraryid + "_filtered.bam", shell=True)
    
    # Index the filtered bam file
    subprocess.call("samtools index " + resultsdir + "/" + libraryid + "_filtered.bam", shell=True)
    
    # Edit the RG of the filtered bam file
    f5 = open(f"{resultsdir}/picard_out.txt", "w")
    subprocess.call("java -Xms8G -Xmx8G -jar " + workingdir + "/reference/picard.jar AddOrReplaceReadGroups I=" + 
        resultsdir + "/" + libraryid + "_filtered.bam O=" + datadir + "/" + libraryid + ".bam RGID=" + libraryid.replace("-", "_") + 
        " RGLB=" + libraryid + " RGPL=illumina RGPU=unit1 RGSM=" + libraryid, shell=True, stdout=f5, stderr=subprocess.STDOUT)
    f5.close()
    
    # Index the resulting bam file
    subprocess.call("samtools index " + datadir + "/" + libraryid + ".bam", shell=True)
    
    # Run MuTect2
    f3 = open(f"{resultsdir}/MuTect2_2_out.txt", "w")
    subprocess.call("gatk --java-options -Xmx4g Mutect2 -R " + reffile + " --mitochondria-mode true -L MT -mbq " + str(minbq) + 
        " --minimum-mapping-quality " + str(minmapq) + " -I " + datadir + "/" + libraryid + ".bam -tumor result" + libraryid.replace("-","_") + 
        " -O " + resultsdir + "/MuTect2_results/" + libraryid + ".bam.vcf.gz", shell=True, stdout=f3, stderr=subprocess.STDOUT)
    f3.close()

    # Run haplogrep2.1
    f4 = open(f"{resultsdir}/haplogrep_out.txt", "w")
    subprocess.call("java -jar " + workingdir + "/reference/haplogrep/haplogrep-2.1.20.jar --in " + resultsdir + 
        "/MuTect2_results/" + libraryid + ".bam.vcf.gz" + " --format vcf --extend-report --out " + resultsdir +
        "/" + libraryid + "_haplogroups.txt", shell=True, stdout=f4, stderr=subprocess.STDOUT)
    f4.close()


def processfillout(libraryid,resultsdir):
    """
    Run the combined mutation estimation on fillout
    Post-processing of the fillout files
    threshold: the critical threshold for calling a cell wild-type
    """
    print("Running the mutation estimation on the fillout..")
    
    # Import the final fillout file
    filloutfile = pd.read_csv(resultsdir + "/" + libraryid + '.fillout', sep='\t')
    
    # Set rownames
    filloutfile.index = [str(filloutfile['Ref'][i]) + ':' + str(int(filloutfile['Start'][i])) + ':' + 
        str(filloutfile['Alt'][i]) for i in range(len(filloutfile))]
    
    # Import haplogrep result
    haplogrepfile = pd.read_csv(os.path.join(resultsdir + "/" + libraryid + '_haplogroups.txt'), sep='\t')
    germlinepos = [x[:-1] for x in haplogrepfile['Found_Polys'][0].split(" ")]

    # Assign variants with >95% VAF as germline if they are used in haplogroup assignment and as homoplasmic otherwise
    filloutfile['somaticstatus'] = 'somatic'
    filloutfile['somaticstatus'].iloc[np.where(np.logical_and((filloutfile['T_AltCount']/filloutfile['T_TotalDepth'] >= 0.95), 
        (filloutfile['Start'].isin(germlinepos))))] = 'germline'
    filloutfile['somaticstatus'].iloc[np.where(np.logical_and((filloutfile['T_AltCount']/filloutfile['T_TotalDepth'] >= 0.95), 
        ~(filloutfile['Start'].isin(germlinepos))))] = 'homoplasmic'
    
    # Output filtered variant file
    filteredvar = filloutfile.loc[:,['Sample','NormalUsed','Chrom','Start','Ref','Alt','VariantClass','Gene','Exon','somaticstatus']]
    filteredvar.to_csv(resultsdir + "/" + libraryid + '_variants.tsv',sep = '\t')
    

def genmaster(libraryid,reffile,resultsdir):
    """
    Run the combined mutation estimation on fillout
    Post-processing of the fillout files
    threshold: the critical threshold for calling a cell wild-type
    """
    print('Generating a master file and a binary matrix of somatic variants for the sample..')
    
    # Import the relevant files
    variantsfile = pd.read_csv(os.path.join(resultsdir + "/" + libraryid + '_variants.tsv'), sep='\t', index_col=0)
    filloutfile = pd.read_csv(os.path.join(resultsdir + "/" + libraryid + '.fillout'), sep='\t')

    # Obtaining all the unique positions
    allpos = np.array([variants[1] for variants in pd.Series(variantsfile.index.values).str.split(':')])
    _, idx = np.unique(allpos, return_index=True)
    uniqpos = allpos[np.sort(idx)]
    
    # Flip the ref and the alt allele for the germline variants
    start = [variants[2] for variants in pd.Series(variantsfile.index.values[variantsfile['somaticstatus'] == 'germline']).str.split(':')]
    pos = [variants[1] for variants in pd.Series(variantsfile.index.values[variantsfile['somaticstatus'] == 'germline']).str.split(':')]
    end = [variants[0] for variants in pd.Series(variantsfile.index.values[variantsfile['somaticstatus'] == 'germline']).str.split(':')]
    variantsfile.rename(index=dict(zip(variantsfile.index.values[variantsfile['somaticstatus'] == 'germline'],[x + ':' + str(y) + ':' + 
        str(z) for x,y,z in zip(start,pos,end)])), inplace=True)
    
    # Flip the ref allele for somatic variants to follow germline variants
    fixthese = []
    for eachpos in uniqpos: # for each unique positions
        curridx = [s for s in variantsfile.index.values if s.split(':')[1] == eachpos] # all the variants at the unique positions
        if eachpos in pos: # if the position is germline position,
            newstart = [[x for x,y,z in zip(start,pos,end)][pos.index(eachpos)] + eachstart[1:] for eachstart in 
                [variants.split(':')[0] + ':' + variants.split(':')[1] + ':' + variants.split(':')[2] for variants in curridx]]
            fixthese.extend(newstart)
        else:
            fixthese.extend(curridx)
    variantsfile.rename(index=dict(zip(variantsfile.index.values,fixthese)), inplace=True)
    
    # Update the row names of other files to be consistent with the variants file
    filloutfile.index=list(variantsfile.index.values)
    
    # New ref and alt alleles
    newref = [variants[0] for variants in pd.Series(fixthese).str.split(':')]
    newalt = [variants[2] for variants in pd.Series(fixthese).str.split(':')]
    
    # Fix the depth matrix to filter variants that are uncertain and order them based on filteredvariants matrix
    masterfile = filloutfile
    
    # Fix the read counts for individual cells for each row accounting for the germline variants
    sampleid = filloutfile['Sample'][0].split('.bam')[0]
    masterfile = pd.DataFrame(index=filloutfile.index.values, columns=[sampleid])
    masterfile[sampleid] = filloutfile['T_AltCount'].astype(int).astype(str).str.cat(filloutfile['T_TotalDepth'].astype(int).astype(str),sep='/')
    
    # Create a variant annotations file based on the fillout file
    variantannot = pd.DataFrame(index=filloutfile.index.values, columns=['Start','Ref','Alt','VariantClass','Gene','T_AltCount','T_RefCount'])
    variantannot = variantannot.fillna(0)
    
    # Include columns for 'Start','Ref','Alt','VariantClass','Gene','T_AltCount','T_RefCount'
    variantannot['Start'] = filloutfile['Start']
    variantannot['Ref'] = newref
    variantannot['Alt'] = newalt
    variantannot['oldRef'] = filloutfile['Ref']
    variantannot['oldAlt'] = filloutfile['Alt']
    variantannot['VariantClass'] = filloutfile['VariantClass']
    variantannot['Gene'] = filloutfile['Gene']
    variantannot['T_AltCount'] = filloutfile['T_AltCount']
    variantannot['T_RefCount'] = filloutfile['T_RefCount']

    # Obtain the mutation signature
    # Initialize the counts and mutation sigature matrix
    motifs_C = ["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT"]
    motifs_T = ["ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"]
    mutsigfile = pd.DataFrame(index=['counts_CA','counts_CG','counts_CT','counts_TA','counts_TC','counts_TG'], columns=range(16))
    mutsigfile = mutsigfile.fillna(0)
    counts_CA = np.zeros(16)
    counts_CG = np.zeros(16)
    counts_CT = np.zeros(16)
    counts_TA = np.zeros(16)
    counts_TC = np.zeros(16)
    counts_TG = np.zeros(16)

    # Import the reference fasta file
    fasta_sequences = SeqIO.parse(open(reffile),'fasta')
    for fasta in fasta_sequences:
        currheader, currsequence = fasta.id, fasta.seq
        if 'MT' in currheader:
            sequence = [base for base in currsequence]

    # Account for germline variants
    for eachone in range(len(pos)):
        sequence[int(pos[eachone])-1] = start[eachone]

    varref = [variants[0] for variants in pd.Series(variantsfile.index.values).str.split(':')]
    varpos = [variants[1] for variants in pd.Series(variantsfile.index.values).str.split(':')]
    varalt = [variants[2] for variants in pd.Series(variantsfile.index.values).str.split(':')]
    mutsigmotifs = []
    for eachone in range(len(varpos)):
        prevpos = int(varpos[eachone])-2
        currpos = int(varpos[eachone])-1
        nextpos = int(varpos[eachone])
        motif = ''.join([sequence[prevpos],sequence[currpos],sequence[nextpos]])
        mutsigmotifs.append(motif)
        if varref[eachone] == 'C':
            if varalt[eachone] == 'A':
                counts_CA[motifs_C.index(motif)] += 1
            elif varalt[eachone] == 'G':
                counts_CG[motifs_C.index(motif)] += 1
            elif varalt[eachone] == 'T':
                counts_CT[motifs_C.index(motif)] += 1
        elif varref[eachone] == 'T':
            if varalt[eachone] == 'A':
                counts_TA[motifs_T.index(motif)] += 1
            elif varalt[eachone] == 'C':
                counts_TC[motifs_T.index(motif)] += 1
            elif varalt[eachone] == 'G':
                counts_TG[motifs_T.index(motif)] += 1
        elif varref[eachone] == 'G':
            motif = str(Seq(motif).complement())
            if varalt[eachone] == 'A':
                counts_CT[motifs_C.index(motif)] += 1
            elif varalt[eachone] == 'C':
                counts_CG[motifs_C.index(motif)] += 1
            elif varalt[eachone] == 'T':
                counts_CA[motifs_C.index(motif)] += 1
        elif varref[eachone] == 'A':
            motif = str(Seq(motif).complement())
            if varalt[eachone] == 'C':
                counts_TG[motifs_T.index(motif)] += 1
            elif varalt[eachone] == 'G':
                counts_TC[motifs_T.index(motif)] += 1
            elif varalt[eachone] == 'T':
                counts_TA[motifs_T.index(motif)] += 1
    mutsigfile.loc['counts_CA'] = counts_CA
    mutsigfile.loc['counts_CG'] = counts_CG
    mutsigfile.loc['counts_CT'] = counts_CT
    mutsigfile.loc['counts_TA'] = counts_TA
    mutsigfile.loc['counts_TC'] = counts_TC
    mutsigfile.loc['counts_TG'] = counts_TG
    
    # store the mutation signature info in variants file
    variantsfile['mutsig'] = mutsigmotifs
    
    # Saving the mutation signature
    mutsigfile.to_csv(resultsdir + "/" + libraryid + '_mutsig.tsv',sep = '\t')
    
    # combine the matrix with the resulting matrix
    resultMT = pd.concat([variantannot,masterfile],axis=1,sort=False) # concatenate everything together
    
    # Saving the final masterfile
    resultMT.to_csv(resultsdir + "/" + libraryid + '_master.tsv',sep = '\t')


if __name__ == "__main__":
    # Parse necessary arguments
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-d", "--datadir",type=str, help="Directory for BAM files")
    parser.add_argument("-r", "--reffile",type=str, help="Reference fasta file")
    parser.add_argument("-g", "--genome",type=str, help="Genome version",default = "GRCh37")    
    parser.add_argument("-q","--mapq",type=int,help="Minimum mapping quality, default = 20",default = 20)
    parser.add_argument("-Q","--baseq",type=int,help="Minimum base quality, default = 20",default = 20)
    parser.add_argument("-s","--strand",type=int,help="Minimum number of reads mapping to forward and reverse strand to call mutation, default=2",default = 2)
    parser.add_argument("-t","--threshold",type=int,help="The critical threshold for calling a cell wild-type, default=0.1",default = 0.1)
    parser.add_argument("-w", "--workingdir", type=str, help="Working directory")
    parser.add_argument("-v", "--vepdir", type=str, help="Directory for vep")
    parser.add_argument("-vc", "--vepcache", type=str, help="Directory for vep cache")
    parser.add_argument("-re", "--resultsdir", type=str, help="Directory for results")
    
    # read in arguments
    args = parser.parse_args()
    datadir = args.datadir
    reffile = args.reffile
    genome = args.genome
    minmapq = args.mapq
    minbq = args.baseq
    minstrand = args.strand
    threshold = args.threshold
    workingdir = args.workingdir
    vepdir = args.vepdir
    vepcache = args.vepcache
    resultsdir = args.resultsdir
    
    # Noting all the parameters
    print("Miminum mapping quality of " + str(minmapq))
    print("Miminum base quality of " + str(minbq))
    print("Miminum number of reads mapping to forward and reverse strand to call mutation of " + str(minstrand))

    try:
        os.makedirs(f"{workingdir}/multibulk_results")
    except OSError as error:
        print(error)
    bamdata = [f for f in os.listdir(datadir) if f.endswith(".bam")]
    print("Number of bam files: ", len(bamdata))
    for bam in bamdata:
        print(f"STARTING PIPELINE FOR {bam}...")
        libraryid = bam[:-4]
        try:
            os.makedirs(f"{resultsdir}")
        except OSError as error:
            print(error)
        variant_calling(datadir,libraryid,reffile,genome,minmapq,minbq,minstrand,workingdir,vepdir,vepcache,resultsdir)
        variant_processing(libraryid, resultsdir)
        runhaplogrep(datadir,libraryid,reffile, workingdir, resultsdir)
        processfillout(libraryid, resultsdir)
        genmaster(libraryid,reffile,resultsdir)
        print(f"DONE WITH PIPELINE FOR {bam}\n")
    print("DONE WITH MULTIBULK")

