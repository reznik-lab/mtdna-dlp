#!/usr/bin/env python3

import os
import argparse
import re
import pandas as pd
import numpy as np
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq


def reference_detect(reffile):
    print("Determining the mtDNA chromosome name...")
    for sequence in SeqIO.parse(open(reffile), "fasta"):
        if re.search('MT', sequence.description.split(" ")[0]):
            return("MT")
        elif re.search('chrM', sequence.description.split(" ")[0]):
            return("chrM")
    raise Exception("Chromosome is neither MT nor chrM")


def mappingquality(reffile, datadir, libraryid):
    print("Converting mapping qualities...")
    subprocess.run(f"java -Xmx5G -Xms5G -jar {workingdir}/reference/GenomeAnalysisTK.jar " +
        f"-T SplitNCigarReads -R {reffile} -I {datadir}/{libraryid}.bam -o {datadir}/{libraryid}.bam " +
        "-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS", shell=True, check=True)
    subprocess.run(f"samtools index {datadir}/{libraryid}.bam", shell=True, check=True)


def variant_calling_normal(resultsdir,datadir,libraryid,reffile,genome,minmapq,minbq,minstrand,workingdir,vepcache,mtchrom,ncbibuild,species,normal,normaldir,mincounts):
    try:
        os.makedirs(f"{resultsdir}/temp_MuTect2_results")
    except OSError:
        pass
    try:
        os.makedirs(f"{resultsdir}/MuTect2_results")
    except OSError:
        pass
    try:
        os.makedirs(f"{resultsdir}/MTvariant_results")
    except OSError:
        pass

    # Running MTvariantpipeline with matched normal
    print("Running MTvariantpipeline with matched normal..")
    subprocess.run(f"python3 {workingdir}/MTvariantpipeline2.py -d {datadir}/ -v {resultsdir}/TEMPMAFfiles/ " +
        f"-o {resultsdir}/MTvariant_results/ -b {libraryid}.bam -n {normal}.bam -nd {normaldir}/ -f {genome} -q {minmapq} " +
        f"-Q {minbq} -s {minstrand} -w {workingdir}/ -vc {vepcache} -f {reffile} -m {mtchrom} -c {mincounts}", shell=True, check=True)

    # MuTect2 mitochondrial mode on tumor
    print("Running MuTect2 on tumor..")
    subprocess.run(f"gatk --java-options -Xmx4g Mutect2 -R {reffile} --mitochondria-mode true -L {mtchrom} " +
        f"-mbq {minbq} --minimum-mapping-quality {minmapq} -I {datadir}/{libraryid}.bam " +
        f"-tumor {libraryid.replace('-','_')} -O {resultsdir}/temp_MuTect2_results/{libraryid}.bam.vcf.gz", shell=True, check=True)

    # MuTect2 mitochondrial mode on normal
    print("Running MuTect2 on normal..")
    subprocess.run(f"gatk --java-options -Xmx4g Mutect2 -R {reffile} --mitochondria-mode true -L {mtchrom} " +
        f"-mbq {minbq} --minimum-mapping-quality {minmapq} -I {normaldir}/{normal}.bam " +
        f"-tumor {normal.replace('-','_')} -O {resultsdir}/temp_MuTect2_results/{normal}.bam.vcf.gz", shell=True, check=True)

    # Left align MuTect2 results (-m - is there for a reason)
    subprocess.run(f"bcftools norm -m - -f {reffile} {resultsdir}/temp_MuTect2_results/{libraryid}.bam.vcf.gz " +
        f"-o {resultsdir}/temp_MuTect2_results/{libraryid}.bam.vcf", shell=True, check=True)
    subprocess.run(f"bcftools norm -m - -f {reffile} {resultsdir}/temp_MuTect2_results/{normal}.bam.vcf.gz " +
        f"-o {resultsdir}/temp_MuTect2_results/{normal}.bam.vcf", shell=True, check=True)
    
    # Convert the MuTect2 result from vcf to maf file
    subprocess.run(f"perl {workingdir}/vcf2maf/vcf2maf.pl --species {species} --vep-data {vepcache} " +
        f"--ncbi-build {ncbibuild} --input-vcf {resultsdir}/temp_MuTect2_results/{libraryid}.bam.vcf " +
        f"--output-maf {resultsdir}/temp_MuTect2_results/{libraryid}.bam.maf --ref-fasta {reffile}", shell=True, check=True)
    subprocess.run(f"perl {workingdir}/vcf2maf/vcf2maf.pl --species {species} --vep-data {vepcache} " +
        f"--ncbi-build {ncbibuild} --input-vcf {resultsdir}/temp_MuTect2_results/{normal}.bam.vcf " + 
        f"--output-maf {resultsdir}/temp_MuTect2_results/{normal}.bam.maf --ref-fasta {reffile}", shell=True, check=True)

    # Run R script to merge tumor and normal mafs
    print("Merging tumor and normal mafs..")
    subprocess.run(f"Rscript {workingdir}/getMAFfromfile.R {resultsdir}/temp_MuTect2_results/{libraryid}.bam.maf " +
        f"{resultsdir}/temp_MuTect2_results/{normal}.bam.maf {resultsdir}/MuTect2_results/{libraryid}.bam.maf", shell=True, check=True)

    subprocess.run(f"rm {resultsdir}/TEMPMAFfiles/*.bam_temp2.maf", shell=True, check=True)


def variant_processing_normal(libraryid,resultsdir):
    """
    Run MTvariantpipeline and MuTect2 on the filtered cells
    MTvariantpipeline: A simple variant calling and annotation pipeline for mitochondrial DNA variants.
    """
    print("Starting variant processing...")

    # Overlap between MuTect and MTvariantpipeline
    # Read in MTvariantpipeline result
    MTvarfile = pd.read_csv(resultsdir + "/MTvariant_results/" + libraryid + ".bam.maf", sep = "\t", comment='#', low_memory=False)

    # Read in MuTect result
    mutectfile = pd.read_csv(resultsdir + "/MuTect2_results/" + libraryid + ".bam.maf", sep = "\t", header=0, low_memory=False)

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

    # convert Hugo_Symbol column to uppercase
    mutectfile["Hugo_Symbol"] = mutectfile["Hugo_Symbol"].str.upper()
    MTvarfile["Hugo_Symbol"] = MTvarfile["Hugo_Symbol"].str.upper()
    
    # Output the overlap as final maf file
    combinedfile = pd.merge(mutectfile, MTvarfile, how='inner', on=['Chromosome','Start_Position',
        'Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','Variant_Type'])
        # 'EXON'])
    
    # Fix INDELs in the same position i.e. A:11866:AC and A:11866:ACC
    aux = combinedfile.loc[combinedfile['Variant_Type'] == 'INS'].groupby('Start_Position').count()['Hugo_Symbol_y'].reset_index()
    positions = list(aux['Start_Position'].loc[aux['Hugo_Symbol_y'] > 1])
    variants = list(combinedfile['ShortVariantID_y'].loc[(combinedfile['Start_Position'].isin(positions)) & (combinedfile['Variant_Type'] == 'INS')])
    if len(positions) != 0:
        dff = combinedfile.loc[combinedfile['ShortVariantID_y'].isin(variants)]

        # Create an auxuliary file only with the last rows to keep: keep unique positions with the highest TumorVAF
        dffaux = dff.sort_values(by='TumorVAF', ascending = False)
        dffaux = dffaux.drop_duplicates('Start_Position', keep = 'first')
        for i in positions:
            vals = dff[['t_alt_count_y', 't_alt_count_x']].loc[dff['Start_Position'] == i].sum(axis = 0).reset_index()
            dvals = dict(zip(list(vals['index']),list(vals[0])))
            dffaux.loc[dffaux['Start_Position'] == i,'t_alt_count_y'] = dvals['t_alt_count_y']
            dffaux.loc[dffaux['Start_Position'] == i,'t_alt_count_x'] = dvals['t_alt_count_x']

        #Remove all variants with duplicated indels
        combinedfile = combinedfile.loc[(~combinedfile['ShortVariantID_y'].isin(variants))]

        # Add unique indel variants with new values
        combinedfile = pd.concat([combinedfile, dffaux])
        combinedfile = combinedfile.sort_values(by='Start_Position', ascending = True)

        # Recalculate TumorVAF
        # combinedfile['TumorVAF_y'] = combinedfile['t_alt_count_y'] / (combinedfile['t_ref_count_y'] + combinedfile['t_alt_count_y'])
        # combinedfile['TumorVAF_x'] = combinedfile['t_alt_count_x'] / (combinedfile['t_ref_count_x'] + combinedfile['t_alt_count_x'])
    
    # Final annotation
    final_result = combinedfile.loc[:,['Tumor_Sample_Barcode_y','Matched_Norm_Sample_Barcode_y','Chromosome',
        'Start_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','Hugo_Symbol_y','EXON_y',
        'n_depth_y',"n_ref_count_y","n_alt_count_y",'t_depth_y','t_ref_count_y','t_alt_count_y',"t_alt_fwd","t_alt_rev"]]
    final_result.columns = ['Sample','NormalUsed','Chrom','Start','Ref','Alt','VariantClass','Gene','Exon',
        'N_TotalDepth',"N_RefCount","N_AltCount",'T_TotalDepth','T_RefCount','T_AltCount',"T_AltFwd","T_AltRev"]
    
    # output the fillout results
    final_result.to_csv(f"{resultsdir}/{libraryid}_master.tsv",sep = '\t',na_rep='NA',index=False)


def variant_calling(resultsdir,datadir,libraryid,reffile,genome,minmapq,minbq,minstrand,workingdir,vepcache,mtchrom,ncbibuild,species,mincounts):
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
    subprocess.run(f"python3 {workingdir}/MTvariantpipeline.py -d {datadir}/ -v {resultsdir}/TEMPMAFfiles/ " +
        f"-o {resultsdir}/MTvariant_results/ -b {libraryid}.bam -g {genome} -q {minmapq} -Q {minbq} " +
        f"-s {minstrand} -w {workingdir}/ -vc {vepcache} -f {reffile} -m {mtchrom} -c {mincounts}", shell=True, check=True)

    # MuTect2 mitochondrial mode
    print("Running MuTect2..")
    subprocess.run(f"gatk --java-options -Xmx4g Mutect2 -R {reffile} --mitochondria-mode true -L {mtchrom} " +
        f"-mbq {minbq} --minimum-mapping-quality {minmapq} -I {datadir}/{libraryid}.bam " +
        f"-tumor {libraryid.replace('-','_')} -O {resultsdir}/MuTect2_results/{libraryid}.bam.vcf.gz", shell=True, check=True)

    # Left align MuTect2 results (-m - is there for a reason)
    subprocess.run(f"bcftools norm -m - -f {reffile} {resultsdir}/MuTect2_results/{libraryid}.bam.vcf.gz " +
        f"-o {resultsdir}/MuTect2_results/{libraryid}.bam.vcf", shell=True, check=True)
    
    # Convert the MuTect2 result from vcf to maf file
    subprocess.run(f"perl {workingdir}/vcf2maf/vcf2maf.pl --species {species} --vep-data {vepcache} " +
        f"--ncbi-build {ncbibuild} --input-vcf {resultsdir}/MuTect2_results/{libraryid}.bam.vcf " + 
        f"--output-maf {resultsdir}/MuTect2_results/{libraryid}.bam.maf --ref-fasta {reffile}", shell=True, check=True)

    subprocess.run(f"rm {resultsdir}/TEMPMAFfiles/*.bam_temp2.maf", shell=True, check=True)


def variant_processing(libraryid,resultsdir):
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

    # convert Hugo_Symbol column to uppercase
    mutectfile["Hugo_Symbol"] = mutectfile["Hugo_Symbol"].str.upper()
    MTvarfile["Hugo_Symbol"] = MTvarfile["Hugo_Symbol"].str.upper()
    
    # Output the overlap as final maf file
    combinedfile = pd.merge(mutectfile, MTvarfile, how='inner', on=['Chromosome','Start_Position','Reference_Allele',
        'Tumor_Seq_Allele2','Variant_Classification',"Variant_Type",'EXON'])
    
    # Fix INDELs in the same position i.e. A:11866:AC and A:11866:ACC
    aux = combinedfile.loc[combinedfile['Variant_Type'] == 'INS'].groupby('Start_Position').count()['Hugo_Symbol_y'].reset_index()
    positions = list(aux['Start_Position'].loc[aux['Hugo_Symbol_y'] > 1])
    variants = list(combinedfile['ShortVariantID'].loc[(combinedfile['Start_Position'].isin(positions)) & (combinedfile['Variant_Type'] == 'INS')])
    if len(positions) != 0:
        dff = combinedfile.loc[combinedfile['ShortVariantID'].isin(variants)]

        # Create an auxuliary file only with the last rows to keep: keep unique positions with the highest TumorVAF
        dffaux = dff.sort_values(by='TumorVAF', ascending = False)
        dffaux = dffaux.drop_duplicates('Start_Position', keep = 'first')
        for i in positions:
            vals = dff[['t_alt_count_y', 't_alt_count_x']].loc[dff['Start_Position'] == i].sum(axis = 0).reset_index()
            dvals = dict(zip(list(vals['index']),list(vals[0])))
            dffaux.loc[dffaux['Start_Position'] == i,'t_alt_count_y'] = dvals['t_alt_count_y']
            dffaux.loc[dffaux['Start_Position'] == i,'t_alt_count_x'] = dvals['t_alt_count_x']

        #Remove all variants with duplicated indels
        combinedfile = combinedfile.loc[(~combinedfile['ShortVariantID'].isin(variants))]

        # Add unique indel variants with new values
        combinedfile = pd.concat([combinedfile, dffaux])
        combinedfile = combinedfile.sort_values(by='Start_Position', ascending = True)

        # Recalculate TumorVAF
        # combinedfile['TumorVAF_y'] = combinedfile['t_alt_count_y'] / (combinedfile['t_ref_count_y'] + combinedfile['t_alt_count_y'])
        # combinedfile['TumorVAF_x'] = combinedfile['t_alt_count_x'] / (combinedfile['t_ref_count_x'] + combinedfile['t_alt_count_x'])
    
    # Final annotation
    final_result = combinedfile.loc[:,['Tumor_Sample_Barcode_y','Matched_Norm_Sample_Barcode_y','Chromosome',
        'Start_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','Hugo_Symbol_y','EXON',
        'n_depth_y',"n_ref_count_y","n_alt_count_y",'t_depth_y','t_ref_count_y','t_alt_count_y',"t_alt_fwd","t_alt_rev"]]
    final_result.columns = ['Sample','NormalUsed','Chrom','Start','Ref','Alt','VariantClass','Gene','Exon',
        'N_TotalDepth',"N_RefCount","N_AltCount",'T_TotalDepth','T_RefCount','T_AltCount',"T_AltFwd","T_AltRev"]
    
    # output the fillout results
    final_result.to_csv(f"{resultsdir}/{libraryid}_master.tsv",sep = '\t',na_rep='NA',index=False)


def runhaplogrep(datadir,libraryid,reffile, workingdir, resultsdir):
    """
    Run haplogrep to obtain the haplogroup information from the bam file
    """
    print("Running haplogrep..")
    
    # Filter the bam file for unmapped reads and mapping quality less than 1
    subprocess.run(f"samtools view -bF 4 -q 1 {datadir}/{libraryid}.bam > {resultsdir}/{libraryid}_filtered.bam", shell=True, check=True)
    
    # Index the filtered bam file
    subprocess.run(f"samtools index {resultsdir}/{libraryid}_filtered.bam", shell=True, check=True)
    
    # Edit the RG of the filtered bam file
    subprocess.run(f"java -Xms8G -Xmx8G -jar {workingdir}/reference/picard.jar AddOrReplaceReadGroups " +
        f"I={resultsdir}/{libraryid}_filtered.bam O={resultsdir}/haplogroup_{libraryid}.bam " +
        f"RGID={libraryid.replace('-','_')} RGLB={libraryid} RGPL=illumina RGPU=unit1 RGSM={libraryid}", shell=True, check=True)
    
    # Index the resulting bam file
    subprocess.run(f"samtools index {resultsdir}/haplogroup_{libraryid}.bam", shell=True, check=True)
    
    # Run MuTect2
    subprocess.run(f"gatk --java-options -Xmx4g Mutect2 -R {reffile} --mitochondria-mode true -L {mtchrom} " +
        f"-mbq {minbq} --minimum-mapping-quality {minmapq} -I haplogroup_{resultsdir}/{libraryid}.bam " +
        f"-tumor result{libraryid.replace('-','_')} -O {resultsdir}/MuTect2_results/haplogroup_{libraryid}.bam.vcf.gz", shell=True, check=True)

    # Run haplogrep2.1
    subprocess.run(f"java -jar {workingdir}/reference/haplogrep/haplogrep-2.1.20.jar " +
        f"--in {resultsdir}/MuTect2_results/haplogroup_{libraryid}.bam.vcf.gz --format vcf --extend-report " +
        f"--out {resultsdir}/{libraryid}_haplogroups.txt", shell=True, check=True)


def processfillout(libraryid, resultsdir, genome, molecule):
    """
    Run the combined mutation estimation on fillout
    Post-processing of the fillout files
    threshold: the critical threshold for calling a cell wild-type
    """
    print("Running the mutation estimation on the fillout..")
    
    # Import the final fillout file
    filloutfile = pd.read_csv(f"{resultsdir}/{libraryid}_master.tsv",sep = '\t')
    
    # Set rownames
    filloutfile.index = [str(filloutfile['Ref'][i]) + ':' + str(int(filloutfile['Start'][i])) + ':' + 
        str(filloutfile['Alt'][i]) for i in range(len(filloutfile))]
    
    # # Assign variants with >95% VAF as germline if they are used in haplogroup assignment and as homoplasmic otherwise
    # filloutfile['ancestral'] = False
    # if (genome == "GRCh38" or genome == "GRCh37") and molecule == "dna":
    #     haplogrepfile = pd.read_csv(os.path.join(resultsdir + "/" + libraryid + '_haplogroups.txt'), sep='\t')
    #     germlinepos = [x[:-1] for x in haplogrepfile['Found_Polys'][0].split(" ")]
    #     filloutfile['ancestral'].iloc[np.where(np.logical_and((filloutfile['T_AltCount']/filloutfile['T_TotalDepth'] >= 0.95), 
    #         (filloutfile['Start'].isin(germlinepos))))] = True

    # Output filtered variant file
    # filteredvar = filloutfile.loc[:,['Sample','NormalUsed','Chrom','Start','Ref','Alt','VariantClass','Gene','Exon']]
    # filteredvar = filloutfile.loc[:,['Sample','NormalUsed','Chrom','Start','Ref','Alt','VariantClass','Gene','Exon','Ancestral']]
    # filteredvar.to_csv(resultsdir + "/" + libraryid + '_variants.tsv',sep = '\t')
    filloutfile.to_csv(f"{resultsdir}/{libraryid}_master.tsv",sep = '\t')    

def genmaster(libraryid,reffile,resultsdir,genome):
    """
    Run the combined mutation estimation on fillout
    Post-processing of the fillout files
    threshold: the critical threshold for calling a cell wild-type
    """
    print('Generating a master file and a binary matrix of somatic variants for the sample..')
    
    # Import the relevant files
    # variantsfile = pd.read_csv(os.path.join(resultsdir + "/" + libraryid + '_variants.tsv'), sep='\t', index_col=0)
    filloutfile = pd.read_csv(f"{resultsdir}/{libraryid}_master.tsv",sep = '\t',index_col=0)
    
    # # Fix the depth matrix to filter variants that are uncertain and order them based on filteredvariants matrix
    # masterfile = filloutfile
    
    # # Fix the read counts for individual cells for each row accounting for the germline variants
    # sampleid = filloutfile['Sample'][0].split('.bam')[0]
    # masterfile = pd.DataFrame(index=filloutfile.index.values, columns=[sampleid])
    # masterfile[sampleid] = filloutfile['T_AltCount'].astype(int).astype(str).str.cat(filloutfile['T_TotalDepth'].astype(int).astype(str),sep='/')

    # # Create a variant annotations file based on the fillout file
    # variantannot = pd.DataFrame(index=filloutfile.index.values, columns=['Start','Ref','Alt','VariantClass','Gene','T_AltCount','T_RefCount'])
    # variantannot = variantannot.fillna(0)
    
    # # Include columns for 'Start','Ref','Alt','VariantClass','Gene','T_AltCount','T_RefCount'
    # variantannot['Ref'] = filloutfile['Ref']
    # variantannot['Alt'] = filloutfile['Alt']
    # variantannot['Gene'] = filloutfile['Gene']
    # variantannot['VariantClass'] = filloutfile['VariantClass']
    # variantannot['T_AltCount'] = filloutfile['T_AltCount']
    # variantannot['T_RefCount'] = filloutfile['T_RefCount']
    # variantannot['Start'] = filloutfile['Start']

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

    varref = [variants[0] for variants in pd.Series(filloutfile.index.values).str.split(':')]
    varpos = [variants[1] for variants in pd.Series(filloutfile.index.values).str.split(':')]
    varalt = [variants[2] for variants in pd.Series(filloutfile.index.values).str.split(':')]
    # varref = [variants[0] for variants in pd.Series(variantsfile.index.values).str.split(':')]
    # varpos = [variants[1] for variants in pd.Series(variantsfile.index.values).str.split(':')]
    # varalt = [variants[2] for variants in pd.Series(variantsfile.index.values).str.split(':')]
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
    filloutfile['mutsig'] = mutsigmotifs
    
    # Saving the mutation signature
    mutsigfile.to_csv(resultsdir + "/" + libraryid + '_mutsig.tsv',sep = '\t')
    
    # # combine the matrix with the resulting matrix
    # resultMT = pd.concat([variantannot,masterfile],axis=1,sort=False) # concatenate everything together
    
    # # Saving the final masterfile
    # resultMT.to_csv(resultsdir + "/" + libraryid + '_master.tsv',sep = '\t')

    # Calculate heteroplasmy
    filloutfile["Heteroplasmy"] = filloutfile['T_AltCount'].astype(int) / filloutfile['T_TotalDepth'].astype(int)

    # Combined matrices together
    filloutfile.to_csv(f"{resultsdir}/{libraryid}_master.tsv",sep = '\t',na_rep='NA',index=False)

if __name__ == "__main__":
    # Parse necessary arguments
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-d", "--datadir",type=str, help="(REQUIRED) Directory for BAM files", required=True)
    parser.add_argument("-l", "--libraryid",type=str, help="(REQUIRED) Library ID", required=True)
    parser.add_argument("-w", "--workingdir", type=str, help="(REQUIRED) Working directory", required=True)
    parser.add_argument("-re", "--resultsdir", type=str, help="(REQUIRED) Directory for results", required=True)
    parser.add_argument("-r", "--reffile",type=str, help="Reference fasta file",default="")
    parser.add_argument("-g", "--genome",type=str, help="Genome version",default = "GRCh37")    
    parser.add_argument("-q","--mapq",type=int,help="Minimum mapping quality, default = 20",default = 20)
    parser.add_argument("-Q","--baseq",type=int,help="Minimum base quality, default = 20",default = 20)
    parser.add_argument("-s","--strand",type=int,help="Minimum number of reads mapping to forward and reverse strand to call mutation, default=2",default = 2)
    parser.add_argument("-t","--threshold",type=int,help="The critical threshold for calling a cell wild-type, default=0.1",default = 0.1)
    parser.add_argument("-vc", "--vepcache", type=str, help="Directory for vep cache", default="$HOME/.vep")
    parser.add_argument("-n", "--normal", type=str, help="matched normal file",default="")
    parser.add_argument("-nd", "--normaldir", type=str, help="directory that contains matched normal file",default="")
    parser.add_argument("-m", "--molecule",type=str, help="Type of molecule (dna or rna), default=dna", default="dna")
    parser.add_argument("-c","--mincounts",type=int,help="Minimum number of read counts for MTvariantpipeline, default = 100", default = 100)

    # read in arguments
    args = parser.parse_args()
    datadir = args.datadir
    reffile = args.reffile
    genome = args.genome
    minmapq = args.mapq
    minbq = args.baseq
    minstrand = args.strand
    threshold = args.threshold
    libraryid = args.libraryid
    workingdir = args.workingdir
    vepcache = args.vepcache
    resultsdir = args.resultsdir
    normal = args.normal
    normaldir = args.normaldir
    molecule = args.molecule
    mincounts = args.mincounts

    # Set the parameters for the genome build
    if genome == 'GRCh37':
        if reffile == "":
            reffile = workingdir + '/reference/b37/b37_MT.fa'
        ncbibuild = 'GRCh37'
        species = "homo_sapiens"
    elif genome == "GRCm38" or genome == "mm10":
        if reffile == "":
            reffile = workingdir + "/reference/mm10/mm10_MT.fa"
        ncbibuild = 'GRCm38'
        species = "mus_musculus"
    elif genome == 'GRCh38':
        if reffile == "":
            reffile = workingdir + '/reference/GRCh38/genome_MT.fa'
        ncbibuild = 'GRCh38'
        species = "homo_sapiens"
    else:
        raise Exception("The genome build you entered is not supported. Supported genomes are GRCh37, GRCh38, GRCm38, and mm10.")
    
    # Run reference_detect to determine mtchrom
    mtchrom = reference_detect(reffile)
    
    # Noting all the parameters
    print("Miminum mapping quality of " + str(minmapq))
    print("Miminum base quality of " + str(minbq))
    print("Miminum number of reads mapping to forward and reverse strand to call mutation of " + str(minstrand))
    print("Reference file: " + reffile)

    # Filtering of cells
    # if molecule == "rna":
    #     mappingquality(reffile,datadir,libraryid)
    if normal != "":
        variant_calling_normal(resultsdir,datadir,libraryid,reffile,genome,minmapq,minbq,minstrand,workingdir,vepcache,mtchrom,ncbibuild,species,normal,normaldir,mincounts)
        variant_processing_normal(libraryid,resultsdir)
    else:
        variant_calling(resultsdir,datadir,libraryid,reffile,genome,minmapq,minbq,minstrand,workingdir,vepcache,mtchrom,ncbibuild,species,mincounts)
        variant_processing(libraryid,resultsdir)
    # if (genome == "GRCh38" or genome == "GRCh37") and molecule == "dna":
    #     runhaplogrep(datadir,libraryid,reffile, workingdir, resultsdir)
    processfillout(libraryid, resultsdir,genome,molecule)
    genmaster(libraryid,reffile,resultsdir,genome)

    print("DONE WITH BULKPIPELINE")