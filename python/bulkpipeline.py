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


def variant_calling_normal(resultsdir,tumordir,tumor_id,reffile,genome,minmapq,minbq,minstrand,workingdir,vepcache,mtchrom,ncbibuild,species,normal_id,normaldir,molecule,mincounts):
    try:
        os.makedirs(f"{resultsdir}/TEMPMAFfiles/tempMuTect2")
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
    subprocess.run(f"python3 {workingdir}/MTvariantpipeline.py -d {tumordir}/ -v {resultsdir}/TEMPMAFfiles/ -w {workingdir}/ " +
        f"-o {resultsdir}/MTvariant_results/ -b {tumor_id}.bam -n {normal_id}.bam -nd {normaldir}/ -g {genome} -q {minmapq} " +
        f"-Q {minbq} -s {minstrand} -vc {vepcache} -f {reffile} -m {mtchrom} -mo {molecule} -c {mincounts}", shell=True, check=True)

    # MuTect2 mitochondrial mode on tumor
    print("Running MuTect2 on tumor..")
    subprocess.run(f"gatk --java-options -Xmx4g Mutect2 -R {reffile} --mitochondria-mode true -L {mtchrom} " +
        f"-mbq {minbq} --minimum-mapping-quality {minmapq} -I {tumordir}/{tumor_id}.bam " +
        f"-tumor {tumor_id.replace('-','_')} -O {resultsdir}/TEMPMAFfiles/tempMuTect2/{tumor_id}.bam.vcf.gz", shell=True, check=True)

    # MuTect2 mitochondrial mode on normal
    print("Running MuTect2 on normal..")
    subprocess.run(f"gatk --java-options -Xmx4g Mutect2 -R {reffile} --mitochondria-mode true -L {mtchrom} " +
        f"-mbq {minbq} --minimum-mapping-quality {minmapq} -I {normaldir}/{normal_id}.bam " +
        f"-tumor {normal_id.replace('-','_')} -O {resultsdir}/TEMPMAFfiles/tempMuTect2/{normal_id}.bam.vcf.gz", shell=True, check=True)

    # Left align MuTect2 results (-m - is there for a reason)
    subprocess.run(f"bcftools norm -m - -f {reffile} {resultsdir}/TEMPMAFfiles/tempMuTect2/{tumor_id}.bam.vcf.gz " +
        f"-o {resultsdir}/TEMPMAFfiles/tempMuTect2/{tumor_id}.bam.vcf", shell=True, check=True)
    subprocess.run(f"bcftools norm -m - -f {reffile} {resultsdir}/TEMPMAFfiles/tempMuTect2/{normal_id}.bam.vcf.gz " +
        f"-o {resultsdir}/TEMPMAFfiles/tempMuTect2/{normal_id}.bam.vcf", shell=True, check=True)
    
    # Convert the MuTect2 result from vcf to maf file
    subprocess.run(f"perl {workingdir}/vcf2maf/vcf2maf.pl --species {species} --vep-data {vepcache} " +
        f"--ncbi-build {ncbibuild} --input-vcf {resultsdir}/TEMPMAFfiles/tempMuTect2/{tumor_id}.bam.vcf " +
        f"--output-maf {resultsdir}/TEMPMAFfiles/tempMuTect2/{tumor_id}.bam.maf --ref-fasta {reffile}", shell=True, check=True)
    subprocess.run(f"perl {workingdir}/vcf2maf/vcf2maf.pl --species {species} --vep-data {vepcache} " +
        f"--ncbi-build {ncbibuild} --input-vcf {resultsdir}/TEMPMAFfiles/tempMuTect2/{normal_id}.bam.vcf " + 
        f"--output-maf {resultsdir}/TEMPMAFfiles/tempMuTect2/{normal_id}.bam.maf --ref-fasta {reffile}", shell=True, check=True)

    print("Merging tumor and normal mafs..")
    # Read in tumor result
    tumorfile = pd.read_csv(resultsdir + "/TEMPMAFfiles/tempMuTect2/" + tumor_id + ".bam.maf", sep = "\t", header=1, low_memory=False)
    # Read in normal result
    normalfile = pd.read_csv(resultsdir + "/TEMPMAFfiles/tempMuTect2/" + normal_id + ".bam.maf", sep = "\t", header=1, low_memory=False)
    
    # Output the overlap as final maf file
    combinedfile = pd.merge(tumorfile.loc[:,['Chromosome','Start_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','Variant_Type']], 
                            normalfile.loc[:,['Chromosome','Start_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','Variant_Type']], 
                            how='inner', on=['Chromosome','Start_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','Variant_Type'])
    # Combined matrices together
    combinedfile.to_csv(f"{resultsdir}/MuTect2_results/{tumor_id}.bam.maf",sep = '\t',na_rep='NA',index=False)
    # Remove temporary files
    subprocess.run(f"rm {resultsdir}/TEMPMAFfiles/*.bam_temp2.maf", shell=True)


def variant_calling(resultsdir,tumordir,tumor_id,reffile,genome,minmapq,minbq,minstrand,workingdir,vepcache,mtchrom,ncbibuild,species,molecule,mincounts):
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
    subprocess.run(f"python3 {workingdir}/MTvariantpipeline.py -d {tumordir}/ -v {resultsdir}/TEMPMAFfiles/ " +
        f"-o {resultsdir}/MTvariant_results/ -b {tumor_id}.bam -g {genome} -q {minmapq} -Q {minbq} -s {minstrand} " +
        f"-w {workingdir}/ -vc {vepcache} -f {reffile} -m {mtchrom} -mo {molecule} -c {mincounts}", shell=True, check=True)

    # MuTect2 mitochondrial mode
    print("Running MuTect2..")
    subprocess.run(f"gatk --java-options -Xmx4g Mutect2 -R {reffile} --mitochondria-mode true -L {mtchrom} " +
        f"-mbq {minbq} --minimum-mapping-quality {minmapq} -I {tumordir}/{tumor_id}.bam " +
        f"-tumor {tumor_id.replace('-','_')} -O {resultsdir}/MuTect2_results/{tumor_id}.bam.vcf.gz", shell=True, check=True)

    # Left align MuTect2 results (-m - is there for a reason)
    subprocess.run(f"bcftools norm -m - -f {reffile} {resultsdir}/MuTect2_results/{tumor_id}.bam.vcf.gz " +
        f"-o {resultsdir}/MuTect2_results/{tumor_id}.bam.vcf", shell=True, check=True)
    
    # Convert the MuTect2 result from vcf to maf file
    subprocess.run(f"perl {workingdir}/vcf2maf/vcf2maf.pl --species {species} --vep-data {vepcache} " +
        f"--ncbi-build {ncbibuild} --input-vcf {resultsdir}/MuTect2_results/{tumor_id}.bam.vcf " + 
        f"--output-maf {resultsdir}/MuTect2_results/{tumor_id}.bam.maf --ref-fasta {reffile}", shell=True, check=True)

    subprocess.run(f"rm {resultsdir}/TEMPMAFfiles/*.bam_temp2.maf", shell=True)


def variant_processing(tumor_id,resultsdir):
    """
    Run MTvariantpipeline and MuTect2 on the filtered cells
    MTvariantpipeline: A simple variant calling and annotation pipeline for mitochondrial DNA variants.
    """
    print("Starting variant processing...")

    # Overlap between MuTect and MTvariantpipeline
    # Read in MTvariantpipeline result
    MTvarfile = pd.read_csv(resultsdir + "/MTvariant_results/" + tumor_id + ".bam.maf", sep = "\t", comment='#', low_memory=False)
    # Read in MuTect result
    mutectfile = pd.read_csv(resultsdir + "/MuTect2_results/" + tumor_id + ".bam.maf", sep = "\t", comment='#', header=0, low_memory=False)

    # Filter out variants falling in the repeat regions of 513-525, and 3105-3109 (black listed regions)
    # Make sure End_Position is also not in the region
    rmregions = list(range(513,524)) + list(range(3105,3109))
    if len(mutectfile['Start_Position'][mutectfile['Start_Position'].isin(rmregions)]) > 0:
        mutectfile = mutectfile[~mutectfile['Start_Position'].isin(rmregions)]
    if len(MTvarfile['Start_Position'][MTvarfile['Start_Position'].isin(rmregions)]) > 0:
        rmthese = MTvarfile['Start_Position'].isin(rmregions)
        MTvarfile = MTvarfile[~rmthese]
    mutectfile.index = range(len(mutectfile.index))
    MTvarfile.index = range(len(MTvarfile.index))
    
    # Output the overlap as final maf file
    combinedfile = pd.merge(mutectfile, MTvarfile, how='inner', on=['Chromosome','Start_Position',
        'Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','Variant_Type'])
    combinedfile = combinedfile.rename(columns={"Hugo_Symbol_y":"Hugo_Symbol",
        "Entrez_Gene_Id_y":"Entrez_Gene_Id","Center_y":"Center","NCBI_Build_y":"NCBI_Build",
        "End_Position_y":"End_Position","Strand_y":"Strand","Tumor_Seq_Allele1_y":"Tumor_Seq_Allele1",
        "dbSNP_RS_y":"dbSNP_RS","dbSNP_Val_Status_y":"dbSNP_Val_Status","Tumor_Sample_Barcode_y":"Tumor_Sample_Barcode",
        "Matched_Norm_Sample_Barcode_y":"Matched_Norm_Sample_Barcode","Match_Norm_Seq_Allele1_y":"Match_Norm_Seq_Allele1",
        "Match_Norm_Seq_Allele2_y":"Match_Norm_Seq_Allele2","Tumor_Validation_Allele1_y":"Tumor_Validation_Allele1",
        "Tumor_Validation_Allele2_y":"Tumor_Validation_Allele2","Match_Norm_Validation_Allele1_y":"Match_Norm_Validation_Allele1",
        "Match_Norm_Validation_Allele2_y":"Match_Norm_Validation_Allele2","Verification_Status_y":"Verification_Status",
        "Validation_Status_y":"Validation_Status","Mutation_Status_y":"Mutation_Status","Sequencing_Phase_y":"Sequencing_Phase",
        "Sequence_Source_y":"Sequence_Source","Validation_Method_y":"Validation_Method","Score_y":"Score","BAM_File_y":"BAM_File",
        "Sequencer_y":"Sequencer","Tumor_Sample_UUID_y":"Tumor_Sample_UUID","Matched_Norm_Sample_UUID_y":"Matched_Norm_Sample_UUID",
        "HGVSc_y":"HGVSc","HGVSp_y":"HGVSp","HGVSp_Short_y":"HGVSp_Short","Exon_Number_y":"Exon_Number","t_depth_y":"t_depth",
        "t_ref_count_y":"t_ref_count","t_alt_count_y":"t_alt_count","n_depth_y":"n_depth","n_ref_count_y":"n_ref_count",
        "n_alt_count_y":"n_alt_count","all_effects_y":"all_effects","Gene_y":"Gene","Feature_y":"Feature",
        "Feature_type_y":"Feature_type","Consequence_y":"Consequence","cDNA_position_y":"cDNA_position",
        "CDS_position_y":"CDS_position","Protein_position_y":"Protein_position","Amino_acids_y":"Amino_acids",
        "Codons_y":"Codons","Existing_variation_y":"Existing_variation","ALLELE_NUM_y":"ALLELE_NUM","DISTANCE_y":"DISTANCE",
        "STRAND_VEP_y":"STRAND_VEP","SYMBOL_y":"SYMBOL","SYMBOL_SOURCE_y":"SYMBOL_SOURCE","HGNC_ID_y":"HGNC_ID",
        "BIOTYPE_y":"BIOTYPE","CANONICAL_y":"CANONICAL","CCDS_y":"CCDS","ENSP_y":"ENSP","SWISSPROT_y":"SWISSPROT",
        "TREMBL_y":"TREMBL","UNIPARC_y":"UNIPARC","RefSeq_y":"RefSeq","SIFT_y":"SIFT","PolyPhen_y":"PolyPhen","EXON_y":"EXON",
        "INTRON_y":"INTRON","DOMAINS_y":"DOMAINS","AF_y":"AF","AFR_AF_y":"AFR_AF","AMR_AF_y":"AMR_AF","ASN_AF_y":"ASN_AF",
        "EAS_AF_y":"EAS_AF","EUR_AF_y":"EUR_AF","SAS_AF_y":"SAS_AF","AA_AF_y":"AA_AF","EA_AF_y":"EA_AF","CLIN_SIG_y":"CLIN_SIG",
        "SOMATIC_y":"SOMATIC","PUBMED_y":"PUBMED","MOTIF_NAME_y":"MOTIF_NAME","MOTIF_POS_y":"MOTIF_POS","HIGH_INF_POS_y":"HIGH_INF_POS",
        "MOTIF_SCORE_CHANGE_y":"MOTIF_SCORE_CHANGE","IMPACT_y":"IMPACT","PICK_y":"PICK","TSL_y":"TSL","HGVS_OFFSET_y":"HGVS_OFFSET",
        "PHENO_y":"PHENO","MINIMISED_y":"MINIMISED","GENE_PHENO_y":"GENE_PHENO","FILTER_y":"FILTER","flanking_bps_y":"flanking_bps",
        "vcf_id_y":"vcf_id","vcf_qual_y":"vcf_qual","gnomAD_AF_y":"gnomAD_AF","gnomAD_AFR_AF_y":"gnomAD_AFR_AF","gnomAD_AMR_AF_y":"gnomAD_AMR_AF",
        "gnomAD_ASJ_AF_y":"gnomAD_ASJ_AF","gnomAD_EAS_AF_y":"gnomAD_EAS_AF","gnomAD_FIN_AF_y":"gnomAD_FIN_AF","gnomAD_NFE_AF_y":"gnomAD_NFE_AF",
        "gnomAD_OTH_AF_y":"gnomAD_OTH_AF","gnomAD_SAS_AF_y":"gnomAD_SAS_AF"})
    
    # Fix INDELs in the same position i.e. A:11866:AC and A:11866:ACC
    aux = combinedfile.loc[combinedfile['Variant_Type'] == 'INS'].groupby('Start_Position').count()['Hugo_Symbol'].reset_index()
    positions = list(aux['Start_Position'].loc[aux['Hugo_Symbol'] > 1])
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
        
    # Final annotation
    filloutfile = combinedfile.loc[:,['Hugo_Symbol','Entrez_Gene_Id','Center','NCBI_Build','Chromosome',
        'Start_Position','End_Position','Strand','Variant_Classification','Variant_Type','Reference_Allele',
        'Tumor_Seq_Allele1','Tumor_Seq_Allele2','dbSNP_RS','dbSNP_Val_Status','Tumor_Sample_Barcode',
        'Matched_Norm_Sample_Barcode','Match_Norm_Seq_Allele1','Match_Norm_Seq_Allele2','Tumor_Validation_Allele1',
        'Tumor_Validation_Allele2','Match_Norm_Validation_Allele1','Match_Norm_Validation_Allele2','Verification_Status',
        'Validation_Status','Mutation_Status','Sequencing_Phase','Sequence_Source','Validation_Method','Score','BAM_File',
        'Sequencer','Tumor_Sample_UUID','Matched_Norm_Sample_UUID','HGVSc','HGVSp','HGVSp_Short','Exon_Number','t_depth',
        't_ref_count','t_alt_count','n_depth','n_ref_count','n_alt_count','t_alt_fwd','t_alt_rev','all_effects','Gene',
        'Feature','Feature_type','Consequence','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons',
        'Existing_variation','ALLELE_NUM','DISTANCE','STRAND_VEP','SYMBOL','SYMBOL_SOURCE','HGNC_ID','BIOTYPE','CANONICAL',
        'CCDS','ENSP','SWISSPROT','TREMBL','UNIPARC','RefSeq','SIFT','PolyPhen','EXON','INTRON','DOMAINS','AF','AFR_AF',
        'AMR_AF','ASN_AF','EAS_AF','EUR_AF','SAS_AF','AA_AF','EA_AF','CLIN_SIG','SOMATIC','PUBMED','MOTIF_NAME',
        'MOTIF_POS','HIGH_INF_POS','MOTIF_SCORE_CHANGE','IMPACT','PICK','TSL','HGVS_OFFSET','PHENO','MINIMISED',
        'GENE_PHENO','FILTER','flanking_bps','vcf_id','vcf_qual','gnomAD_AF','gnomAD_AFR_AF','gnomAD_AMR_AF',
        'gnomAD_ASJ_AF','gnomAD_EAS_AF','gnomAD_FIN_AF','gnomAD_NFE_AF','gnomAD_OTH_AF','gnomAD_SAS_AF']]
    filloutfile.columns = ['Hugo_Symbol','Entrez_Gene_Id','Center','NCBI_Build','Chrom','Start','End_Position','Strand',
        'VariantClass','Variant_Type','Ref','Tumor_Seq_Allele1','Alt','dbSNP_RS','dbSNP_Val_Status','Sample','NormalUsed',
        'Match_Norm_Seq_Allele1','Match_Norm_Seq_Allele2','Tumor_Validation_Allele1','Tumor_Validation_Allele2',
        'Match_Norm_Validation_Allele1','Match_Norm_Validation_Allele2','Verification_Status','Validation_Status',
        'Mutation_Status','Sequencing_Phase','Sequence_Source','Validation_Method','Score','BAM_File','Sequencer',
        'Tumor_Sample_UUID','Matched_Norm_Sample_UUID','HGVSc','HGVSp','HGVSp_Short','Exon_Number','T_TotalDepth',
        'T_RefCount','T_AltCount','N_TotalDepth','N_RefCount','N_AltCount','T_AltFwd','T_AltRev','all_effects','Gene',
        'Feature','Feature_type','Consequence','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons',
        'Existing_variation','ALLELE_NUM','DISTANCE','STRAND_VEP','SYMBOL','SYMBOL_SOURCE','HGNC_ID','BIOTYPE','CANONICAL',
        'CCDS','ENSP','SWISSPROT','TREMBL','UNIPARC','RefSeq','SIFT','PolyPhen','Exon','INTRON','DOMAINS','AF','AFR_AF',
        'AMR_AF','ASN_AF','EAS_AF','EUR_AF','SAS_AF','AA_AF','EA_AF','CLIN_SIG','SOMATIC','PUBMED','MOTIF_NAME','MOTIF_POS',
        'HIGH_INF_POS','MOTIF_SCORE_CHANGE','IMPACT','PICK','TSL','HGVS_OFFSET','PHENO','MINIMISED','GENE_PHENO','FILTER',
        'flanking_bps','vcf_id','vcf_qual','gnomAD_AF','gnomAD_AFR_AF','gnomAD_AMR_AF','gnomAD_ASJ_AF','gnomAD_EAS_AF',
        'gnomAD_FIN_AF','gnomAD_NFE_AF','gnomAD_OTH_AF','gnomAD_SAS_AF']
    filloutfile.index = [str(filloutfile['Ref'][i]) + ':' + str(int(filloutfile['Start'][i])) + ':' + 
                         str(filloutfile['Alt'][i]) for i in range(len(filloutfile))]
    
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
        if 'chrM' in currheader:
            sequence = [base for base in currsequence]
    varref = [variants[0] for variants in pd.Series(filloutfile.index.values).str.split(':')]
    varpos = [variants[1] for variants in pd.Series(filloutfile.index.values).str.split(':')]
    varalt = [variants[2] for variants in pd.Series(filloutfile.index.values).str.split(':')]
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
    mutsigfile.to_csv(resultsdir + "/" + tumor_id + '_mutsig.tsv',sep = '\t')
    # Calculate heteroplasmy
    filloutfile["Heteroplasmy"] = filloutfile['T_AltCount'].astype(int) / filloutfile['T_TotalDepth'].astype(int)
    # Combined matrices together
    filloutfile.to_csv(f"{resultsdir}/{tumor_id}.bam.maf",sep = '\t',na_rep='',index=False)


if __name__ == "__main__":
    # Parse necessary arguments
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-t", "--tumor_id",type=str, help="(REQUIRED) Path of the tumor sample", required=True)
    parser.add_argument("-w", "--workingdir", type=str, help="(REQUIRED) Working directory with scripts", required=True)
    parser.add_argument("-re", "--resultsdir", type=str, help="(REQUIRED) Directory for results", required=True)
    parser.add_argument("-r", "--reffile",type=str, help="Reference fasta file",default="")
    parser.add_argument("-g", "--genome",type=str, help="Genome version, default=GRCh37",default = "GRCh37")    
    parser.add_argument("-q","--mapq",type=int,help="Minimum mapping quality, default = 20",default = 20)
    parser.add_argument("-Q","--baseq",type=int,help="Minimum base quality, default = 20",default = 20)
    parser.add_argument("-s","--strand",type=int,help="Minimum number of reads mapping to forward and reverse strand to call mutation, default=2",default = 2)
    parser.add_argument("-th","--threshold",type=int,help="The critical threshold for calling a cell wild-type, default=0.1",default = 0.1)
    parser.add_argument("-vc", "--vepcache", type=str, help="Directory for vep cache", default="$HOME/.vep")
    parser.add_argument("-n", "--normal_id", type=str, help="Path of the normal sample",default="")
    parser.add_argument("-m", "--molecule",type=str, help="Type of molecule (dna or rna), default=dna", default="dna")
    parser.add_argument("-c","--mincounts",type=int,help="Minimum number of read counts for MTvariantpipeline, default = 100", default = 100)

    # read in arguments
    args = parser.parse_args()
    tumor_id = args.tumor_id
    reffile = args.reffile
    genome = args.genome
    minmapq = args.mapq
    minbq = args.baseq
    minstrand = args.strand
    threshold = args.threshold
    workingdir = args.workingdir
    vepcache = args.vepcache
    resultsdir = args.resultsdir
    normal_id = args.normal_id
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
    
    # Fix the Tumor ID
    if tumor_id.find('/') != -1:
        # It is a path
        tumordir = os.path.dirname(tumor_id)
    else:
        # It is a file name
        tumordir = os.getcwd()
    # Get the basename for Tumor ID
    tumor_id = os.path.basename(tumor_id)
    # Remove the bam suffix if there is any
    if tumor_id.endswith('.bam'):
        tumor_id = tumor_id[:-4]
    
    # Run reference_detect to determine mtchrom
    mtchrom = reference_detect(reffile)
    
    # Noting all the parameters
    print("Miminum mapping quality of " + str(minmapq))
    print("Miminum base quality of " + str(minbq))
    print("Miminum number of reads mapping to forward and reverse strand to call mutation of " + str(minstrand))
    print("Reference file: " + reffile)

    # Filtering of cells
    if normal_id != "":
        # Fix the Normal ID
        if normal_id.find('/') != -1:
            # It is a path
            normaldir = os.path.dirname(normal_id)
        else:
            # It is a file name
            normaldir = os.getcwd()
        # Get the basename for Normal ID
        normal_id = os.path.basename(normal_id)
        if normal_id.endswith('.bam'):
            normal_id = normal_id[:-4]
        variant_calling_normal(resultsdir,tumordir,tumor_id,reffile,genome,minmapq,minbq,minstrand,workingdir,
                               vepcache,mtchrom,ncbibuild,species,normal_id,normaldir,molecule,mincounts)
    else:
        variant_calling(resultsdir,tumordir,tumor_id,reffile,genome,minmapq,minbq,minstrand,workingdir,
                        vepcache,mtchrom,ncbibuild,species,molecule,mincounts)
    variant_processing(tumor_id,resultsdir)
    print("DONE WITH BULKPIPELINE")