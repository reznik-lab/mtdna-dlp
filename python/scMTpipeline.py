#!/usr/bin/env python3

import os
import argparse
import re
from re import T, sub
from tkinter.tix import Tree
import pandas as pd
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import numpy as np
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
import glob
from pybedtools import BedTool



def reference_detect(reffile):
    print("Determining the mtDNA chromosome name..")
    for sequence in SeqIO.parse(open(reffile), "fasta"):
        if re.search('MT', sequence.description):
            mtchrom = 'MT'
        elif re.search('chrM', sequence.description):
            mtchrom = 'chrM'

def merging_bams(datadir,libraryid,resultsdir):
    print("Merging the cells..")
    try:
        os.makedirs(f"{resultsdir}/merged")
    except OSError:
        pass
    # Merging filtered cells into pseudobulk and indexing the merged file
    subprocess.call("samtools merge " + resultsdir + "/merged/" + libraryid + "-merged.bam " + datadir + "/" + "*.bam", shell=True) 
    subprocess.call("samtools index " + resultsdir + "/merged/" + libraryid + "-merged.bam", shell=True)
    

def preproccess_bams(datadir, reffile, workingdir, vepcache, resultsdir, genome, species, ncbibuild):
    print("Preproccessing bams...")
    try:
        os.makedirs(f"{resultsdir}/MuTect2_results")
    except OSError:
        pass
    try:
        os.makedirs(f"{resultsdir}/MTvariant_results")
    except OSError:
        pass
    try:
        os.makedirs(f"{resultsdir}/filteredfiles")
    except OSError:
        pass
    
    # Run reference_detect to determine mtchrom
    reference_detect(reffile)
    
    for file in os.listdir(datadir):
        if file.endswith(".bam"):
            libraryid = file[:-4]
            print("Running MTvariantpipeline..")
            subprocess.call("python3 " + workingdir + "/MTvariantpipeline.py -d " + datadir + "/ -v " + resultsdir + "/TEMPMAFfiles/ -o " + 
                resultsdir + "/MTvariant_results/ -b " + libraryid + ".bam -g " + genome + " -q " + str(minmapq) + " -Q " + str(minbq) + 
                " -s " + str(minstrand) + " -w " + workingdir + "/ -vc " + vepcache + " -f " + reffile, shell=True)
            # MuTect2 mitochondrial mode
            print("Running MuTect2..")
            subprocess.call("gatk --java-options -Xmx4g Mutect2 -R " + reffile + " --mitochondria-mode true -L " + mtchrom + " -mbq " + str(minbq) + 
                " --minimum-mapping-quality " + str(minmapq) + " -I " + datadir + "/" + libraryid + ".bam -tumor " + libraryid.replace("-","_") + 
                " -O " + resultsdir + "/MuTect2_results/" + libraryid + ".bam.vcf.gz", shell=True)
            # Left align MuTect2 results
            subprocess.call("bcftools norm -m - -f " + reffile + " " + resultsdir + "/MuTect2_results/" + libraryid + ".bam.vcf.gz" + 
                " -o " + resultsdir + "/MuTect2_results/" + libraryid + ".bam.vcf", shell=True)
            # Convert the MuTect2 result from vcf to maf file
            subprocess.call("perl " + workingdir + "/vcf2maf/vcf2maf.pl --species " + species + " --vep-data " + vepcache + " --input-vcf " + 
                resultsdir + "/MuTect2_results/" + libraryid + ".bam.vcf" + " --output-maf " + resultsdir + "/MuTect2_results/" + libraryid + 
                ".bam.maf" + " --ncbi-build " + ncbibuild + ' --ref-fasta ' + reffile, shell=True)
            # Create filtered files
            subprocess.call(f"samtools view -h {datadir}/{file} | grep -h 'X0\|@' > {resultsdir}/filteredfiles/{file}.sam", shell=True)
            subprocess.call(f"samtools view {datadir}/{file} >> {resultsdir}/filteredfiles/{file}.sam", shell=True)
            subprocess.call(f"samtools view -bSq 20 {resultsdir}/filteredfiles/{file}.sam > {resultsdir}/filteredfiles/filtered{file}", shell=True)
            # subprocess.call(f"rm {resultsdir}/filteredfiles/{file}.sam")

    
def variant_calling(datadir,libraryid,reffile,genome,minmapq,minbq,minstrand,workingdir,vepcache,resultsdir,mtchrom,species,ncbibuild):
    try:
        os.makedirs(f"{resultsdir}/MuTect2_results")
    except OSError:
        pass
    try:
        os.makedirs(f"{resultsdir}/MTvariant_results")
    except OSError:
        pass

    # MTvariantpipeline without matched normal
    print("Running MTvariantpipeline..")

    # Running MTvariantpipeline
    subprocess.call("python3 " + workingdir + "/MTvariantpipeline.py -d " + resultsdir + "/merged/ -v " + resultsdir + "/mergedTEMPMAFfiles/ -o " + 
        resultsdir + "/MTvariant_results/ -b " + libraryid + "-merged.bam -g " + genome + " -q " + str(minmapq) + " -Q " + 
        str(minbq) + " -s " + str(minstrand) + " -w " + workingdir + "/ -vc " + vepcache + " -f " + reffile + " -m " + mtchrom, shell=True)

    # MuTect2 mitochondrial mode
    print("Running MuTect2..")
    subprocess.call("gatk --java-options -Xmx4g Mutect2 -R " + reffile + " --mitochondria-mode true -L " + mtchrom + " -mbq " + str(minbq) + 
        " --minimum-mapping-quality " + str(minmapq) + " -I " + resultsdir + "/merged/" + libraryid + "-merged.bam -tumor " + 
        libraryid.replace("-","_") + " -O " + resultsdir + "/MuTect2_results/" + libraryid + "-merged.bam.vcf.gz", shell=True)
    
    # Left align MuTect2 results
    subprocess.call("bcftools norm -m - -f " + reffile + " " + resultsdir + "/MuTect2_results/" + libraryid + "-merged.bam.vcf.gz" + 
        " -o " + resultsdir + "/MuTect2_results/" + libraryid + "-merged.bam.vcf", shell=True)

    # Convert the MuTect2 result from vcf to maf file
    subprocess.call("perl " + workingdir + "/vcf2maf/vcf2maf.pl --species " + species + " --vep-data " + vepcache + " --input-vcf " + 
        resultsdir + "/MuTect2_results/" + libraryid + "-merged.bam.vcf" + " --output-maf " + resultsdir + 
        "/MuTect2_results/" + libraryid + "-merged.bam.maf" + " --ncbi-build " + ncbibuild + ' --ref-fasta ' + reffile, shell=True)


def variant_processing(datadir,libraryid,reffile,patternlist,resultsdir):
    """
    Run MTvariantpipeline and MuTect2 on the filtered cells
    MTvariantpipeline: A simple variant calling and annotation pipeline for mitochondrial DNA variants.
    """
    # Overlap between MuTect and MTvariantpipeline
    # Read in MTvariantpipeline result
    MTvarfile = pd.read_csv(resultsdir + "/MTvariant_results/" + libraryid + "-merged.bam.maf", sep = "\t", comment='#', low_memory=False)
    
    # Read in MuTect result
    mutectfile = pd.read_csv(resultsdir + "/MuTect2_results/" + libraryid + "-merged.bam.maf", sep = "\t", header=1, low_memory=False)

    # Filter out variants falling in the repeat regions of 302-315, 513-525, and 3105-3109 (black listed regions)
    rmregions = list(range(301,314)) + list(range(513,524)) + list(range(3105,3109))
    if len(mutectfile['Start_Position'][mutectfile['Start_Position'].isin(rmregions)]) > 0:
        mutectfile = mutectfile[~mutectfile['Start_Position'].isin(rmregions)]
    if len(MTvarfile['Start_Position'][MTvarfile['Start_Position'].isin(rmregions)]) > 0:
        rmthese = MTvarfile['Start_Position'].isin(rmregions)
        MTvarfile = MTvarfile[~rmthese]
    mutectfile.index = range(len(mutectfile.index))
    MTvarfile.index = range(len(MTvarfile.index))

    # combine pseudo-bulk results
    MTvarfile = pd.merge(mutectfile, MTvarfile, how='inner', on=['Chromosome','Start_Position','Reference_Allele',
        'Tumor_Seq_Allele2','Variant_Classification',"Variant_Type",'EXON'])

    # Import the reference fasta file
    fasta_sequences = SeqIO.parse(open(reffile),'fasta')
    for fasta in fasta_sequences:
        currheader, currsequence = fasta.id, fasta.seq
        if 'MT' in currheader:
            sequence = [base for base in currsequence]
    
    # Fix the INDEL positions and alleles and find homopolymers for re-calculating the reference read counts
    indels = pd.DataFrame(index=range(len([i for i,x in enumerate((MTvarfile['Variant_Type'] == 'INS') | (MTvarfile['Variant_Type'] == 'DEL')) if x])), columns = [0,1])
    indels = indels.fillna(0)
    eachindel = 0
    # Iterate through deletions
    for deletion in [ i for i,x in enumerate(MTvarfile['Variant_Type'] == 'DEL') if x]:
        # Fix the position to match the mutect position
        MTvarfile.loc[MTvarfile.index.values[deletion],'Start_Position'] = MTvarfile.loc[MTvarfile.index.values[deletion],'Start_Position'] - 1
        # Fix the ref allele and the Tumor_Seq_Allele1 using the already changed start position
        MTvarfile.loc[MTvarfile.index.values[deletion],'Reference_Allele'] = sequence[int(MTvarfile.loc[MTvarfile.index.values[deletion],'Start_Position'])-1] + MTvarfile.loc[MTvarfile.index.values[deletion],'Reference_Allele']
        MTvarfile.loc[MTvarfile.index.values[deletion],'Tumor_Seq_Allele1'] = MTvarfile.loc[MTvarfile.index.values[deletion],'Reference_Allele']
        # current allele
        currallele = sequence[int(MTvarfile.loc[MTvarfile.index.values[deletion],'Start_Position']) - 1]
        i = 0
        while (sequence[int(MTvarfile.loc[MTvarfile.index.values[deletion],'Start_Position']) + i - 1] == currallele):
            i = i + 1
        indels.iloc[eachindel,:] = [int(MTvarfile.loc[MTvarfile.index.values[deletion],'Start_Position']),int(MTvarfile.loc[MTvarfile.index.values[deletion],'Start_Position'] + i + 1)]
        eachindel = eachindel + 1
        # Fix the alt allele
        MTvarfile.loc[MTvarfile.index.values[deletion],'Tumor_Seq_Allele2'] = sequence[int(MTvarfile.loc[MTvarfile.index.values[deletion],'Start_Position']) - 1]
    # Iterate through insertions
    for insertion in [ i for i,x in enumerate(MTvarfile['Variant_Type'] == 'INS') if x]:
        # Fix the ref allele and the Tumor_Seq_Allele1
        MTvarfile.loc[MTvarfile.index.values[insertion],'Reference_Allele'] = sequence[int(MTvarfile.loc[MTvarfile.index.values[insertion],'Start_Position']) - 1]
        MTvarfile.loc[MTvarfile.index.values[insertion],'Tumor_Seq_Allele1'] = sequence[int(MTvarfile.loc[MTvarfile.index.values[insertion],'Start_Position']) - 1]
        # current allele
        currallele = MTvarfile.loc[MTvarfile.index.values[insertion],'Reference_Allele']
        i = 0
        while (sequence[int(MTvarfile.loc[MTvarfile.index.values[insertion],'Start_Position']) + i - 1] == currallele):
            i = i + 1
        indels.iloc[eachindel,:] = [int(MTvarfile.loc[MTvarfile.index.values[insertion],'Start_Position']),int(MTvarfile.loc[MTvarfile.index.values[insertion],'Start_Position'])+i+1]
        eachindel = eachindel + 1
        # Fix the alt allele
        MTvarfile.loc[MTvarfile.index.values[insertion],'Tumor_Seq_Allele2'] = sequence[int(MTvarfile.loc[MTvarfile.index.values[insertion],'Start_Position']) - 1] + MTvarfile.loc[MTvarfile.index.values[insertion],'Tumor_Seq_Allele2']

    # Save the indel annotations in a BedTool format
    indelfile = ''
    for eachindels in range(len(indels.index)):
        indelfile = indelfile + 'MT\t'+ str(indels.iloc[eachindels,0]) + '\t' + str(indels.iloc[eachindels,1]) + '\t.\n'
    indels_res = BedTool(indelfile, from_string=True).sort().merge()
    indels = pd.read_csv(indels_res.fn, names=['0', '1'],sep='\t')
    indivcol = pd.DataFrame(index=MTvarfile.index.values)
    indivcol = indivcol.fillna(0)
    
    # If the list exists, then read the list
    if patternlist != "":
        processbams = pd.read_csv(resultsdir + "/" + patternlist, header = None, low_memory=False)[0]
        saveasthis = resultsdir + "/" + libraryid + "-merged_" + patternlist + ".fillout"
    else:
        processbams = glob.glob(resultsdir + '/TEMPMAFfiles/*.bam_temp.maf')
        saveasthis = resultsdir + "/" + libraryid + "-merged.fillout"
    for filepath in processbams:
        file = os.path.basename(filepath)
        if os.stat(os.path.join(resultsdir + '/TEMPMAFfiles/' + file)).st_size != 0:
            print("Processing " + file)
            indivfile = pd.read_csv(os.path.join(resultsdir + '/TEMPMAFfiles/' + file), header = None, sep = "\t")
            # for each bam
            for eachrow in MTvarfile.index.values:
                # if the position is present
                if any(indivfile.iloc[:,1] == int(MTvarfile.loc[eachrow,'Start_Position'])):
                    # if the position is present and it is an SNV
                    if MTvarfile.loc[eachrow,'Variant_Type'] == 'SNP':
                        posmatch = indivfile.index[indivfile.iloc[:,1] == int(MTvarfile.loc[eachrow,'Start_Position'])]
                        # if the position is present and it is an SNV and the alt read matches
                        if MTvarfile.loc[eachrow,'Tumor_Seq_Allele2'] in indivfile.iloc[indivfile.index[indivfile.iloc[:,1] == int(MTvarfile.loc[eachrow,'Start_Position'])],3].values:
                            corrrow = posmatch[indivfile.iloc[posmatch,3] == MTvarfile.loc[eachrow,'Tumor_Seq_Allele2']]
                            dp = int(indivfile.iloc[corrrow][4].values[0].split(',')[0]) + int(indivfile.iloc[corrrow][4].values[0].split(',')[1]) # int(indivfile.iloc[corrrow][5]) # 
                            rd = int(indivfile.iloc[corrrow][4].values[0].split(',')[0])
                            ad = int(indivfile.iloc[corrrow][4].values[0].split(',')[1])
                            indivcol.loc[eachrow,str(file.split('_')[0])] = "DP=" + str(dp) + ";" + "RD=" + str(rd) + ";" + "AD="  + str(ad)
                        # if the position is present and it is an SNV but there are no alt reads
                        elif "." in indivfile.iloc[indivfile.index[indivfile.iloc[:,1] == int(MTvarfile.loc[eachrow,'Start_Position'])],3].values:
                            corrrow = posmatch[indivfile.iloc[posmatch,3] == "."]
                            dp = int(indivfile.iloc[indivfile.index[(indivfile.iloc[:,1] == int(MTvarfile.loc[eachrow,'Start_Position'])) & (indivfile.iloc[:,3] == ".")].tolist()][4])
                            rd = int(indivfile.iloc[indivfile.index[(indivfile.iloc[:,1] == int(MTvarfile.loc[eachrow,'Start_Position'])) & (indivfile.iloc[:,3] == ".")].tolist()][4])
                            ad = 0
                            indivcol.loc[eachrow,str(file.split('_')[0])] = "DP=" + str(dp) + ";" + "RD=" + str(rd) + ";" + "AD="  + str(ad)
                        # if the position is present and it is an SNV but there are only other multi-allelic reads
                        else:
                            corrrow = indivfile.iloc[posmatch,3].index.values
                            dp = int(indivfile.iloc[max(corrrow),4].split(',')[0])
                            rd = int(indivfile.iloc[max(corrrow),4].split(',')[0])
                            ad = 0
                            indivcol.loc[eachrow,str(file.split('_')[0])] = "DP=" + str(dp) + ";" + "RD=" + str(rd) + ";" + "AD="  + str(ad)
                    # if the position is present and it is an INS, check the alt reads
                    elif MTvarfile.loc[eachrow,'Variant_Type'] == 'INS':
                        posmatch = indivfile.index[indivfile.iloc[:,1] == int(MTvarfile.loc[eachrow,'Start_Position'])]
                        # if the position is present and it is an INDEL and the alt allele matches
                        if sum(indivfile.iloc[posmatch,3].apply(len) == len(MTvarfile.loc[eachrow,'Tumor_Seq_Allele2'])) > 0:
                            # all the alt alleles for INS
                            corrrow = indivfile.iloc[posmatch,3].index[indivfile.iloc[posmatch,3].apply(len) == len(MTvarfile.loc[eachrow,'Tumor_Seq_Allele2'])]
                            # take the last values if accounting for multiallelic
                            dp = int(indivfile.iloc[max(corrrow),4].split(',')[0]) + int(indivfile.iloc[max(corrrow),4].split(',')[1])
                            rd = int(indivfile.iloc[max(corrrow),4].split(',')[0])
                            ad = int(indivfile.iloc[max(corrrow),4].split(',')[1])
                            indivcol.loc[eachrow,str(file.split('_')[0])] = "DP=" + str(dp) + ";" + "RD=" + str(rd) + ";" + "AD="  + str(ad)
                        # if the position is present and it is an INS but there are no alt reads
                        elif "." in indivfile.iloc[indivfile.index[indivfile.iloc[:,1] == int(MTvarfile.loc[eachrow,'Start_Position'])],3].values:
                            for eachindel in range(len(indels.index)):
                                if len(set(range(indels.iloc[eachindel][0],indels.iloc[eachindel][1]+1)) & set([int(MTvarfile.loc[eachrow,'Start_Position'])])) > 0:
                                    startpos = indels.iloc[eachindel][0]
                                    endpos = indels.iloc[eachindel][1]
                            bedfile = BedTool('MT\t'+ str(startpos-1) + '\t' + str(endpos+1) + '\t.', from_string=True)
                            x = 0
                            # for each read
                            for i in BedTool(resultsdir + "/filteredfiles/filtered" + file.replace('_temp.maf','')).intersect(bedfile,F=1):
                                # for each indel position
                                x = x + 1
                            dp = x
                            rd = x
                            ad = 0
                            indivcol.loc[eachrow,str(file.split('_')[0])] = "DP=" + str(dp) + ";" + "RD=" + str(rd) + ";" + "AD="  + str(ad)
                        # if the position is present and it is an INS, but there are only other multi-allelic reads
                        else:
                            corrrow = indivfile.iloc[posmatch,3].index.values
                            dp = int(indivfile.iloc[max(corrrow),4].split(',')[0])
                            rd = int(indivfile.iloc[max(corrrow),4].split(',')[0])
                            ad = 0
                            indivcol.loc[eachrow,str(file.split('_')[0])] = "DP=" + str(dp) + ";" + "RD=" + str(rd) + ";" + "AD="  + str(ad)
                    # if the position is present and it is a DEL, check the alt reads
                    elif MTvarfile.loc[eachrow,'Variant_Type'] == 'DEL':
                        posmatch = indivfile.index[indivfile.iloc[:,1] == int(MTvarfile.loc[eachrow,'Start_Position'])]
                        if sum(indivfile.iloc[posmatch,2].apply(len) == len(MTvarfile.loc[eachrow,'Reference_Allele'])) > 0:
                            corrrow = indivfile.iloc[posmatch,2].index[indivfile.iloc[posmatch,2].apply(len) == len(MTvarfile.loc[eachrow,'Reference_Allele'])]
                            dp = int(indivfile.iloc[max(corrrow),4].split(',')[0]) + int(indivfile.iloc[max(corrrow),4].split(',')[1]) # int(indivfile.iloc[max(corrrow),5])
                            rd = int(indivfile.iloc[max(corrrow),4].split(',')[0])
                            ad = int(indivfile.iloc[max(corrrow),4].split(',')[1])
                            indivcol.loc[eachrow,str(file.split('_')[0])] = "DP=" + str(dp) + ";" + "RD=" + str(rd) + ";" + "AD="  + str(ad)
                        # the position is present but is not a DEL
                        elif "." in indivfile.iloc[indivfile.index[indivfile.iloc[:,1] == int(MTvarfile.loc[eachrow,'Start_Position'])],3].values:
                            for eachindel in range(len(indels.index)):
                                if len(set(range(indels.iloc[eachindel][0],indels.iloc[eachindel][1]+1)) & set([int(MTvarfile.loc[eachrow,'Start_Position'])])) > 0:
                                    ######
                                    startpos = indels.iloc[eachindel][0]
                                    endpos = indels.iloc[eachindel][1]
                            bedfile = BedTool('MT\t'+ str(startpos-1) + '\t' + str(endpos+1) + '\t.', from_string=True)
                            x = 0
                            # for each read
                            for i in BedTool(resultsdir + "/filteredfiles/filtered" + file.replace('_temp.maf','')).intersect(bedfile,F=1):
                                # for each indel position
                                x = x + 1
                            dp = x
                            rd = x
                            ad = 0
                            indivcol.loc[eachrow,str(file.split('_')[0])] = "DP=" + str(dp) + ";" + "RD=" + str(rd) + ";" + "AD="  + str(ad)
                        # if the position is present and it is a DEL, but there are only other multi-allelic reads
                        else:
                            corrrow = indivfile.iloc[posmatch,3].index.values
                            dp = int(indivfile.iloc[max(corrrow),4].split(',')[0])
                            rd = int(indivfile.iloc[max(corrrow),4].split(',')[0])
                            ad = 0
                            indivcol.loc[eachrow,str(file.split('_')[0])] = "DP=" + str(dp) + ";" + "RD=" + str(rd) + ";" + "AD="  + str(ad)
                # if the position is not present
                else:
                    dp = 0
                    rd = 0
                    ad = 0
                    indivcol.loc[eachrow,str(file.split('_')[0])] = "DP=" + str(dp) + ";" + "RD=" + str(rd) + ";" + "AD="  + str(ad)
        else:
            # for bam without any TEMP file
            for eachrow in range(len(MTvarfile.index)):
                dp = 0
                rd = 0
                ad = 0
                indivcol.loc[eachrow,str(file.split('_')[0])] = "DP=" + str(dp) + ";" + "RD=" + str(rd) + ";" + "AD="  + str(ad)

    # Fix rd of the multi-allelic positions for INDELs
    for eachmultpos in MTvarfile['Start_Position'][(MTvarfile['Variant_Type'] == "INS") | (MTvarfile['Variant_Type'] == "DEL")].duplicated().index[MTvarfile['Start_Position'][(MTvarfile['Variant_Type'] == "INS") | (MTvarfile['Variant_Type'] == "DEL")].duplicated()]:
        prevvalues = []
        # Obtaining the values from the first line
        if eachmultpos-1 not in MTvarfile['Start_Position'][(MTvarfile['Variant_Type'] == "INS") | (MTvarfile['Variant_Type'] == "DEL")].duplicated().index[MTvarfile['Start_Position'][(MTvarfile['Variant_Type'] == "INS") | (MTvarfile['Variant_Type'] == "DEL")].duplicated()]:
            newdpvalues = [] # if it is first line, start with an empty new dp values
            # Iterate through multi-alleleic positions in the fillout file
            for eachcell in range(len(indivcol.loc[int(MTvarfile.index.values[[i-1 for i,x in enumerate(MTvarfile.index.values == eachmultpos) if x]]),:].values)):
                start = int(indivcol.loc[int(MTvarfile.index.values[[i-1 for i,x in enumerate(MTvarfile.index.values == eachmultpos) if x]]),:].str.find("RD=").values[eachcell]+3)
                end = int(indivcol.loc[int(MTvarfile.index.values[[i-1 for i,x in enumerate(MTvarfile.index.values == eachmultpos) if x]]),:].str.find(";A").values[eachcell])
                prevvalues = prevvalues + [indivcol.loc[int(MTvarfile.index.values[[i-1 for i,x in enumerate(MTvarfile.index.values == eachmultpos) if x]]),:][eachcell][start:end]]
                testmultpos = eachmultpos
                newdpvalue = int(indivcol.loc[int(MTvarfile.index.values[[i-1 for i,x in enumerate(MTvarfile.index.values == eachmultpos) if x]]),:][eachcell][start:end]) + int(indivcol.loc[int(MTvarfile.index.values[[i-1 for i,x in enumerate(MTvarfile.index.values == eachmultpos) if x]]),:][eachcell][(end+4):]) # int(MTvarfile[MTvarfile.index == eachmultpos]['Start_Position'])-1
                while testmultpos in MTvarfile['Start_Position'][(MTvarfile['Variant_Type'] == "INS") | (MTvarfile['Variant_Type'] == "DEL")].duplicated().index[MTvarfile['Start_Position'][(MTvarfile['Variant_Type'] == "INS") | (MTvarfile['Variant_Type'] == "DEL")].duplicated()]:
                    newdpvalue = newdpvalue + int(indivcol.loc[testmultpos,:][eachcell][int(indivcol.loc[testmultpos,:].str.find(";A").values[eachcell]+4):])
                    testmultpos = testmultpos + 1
                newdpvalues = newdpvalues + [newdpvalue]
            # Updating the values from the subsequent lines based on the first line
            if eachmultpos-1 not in MTvarfile['Start_Position'][(MTvarfile['Variant_Type'] == "INS") | (MTvarfile['Variant_Type'] == "DEL")].duplicated().index[MTvarfile['Start_Position'][(MTvarfile['Variant_Type'] == "INS") | (MTvarfile['Variant_Type'] == "DEL")].duplicated()]:
                for eachcell in range(len(indivcol.loc[int(MTvarfile.index.values[[i-1 for i,x in enumerate(MTvarfile.index.values == eachmultpos) if x]]),:].values)):
                    start = int(indivcol.loc[int(MTvarfile.index.values[[i-1 for i,x in enumerate(MTvarfile.index.values == eachmultpos) if x]]),:].str.find("RD=").values[eachcell]+3)
                    end = int(indivcol.loc[int(MTvarfile.index.values[[i-1 for i,x in enumerate(MTvarfile.index.values == eachmultpos) if x]]),:].str.find(";A").values[eachcell])
                    indivcol.loc[int(MTvarfile.index.values[[i-1 for i,x in enumerate(MTvarfile.index.values == eachmultpos) if x]]),:][eachcell] = indivcol.loc[int(MTvarfile.index.values[[i-1 for i,x in enumerate(MTvarfile.index.values == eachmultpos) if x]]),:][eachcell].replace('DP=' + indivcol.loc[int(MTvarfile.index.values[[i-1 for i,x in enumerate(MTvarfile.index.values == eachmultpos) if x]]),:][eachcell][3:indivcol.loc[int(MTvarfile.index.values[[i-1 for i,x in enumerate(MTvarfile.index.values == eachmultpos) if x]]),:][eachcell].find(";R")], 'DP=' + str(newdpvalues[eachcell]))
            for eachcell in range(len(indivcol.loc[eachmultpos,:].values)):
                start = int(indivcol.loc[eachmultpos,:].str.find("RD=").values[eachcell]+3)
                end = int(indivcol.loc[eachmultpos,:].str.find(";A").values[eachcell])
                tempres = indivcol.loc[eachmultpos,:][eachcell].replace(';RD=' + str(indivcol.loc[eachmultpos,:][eachcell][start:end]), ';RD=' + prevvalues[eachcell])
                indivcol.loc[eachmultpos,:][eachcell] = tempres.replace('DP=' + tempres[3:tempres.find(";R")], 'DP=' + str(newdpvalues[eachcell]))
    
    indivsumcol = pd.DataFrame(index=MTvarfile.index.values)
    indivsumcol = indivsumcol.fillna(0)
    for eachpos in indivsumcol.index.values:
        sumDP = 0
        sumRD = 0
        sumAD = 0
        for eachcell in range(len(indivcol.columns)):
            sumDP = sumDP + int(indivcol.loc[eachpos][eachcell][3:int(indivcol.loc[eachpos][eachcell].find(';RD='))])
            sumRD = sumRD + int(indivcol.loc[eachpos][eachcell][int(indivcol.loc[eachpos][eachcell].find(';RD=')+4):int(indivcol.loc[eachpos][eachcell].find(';AD='))])
            sumAD = sumAD + int(indivcol.loc[eachpos][eachcell][int(indivcol.loc[eachpos][eachcell].find(';AD=')+4):])
        indivsumcol.loc[eachpos,'S_TotalDepth'] = sumDP
        indivsumcol.loc[eachpos,'S_RefCount'] = sumRD
        indivsumcol.loc[eachpos,'S_AltCount'] = sumAD
    
    # MTvarfile.to_csv("/home/parkt/sc_combined.tsv",sep='\t',index=False)

    # Final annotation
    # final_result = MTvarfile.loc[:,['Tumor_Sample_Barcode','Matched_Norm_Sample_Barcode','Chromosome',
    #     'Start_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','Hugo_Symbol','EXON',
    #     'n_depth','t_depth','t_ref_count','t_alt_count']]
    # final_result.columns = ['Sample','NormalUsed','Chrom','Start','Ref','Alt','VariantClass','Gene','Exon',
    #     'N_TotalDepth','T_TotalDepth','T_RefCount','T_AltCount']
    final_result = MTvarfile.loc[:,['Tumor_Sample_Barcode_y','Matched_Norm_Sample_Barcode_y','Chromosome',
        'Start_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','Hugo_Symbol_y','EXON',
        'n_depth_y',"n_ref_count_y","n_alt_count_y",'t_depth_y','t_ref_count_y','t_alt_count_y',"t_alt_fwd","t_alt_rev"]]
    final_result.columns = ['Sample','NormalUsed','Chrom','Start','Ref','Alt','VariantClass','Gene','Exon',
        'N_TotalDepth',"N_RefCount","N_AltCount",'T_TotalDepth','T_RefCount','T_AltCount',"T_AltFwd","T_AltRev"]
    indivsumcol.index = indivcol.index.values
    final_result = pd.concat([final_result, indivsumcol, indivcol], axis = 1)
    
    # output the fillout results
    final_result.to_csv(saveasthis,sep = '\t',na_rep='NA',index=False)
    

def runhaplogrep(datadir,libraryid,reffile,workingdir,resultsdir):
    """
    Run haplogrep to obtain the haplogroup information from the merged bam file
    """
    print("Preparing haplogrep..")
    
    # Filter the bam file for unmapped reads and mapping quality less than 1
    subprocess.call("samtools view -bF 4 -q 1 " + resultsdir + "/merged/" + libraryid + "-merged.bam > " + 
        resultsdir + "/filtered" + libraryid + "-merged.bam", shell=True)

    # Index the filtered bam file
    subprocess.call("samtools index " + resultsdir + "/filtered" + libraryid + "-merged.bam", shell=True)
    
    # Edit the RG of the filtered bam file
    subprocess.call("java -Xms8G -Xmx8G -jar " + workingdir + "/reference/picard.jar AddOrReplaceReadGroups I=" + 
        resultsdir + "/filtered" + libraryid + "-merged.bam O=" + resultsdir + "/merged/result" + libraryid + 
        "-merged.bam RGID=" + libraryid.replace("-", "_") + " RGLB=" + libraryid + 
        " RGPL=illumina RGPU=unit1 RGSM=" + libraryid, shell=True)

    # Index the resulting bam file
    subprocess.call("samtools index " + resultsdir + "/merged/result" + libraryid + "-merged.bam", shell=True)
    
    # Run MuTect2
    subprocess.call("gatk --java-options -Xmx4g Mutect2 -R " + reffile + " --mitochondria-mode true -L MT -mbq " + str(minbq) + 
        " --minimum-mapping-quality " + str(minmapq) + " -I " + resultsdir + "/merged/result"  + libraryid + "-merged.bam -tumor result" + 
        libraryid.replace("-","_") + " -O " + resultsdir + "/MuTect2_results/result" + libraryid + "-merged.bam.vcf.gz", shell=True)

    # Run haplogrep2.1
    subprocess.call("java -jar " + workingdir + "/reference/haplogrep/haplogrep-2.1.20.jar --in " + resultsdir + 
        "/MuTect2_results/result" + libraryid + "-merged.bam.vcf.gz" + " --format vcf --extend-report --out " + 
        resultsdir + "/" + libraryid + "_haplogroups.txt", shell=True)


def calcprob(depth,hetprob):
    """
    Calculates probability of mutation given read depth and bulk heteroplasmy. done as a matrix operation
    depth = number of reads
    hetprob = mean heteroplasmy of mutation across clearly mutant cells
    """
    wtprob = 1-hetprob
    mutprobs = (wtprob)**depth
    return mutprobs


def splitfout(str,flip):
    """
    splits the fillout data
    """
    splits = str.split(';')
    dp = int(float(splits[0].split('=')[1]))
    rd = int(float(splits[1].split('=')[1]))
    ad = int(float(splits[2].split('=')[1]))
    if dp == 0:
        vaf = float('nan')
    elif flip:
        vaf = rd/dp
    else:
        vaf = ad/dp
    return [{'vaf': float(vaf),'dp' : dp}]


def makeMTdf(filloutfile):
    """
    Making the data frame of VAF and depth
    """
    fout = pd.read_csv(filloutfile, sep='\t')
    
    # Set rownames
    fout.index = [str(fout['Ref'][i]) + ':' + str(int(fout['Start'][i])) + ':' + str(fout['Alt'][i]) for i in range(len(fout))]
    
    # Get cells
    if 'S_AltCount' in fout.columns:
        startix = fout.columns.get_loc('S_AltCount') + 1
        endix = fout.shape[1]
        cids = fout.columns[startix:endix]
    else:
        startix = fout.columns.get_loc('Exon') + 1
        endix = fout.shape[1]
        cids = fout.columns[startix:endix]
    
    # Make several dataframes
    vaf = np.zeros((len(fout.index),len(cids)))
    vaf = pd.DataFrame(vaf)
    vaf.columns = cids
    depth = np.zeros((len(fout.index),len(cids)))
    depth = pd.DataFrame(depth)
    depth.columns = cids
    
    # Iterate through the variants and calculate the VAF and read depth
    for varid in range(len(fout.index)):
        # print(fout.iloc[varid,startix:endix])
        tempres = fout.iloc[varid,startix:endix].apply(splitfout,args=(False,))
        vaf.iloc[varid] = [cell[0]['vaf'] for cell in tempres]
        depth.iloc[varid] = [cell[0]['dp'] for cell in tempres]
    
    # Store the row names and column names
    vaf.index = fout.index.values
    depth.index = fout.index.values
    
    # Return VAF and read depth info
    return [{'vaf': vaf,'depth' : depth}]


def processfillout(libraryid,threshold,resultsdir,genome):
    """
    Run the combined mutation estimation on fillout
    Post-processing of the fillout files
    threshold: the critical threshold for calling a cell wild-type
    """
    print("Running the mutation estimation on the fillout..")
    
    # Import the final fillout file
    filloutfile = os.path.join(resultsdir + "/" + libraryid + '-merged.fillout')
    res = makeMTdf(filloutfile)
    
    # Import haplogrep result
    if genome == "GRCh38" or genome == "GRCh37":
        haplogrepfile = pd.read_csv(os.path.join(resultsdir + "/" + libraryid + '_haplogroups.txt'), sep='\t')
        germlinepos = [x[:-1] for x in haplogrepfile['Found_Polys'][0].split(" ")]

    # Remove variants that are not present in all cells
    removethese = [i for i, x in enumerate(res[0]['vaf'].sum(axis = 1, skipna = True) == 0) if x]
    if len(removethese) > 0:
        res[0]['vaf'].drop(res[0]['vaf'].index[removethese])
        res[0]['depth'].drop(res[0]['depth'].index[removethese])
    
    # Output VAF and read depth files
    res[0]['vaf'].to_csv(resultsdir + "/" + libraryid + '_vaf.tsv',sep = '\t',na_rep='NA')
    res[0]['depth'].to_csv(resultsdir + "/" + libraryid + '_depth.tsv',sep = '\t')
    
    # Initialize mutation probability matrix
    mutprob = np.zeros((len(res[0]['vaf'].index),len(res[0]['vaf'].columns)))
    mutprob = pd.DataFrame(mutprob)
    mutprob.columns = res[0]['vaf'].columns
    mutprob.index = res[0]['vaf'].index.values
    
    # Iterate through variants and calculate the heteroplasmy for the cells with an evident mutation
    for varid in mutprob.index.values:
        # Identify mutant cells
        mutcells = res[0]['vaf'].columns[res[0]['vaf'].loc[varid] > 0]
        # Skip the variant if no mutant cells
        if len(mutcells) == 0:
            print('Skipping ' + varid)
        # Identify wildtype cells
        wtcells = res[0]['vaf'].columns[res[0]['vaf'].loc[varid] == 0]
        # Calculate the heteroplasmy for the mutant cells only
        het = res[0]['vaf'].loc[varid,mutcells].mean(skipna = True)
        tempres = calcprob(res[0]['depth'].loc[varid,wtcells], het)
        mutprob.loc[varid,mutcells] = 1 # Mutant cells have probability of 1
        # Calculate the probability by calling calcprob function
        mutprob.loc[varid,wtcells] = tempres.loc[wtcells]

    # Output mutation probability matrix
    mutprob.to_csv(resultsdir + "/" + libraryid + '_mutprob.tsv',sep = '\t')
    
    # Calculate the average read depth for each variant
    ad = round(res[0]['vaf']*res[0]['depth'])
    filteredvar = np.zeros((0,0))
    filteredvar = pd.DataFrame(filteredvar)
    
    # Calculate bulk heteroplasmy and calculate proportion of mutant, WT cells for each variant
    for varid in mutprob.index.values:
        filteredvar.loc[varid,'bulk'] = ad.loc[varid,].sum(skipna = True) / res[0]['depth'].loc[varid,].sum(skipna = True)
        filteredvar.loc[varid,'mutant'] = sum(mutprob.loc[varid] == 1)
        filteredvar.loc[varid,'wt'] = sum(mutprob.loc[varid] < threshold)
        filteredvar.loc[varid,'avgindiv'] = res[0]['vaf'].loc[varid,res[0]['vaf'].columns[res[0]['vaf'].loc[varid] != 0]].mean(skipna = True)
        filteredvar.loc[varid,'stdindiv'] = res[0]['vaf'].loc[varid,res[0]['vaf'].columns[res[0]['vaf'].loc[varid] != 0]].std(skipna = True)
    filteredvar['mutprop'] = filteredvar['mutant'] / (filteredvar['mutant'] + filteredvar['wt'])
    filteredvar['numcells'] = filteredvar['mutant'] + filteredvar['wt']
    filteredvar['detectprop'] = filteredvar['numcells']/res[0]['vaf'].shape[1]
    prevpos = 0

    # Assign variants with >95% VAF as germline if they are used in haplogroup assignment and as homoplasmic otherwise
    for varid in mutprob.index.values:
        # print(varid.split(":")[1])
        if int(varid.split(":")[1]) != prevpos and filteredvar.loc[varid,'bulk'] >= 0.95:
            if (genome == "GRCh38" or genome == "GRCh37") and varid.split(":")[1] in germlinepos:
                filteredvar.loc[varid,'somaticstatus'] = 'germline'
            else:
                filteredvar.loc[varid,'somaticstatus'] = 'homoplasmic'
        else:
            filteredvar.loc[varid,'somaticstatus'] = 'somatic'
        prevpos = int(varid.split(":")[1])

    # Output filtered variant file
    filteredvar.to_csv(resultsdir + "/" + libraryid + '_variants.tsv',sep = '\t')
    
    # Generate a scatterplot with y = x line
    print("Plotting..")
    plt.plot()
    plt.scatter(filteredvar['mutprop'],filteredvar['bulk'],s=filteredvar['numcells']) # , c='blue', alpha=0.05, edgecolors='none')
    plt.plot([0,1], [0,1])
    plt.axvline(x=0.95,linestyle='dashed')
    plt.axhline(y=0.95,linestyle='dashed')
    plt.xlabel("Proportion of Cells with Mutation")
    plt.ylabel("Bulk Heteroplasmy")
    plt.savefig(resultsdir + "/" + libraryid + '_heteroplasmy.pdf')
    plt.clf()
    

def genmaster(libraryid,reffile,resultsdir,genome):
    """
    Run the combined mutation estimation on fillout
    Post-processing of the fillout files
    threshold: the critical threshold for calling a cell wild-type
    """
    print('Generating a master file and a binary matrix of somatic variants for the sample..')
    
    # Import the relevant files
    # allcellsfile = pd.read_csv(os.path.join(allcells), sep='\t')
    variantsfile = pd.read_csv(os.path.join(resultsdir + "/" + libraryid + '_variants.tsv'), sep='\t', index_col=0)
    filloutfile = pd.read_csv(os.path.join(resultsdir + "/" + libraryid + '-merged.fillout'), sep='\t')
    vaffile = pd.read_csv(os.path.join(resultsdir + "/" + libraryid + '_vaf.tsv'), sep='\t', index_col=0)
    depthfile = pd.read_csv(os.path.join(resultsdir + "/" + libraryid + '_depth.tsv'), sep='\t', index_col=0)
    vaffile.index=list(depthfile.index.values)
    varfile = vaffile.mul(depthfile, fill_value=0)
    
    if genome == "GRCh38" or genome == "GRCh37":
        # Obtaining all the unique positions
        allpos = np.array([variants[1] for variants in pd.Series(variantsfile.index.values).str.split(':')])
        _, idx = np.unique(allpos, return_index=True)
        uniqpos = allpos[np.sort(idx)]

        # Flip the ref and the alt allele for the germline variants
        start = [variants[2] for variants in pd.Series(variantsfile.index.values[variantsfile['somaticstatus'] == 'germline']).str.split(':')]
        pos = [variants[1] for variants in pd.Series(variantsfile.index.values[variantsfile['somaticstatus'] == 'germline']).str.split(':')]
        end = [variants[0] for variants in pd.Series(variantsfile.index.values[variantsfile['somaticstatus'] == 'germline']).str.split(':')]
        variantsfile.rename(index=dict(zip(variantsfile.index.values[variantsfile['somaticstatus'] == 'germline'],[x+':'+str(y)+':'+str(z) for x,y,z in zip(start,pos,end)])), inplace=True)
        
        # Flip the ref allele for somatic variants to follow germline variants
        fixthese = []
        for eachpos in uniqpos: # for each unique positions
            curridx = [s for s in variantsfile.index.values if s.split(':')[1] == eachpos] # all the variants at the unique positions
            if eachpos in pos: # if the position is germline position,
                newstart = [[x for x,y,z in zip(start,pos,end)][pos.index(eachpos)] + eachstart[1:] for eachstart in [variants.split(':')[0] + ':' + variants.split(':')[1] + ':' + variants.split(':')[2] for variants in curridx]]
                fixthese.extend(newstart)
            else:
                fixthese.extend(curridx)
        variantsfile.rename(index=dict(zip(variantsfile.index.values,fixthese)), inplace=True)
        
        # Update the row names of other files to be consistent with the variants file
        varfile.index=list(variantsfile.index.values)
        depthfile.index=list(variantsfile.index.values)
        filloutfile.index=list(variantsfile.index.values)
        vaffile.index=list(variantsfile.index.values)
        
        # New ref and alt alleles
        newref = [variants[0] for variants in pd.Series(fixthese).str.split(':')]
        newalt = [variants[2] for variants in pd.Series(fixthese).str.split(':')]
    
    # Interpreting Picard results that are in metrics.txt file
    resultMTcoverage = pd.DataFrame(index=[], columns=['sampleid','MTreadcounts'])
    totfiles = filter(lambda fname: 'temp.maf' in fname, os.listdir(os.path.join(resultsdir + "/TEMPMAFfiles/")))
    for file in totfiles:
        sampleid = file.replace('_MT.bam_temp.maf','')
        currrow = pd.DataFrame([sampleid]).T
        currrow.columns = ['sampleid']
        resultMTcoverage = pd.concat([resultMTcoverage,currrow],axis=0)
    
    # Fix the depth matrix to filter variants that are uncertain and order them based on filteredvariants matrix
    masterfile = pd.DataFrame(index=variantsfile.index.values, columns=depthfile.columns)
    masterfile = masterfile.fillna(0)
    
    # Fix the read counts for individual cells for each row accounting for the germline variants
    masterfile = varfile.astype(int).astype(str) + '/' + depthfile.astype(int).astype(str)
    
    # Create a variant annotations file based on the fillout file
    variantannot = pd.DataFrame(index=filloutfile.index.values, columns=filloutfile[['Start','Ref','Alt',
        'VariantClass','Gene','T_AltCount','T_RefCount','S_AltCount','S_RefCount']].columns)
    variantannot = variantannot.fillna(0)
    
    # Include columns for 'Start','Ref','Alt','VariantClass','Gene','T_AltCount','T_RefCount'
    if genome == "GRCh38" or genome == "GRCh37":
        variantannot['Ref'] = newref
        variantannot['Alt'] = newalt
        variantannot['oldRef'] = filloutfile['Ref']
        variantannot['oldAlt'] = filloutfile['Alt']
        variantannot['Gene'] = filloutfile['Gene']
    elif genome == "GRCm38" or genome == "mm10":
        variantannot['Ref'] = filloutfile['Ref']
        variantannot['Alt'] = filloutfile['Alt']
    variantannot['VariantClass'] = filloutfile['VariantClass']
    variantannot['T_AltCount'] = filloutfile['T_AltCount']
    variantannot['T_RefCount'] = filloutfile['T_RefCount']
    variantannot['Start'] = filloutfile['Start']
    
    # Recalculating the read counts
    if patternlist != "":
        for eachrow in variantannot.index.values:
            if str(variantannot.loc[eachrow,'Start']) in pos: # it's a germline variant position
                # sum the numerators of the masterfile for the row
                variantannot.loc[eachrow,'S_AltCount'] = sum(pd.to_numeric(masterfile.loc[eachrow,:].apply(lambda x : x.split('/')[1]).values)) - sum(pd.to_numeric(masterfile.loc[eachrow,:].apply(lambda x : x.split('/')[0]).values))
                variantannot.loc[eachrow,'S_RefCount'] = sum(pd.to_numeric(masterfile.loc[eachrow,:].apply(lambda x : x.split('/')[0]).values))
            variantannot.loc[eachrow,'S_AltCount'] = sum(pd.to_numeric(masterfile.loc[eachrow,:].apply(lambda x : x.split('/')[0]).values))
            variantannot.loc[eachrow,'S_RefCount'] = sum(pd.to_numeric(masterfile.loc[eachrow,:].apply(lambda x : x.split('/')[1]).values)) - sum(pd.to_numeric(masterfile.loc[eachrow,:].apply(lambda x : x.split('/')[0]).values))
        variantannot['VAF_total'] = [str(variantannot['S_AltCount'][i]) + '/' + str(variantannot['S_AltCount'][i] + variantannot['S_RefCount'][i]) for i in range(len(variantannot))]
    else:
        variantannot['S_AltCount'] = filloutfile["S_AltCount"]
        variantannot['S_RefCount'] = filloutfile["S_RefCount"]
        variantannot['VAF_total'] = filloutfile["S_AltCount"].astype(str) + '/' + (filloutfile["S_AltCount"].values + filloutfile["S_RefCount"].values).astype(str)
    
    # Iterate through the variants to recalculate the bulk proportions
    for eachrow in variantannot.index.values:
        variantannot.loc[eachrow,'bulkprop'] = float(variantannot.loc[eachrow,'S_AltCount']/(variantannot.loc[eachrow,'S_AltCount'] + variantannot.loc[eachrow,'S_RefCount']))
    
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
    if genome == "GRCh38" or genome == "GRCh37":
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
    
    # prepare the coverage and copy number information for each cell and combine the matrix with the resulting matrix
    resultMT1 = pd.concat([variantannot, variantsfile],axis=1,sort=False)
    resultMT = pd.concat([resultMT1,masterfile],axis=1,sort=False)
    
    # Re-calculate mutant cells and numcells for each variant
    for varid in resultMT.index.values[5:len(resultMT.index)]:
        if resultMT.loc[varid,'somaticstatus'] == 'germline':
            resultMT.loc[varid,'mutant'] = sum((resultMT.loc[varid,masterfile.columns[masterfile.columns.isin(resultMTcoverage.columns)]].apply(lambda x : int(x.split('/')[1])) - resultMT.loc[varid,masterfile.columns[masterfile.columns.isin(resultMTcoverage.columns)]].apply(lambda x : int(x.split('/')[0])))/resultMT.loc[varid,masterfile.columns[masterfile.columns.isin(resultMTcoverage.columns)]].apply(lambda x : int(x.split('/')[1])) > 0)
        else:
            resultMT.loc[varid,'mutant'] = sum(resultMT.loc[varid,masterfile.columns[masterfile.columns.isin(resultMTcoverage.columns)]].apply(lambda x : int(x.split('/')[0]))/resultMT.loc[varid,masterfile.columns[masterfile.columns.isin(resultMTcoverage.columns)]].apply(lambda x : int(x.split('/')[1])) > 0)
        resultMT.loc[varid,'numcells'] = resultMT.loc[varid,'mutant'] + resultMT.loc[varid,'wt']
    
    # Saving the final masterfile
    resultMT.to_csv(resultsdir + "/" + libraryid + '_master.tsv',sep = '\t')
    
    # Process the master file to generate binary and vaf matrix of filtered variants
    if genome == "GRCh38" or genome == "GRCh37":
        # Iterate through the germline variants to flip the vaf for the individual cells
        for eachrow in vaffile.index.values:
            if eachrow.split(':')[1] in pos: # it's a germline variant position
                vaffile.loc[eachrow] = 1 - vaffile.loc[eachrow]
    # Saving the final vaf file
    vaffile.to_csv(resultsdir + "/" + libraryid + '_vaf.tsv',sep = '\t',na_rep='NA')
    
    # Make a binary file from the vaf file
    for eachrow in vaffile.index.values:
        vaffile.loc[eachrow][vaffile.loc[eachrow] > 0] = 1
    vaffile.to_csv(resultsdir + "/" + libraryid + '_binary.tsv',sep = '\t',na_rep='NA')


if __name__ == "__main__":
    # Parse necessary arguments
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("-d", "--datadir",type=str, help="(REQUIRED) Directory for BAM files", required=True)
    parser.add_argument("-l", "--libraryid",type=str, help="(REQUIRED) Library ID", required=True)
    parser.add_argument("-w", "--workingdir", type=str, help="(REQUIRED) Working directory", required=True)
    parser.add_argument("-re", "--resultsdir", type=str, help="(REQUIRED) Directory for results", required=True)
    parser.add_argument("-q","--mapq",type=int,help="Minimum mapping quality, default = 20", default = 20)
    parser.add_argument("-Q","--baseq",type=int,help="Minimum base quality, default = 20", default = 20)
    parser.add_argument("-s","--strand",type=int,help="Minimum number of reads mapping to forward and reverse strand to call mutation, default=2",default = 2)
    parser.add_argument("-p", "--patternlist",type=str, help="File containing a list of filenames to process at a time", default = "")
    parser.add_argument("-t","--threshold",type=int,help="The critical threshold for calling a cell wild-type, default=0.1", default = 0.1)
    parser.add_argument("-vc", "--vepcache", type=str, help="Directory for vep cache", default="$HOME/.vep")
    parser.add_argument("-g", "--genome",type=str, help="Genome version",default = "GRCh37") 
    parser.add_argument("-r", "--reffile",type=str, help="Reference fasta file", default="")
    
    # read in arguments    
    args = parser.parse_args()
    datadir = args.datadir
    reffile = args.reffile
    minmapq = args.mapq
    minbq = args.baseq
    minstrand = args.strand
    threshold = args.threshold
    libraryid = args.libraryid
    patternlist = args.patternlist
    workingdir = args.workingdir
    vepcache = args.vepcache
    resultsdir = args.resultsdir
    genome = args.genome

    # Set the parameters for the genome build
    if genome == 'GRCh37':
        if reffile == "":
            reffile = workingdir + '/reference/b37/b37_MT.fa'
        # mtchrom = 'MT'
        ncbibuild = 'GRCh37'
        species = "homo_sapiens"
    elif genome == "GRCm38" or genome == "mm10":
        if reffile == "":
            reffile = workingdir + "/reference/mm10/mm10_MT.fa"
        # mtchrom = 'chrM'
        # mtchrom = "MT"
        ncbibuild = 'GRCm38'
        species = "mus_musculus"
    elif genome == 'GRCh38':
        if reffile == "":
            reffile = workingdir + '/reference/GRCh38/genome_MT.fa'
        # mtchrom = 'MT'
        ncbibuild = 'GRCh38'
        species = "homo_sapiens"
    else:
        raise Exception("The genome build you entered is not supported. Supported genomes are GRCh37, GRCh38, GRCm38, and mm10.")

    # Noting all the parameters
    print("Miminum mapping quality of " + str(minmapq))
    print("Miminum base quality of " + str(minbq))
    print("Miminum number of reads mapping to forward and reverse strand to call mutation of " + str(minstrand))

    # Filtering of cells
    merging_bams(datadir,libraryid,resultsdir)
    preproccess_bams(datadir,reffile,workingdir,vepcache,resultsdir,genome,species,ncbibuild)
    variant_calling(datadir,libraryid,reffile,genome,minmapq,minbq,minstrand,workingdir,vepcache,resultsdir,species,ncbibuild)
    variant_processing(datadir,libraryid,reffile,patternlist,resultsdir)
    if genome == "GRCh38" or genome == "GRCh37":
        runhaplogrep(datadir,libraryid,reffile,workingdir,resultsdir)
    processfillout(libraryid,threshold,resultsdir,genome)
    genmaster(libraryid,reffile,resultsdir,genome)
    
