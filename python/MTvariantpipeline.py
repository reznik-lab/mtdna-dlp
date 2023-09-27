# imports
import os
import sys
import numpy as np
import pandas as pd
import argparse
import pysam


def variant_calling(datadir, tumorbam, normalbam, normaldir, vcfdir, outdir, workingdir, vepcache, fasta, minmapq, minbq, minstrand, genome, mtchrom, mincounts):
    '''
    Variant calling
    '''
    print("Starting variant calling...")

    # set parameters for genome build
    if genome == 'GRCh37':
        if fasta == "":
            fasta = workingdir + '/reference/b37/b37_MT.fa'
        mtchrom = 'MT'
        ncbibuild = 'GRCh37'
        maf2maf_fasta = fasta
        bcfploidy_genome = 'GRCh37'
    elif genome == "GRCm38" or genome == "mm10":
        if fasta == "":
            fasta = workingdir + "/reference/mm10/mm10_MT.fa"
        # mtchrom = 'chrM'
        mtchrom = 'MT'
        ncbibuild = 'mm10'
        maf2maf_fasta = fasta
        bcfploidy_genome = 'mm10'
    elif genome == 'GRCh38':
        if fasta == "":
            fasta = workingdir + '/reference/GRCh38/genome_MT.fa'
        mtchrom = 'MT'
        ncbibuild = 'GRCh38'
        maf2maf_fasta = fasta
        bcfploidy_genome = 'GRCh38'
    else:
        print('The genome build you entered is not supported.')
        sys.exit(0)

    # create lists for column names
    mafnamedict = {4:['t_ref_count','t_alt_count'], 6:['t_ref_fwd','t_alt_fwd'], 
        7:['t_ref_rev','t_alt_rev'], 8:['n_ref_count','n_alt_count'], 
        10:['n_ref_fwd','n_alt_fwd'], 11:['n_ref_rev','n_alt_rev']}
    retaincols = ','.join(['Tumor_Sample_Barcode','Matched_Norm_Sample_Barcode',
        't_ref_count','t_alt_count','t_ref_fwd','t_alt_fwd','t_ref_rev','t_alt_rev',
        'n_ref_count','n_alt_count','n_ref_fwd','n_alt_fwd','n_ref_rev','n_alt_rev'])

    # Create a dataframe that stores the depth of MT coverage for each sample
    mtcounts = pd.DataFrame( columns = ['MTCounts','MTCountsNormal'] )

    # Try to get the readcounts
    try:
        mt = pysam.view('-c',datadir + "/" + tumorbam,'-q 10','-F 1536',mtchrom)
        mtcounts.at[tumorbam,'MTCounts'] = int(mt)
    except:
        raise Exception('Error in getting read counts for tumor bam ' + tumorbam)
    if normalbam != "":
        try:
            mtnorm = pysam.view('-c',normaldir + "/" + normalbam,'-q 10','-F 1536',mtchrom)
            mtcounts.at[tumorbam,'MTCountsNormal'] = mtnorm
        except:
            raise Exception('Error in getting read counts for normal bam ' + normalbam)

    # check to make sure number of counts is big enough
    if int(mt) < mincounts:
        raise Exception('Not enough MT reads')

    # countcall
    if normalbam != "":
        print('We have a matched normal bam file for ' + tumorbam)
        if bcfploidy_genome == 'mm10':
            countcall = f"samtools mpileup --region {mtchrom} --count-orphans --no-BAQ --min-MQ {minmapq} --min-BQ {minbq} " \
                + "--ignore-RG --excl-flags UNMAP,SECONDARY,QCFAIL,DUP --BCF --output-tags DP,AD,ADF,ADR --gap-frac 0.005 " \
                + f"--tandem-qual 80 -L 1000000 -d 1000000 --open-prob 30 --fasta-ref {fasta} {datadir}/{tumorbam} {normaldir}/{normalbam} " \
                + f"| bcftools call --multiallelic-caller --ploidy-file {workingdir}/reference/chrM_ploidy --keep-alts " \
                + f"| bcftools norm --multiallelics -any --do-not-normalize | vt normalize -r {fasta} - 2>/dev/null " \
                + f"| bcftools query --format '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%DP\t%ADF\t%ADR]\n' > {vcfdir}/{tumorbam}_temp.maf"
            mafcall = f"perl {workingdir}/vcf2maf/maf2maf.pl --vep-data {vepcache}/ --species mus_musculus --input-maf {vcfdir}/{tumorbam}_temp2.maf " \
                + f"--output-maf {outdir}/{tumorbam}.maf --retain-cols {retaincols} --ncbi-build GRCm38 --ref-fasta {fasta}"
        else:
            countcall = f"samtools mpileup --region {mtchrom} --count-orphans --no-BAQ --min-MQ {minmapq} --min-BQ {minbq} " \
                + "--ignore-RG --excl-flags UNMAP,SECONDARY,QCFAIL,DUP --BCF --output-tags DP,AD,ADF,ADR --gap-frac 0.005 " \
                + f"--tandem-qual 80 -L 1000000 -d 1000000 --open-prob 30 --fasta-ref {fasta} {datadir}/{tumorbam} {normaldir}/{normalbam} " \
                + f"| bcftools call --multiallelic-caller --ploidy {bcfploidy_genome} --keep-alts | bcftools norm " \
                + f"--multiallelics -any --do-not-normalize | vt normalize -r {fasta} - 2>/dev/null | bcftools query " \
                + f"--format '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%DP\t%ADF\t%ADR]\n' > {vcfdir}/{tumorbam}_temp.maf"
            mafcall = f"perl {workingdir}/vcf2maf/maf2maf.pl --vep-data {vepcache}/ --input-maf {vcfdir}/{tumorbam}_temp2.maf " \
                + f"--output-maf {outdir}/{tumorbam}.maf --retain-cols {retaincols} --ncbi-build {ncbibuild} --ref-fasta {fasta}"
    else:
        print('We do not have a normal bam file for ' + tumorbam)
        if bcfploidy_genome == 'mm10':
            countcall = f"samtools mpileup --region {mtchrom} --count-orphans --no-BAQ --min-MQ {minmapq} --min-BQ {minbq} " \
                + "--ignore-RG --excl-flags UNMAP,SECONDARY,QCFAIL,DUP --BCF --output-tags DP,AD,ADF,ADR --gap-frac 0.005 " \
                + f"--tandem-qual 80 -L 1000000000 -d 1000000000 --open-prob 30 --fasta-ref {fasta} {datadir}/{tumorbam} " \
                + f"| bcftools call --multiallelic-caller --ploidy-file {workingdir}/reference/chrM_ploidy --keep-alts " \
                + f"| bcftools norm --multiallelics -any --do-not-normalize | vt normalize -r {fasta} - 2>/dev/null " \
                + f"| bcftools query --format '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%DP\t%ADF\t%ADR]\n' > {vcfdir}/{tumorbam}_temp.maf"
            mafcall = f"perl {workingdir}/vcf2maf/maf2maf.pl --vep-data {vepcache}/ --species mus_musculus --input-maf {vcfdir}/{tumorbam}_temp2.maf " \
                + f"--output-maf {outdir}/{tumorbam}.maf --retain-cols {retaincols} --ncbi-build GRCm38 --ref-fasta {fasta}"
        else:
            countcall = f"samtools mpileup --region {mtchrom} --count-orphans --no-BAQ --min-MQ {minmapq} --min-BQ {minbq} " \
                + "--ignore-RG --excl-flags UNMAP,SECONDARY,QCFAIL,DUP --BCF --output-tags DP,AD,ADF,ADR --gap-frac 0.005 " \
                + f"--tandem-qual 80 -L 1000000000 -d 1000000000 --open-prob 30 --fasta-ref {fasta} {datadir}/{tumorbam} " \
                + f"| bcftools call --multiallelic-caller --ploidy {bcfploidy_genome} --keep-alts | bcftools norm " \
                + f"--multiallelics -any --do-not-normalize | vt normalize -r {fasta} - 2>/dev/null | bcftools query " \
                + f"--format '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%DP\t%ADF\t%ADR]\n' > {vcfdir}/{tumorbam}_temp.maf"
            mafcall = f"perl {workingdir}/vcf2maf/maf2maf.pl --vep-data {vepcache}/ --input-maf {vcfdir}/{tumorbam}_temp2.maf " \
                + f"--output-maf {outdir}/{tumorbam}.maf --retain-cols {retaincols} --ncbi-build {ncbibuild} --ref-fasta {fasta}"
    print("COUNTCALL: ", countcall)
    os.system(countcall)
    print("DONE WITH COUNTCALL")

    # Process prelim MAF file and remove any rows that have 0 non-ref reads.
    tempmaf = pd.read_csv(vcfdir + tumorbam + "_temp.maf",header = None,sep = '\t')
    tempmaf = tempmaf[ tempmaf[3] != '.' ]
    if normalbam != "":
        tempmaf.columns = ['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',4,'t_depth',6,7,8,'n_depth',10,11]
        tempmaf['Tumor_Sample_Barcode'] = tumorbam
        tempmaf['Matched_Norm_Sample_Barcode'] = normalbam
        cols2use = [4,6,7,8,10,11]
    else:
        tempmaf.columns = ['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',4,'t_depth',6,7]
        tempmaf['Tumor_Sample_Barcode'] = tumorbam
        tempmaf['Matched_Norm_Sample_Barcode'] = ''
        for nacol in ['n_depth','n_ref_count', 'n_alt_count', 'n_ref_fwd', 'n_alt_fwd', 'n_ref_rev', 'n_alt_rev']:
            tempmaf[nacol] = np.nan
        cols2use = [4,6,7]
    for col in cols2use:
        cname1 = mafnamedict[col][0]
        cname2 = mafnamedict[col][1]
        tempmaf[cname1] = [item.split(',')[0] for item in tempmaf[col]]
        tempmaf[cname2] = [item.split(',')[1] for item in tempmaf[col]]

    # Drop unnecessary columns
    tempmaf = tempmaf.drop(cols2use, 1)

    # Make sure that any variants we keep have at least minimum reads, with at least one alternate read in both directions
    tempmaf = tempmaf[ tempmaf['t_alt_fwd'].map(int) >= minstrand ]
    tempmaf = tempmaf[ tempmaf['t_alt_rev'].map(int) >= minstrand ]

    # remove position 3106
    tempmaf = tempmaf[~tempmaf['Start_Position'].isin(list(range(3106, 3107)))]

    # Write out to a second temporary MAF file, and then call maf2maf
    tempmaf.to_csv(vcfdir + tumorbam + "_temp2.maf",index = None,sep = '\t')
    if tempmaf.shape[0] == 0:
        raise Exception("No MT variants for " + tumorbam)
    print("MAFCALL: ", mafcall)
    os.system(mafcall)
    print("DONE WITH MAFCALL")

    # # write out final counts file
    # mtcounts.to_csv(outdir + "/" + args.bamfiles.split('/')[-1] + 'Counts.txt')


def final_processing(outdir, workingdir, tumorbam, normalbam):
    '''
    Annotate SNPs with tRNA and mitimpact data.
    Assume that all frameshifts/nonsense are potentially pathogenic.
    Add information on supporting forward and reverse reads.
    Assign conditions for pathogenicity.
    Modify the gene names to include rRNA and tRNA.
    Write out final maf file.
    '''
    print("Starting final processing...")

    # Read in the annotations
    trna = pd.read_csv(workingdir + '/reference/MitoTIP_August2017.txt',header = 0,sep = '\t')
    mitimpact = pd.read_csv(workingdir + '/reference/MitImpact_db_2.7.txt',header = 0,sep = '\t',decimal = ',') # note that the decimal point here is indicated as a comma, mitimpact is funnypontrna

    # Make the indices searchable for annotation for tRNA data
    trna.index = [trna.at[item,'rCRS base'] + str(trna.at[item,'Position']) + 
        trna.at[item,'Change'] if trna.at[item,'Change'] != 'del' else trna.at[item,'rCRS base'] + 
        str(trna.at[item,'Position']) + trna.at[item,'Change'] + trna.at[item,'rCRS base'] for item in trna.index]

    # For Mitimpact, consider only snps, and make data searchable
    mitimpact = mitimpact[ mitimpact['Start'] == mitimpact['End'] ]
    mitimpact.index = [mitimpact.at[item,'Ref'] + str(mitimpact.at[item,'Start']) + 
        mitimpact.at[item,'Alt'] for item in mitimpact.index]

    # Make sure there are no duplicate indices for the annotation
    trna = trna[~trna.index.duplicated(keep='first')]
    mitimpact = mitimpact[~mitimpact.index.duplicated(keep='first')]

    # Indicate the columns to keep for trna and mitimpact when annotating
    trna_cols = ['Predictive score']
    mitimpact_cols = ['APOGEE_boost_mean_prob', 'Mitomap_Dec2016_Status', 'Mitomap_Dec2016_Disease']

    # These are the "bad" mutations we automatically call pathogenic
    badmuts = ['Nonsense_Mutation','Nonstop_Mutation','Frame_Shift_Del','Frame_Shift_Ins']

    # Read in the gene positions
    genepos = pd.read_csv(workingdir + '/reference/GenePositions_imported.csv',header = 0,index_col = 0)
    
    # Read in the MAF file
    maf = pd.read_csv(outdir + "/" + tumorbam + '.maf',header = 0,sep = '\t',comment = '#')

    # Make a short name for each variant. 
    maf['ShortVariantID'] = maf['Reference_Allele'] + maf['Start_Position'].map(str) + maf['Tumor_Seq_Allele2']

    # Add the annotations
    for col in trna_cols + mitimpact_cols:
        colorder = maf.columns.tolist()
        maf = pd.concat( [maf,pd.DataFrame( columns = [col] )] , sort = False )
        
        # Keep column corder
        maf = maf[colorder + [col]]
        
    if len(np.intersect1d(trna.index.values, maf[ 'ShortVariantID' ])) > 0:
        maf[trna_cols] = trna.reindex(columns = trna_cols, index = maf['ShortVariantID'].values)
    else:
        maf[trna_cols] = np.NaN
    if sum(maf['ShortVariantID'].isin(mitimpact.index.values)) > 0:
        maf[mitimpact_cols] = mitimpact.reindex(columns = mitimpact_cols, index = maf['ShortVariantID'].values)
    else:
        maf[mitimpact_cols] = 0
    maf['TumorVAF'] = maf['t_alt_count']/maf['t_depth']
    if normalbam != "":
        maf['NormalVAF'] = maf['n_alt_count']/maf['n_depth']
    else:
        maf['NormalVAF'] = 'NA'

    # create columns for pathogenicity
    maf['Pathogenic_Reason'] = ''
    maf['Pathogenic_mtDNA_Variant'] = False
    
    # Anything with APOGEE score greater than 0.9
    maf.loc[maf['APOGEE_boost_mean_prob'] > 0.9,'Pathogenic_mtDNA_Variant'] = True
    maf.loc[maf['APOGEE_boost_mean_prob'] > 0.9,'Pathogenic_Reason'] = 'APOGEE'
    
    # Anything confirmed in MITOMAP
    maf.loc[ maf['Mitomap_Dec2016_Status'].isin(['Confirmed','Cfrm']), 'Pathogenic_mtDNA_Variant'] = True
    maf.loc[ maf['Mitomap_Dec2016_Status'].isin(['Confirmed','Cfrm']), 'Pathogenic_Reason'] = 'MITOMAP'
    
    # Anything with tRNA score greater than 18
    maf.loc[maf['Predictive score'] > 16.2,'Pathogenic_mtDNA_Variant'] = True
    maf.loc[maf['Predictive score'] > 16.2,'Pathogenic_Reason'] = 'tRNA Predictive Score'
    
    # Any frameshift/nonsense variants
    maf.loc[maf['Variant_Classification'].isin(badmuts),'Pathogenic_mtDNA_Variant'] = True
    maf.loc[maf['Variant_Classification'].isin(badmuts),'Pathogenic_Reason'] = 'Frameshift/Nonsense'

    # modify the gene names to include rRNA and tRNA and change the control region symbols
    maf['Hugo_Symbol'] = genepos.loc[maf['Start_Position'],'Gene'].reset_index(drop = True)
    maf.loc[(maf['Start_Position'].map(int) <= 576) | (maf['Start_Position'].map(int) >= 16024),'Hugo_Symbol'] = 'ControlRegion'

    # write out final files
    maf.to_csv(outdir + "/" + tumorbam + '.maf',index = None,sep = '\t')


if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("-d", "--datadir",type=str, help="directory for input BAM files", required=True)
    parser.add_argument("-v", "--vcfdir", type=str, help="directory for intermediate VCF files", required=True)
    parser.add_argument("-o","--outdir",type=str,help="directory for output MAF files", required=True)
    parser.add_argument("-w", "--workingdir", type=str, help="Working directory", required=True)
    parser.add_argument("-vc", "--vepcache", type=str, help="Directory for vep cache", required=True)
    parser.add_argument("-b", "--tumorbam", type=str, help="name of tumor bam, should not include full path", required=True)
    parser.add_argument("-n", "--normalbam", type=str, help="name of normal bam, should not include full path", default="")
    parser.add_argument("-nd", "--normaldir", type=str, help="directory that contains matched normal file",default="")
    parser.add_argument("-g","--genome",type=str,help="Genome build to use, default = GRCh37", default = "GRCh37")
    parser.add_argument("-f", "--fasta", type=str, help="path to fasta file", default="")
    parser.add_argument("-q","--mapq",type=int,help="minimum mapping quality, default = 10",default = 10)
    parser.add_argument("-Q","--baseq",type=int,help="minimum base quality, default = 10",default = 10)
    parser.add_argument("-s","--strand",type=int,help="Minimum number of reads mapping to forward and reverse strand to call mutation, default=2",default = 2)
    parser.add_argument("-m", "--mtchrom",type=str, help="Chromosome type", default="MT")
    parser.add_argument("-c", "--mincounts",type=int, help="Minimum number of read counts, default = 100", default=100)
    
    # read arguments
    args = parser.parse_args()
    datadir = args.datadir
    vcfdir = args.vcfdir
    outdir = args.outdir
    workingdir = args.workingdir
    vepcache = args.vepcache
    tumorbam = args.tumorbam
    normalbam = args.normalbam
    normaldir = args.normaldir
    genome = args.genome
    fasta = args.fasta
    minmapq = args.mapq
    minbq = args.baseq
    minstrand = args.strand
    mtchrom = args.mtchrom
    mincounts = args.mincounts
    
    # make output directories if they don't exist
    if not os.path.exists(vcfdir):
        os.makedirs(vcfdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    variant_calling(datadir, tumorbam, normalbam, normaldir, vcfdir, outdir, workingdir, vepcache, fasta, minmapq, minbq, minstrand, genome, mtchrom, mincounts)
    final_processing(outdir, workingdir, tumorbam, normalbam)

    print("DONE WITH MT VARIANT PIPELINE")