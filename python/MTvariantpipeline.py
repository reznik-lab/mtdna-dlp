# Function to call vcf2maf for mitochondrial variants

import os, sys, numpy as np, pandas as pd, argparse, pysam
#pdb
#scipy as sp
#time

# Parse necessary arguments
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("-d", "--datadir",type=str, help="directory for BAM files")
parser.add_argument("-v", "--vcfdir", type=str, help="directory for intermediate VCF files")
parser.add_argument("-o","--outdir",type=str,help="directory for MAF files")
parser.add_argument("-b","--bamfiles",type=str,help="path to tab-delimited, no header file of tumor/normal bams. BAM names should not include the full path, just the name of the files in datadir to call variants on. If not a .csv file, just include the single name of the tumor bam file you would like to call.")
parser.add_argument("-q","--mapq",type=int,help="minimum mapping quality, default = 10",default = 10)
parser.add_argument("-Q","--baseq",type=int,help="minimum base quality, default = 10",default = 10)
parser.add_argument("-n","--normalflag",type=str,help="matched normal available, default = False",default = "False")
parser.add_argument("-g","--genome",type=str,help="Genome build to use, default = GRCh37, for gdc use GRCh38",default = "GRCh37")
parser.add_argument("-s","--strand",type=int,help="Minimum number of reads mapping to forward and reverse strand to call mutation, default=2",default = 2)
parser.add_argument("-h","--help",action='help', default=argparse.SUPPRESS,
                    help='A simple variant calling and annotation pipeline for mitochondrial DNA variants. Accepts both individual BAM files and paired tumor/normal BAM files. Because of the hairiness of multiallelic calling, we keep all positions in the VCF and then filter to remove any positions without at least 10 reads supporting the putative variant. FIX TO DO THIS FOR BOTH NORMAL AND TUMOR SAMPLES! To send a call to bsub, try bsub -R "rusage[mem=16]" -M 32 -We 120 -W 4800 -e $HOME/work/mtimpact/scratch/ -o /home/reznik/work/mtimpact/scratch/ python MTvariantpipeline.py ... with the suitable options specified. Note that we use the CMO version of b37 (which seems to use rCRS), but is named HG19.')
parser.add_argument("-hd", "--homedir", type=str, help="Home directory for a lot of stuff")
parser.add_argument("-vd", "--vepdir", type=str, help="Directory for vep")
parser.add_argument("-vc", "--vepcache", type=str, help="Directory for vep cache")

args = parser.parse_args()

# Read in the arguments
datadir = args.datadir
vcfdir = args.vcfdir
outdir = args.outdir
genome = args.genome
minstrand = args.strand
normalflag = args.normalflag
homedir = args.homedir
vepdir = args.vepdir
vepcache = args.vepcache

# Set key parameters
minmapq = args.mapq
minbq = args.baseq

# Make sure the output directories are created
if not os.path.exists(vcfdir):
    os.makedirs(vcfdir)
if not os.path.exists(outdir):
    os.makedirs(outdir)
    
# Set the parameters for the genome build
if genome == 'GRCh37':
    fasta = homedir + 'reference/b37/b37_MT.fa'
    mtchrom = 'MT'
    ncbibuild = 'GRCh37'
    maf2maf_fasta = fasta
    bcfploidy_genome = 'GRCh37'
elif genome == 'GRCh38':
    fasta = homedir + 'reference/GRCh38/genome_MT.fa'
    mtchrom = 'MT'
    ncbibuild = 'GRCh38'
    maf2maf_fasta = fasta
    bcfploidy_genome = 'GRCh38'
# elif genome == 'mm10':
#     fasta = homedir + 'reference/mm10_genome.fa'
#     mtchrom = 'chrM'
#     ncbibuild = 'mm10'
#     maf2maf_fasta = fasta
#     bcfploidy_genome = 'mm10'
# elif genome == 'hg19':
#     # this is the bona fide hg19, with chrM of length 16571
#     fasta = '/ifs/depot/pi/resources/genomes/hg19/fasta/hg19.fasta' #where to get this file?
#     mtchrom = 'chrM'
#     ncbibuild = 'GRCh38'
    
    # we need the mapping from yoruba to rcrs
    rcrsmapping = pd.read_csv(homedir + 'reference/Yoruba2rCRS.txt',header = 0,index_col = 0)
    
    # we also need to make sure maf2maf uses the correct hg38 fasta
    maf2maf_fasta = '/ifs/depot/pi/resources/genomes/GRCh38/fasta/GRCh38.d1.vd1.fa' #where to get this file?
    
    # ploidy for bcftools same as GRCh37
    bcfploidy_genome = 'GRCh37'
else:
    print('The genome build you entered is not supported.')
    sys.exit(0)

# List of what to map temporary MAF columns to     
mafnamedict = {4:['t_ref_count','t_alt_count'], 6:['t_ref_fwd','t_alt_fwd'], 7:['t_ref_rev','t_alt_rev'], 8:['n_ref_count','n_alt_count'], 10:['n_ref_fwd','n_alt_fwd'], 11:['n_ref_rev','n_alt_rev']}

retaincols = ','.join( ['Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 't_ref_count','t_alt_count','t_ref_fwd','t_alt_fwd','t_ref_rev','t_alt_rev', 'n_ref_count','n_alt_count','n_ref_fwd','n_alt_fwd','n_ref_rev','n_alt_rev'] )
    
# Read in the annotations
trna = pd.read_csv(homedir + 'reference/MitoTIP_August2017.txt',header = 0,sep = '\t')
mitimpact = pd.read_csv(homedir + 'reference/MitImpact_db_2.7.txt',header = 0,sep = '\t',decimal = ',') # note that the decimal point here is indicated as a comma, mitimpact is funnypontrna

# Make the indices searchable for annotation for tRNA data
trna.index = [trna.at[item,'rCRS base'] + str(trna.at[item,'Position']) + trna.at[item,'Change'] if trna.at[item,'Change'] != 'del' else trna.at[item,'rCRS base'] + str(trna.at[item,'Position']) + trna.at[item,'Change'] + trna.at[item,'rCRS base'] for item in trna.index]

# For Mitimpact, consider only snps, and make data searchable
mitimpact = mitimpact[ mitimpact['Start'] == mitimpact['End'] ]
mitimpact.index = [mitimpact.at[item,'Ref'] + str(mitimpact.at[item,'Start']) + mitimpact.at[item,'Alt'] for item in mitimpact.index]

# Make sure there are no duplicate indices for the annotation
trna = trna[~trna.index.duplicated(keep='first')]
mitimpact = mitimpact[~mitimpact.index.duplicated(keep='first')]

# Indicate the columns to keep for trna and mitimpact when annotating
trna_cols = ['Predictive score']
mitimpact_cols = ['APOGEE_boost_mean_prob', 'Mitomap_Dec2016_Status', 'Mitomap_Dec2016_Disease']

# These are the "bad" mutations we automatically call pathogenic
badmuts = ['Nonsense_Mutation','Nonstop_Mutation','Frame_Shift_Del','Frame_Shift_Ins']

# Read in the gene positions
genepos = pd.read_csv(homedir + 'reference/GenePositions_imported.csv',header = 0,index_col = 0)

# Read in the paired BAM files
if args.bamfiles.endswith('.csv') or args.bamfiles.endswith('.txt'):
    bamfiles = pd.read_csv(args.bamfiles,header = None,sep = '\t')
    fs = bamfiles.iloc[:,0]
else:
    # We just have a single bamfile, call from that
    bamfiles = pd.DataFrame(columns = [0,1])
    bamfiles.at[0,0] = args.bamfiles

# Also create a dataframe that stores the depth of MT coverage for each sample
mtcounts = pd.DataFrame( columns = ['MTCounts','MTCountsNormal'] )

for ii in range(bamfiles.shape[0]):
    
    f = bamfiles.at[ii,0].strip() # remove leading and trailing whitespace
    
    print('Working on ' + f + '...')
    
    # Check if we have a normal file
    if pd.isnull(bamfiles.at[ii,1]):
        normalflag = False
    else:
        normalbam = bamfiles.at[ii,1].strip() # remove leading and trailing whitespace
        normalflag = True
    
    # Try to get the readcounts for the file
    try:
        mt = pysam.view('-c',datadir + f,'-q 10','-F 1536',mtchrom)
        mtcounts.at[f,'MTCounts'] = int(mt)
    except:
        print('Error in getting read counts for ' + f + ', moving on...')
        continue
    
    if normalflag:
        try:
            mtnorm = pysam.view('-c',datadir + normalbam,'-q 10','-F 1536',mtchrom)
            mtcounts.at[f,'MTCountsNormal'] = mtnorm
        except:
            print('Error in getting read counts for ' + normalbam + ',skipping this part.')
    
    # If the number of counts is small, just move on
    if int(mt) < 100:
        print('Skipping ' + f + '...no MT reads to call variants with.')
        continue
    
    if normalflag:
        # We have a matched normal bam
        print('We have a matched normal bam file for ' + f)
        
        if bcfploidy_genome == 'mm10':
            countcall = ' '.join(["samtools mpileup --region", mtchrom, "--count-orphans --no-BAQ --min-MQ",str(minmapq), "--min-BQ", str(minbq), "--ignore-RG --excl-flags UNMAP,SECONDARY,QCFAIL,DUP --BCF --output-tags DP,AD,ADF,ADR --gap-frac 0.005 --tandem-qual 80 -L 1000000 -d 1000000 --open-prob 30 --fasta-ref", fasta, datadir  + f , datadir + normalbam + "| bcftools call --multiallelic-caller --ploidy-file " + homedir + "chrM_ploidy --keep-alts | bcftools norm --multiallelics -any --do-not-normalize | " + homedir + "vt normalize -r " + fasta + " - 2>/dev/null | bcftools query --format '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%DP\t%ADF\t%ADR]\n'", ">", vcfdir + f + "_temp.maf"])
            mafcall = ' '.join( ["perl " + homedir + "vcf2maf/maf2maf.pl --vep-data " + vepcache + " --vep-path " + vepdir + " --species mus_musculus --input-maf", vcfdir + f + "_temp2.maf","--output-maf", outdir + f + ".maf","--retain-cols", retaincols, "--ncbi-build GRCm38 --ref-fasta",fasta, "--vep-data /opt/common/CentOS_6-dev/vep/cache"] )

        else:
            countcall = ' '.join(["samtools mpileup --region", mtchrom, "--count-orphans --no-BAQ --min-MQ",str(minmapq), "--min-BQ", str(minbq), "--ignore-RG --excl-flags UNMAP,SECONDARY,QCFAIL,DUP --BCF --output-tags DP,AD,ADF,ADR --gap-frac 0.005 --tandem-qual 80 -L 1000000 -d 1000000 --open-prob 30 --fasta-ref", fasta, datadir  + f, datadir + normalbam + "| bcftools call --multiallelic-caller --ploidy", bcfploidy_genome, "--keep-alts | bcftools norm --multiallelics -any --do-not-normalize | " + homedir + "vt normalize -r " + fasta + " - 2>/dev/null | bcftools query --format '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%DP\t%ADF\t%ADR]\n'", ">", vcfdir + f + "_temp.maf"])
            mafcall = ' '.join( ["perl " + homedir + "vcf2maf/maf2maf.pl --vep-data " + vepcache + " --vep-path " + vepdir + " --input-maf", vcfdir + f + "_temp2.maf","--output-maf", outdir + f + ".maf","--retain-cols", retaincols, "--ncbi-build", ncbibuild, '--ref-fasta',maf2maf_fasta, "--vep-data /opt/common/CentOS_6-dev/vep/cache"])
        
    if not normalflag:
        # We don't have a normal bam
        print('We do not have a normal bam file for ' + f)
        
        if bcfploidy_genome == 'mm10':
            countcall = ' '.join(["samtools mpileup --region", mtchrom, "--count-orphans --no-BAQ --min-MQ",str(minmapq), "--min-BQ", str(minbq), "--ignore-RG --excl-flags UNMAP,SECONDARY,QCFAIL,DUP --BCF --output-tags DP,AD,ADF,ADR --gap-frac 0.005 --tandem-qual 80 -L 1000000000 -d 1000000000 --open-prob 30 --fasta-ref", fasta, datadir  + f + "| bcftools call --multiallelic-caller --ploidy-file " + homedir + "chrM_ploidy --keep-alts | bcftools norm --multiallelics -any --do-not-normalize | " + homedir + "vt normalize -r " + fasta + " - 2>/dev/null | bcftools query --format '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%DP\t%ADF\t%ADR]\n'", ">", vcfdir + f + "_temp.maf"])
            mafcall = ' '.join( ["perl " + homedir + "vcf2maf/maf2maf.pl --vep-data " + vepcache + " --vep-path " + vepdir + " --species mus_musculus --input-maf", vcfdir + f + "_temp2.maf","--output-maf", outdir + f + ".maf","--retain-cols", retaincols, "--ncbi-build GRCm38 --ref-fasta",fasta, "--vep-data /opt/common/CentOS_6-dev/vep/cache"])

        else:
            countcall = ' '.join(["samtools mpileup --region", mtchrom, "--count-orphans --no-BAQ --min-MQ",str(minmapq), "--min-BQ", str(minbq), "--ignore-RG --excl-flags UNMAP,SECONDARY,QCFAIL,DUP --BCF --output-tags DP,AD,ADF,ADR --gap-frac 0.005 --tandem-qual 80 -L 1000000000 -d 1000000000 --open-prob 30 --fasta-ref", fasta, datadir  + f + "| bcftools call --multiallelic-caller --ploidy", bcfploidy_genome, "--keep-alts | bcftools norm --multiallelics -any --do-not-normalize | bcftools query --format '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%DP\t%ADF\t%ADR]\n'", ">", vcfdir + f + "_temp.maf"])
            mafcall = ' '.join( ["perl " + homedir + "vcf2maf/maf2maf.pl --vep-data " + vepcache + " --vep-path " + vepdir + " --input-maf", vcfdir + f + "_temp2.maf","--output-maf", outdir + f + ".maf","--retain-cols", retaincols, "--ncbi-build", ncbibuild, '--ref-fasta',fasta])        

    print("COUNTCALL: ", countcall)
    os.system(countcall)
    print("DONE WITH COUNTCALL")
    
    # Read in the prelim MAF file, and remove any rows that have 0 non-ref reads.
    tempmaf = pd.read_csv(vcfdir + f + "_temp.maf",header = None,sep = '\t')
    tempmaf = tempmaf[ tempmaf[3] != '.' ]
        
    if normalflag:
        tempmaf.columns = ['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',4,'t_depth',6,7,8,'n_depth',10,11]
        
        tempmaf['Tumor_Sample_Barcode'] = f
        tempmaf['Matched_Norm_Sample_Barcode'] = normalbam
        
        cols2use = [4,6,7,8,10,11]

    if not normalflag:
        
        tempmaf.columns = ['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',4,'t_depth',6,7]
        
        tempmaf['Tumor_Sample_Barcode'] = f
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
    tempmaf = tempmaf.drop( cols2use, 1 )
    
    # Make sure that any variants we keep have at least minimum reads, with at least one alternate read in both directions
    tempmaf = tempmaf[ tempmaf['t_alt_fwd'].map(int) >= minstrand ]
    tempmaf = tempmaf[ tempmaf['t_alt_rev'].map(int) >= minstrand ]
    
    # If we are using hg19 with the 16571 long MT chromosome, change the mapping so that we can call variants correctly. 
    if genome == 'hg19':
        tempmaf['YorubaPosition'] = tempmaf['Start_Position'].copy()
        newstarts = rcrsmapping.ix[ 'Position:' + tempmaf['YorubaPosition'].map(int).map(str),0 ].reset_index(drop = True)

        tempmaf['Start_Position'] = newstarts.tolist()
        tempmaf = tempmaf[ tempmaf['Start_Position'].notnull() ]

        tempmaf['Start_Position'] = tempmaf['Start_Position'].map(int)
    
    # Write out to a second temporary MAF file, and then call maf2maf
    tempmaf.to_csv(vcfdir + f + "_temp2.maf",index = None,sep = '\t')
    if tempmaf.shape[0] == 0:
        print('Skipping ' + f + '...no MT variants.')
        continue
    
    print("MAFCALL: ", mafcall)
    os.system(mafcall)
    
    ####################################################################################
    ####################################################################################
    
    # Part 2: Variant annotation. The general workflow is to annotate SNPs with tRNA and mitimpact data, and to assume that all frameshifts/nonsense are potentially pathogenic. We also add information on supporting forward and reverse reads.
    
    # Read in the MAF file
    print("reading in MAF")
    maf = pd.read_csv(outdir + f + '.maf',header = 0,sep = '\t',comment = '#')
    
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
    if normalflag:
        maf['NormalVAF'] = maf['n_alt_count']/maf['n_depth']
    else:
        maf['NormalVAF'] = 'NA'

    ####################################################################################
    ####################################################################################
    #pdb.set_trace()
    # Part 3: Assign conditions for pathogenicity
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
    
    ####################################################################################
    ####################################################################################
    
    # Part 4: Modify the gene names to include rRNA and tRNA
    maf['Hugo_Symbol'] = genepos.loc[maf['Start_Position'],'Gene'].reset_index(drop = True)
    
    # Also change the control region symbols
    maf.loc[(maf['Start_Position'].map(int) <= 576) | (maf['Start_Position'].map(int) >= 16024),'Hugo_Symbol'] = 'ControlRegion'
    
    # Part 5: Write out file
    maf.to_csv(outdir + f + '.maf',index = None,sep = '\t')

# Write out the counts file as well
mtcounts.to_csv(outdir + args.bamfiles.split('/')[-1] + 'Counts.txt')

# Write out the options that we used
optsfile = open(outdir + args.bamfiles.split('/')[-1] + 'Options.txt','w')
optsfile.write('datadir: ' + datadir + '\n')
optsfile.write('vcfdir: ' + vcfdir + '\n')
optsfile.write('outdir: ' + outdir + '\n')
optsfile.write('genome: ' + genome + '\n')
optsfile.write('minstrand: ' + str(minstrand) + '\n')
optsfile.write('minmapq: ' + str(minmapq) + '\n')
optsfile.write('minbq:' + str(minbq) + '\n')
optsfile.close()

print("DONE WITH MT VARIANT PIPELINE")
