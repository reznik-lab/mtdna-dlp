import os
import pysam
import argparse
import warnings
import pandas as pd

class SplitBam(object):
    """
    calculate reads per bin from the input bam file
    """
    def __init__(self, bam_file, odir, chromosomes, barcode_csv):
        self.bam = self.__get_bam_reader(bam_file)
        self.odir = odir
        if chromosomes:
            self.chromosomes = chromosomes
        else:
            self.chromosomes = self.__get_chr_names()
        self.chr_lengths = self.__get_chr_lengths()
        self.barcodes = self.__get_barcodes(barcode_csv)

    def __get_barcodes(self, barcode_file):
        df = pd.read_csv(barcode_file, header=None, sep="\t")
        return set(df.iloc[:, 0].values.tolist())

    def __enter__(self):
        return self

    def __exit__(self):
        # clean up output if there are any exceptions
        self.bam.close()

    def __get_chr_lengths(self):
        """ returns dict with chromosome names and lengths
        :returns dictionary with chromosome name (str) and lengths(int)
        :rtype dictionary
        """
        names = self.bam.references
        lengths = self.bam.lengths
        return {name: length for name, length in zip(names, lengths)}

    def __get_bam_reader(self, bam_file):
        """returns pysam bam object
        :returns pysam bam object
        """
        return pysam.AlignmentFile(bam_file, 'rb')

    def __get_chr_names(self):
        """extracts chromosome names from the bam file  CCTACACTCGATCCGG-1_1.fq.gz
        :returns list of chromosome names
        :rtype list
        """
        return self.bam.references

    def __get_file_for_tag(self, tag_to_file_map, tag_value):
        # TODO make a class for this functionality
        """
        Get a file handle to a file in the specified output dir,
        for writing alignments with the specified tag_value.

        :param tag_to_file_map: map from tag value to file handler
        :param tag_value: tag value
        :param odir: directory where the file should be written
        :param alignment_file: handle to the alignment file
        :return: File handle for the specified tag value.
        """
        if tag_value not in tag_to_file_map.keys():
            tag_file_name = os.path.join(self.odir, tag_value + ".bam")
            file_handle = pysam.AlignmentFile(tag_file_name, "wb", template=self.bam)
            tag_to_file_map[tag_value] = file_handle
        return tag_to_file_map[tag_value]

    def main(self):
        read_count = 0
        tag_to_file_map = {}
        for read in self.bam.fetch():
            read_count += 1
            tags = dict(read.tags)
            if 'CB' in tags.keys():
                if tags['CB'] in self.barcodes:  # filter by barcodes
                    outfile = self.__get_file_for_tag(tag_to_file_map, tags['CB'])
                    outfile.write(read)
            else:
                outfile = self.__get_file_for_tag(tag_to_file_map, "undetermined")
                outfile.write(read)
        print("closing files: ")
        for file in tag_to_file_map.values():
            try:
                file.close()
            except:
                warnings('error closing file ' + file)
        print("num reads = "+str(read_count))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help='specify the path to the input bam file')
    parser.add_argument('odir', help='specify path to the output dir')
    parser.add_argument('--chromosomes', nargs='*', default=list(map(str, range(1, 23))) + ['chrM'], help='specify target chromosomes')
    parser.add_argument('--barcode_csv', help='threshold for the mapping quality, reads with quality lower than threshold will be ignored')
    args = parser.parse_args()

    with SplitBam(args.bam, args.odir, args.chromosomes, args.barcode_csv) as splitbam:
        splitbam.main()


# Example:
# python3.7 split_bam.py [post sorted MT subset bamfile] [splits dir] 
#     --barcode_csv [barcode file]