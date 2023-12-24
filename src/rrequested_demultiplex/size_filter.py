import sys
from Bio import SeqIO
import gzip
import os
from rrequested_pkg.file_system import FileSystem

class Size_Filter:
    def __init__(self, infile, minsize, maxsize):
        """
        Initializes the Size_Filter class with input file path, minimum size, and maximum size.

        Parameters:
        - infile (str): Path to the input file.
        - minsize (int): Minimum size for filtering sequences.
        - maxsize (int): Maximum size for filtering sequences.
        """
        self.infile = infile
        self.minsize = minsize
        self.maxsize = maxsize

    def load_data(self):
        """
        Load data from the input file based on the file format (fastq or fasta).
        Returns a dictionary or a list of sequences depending on the file format.

        Returns:
        - dict or list: A dictionary or a list containing read headers as keys and corresponding sequences.
        """
        print("Reading data from input file...", file=sys.stderr)
        if self.infile.endswith(".gz"):  # If file is gzipped, unzip it
            y = gzip.open(self.infile, "rt", encoding="latin-1")
            if self.infile.endswith((".fastq.gz", ".fasta.gz", ".fna.gz", ".fas.gz", ".fa.gz")):  # Read file as fastq or fasta
                if self.infile.endswith((".fastq.gz")):
                    records = SeqIO.parse(y, "fastq")
                    seq_dict = {}  # Create a dictionary to store data from the file
                    for record in records:
                        seq_dict.update({record.id: [str(record.seq), record.format("fastq").split("\n")[3]]})
                elif self.infile.endswith((".fasta.gz", ".fna.gz", ".fas.gz", ".fa.gz")):
                    records = SeqIO.parse(y, "fasta")
                    seq_dict = {}
                    for record in records:
                        seq_dict.update({str(record.id): str(record.seq)})
                y.close()
                return seq_dict
        elif self.infile.endswith((".fastq", ".fasta", ".fna", ".fas", ".fa")):
            with open(self.infile, "r") as y:
                if self.infile.endswith((".fastq")):
                    records = SeqIO.parse(y, "fastq")
                    seq_dict = {}  # Create a dictionary to store data from the file
                    for record in records:
                        seq_dict.update({record.id: [str(record.seq), record.format("fastq").split("\n")[3]]})
                elif self.infile.endswith((".fasta", ".fna", ".fas", ".fa")):
                    records = SeqIO.parse(y, "fasta")
                    seq_dict = {}
                    for record in records:
                        seq_dict.update({str(record.id): str(record.seq)})
                y.close()
                return seq_dict
        else:
            raise ValueError("File is the wrong format")
        print("Done", file=sys.stderr)

    def size_filter(self):
        """
        Filters sequences based on size criteria and writes the selected sequences to a new file.
        Prints a summary of the size selection process.
        """
        print("Size selection process of sequences with length over " + str(round(self.minsize*0.9, 0)) +
              " bp and under " + str(round(self.maxsize*1.1, 0)) + " has been started...", file=sys.stderr)
        
        with open(self.infile, "r") as ifp:
            ls = ifp.readlines()
        
        ifp.close()
        ext = os.path.splitext(self.infile)[1]
        
        if ext == ".fasta":
            outfile = os.path.join(FileSystem.get_base_dir(self.infile)[0], FileSystem.get_base_dir(self.infile)[1].split(".")[0] + "_size_selected.fasta")
            inf = Size_Filter(self.infile, self.minsize, self.maxsize)
            readdict = inf.load_data()
            
            with open(outfile, "w") as ofp:
                print("Writing size-selected sequences from %s to %s" % (self.infile, outfile), file=sys.stderr)
                total = int(len(list(readdict.keys())))
                ass = 0
                
                for el in list(readdict.keys()):
                    if self.minsize * 0.9 <= len(readdict[el]) <= self.maxsize * 1.1:
                        ofp.write(">" + el + "\n")
                        ofp.write(readdict[el] + "\n")
                    else:
                        ass += 1
                        pass
            print("Finished writing", file=sys.stderr)
            ofp.close()
            
            print("Size selection process ended: %d sequences were discarded out of %d total sequences (%g percent)" %
                  (ass, total, 100 - (round(1 - ass / total, 4)) * 100), file=sys.stderr)
        
        if ext == ".fastq":
            outfile = os.path.join(FileSystem.get_base_dir(self.infile)[0], FileSystem.get_base_dir(self.infile)[1].split(".")[0] + "_size_selected.fastq")
            inf = Size_Filter(self.infile, self.minsize, self.maxsize)
            readdict = inf.load_data()
            total = int(len(list(readdict.keys())))
            ass = 0
            
            with open(outfile, "w") as ofp:
                print("Writing size-selected sequences from %s to %s" % (self.infile, outfile), file=sys.stderr)
                for el in list(readdict.keys()):
                    if self.minsize * 0.9 <= len(readdict[el][0]) <= self.maxsize * 1.1:
                        ofp.write("@" + el + "\n")
                        ofp.write(readdict[el][0] + "\n")
                        ofp.write("+" + el + "\n")
                        ofp.write(readdict[el][1] + "\n")
                    else:
                        ass += 1
                        pass
            print("Finished writing", file=sys.stderr)
            ofp.close()
            
            print("Size selection process ended: %d sequences were discarded out of %d total sequences (%g percent)" %
                  (ass, total, 100 - (round(1 - ass / total, 4)) * 100), file=sys.stderr)
