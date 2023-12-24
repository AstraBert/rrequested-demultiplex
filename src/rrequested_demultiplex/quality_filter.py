import sys
from Bio import SeqIO
import gzip

class Quality_Filter:
    def __init__(self, infile, quality):
        """
        Initializes the Quality_Filter class with input file path and quality threshold.

        Parameters:
        - infile (str): Path to the input file.
        - quality (float): Quality threshold for filtering reads.
        """
        self.infile = infile
        self.quality = quality

    def load_data(self):
        """
        Load data from the input file if it is in fastq format (unzips if it is zipped).

        Returns:
        - dict: A dictionary containing read headers as keys and corresponding sequences and quality strings as values.
        """
        print("Reading data from input file...", file=sys.stderr)
        if self.infile.endswith(".gz"):  # If file is gzipped, unzip it
            y = gzip.open(self.infile, "rt", encoding="latin-1")
            if self.infile.endswith(".fastq.gz"):  # Read file as fastq if it is fastq
                records = SeqIO.parse(y, "fastq")
                seq_dict = {}  # Create a dictionary to store data from the file
                for record in records:
                    seq_dict.update({record.id: [str(record.seq), record.format("fastq").split("\n")[3]]})
                y.close()
                return seq_dict
        elif self.infile.endswith(".fastq"):
            with open(self.infile, "r") as y:
                records = SeqIO.parse(y, "fastq")
                seq_dict = {}  # Create a dictionary to store data from the file
                for record in records:
                    seq_dict.update({record.id: [str(record.seq), record.format("fastq").split("\n")[3]]})
                y.close()
                return seq_dict
        else:
            raise ValueError("File is the wrong format")
        print("Done", file=sys.stderr)

    @staticmethod
    def ascii_conv_and_mean(line):
        """
        Converts ASCII characters to Phred quality scores and calculates the mean.

        Parameters:
        - line (str): Phred quality string.

        Returns:
        - float: Mean Phred quality score.
        """
        phred_quality_dict = {'!': 0, '"': 1, '#': 2, '$': 3, '%' : 4, '&' : 5, "'" : 6, '(' : 7, ')' : 8, '*' : 9, '+' : 10,
                              ',' : 11, '-' : 12, '.' : 13, '/': 14, '0': 15, '1': 16, '2': 17, '3': 18, '4': 19, '5': 20,
                              '6': 21, '7': 22, '8': 23, '9': 24, ':': 25, ';': 26, '<': 27, '=': 28, '>': 29, '?': 30,
                              '@': 31, 'A': 32, 'B': 33, 'C': 34, 'D': 35, 'E': 36, 'F': 37, 'G': 38, 'H': 39, 'I': 40,
                              'J': 41, 'K': 42, 'L': 43, 'M': 44, 'N': 45, 'O': 46, 'P': 47, 'Q': 48, 'R': 49, 'S': 50,
                              'T': 51, 'U': 52, 'V': 53, 'W': 54, 'X': 55, 'Y': 56, 'Z': 57, '[': 58, '\\': 59, ']': 60,
                              '^': 61, '_': 62, '`': 63, 'a': 64, 'b': 65, 'c': 66, 'd': 67, 'e': 68, 'f': 69, 'g': 70,
                              'h': 71, 'i': 72, 'j': 73, 'k': 74, 'l': 75, 'm': 76, 'n': 77, 'o': 78, 'p': 79, 'q': 80,
                              'r': 81, 's': 82, 't': 83, 'u': 84, 'v': 85, 'w': 86, 'x': 87, 'y': 88, 'z': 89, '{': 90,
                              '|': 91, '}': 92, '~': 93}
        mean_list = [phred_quality_dict[i] for i in line if i != '\n']
        return round(sum(mean_list) / len(mean_list), 3)

    def filter(self):
        """
        Filters reads based on the quality threshold and updates the input file with filtered reads.
        Prints the summary of the filtering process.
        """
        print("Filtering of reads over quality threshold " + str(self.quality) + " has been started", file=sys.stderr)
        fastq = Quality_Filter(self.infile, self.quality)
        readdict = fastq.load_data()
        with open(self.infile, "r+") as f:
            f.truncate()
        disc = 0
        total = int(len(list(readdict.keys())))
        with open(self.infile, "w") as fp:
            for i in list(readdict.keys()):
                if Quality_Filter.ascii_conv_and_mean(readdict[i][1]) < self.quality:
                    disc += 1
                else:
                    fp.write("@" + i + "\n")
                    fp.write(readdict[i][0] + "\n")
                    fp.write("+" + i + "\n")
                    fp.write(readdict[i][1] + "\n")
        print("Filtering finished: %d reads were discarded out of %d (%g percent), now re-compiling the file with filtered reads..." %
              (disc, total, 100 * (round(disc / total, 4))), file=sys.stderr)
