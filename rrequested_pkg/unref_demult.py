import sys
from Bio.Seq import Seq
import os
import edlib
from rrequested_pkg.file_system import FileSystem
from rrequested_pkg.size_filter import Size_Filter

class Unreferenced_Demultiplexing:
    def __init__(self, infile):
        """
        Initializes the Unreferenced_Demultiplexing class with the input file path.

        Parameters:
        - infile (str): Path to the input file.
        """
        self.infile = infile

    @staticmethod
    def find_n_remove(l1, l2):
        """
        Finds and removes elements from l2 that are present in l1.

        Parameters:
        - l1 (list): List containing elements to be found and removed.
        - l2 (list): List from which elements are to be removed.

        Returns:
        - list: Modified list l2 after removal of elements.
        """
        for j in range(len(l1)):
            for k in range(len(l2)):
                if k < len(l2) and l1[j] == l2[k]:
                    l2.remove(l2[k])
                else:
                    continue
        return l2

    @staticmethod
    def reverse_complement(seq: str):
        """
        Returns the reverse complement of a DNA sequence (also degenerate).

        Parameters:
        - seq (str): Input DNA sequence.

        Returns:
        - str: Reverse complement of the input DNA sequence.
        """
        dna = Seq(seq)
        return str(dna.reverse_complement())

    def find_the_num(self):
        """
        Finds highly divergent sequences (>50% different from any other) from the input file.

        Returns:
        - list: List of highly divergent sequences.
        """
        print("Started the search for highly diverging sequences (>50% apart from any other)", file=sys.stderr)
        the_num = []
        inf = Size_Filter(self.infile, 0, 0)
        readdict = inf.load_data()
        ks = list(readdict.keys())
        if type(readdict[ks[0]]) == list:
            lines = [val[0] for val in list(readdict.values())]
        if type(readdict[ks[0]]) == str:
            lines = [val for val in list(readdict.values())]
        try:
            for i in range(len(lines)):
                if len(the_num) == 0:
                    the_num.append(lines[i])
                else:
                    c = 0
                    for j in the_num:
                        alignment1 = edlib.align(lines[i], j, k=len(j) * 0.5, task="distance", mode="NW")
                        if alignment1["editDistance"] == -1:
                            alignment2 = edlib.align(lines[i], Unreferenced_Demultiplexing.reverse_complement(j),
                                                     k=len(j) * 0.5, task="distance", mode="NW")
                            if alignment2["editDistance"] == -1:
                                c += 1
                            else:
                                break
                        else:
                            break
                        if c == len(the_num):
                            the_num.append(lines[i])
                        else:
                            continue
            print("Search finished:\nNumber of analyzed sequences: %d\nNumber of highly divergent sequences: %d" %
                  (len(lines), len(the_num)), file=sys.stderr)
            return the_num
        except KeyboardInterrupt:
            sys.exit()

    @staticmethod
    def find_the_num_list(lines):
        """
        Finds highly divergent sequences (>50% different from any other) from a given list of sequences.

        Parameters:
        - lines (list): List of sequences.

        Returns:
        - list: List of highly divergent sequences.
        """
        print("Started the search for highly diverging sequences (>50% apart from any other)", file=sys.stderr)
        the_num = []
        try:
            seqs = 0
            for i in range(len(lines)):
                seqs += 1
                if len(the_num) == 0:
                    the_num.append(lines[i])
                else:
                    c = 0
                    for j in the_num:
                        alignment1 = edlib.align(lines[i], j, k=len(j) * 0.5, task="distance", mode="NW")
                        if alignment1["editDistance"] == -1:
                            alignment2 = edlib.align(lines[i], Unreferenced_Demultiplexing.reverse_complement(j),
                                                     k=len(j) * 0.5, task="distance", mode="NW")
                            if alignment2["editDistance"] == -1:
                                c += 1
                            else:
                                break
                        else:
                            break
                        if c == len(the_num):
                            the_num.append(lines[i])
                        else:
                            continue
            print("Search finished:\nNumber of analyzed sequences: %d\nNumber of highly divergent sequences: %d" %
                  (seqs, len(the_num)), file=sys.stderr)
            return the_num
        except KeyboardInterrupt:
            sys.exit()

    def demultiplex(self):
        """
        Demultiplexes sequences from the input file into separate groups based on their similarity.
        Writes the demultiplexed sequences into individual files.
        """
        try:
            inf = Unreferenced_Demultiplexing(self.infile)
            value = False
            refseqs = inf.find_the_num()
            groups = []
            ind = -1
            with open(self.infile, "r") as f:
                lines = f.readlines()
            f.close()
            nf = Size_Filter(self.infile, 0, 0)
            readdict = nf.load_data()
            ks = list(readdict.keys())
            if type(readdict[ks[0]]) == list:
                seq = [val[0] for val in list(readdict.values())]
            if type(readdict[ks[0]]) == str:
                seq = [val for val in list(readdict.values())]
            tot = len(seq)
            print("Pairwise alignment assigning of raw reads to the highly divergent sequences started", file=sys.stderr)
            processed = []
            for j in refseqs:
                groups.append([])
                ind += 1
                for i in range(len(seq)):
                    alignment1 = edlib.align(seq[i], j, k=len(j) * 0.3, task="distance", mode="NW")
                    if alignment1["editDistance"] == -1:
                        alignment2 = edlib.align(seq[i], Unreferenced_Demultiplexing.reverse_complement(j),
                                                 k=len(j) * 0.3, task="distance", mode="NW")
                        if alignment2["editDistance"] == -1:
                            continue
                        else:
                            groups[ind].append(seq[i])
                            processed.append(seq[i])
                    else:
                        groups[ind].append(seq[i])
                        processed.append(seq[i])
                print("Assigned %d sequences to group %d (%g perc)" % (len(groups[ind]), ind,
                                                                       round(((len(groups[ind]) / len(seq)) * 100), 3)),
                      file=sys.stderr)
            no_group = Unreferenced_Demultiplexing.find_n_remove(processed, seq)
            run = 0
            if len(no_group) > 0:
                print("Non-grouped sequences are %d, starting reassignment..." % (len(no_group)), file=sys.stderr)
                while len(no_group) > 0 and run < 5:
                    print("Started round %d" % (run + 1), file=sys.stderr)
                    refs = Unreferenced_Demultiplexing.find_the_num_list(no_group)
                    ranks = []
                    proc = []
                    indx = -1
                    for j in refs:
                        ranks.append([])
                        indx += 1
                        for i in range(len(no_group)):
                            alignment1 = edlib.align(no_group[i], j, k=len(j) * 0.3, task="distance", mode="NW")
                            if alignment1["editDistance"] == -1:
                                alignment2 = edlib.align(no_group[i], Unreferenced_Demultiplexing.reverse_complement(j),
                                                         k=len(j) * 0.3, task="distance", mode="NW")
                                if alignment2["editDistance"] == -1:
                                    continue
                                else:
                                    ranks[indx].append(no_group[i])
                                    proc.append(no_group[i])
                            else:
                                ranks[indx].append(no_group[i])
                                proc.append(no_group[i])
                        print("Assigned %d sequences to no_group %d (%g perc)" % (len(ranks[indx]), indx,
                                                                                   round(((len(ranks[indx]) / len(no_group)) * 100), 3)),
                              file=sys.stderr)
                    for r in ranks:
                        groups.append(r)
                    no_group = Unreferenced_Demultiplexing.find_n_remove(proc, no_group)
                    run += 1
                if run >= 4 and len(no_group) > 0:
                    print("Non-grouped sequences are still %d (%g perc), they will be grouped together" %
                          (len(no_group), round((len(no_group) / tot) * 100, 3)), file=sys.stderr)
                    groups.append(no_group)
                    value = True
            else:
                print("All sequences have been demultiplexed at the first attempt, proceeding with the analysis")
            path = FileSystem.makedir_orchange(
                os.path.join(FileSystem.get_base_dir(self.infile)[0],
                             FileSystem.get_base_dir(self.infile)[1] + "-" + os.path.splitext(self.infile)[1][1] +
                             os.path.splitext(self.infile)[1][len(os.path.splitext(self.infile)[1]) - 1] + "-demultiplexed"))
            print("Creating a new folder at " + path + " to write demultiplexed files in there", file=sys.stderr)
            print("Writing demultiplexed sequences in separated files", file=sys.stderr)
            count = 0
            for g in groups:
                if len(g) >= (tot) * 0.01 and groups.index(g) != len(groups) - 1:
                    count += 1
                    with open(os.path.join(path, str(count) + ".fasta"), "w") as fp:
                        c = 1
                        for s in g:
                            fp.write(">" + str(c) + "\n")
                            fp.write(s + "\n")
                            c += 1
                if len(g) >= (tot) * 0.01 and groups.index(g) == len(groups) - 1 and value:
                    with open(os.path.join(path, "nogroup.fasta"), "w") as fp:
                        c = 1
                        for s in g:
                            fp.write(">" + str(c) + "\n")
                            fp.write(s + "\n")
                            c += 1
                if len(g) >= (tot) * 0.01 and groups.index(g) == len(groups) - 1 and value == False:
                    count += 1
                    with open(os.path.join(path, str(count) + ".fasta"), "w") as fp:
                        c = 1
                        for s in g:
                            fp.write(">" + str(c) + "\n")
                            fp.write(s + "\n")
                            c += 1
                else:
                    continue
            print("The program ended its execution successfully", file=sys.stderr)
        except KeyboardInterrupt:
            sys.exit()
