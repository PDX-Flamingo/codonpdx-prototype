# !/usr/bin/env python
""" Codon Usage Analysis"""

from clint.textui import puts, colored
from collections import deque
from random import randint


class SingleTSeqRecord(object):
    organism = None
    definition = None
    ncbi_taxid = None  # Taxonomy id
    ncbi_geninfo = None  # Unique To Every New and Updated Sequence
    ncbi_assession = None  # Remains Constant on Update
    ncbi_seqlength = None


    def __init__(self):
        """

        :rtype :
        """
        self.data = []
        self.codons = {'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0,
                       'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0,
                       'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0,
                       'GTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0,
                       'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0,
                       'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
                       'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0,
                       'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
                       'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0,
                       'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0,
                       'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0,
                       'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
                       'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}
        self.__token_queue__ = deque()


    """Codon Dictionary Utility Functions"""

    def __clear_dict__(self):
        for key, value in self.codons.iteritems():
            self.codons[key] = 0
            return

    def is_dict_empty(self):
        result = True
        for key, value in self.codons.iteritems():
            if value != 0:
                result = False
        return result

    def randomize_dict(self):
        for key, value in self.codons.iteritems():
            self.codons[key] = randint(0, 50)
        return

    def print_dictionary(self):
        for key, value in self.codons.iteritems():
            puts(colored.green(key) + " : " + colored.cyan(value) + ", ",newline=False)
        puts()
        return

    """TOKENIZE (Optimized)"""

    def __tokenize__(self, sequence):
        indexes = [(3*x, (3*x)+3) for x in range(len(sequence)/3)]
        for index in indexes:
            self.__token_queue__.append(sequence[index[0]:index[1]])
        return

    """Load"""

    def load(self, seq):
        self.organism = seq.tseq_orgname.string
        self.ncbi_geninfo = seq.tseq_gi.string
        self.ncbi_assession = seq.tseq_accver.string
        self.ncbi_taxid = seq.tseq_taxid.string
        self.definition = seq.tseq_defline.string
        self.ncbi_seqlength = seq.tseq_length.string
        self.__tokenize__(str(seq.tseq_sequence.string))

        ### Process Tokens
        while self.__token_queue__:
            token = self.__token_queue__.popleft()
            if self.codons.has_key(token):
                self.codons[token] += 1


class SequenceAnalysis(SingleTSeqRecord):

    def __init__(self):
        SingleTSeqRecord.__init__(self)
        self.data = []
        self.seq_stack = []
        self.total_seq_length = 0




    def add(self, batch_sequence):
        for seq in batch_sequence:
            req = SingleTSeqRecord()
            req.load(seq)
            self.total_seq_length += self.num(seq.tseq_length.string)
            self.seq_stack.append(req)
        puts(colored.yellow("Total Length of Sequences: " + str(self.total_seq_length)))
        return self


    @property
    def get_seq_stack(self):
        return self.seq_stack
    
    @property
    def pop(self):
        return self.seq_stack.pop()

    @staticmethod
    def num(s):
        try:
            return int(s)
        except ValueError:
            return float(s)

    @staticmethod
    def main():
        puts(colored.yellow('** MAIN CLASS **'))
        return


if __name__ == "__main__":
    SequenceAnalysis.main()