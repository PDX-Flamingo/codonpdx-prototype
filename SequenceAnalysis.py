__author__ = 'ptseng'

from clint.textui import puts, colored
from random import randint


class SequenceAnalysis(object):

    Organism = ''
    Definition = ''
    NCBI_GenInfo = ''
    NCBI_Assession = ''

    CodonsDict = {'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0,
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

    def __init__(self):
        self.data = []
        self.__clear_dict__()

    """Codon Dictionary Utility Functions"""
    def __clear_dict__(self):
        for key, value in self.CodonsDict.iteritems():
            self.CodonsDict[key] = 0
            return

    def is_dict_empty(self):
        result = True
        for key, value in self.CodonsDict.iteritems():
            if value != 0:
                result = False
        return result

    def randomize_dict(self):
        for key, value in self.CodonsDict.iteritems():
            self.CodonsDict[key] = randint(0, 50)
        return

    def print_dictionary(self):
        for key, value in self.CodonsDict.iteritems():
            puts(colored.green(key) + " : " + colored.cyan(value))
        return

    @staticmethod
    def main():
        puts(colored.yellow('** MAIN CLASS **'))
        return


if __name__ == "__main__":
    SequenceAnalysis.main()