""" Unit Tests for SequenceAnalysis """

import os.path
import re
import time

from NCBIQuery import NCBIQuery
from nose.tools import ok_, eq_, istest


class TestNCBIQuery():
    term_all = "All[Filter]"
    term1 = "KJ408799[Accession]"
    term1_seq = ("ACCAGCGCACTTCGGCAGCGGCAGCACCTCGGCAGCGTC"
                 "AGTGAAAATGCCAAGCAAGAAAAGCGGCCCGCAACCCCA"
                 "TAAGAGGTGGGTGTTCACCCTTAATAATCCTTCCGAGGA"
                 "GGAGAAAAACAAAATACGGGAGCTTCCAATCTCCCTTTT"
                 "TGATTATTTTGTTTGCGGAGAGGAAGGTTTGGAAGAGGG"
                 "TAGAACTCCTCACCTCCAGGGGTTTGCGAATTTTGCTAA"
                 "GAAGCAGACTTTTAACAAGGTGAAGTGGTATTTTGGTGC"
                 "CCGCTGCCACATCGAGAAAGCGAAAGGAACCGACCAGCA"
                 "GAATAAAGAATACTGCAGTAAAGAAGGCCACATACTTAT"
                 "CGAGTGTGGAGCTCCGCGGAACCAGGGGAAGCGCAGCGA"
                 "CCTGTCTACTGCTGTGAGTACCCTTTTGGAGACGGGGTC"
                 "TTTGGTGACTGTAGCCGAGCAGTTCCCTGTAACGTATGT"
                 "GAGAAATTTCCGCGGGCTGGCTGAACTTTTGAAAGTGAG"
                 "CGGGAAGATGCAGCAGCGTGATTGGAAGACAGCTGTACA"
                 "CGTCATAGTGGGCCCGCCCGGTTGTGGGAAGAGCCAGTG"
                 "GGCCCGTAATTTTGCTGAGCCTAGCGACACCTACTGGAA"
                 "GCCTAGTAGAAATAAGTGGTGGGATGGATATCATGGAGA"
                 "AGAAGTTGTTGTTTTGGATGATTTTTATGGCTGGTTACC"
                 "TTGGGATGATCTACTGAGACTGTGTGACCGGTATCCATT"
                 "GACTGTAGAGACTAAAGGCGGTACTGTTCCTTTTTTGGC"
                 "CCGCAGTATTTTGATTACCAGCAATCAGGCCCCCCAGGA"
                 "ATGGTACTCCTCAACTGCTGTCCCAGCTGTAGAAGCTCT"
                 "CTATCGGAGGATTACTACTTTGCAATTTTGGAAGACTGC"
                 "TGGAGAACAATCCACGGAGGTACCCGAAGGCCGATTTGA"
                 "AGCAGTGGACCCACCCTGTGCCCTTTTCCCATATAAAAT"
                 "AAATTACTGAGTCTTTTTTGTTATCACATCGTAATGGTT"
                 "TTTATTTTTATTTATTTAGAGGGTCTTTTAGGATAAATT"
                 "CTCTGAATTGTACATAAATAGTCAGCCTTACCACATAAT"
                 "TTTGGGCTGTGGCTGCATTTTGGAGCGCATAGCCGAGGC"
                 "CTGTGTGCTCGACATTGGTGTGGGTATTTAAATGGAGCC"
                 "ACAGCTGGTTTCTTTTATTATTTGGTTGGAACCAATCAA"
                 "TTGTTTGGTCCAGCTCAGGTTTGGGGGTGAAGTACCTGG"
                 "AGTGGTAGGTAAAGGGCTGCCTTATGGTGTGGCGGGAGG"
                 "AGTAGTTAATATAGGGGTCATAGGCCAAGTTGGTGGAGG"
                 "GGGTTACAAAGTTGGCATCCAAGATAACAACAGTGGACC"
                 "CAACACCTCTTTCATTAGAGGTGATGGGGTCTCTGGGGT"
                 "AAAATTCATATTTAGCCTTTCTAATACGGTAGTATTGGA"
                 "AAGGTAGGGGTAGGGGGTTGGTGCCGCCTGAGGGGGGGA"
                 "GGAACTGGCCGATGTTGAATTTGAGGTGGTTAACATGCC"
                 "AAGATGGCTGCGAGTATCCTCCTTTTATGGTGAGTACAA"
                 "ATTCTGTAGAAAGGCGGGAATTGAAGATACCCGTCTTTC"
                 "GGCGCCATCTGTAACGGTTTCTGAAGGCGGGGTGTGCCA"
                 "AATATGGTCTTCTCCGGAGGATGTTTCCAAGATGGCTGC"
                 "GGGGGCGGGTCCTTCGTCTGCGGTAACGCCTCCTTGGCC"
                 "ACGTCATCCTATAAAAGTGAAAGAAGTGCGCTGCTGTAG"
                 "TATT"
                 )

    def __init__(self):
        self.data = []

    def setup(self):
        pass

    def teardown(self):
        tinyseqxml = self.term1 + '.tseqxml'
        if os.path.isfile(tinyseqxml):
            os.remove(tinyseqxml)

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    @istest
    def test_server_info(self):
        obj = NCBIQuery()
        info = obj.server_info()
        eq_(info["DbInfo"]["Description"], "Core Nucleotide db", str(info["DbInfo"]["Description"]))

        try:
            time.strptime(info["DbInfo"]["LastUpdate"], '%Y/%m/%d  %H:%M')
        except ValueError:
            ok_(False, info["DbInfo"]["LastUpdate"] + " Failed time.strptime format")
        ok_(True, info["DbInfo"]["LastUpdate"])

        ok_(info["DbInfo"]["Count"] >= 132454040)

    @istest
    def test_search_info(self):
        results = NCBIQuery().search_info(self.term1)
        eq_(results["Count"], 1)
        eq_(results["IdList"], ['594139292'])
        eq_(results["ReturnMax"], 1)
        eq_(results["QueryKey"], "1")
        ok_(re.match('NCID.+', results["WebEnv"]), results["WebEnv"])

    @istest
    def test_search_info_all(self):
        results = NCBIQuery().search_info(self.term_all)
        ok_(results["Count"] >= 132454040)
        eq_(results["ReturnMax"], 20)
        eq_(results["QueryKey"], "1")
        ok_(re.match('NCID.+', results["WebEnv"]), results["WebEnv"])

    @istest
    def test_search(self):
        results = NCBIQuery().search(self.term1)
        eq_(results[0].tseq_gi.string, "594139292")
        eq_(results[0].tseq_accver.string, "KJ408799.1")
        eq_(results[0].tseq_taxid.string, "133704")
        eq_(results[0].tseq_orgname.string, "Porcine circovirus 1")
        eq_(results[0].tseq_defline.string, "Porcine circovirus 1 isolate PCV1-Hun, complete genome")
        eq_(results[0].tseq_length.string, "1759")
        eq_(results[0].tseq_sequence.string, self.term1_seq)