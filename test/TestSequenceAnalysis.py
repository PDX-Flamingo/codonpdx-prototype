__author__ = 'ptseng'
""" Unit Tests for SequenceAnalysis """

import SequenceAnalysis

from nose.tools import eq_, istest


class TestSequenceAnalysis():

    def __init__(self):
        pass

    def setup(self):
        print ("TestUM:setup() before each test method")

    def teardown(self):
        print ("TestUM:teardown() after each test method")

    @classmethod
    def setup_class(cls):
        print ("setup_class() before any methods in this class")

    @classmethod
    def teardown_class(cls):
        print ("teardown_class() after any methods in this class")

    @istest
    def test_initial_codon_dictionary_is_empty(self):
        seqana = SequenceAnalysis.SequenceAnalysis()
        eq_(64, len(seqana.CodonsDict))
        eq_(seqana.is_dict_empty(), True)
        seqana.randomize_dict()
        eq_(seqana.is_dict_empty(), False)

    @istest
    def test_print_codon_dictionary(self):
        seqana = SequenceAnalysis.SequenceAnalysis()
        seqana.randomize_dict()
        seqana.print_dictionary()
