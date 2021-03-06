from __future__ import with_statement

""" Unit Tests for SequenceAnalysis """
import os

from StringIO import StringIO
from TokenizeMRJob import TokenizeMRJob
from collections import deque
from NCBIQuery import NCBIQuery
from SequenceAnalysis import SequenceAnalysis, SingleTSeqRecord
from nose.tools import eq_, istest


class TestSequenceAnalysis():

    term1 = "KJ408799[Accession]"
    seq = None

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
    def test_initial_codon_dictionary_is_empty(self):
        seqana = SequenceAnalysis()
        eq_(64, len(seqana.codons))
        eq_(seqana.is_dict_empty(), True)
        seqana.randomize_dict()
        eq_(seqana.is_dict_empty(), False)

    @istest
    def test_tokenize(self):
        seq = SequenceAnalysis()
        seq.__tokenize__("AAATTTCCCGGG")
        eq_(seq.__token_queue__, deque(['AAA', 'TTT', 'CCC', 'GGG']))

    @istest
    def test_load_single_sequence(self):
        tseq = NCBIQuery().search(self.term1, progress=False)
        seq = SequenceAnalysis()
        seq.add(tseq)

        for plate in seq.seq_stack:
            assert isinstance(plate, SingleTSeqRecord)
            eq_(plate.organism, 'Porcine circovirus 1')
            eq_(plate.definition, 'Porcine circovirus 1 isolate PCV1-Hun, complete genome')
            eq_(plate.ncbi_taxid, '133704')
            eq_(plate.ncbi_geninfo, '594139292')
            eq_(plate.ncbi_assession, 'KJ408799.1')
            eq_(plate.ncbi_seqlength, '1759')
            eq_(len(plate.__token_queue__), 0)  # Should be processed


    @istest
    def test_tokenize_mapreduce(self):
        mr_job = TokenizeMRJob(['-r', 'local', '--no-conf', '-'])
        mr_job.sandbox(stdin=StringIO("AAATTTCCCGGG"))

        tokens = []
        counts = []
        with mr_job.make_runner() as runner:
            runner.run()
            for line in runner.stream_output():
                key, value = mr_job.parse_output_line(line)
                tokens.append(key)
                counts.append(value)

        eq_(tokens, ['AAA', 'CCC', 'GGG', 'TTT'])
        eq_(counts, [1, 1, 1, 1])
