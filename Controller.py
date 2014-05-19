#!/usr/bin/env python
""" CONTROLLER """

import datetime

from clint.textui import colored, puts

from NCBIQuery import NCBIQuery
from SequenceAnalysis import SequenceAnalysis, SingleTSeqRecord

""" Utility Functions """
def return_date_from_now(offset):
    """
    :param offset: Number of days to roll back the clock
    :return: Date from today after applying offset
    """
    now = datetime.date.today()
    delta = datetime.timedelta(days=offset)
    now = now - delta
    return now.strftime("%Y/%m/%d")

class Controller(object):

    @staticmethod
    def main():
        query = NCBIQuery()
        record = query.server_info()
        puts(colored.cyan(record["DbInfo"]["Description"]))
        puts(colored.blue("Last Updated: " + record["DbInfo"]["LastUpdate"]))
        puts(colored.blue("Total Records: " + record["DbInfo"]["Count"]))

        puts()
        puts(colored.green("** CODON Prototype **"))
        puts()
        puts("1. Galdieria sulphuraria (5 minute download)")
        puts("2. Opuntia[orgn] and rpl16 (5 second download)")
        puts("3. New and Updated Genomes Within Last 2 Days")
        puts()

        choice = raw_input(colored.cyan("Please enter a species: "))

        if int(choice) == 1:
            term = 'Galdieria sulphuraria'
        elif int(choice) == 2:
            term = 'Opuntia[orgn] and rpl16'
        elif int(choice) == 3:
            term = ("all[FILT] AND (" + return_date_from_now(2) +
                    "[PDAT]:" + return_date_from_now(0) +
                    "[PDAT] OR " + return_date_from_now(2) +
                    "[MDAT]:" + return_date_from_now(0) +
                    "[MDAT])"
                    )
        else:
            term = choice

        print "Searching for: ", term
        print

        results = query.search(term)

        for tseq in results:
            table = {'NCBI Assession Number':tseq.tseq_accver.string,
                     'NCBI GenInfo Number': tseq.tseq_gi.string,
                     'Definition': tseq.tseq_defline.string,
                     'Organism': tseq.tseq_orgname.string,
                     'Sequence Type': tseq.tseq_seqtype['value']}
            #'Sequence': tseq.tseq_sequence.string}
            assert str(tseq.tseq_seqtype['value']) == 'nucleotide'
            print
            for name, item in table.items():
                print '{0:10} ==> {1}'.format(name, item)

            assert int(tseq.tseq_length.string) == len(tseq.tseq_sequence.string)

        codonanalysis = SequenceAnalysis().add(results)
        seqs = codonanalysis.get_seq_stack

        for seq in seqs:
            assert isinstance(seq, SingleTSeqRecord)
            seq.print_dictionary()


if __name__ == "__main__":
    Controller.main()
