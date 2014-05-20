#!/usr/bin/env python
""" CONTROLLER """

import datetime
import time

from clint.textui import colored, puts

from NCBIQuery import NCBIQuery
from SequenceAnalysis import SequenceAnalysis, SingleTSeqRecord

""" Utility Functions """
CURRENT_MILLI_TIME = lambda: int(round(time.clock() * 1000))

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
        puts("1. Galdieria sulphuraria")
        puts("2. Opuntia")
        puts("3. New and Updated Genomes Within Last 4 Days")
        puts("4. Mycoplasma genitalium")
        puts("5. Pseudomonas aeruginosa")
        puts("6. Human Chromosome 2")
        puts()
        puts()

        choice = raw_input(colored.cyan("Please enter a species: "))

        if int(choice) == 1:
            term = 'Galdieria sulphuraria'
            option = '[Organism] AND srcdb_genbank[PROP]'
        elif int(choice) == 2:
            term = 'Opuntia'
            option = '[Organism] AND srcdb_genbank[PROP]'
        elif int(choice) == 3:
            term = 'all'
            option = ("[FILT] AND (" + return_date_from_now(4) +
                    "[PDAT]:" + return_date_from_now(0) +
                    "[PDAT] OR " + return_date_from_now(4) +
                    "[MDAT]:" + return_date_from_now(0) +
                    "[MDAT]) AND srcdb_genbank[PROP]"
                      )
        elif int(choice) == 4:
            term = 'Mycoplasma genitalium'
            option = '[Organism] AND srcdb_genbank[PROP]'
        elif int(choice) == 5:
            term = 'BAUN01000001[Accession]'
            option = ''
        elif int(choice) == 6:
            term = 'KE150089.1[Accession]'
            option = ''
        else:
            term = choice
            option = ''

        print "Searching for: ", term
        print

        info = query.search_info(term, option)
        puts(colored.yellow("Fetching " + str(info["Count"]) + " Results"))
        results = query.search(term, option)

        '''
        for tseq in results:
            table = {'NCBI Assession Number':tseq.tseq_accver.string,
                     'NCBI GenInfo Number': tseq.tseq_gi.string,
                     'Definition': tseq.tseq_defline.string,
                     'Length': tseq.tseq_length.string,
                     #'Organism': tseq.tseq_orgname.string,
                     'Sequence Type': tseq.tseq_seqtype['value']}
            #'Sequence': tseq.tseq_sequence.string}
            assert str(tseq.tseq_seqtype['value']) == 'nucleotide'
            print
            for name, item in table.items():
                print '{0:10} ==> {1}'.format(name, item)

            assert int(tseq.tseq_length.string) == len(tseq.tseq_sequence.string)
        '''

        start = CURRENT_MILLI_TIME()
        codonanalysis = SequenceAnalysis().add(results)
        end = CURRENT_MILLI_TIME()
        print "Tokenize and Analyze in " + str(end - start) + " ms"
        seqs = codonanalysis.get_seq_stack

        for seq in seqs:
            assert isinstance(seq, SingleTSeqRecord)
            #seq.print_dictionary()


if __name__ == "__main__":
    Controller.main()
