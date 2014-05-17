#!/usr/bin/env python
""" This Script Queries Genbank """

import datetime
import os
import os.path
import sys
import time

from Bio_Eutils import Entrez
from bs4 import BeautifulSoup
from clint.textui import colored, puts
from progressbar import Bar, ETA, \
    FileTransferSpeed, Percentage, \
    ProgressBar, RotatingMarker

### Utility Methods
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


### MAIN
def main():
    """
    :return: None
    """
    sys.setrecursionlimit(20000)
    rows, columns = os.popen('stty size', 'r').read().split()

    # ## Entrez Identification
    Entrez.email = "ptseng@pdx.edu"

    handle = Entrez.einfo(db="nucleotide")
    record = Entrez.read(handle)
    puts(colored.blue(record["DbInfo"]["Description"] +
                      ", Last Updated: " +
                      record["DbInfo"]["LastUpdate"]))
    puts(colored.blue("Total Records: " + record["DbInfo"]["Count"]))
    # keys = record["DbInfo"].keys()
    # print keys

    # for field in record["DbInfo"]["FieldList"]:
    #    print("%(Name)s, %(FullName)s, %(Description)s" % field)

    #print return_date_from_now(1)

    handle.close()

    puts()
    puts(colored.green("** Genbank Prototype **"))
    puts()
    puts("1. Galdieria sulphuraria (5 minute download)")
    puts("2. Opuntia[orgn] and rpl16 (5 second download)")
    puts("3. New and Updated Genomes Within Last 2 Days")
    puts()

    choice = raw_input(colored.cyan("Please enter a species: "))

    if int(choice) == 1:
        var = 'Galdieria sulphuraria'
        search_handle = Entrez.esearch(db="nucleotide",
                                       term=var,
                                       usehistory="y")
    elif int(choice) == 2:
        var = 'Opuntia[orgn] and rpl16'
        search_handle = Entrez.esearch(db="nucleotide",
                                       term=var,
                                       usehistory="y")
    elif int(choice) == 3:
        var = 'updatednew_last_2_days'
        search_handle = Entrez.esearch(db="nucleotide",
                                       term=
                                       "all[FILT] AND (" +
                                       return_date_from_now(2) +
                                       "[PDAT]:" + return_date_from_now(0) +
                                       "[PDAT] OR " + return_date_from_now(2) +
                                       "[MDAT]:" + return_date_from_now(0) +
                                       "[MDAT])",
                                       usehistory="y", )
    else:
        var = ''
        search_handle = None
    print "Searching for: ", var
    print

    search_results = Entrez.read(search_handle)
    search_handle.close()

    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]

    gi_list = search_results["IdList"]
    retmax = int(search_results["RetMax"])
    count = int(search_results["Count"])
    length_gilist = len(gi_list)

    assert retmax == length_gilist

    table = {'Total Records': count, 'Records Returned': length_gilist}
    for name, number in table.items():
        print '{0:20} ==> {1:10d}'.format(name, number)

    print

    ### Download data referred to in the previous search in batches
    tinyseqxml = var + '.tseqxml'
    if not os.path.isfile(tinyseqxml):
        batch_size = 50
        widgets = [tinyseqxml + "\t", Percentage(), ' ',
                   Bar(marker=RotatingMarker()),
                   ' ', ETA(), ' ', FileTransferSpeed(unit='r')]
        out_handle = open(tinyseqxml, "w")
        pbar = ProgressBar(widgets=widgets,
                           maxval=count,
                           term_width=int(columns) - 20)
        for start in range(0, count, batch_size):
            end = min(count, start + batch_size)
            #print("Going to download record %i to %i" % (start+1, end))
            fetch_handle = Entrez.efetch(db="nucleotide",
                                         rettype="fasta",
                                         retmode="xml",
                                         retstart=start,
                                         retmax=batch_size,
                                         webenv=webenv,
                                         query_key=query_key)
            data = fetch_handle.read()
            fetch_handle.close()
            out_handle.write(data)
            status = (start + 1) / float(count)
            if pbar.seconds_elapsed == 0:
                pbar.start()
            if status > 0.001:
                pbar.term_width = int(columns) - 20
                pbar.update(status * count)
        pbar.finish()
        out_handle.close()

    ## At this point, file is stream of multiple XML documents.
    ### Do not parse as one single XML document.

    print
    print "Parsing XML"
    ### Beautiful Soup PUNCH! (Think Donkey Kong)
    start = CURRENT_MILLI_TIME()
    soup = BeautifulSoup(open(tinyseqxml, "r"), "lxml")
    end = CURRENT_MILLI_TIME()
    print "Finished Parsing in " + str(end - start) + " ms"

    sequences = soup.find_all('tseq')
    seq_build = {}

    for tseq in sequences:
        table = {'NCBI Assession Number': tseq.tseq_accver.string,
                 'NCBI GenInfo Number': tseq.tseq_gi.string,
                 'Definition': tseq.tseq_defline.string}
                 #'Sequence': tseq.tseq_sequence.string}
        print
        for name, item in table.items():
            print '{0:10} ==> {1}'.format(name, item)
        assert int(tseq.tseq_length.string) == len(tseq.tseq_sequence.string)
        seq_build[tseq.tseq_gi.string] = tseq.tseq_sequence.string

    assert len(seq_build) == count

    print seq_build


if __name__ == "__main__":
    main()


