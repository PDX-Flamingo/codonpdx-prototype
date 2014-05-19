#!/usr/bin/env python
""" This Script Queries Genbank """

import sys
import time

import os
import os.path

from Bio_Eutils import Entrez
from bs4 import BeautifulSoup
from progressbar import Bar, ETA, \
    FileTransferSpeed, Percentage, \
    ProgressBar, RotatingMarker

### Utility Methods
CURRENT_MILLI_TIME = lambda: int(round(time.clock() * 1000))


class NCBIQuery():

    def __init__(self, email="ptseng@pdx.edu", recursionlimit=20000):
        self.data = []
        sys.setrecursionlimit(recursionlimit)
        #self.rows, self.columns = os.popen('stty size', 'r').read().split()
        self.columns = 125
        self.email = email

    def server_info(self):
        Entrez.email = self.email
        handle = Entrez.einfo(db="nucleotide")
        record = Entrez.read(handle)
        handle.close()
        return record

    def search_info(self, term, options=''):
        Entrez.email = self.email
        search_handle = Entrez.esearch(db="nucleotide",
                                       term=term+options,
                                       usehistory="y")
        search_results = Entrez.read(search_handle)
        search_handle.close()

        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]

        gi_list = search_results["IdList"]
        retmax = int(search_results["RetMax"])
        count = int(search_results["Count"])

        assert retmax == len(gi_list)

        dict = {"Count": count, "IdList": gi_list, "ReturnMax": retmax,
                "WebEnv": webenv, "QueryKey": query_key}
        return dict

    def search(self, organism='', options='', batch_size=50, progress=True):
        Entrez.email = self.email
        search_handle = Entrez.esearch(db="nucleotide",
                                       term=organism+options,
                                       usehistory="y")
        search_results = Entrez.read(search_handle)
        search_handle.close()

        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]

        gi_list = search_results["IdList"]
        retmax = int(search_results["RetMax"])
        count = int(search_results["Count"])

        if count == 0:
            return

        assert retmax == len(gi_list)
        ### Download data referred to in the previous search in batches
        tinyseqxml = organism + '.tseqxml'
        if not os.path.isfile(tinyseqxml):
            widgets = [tinyseqxml + "\t", Percentage(), ' ',
                       Bar(marker=RotatingMarker()),
                       ' ', ETA(), ' ', FileTransferSpeed(unit='r')]
            out_handle = open(tinyseqxml, "w")
            pbar = ProgressBar(widgets=widgets,
                               maxval=count,
                               term_width=int(self.columns) - 20)
            for start in range(0, count, batch_size):
                #end = min(count, start + batch_size)
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
                if progress and pbar.seconds_elapsed == 0:
                    pbar.start()
                if status > 0.001 and progress:
                    pbar.term_width = int(self.columns) - 20
                    pbar.update(status * count)
            if progress:
                pbar.finish()
            out_handle.close()

        ## At this point, file may be a stream of multiple XML documents.
        ### Do not parse as one single XML document.

        ### Beautiful Soup PUNCH! (Think Donkey Kong)
        start = CURRENT_MILLI_TIME()
        soup = BeautifulSoup(open(tinyseqxml, "r"), "lxml")
        end = CURRENT_MILLI_TIME()
        print "Finished Parsing XML in " + str(end - start) + " ms"

        return soup.find_all('tseq')

