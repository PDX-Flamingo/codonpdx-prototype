from collections import deque
from mrjob.job import MRJob


class TokenizeMRJob(MRJob):

    def mapper(self, _, line):
        indexes = deque((3*x, (3*x)+3) for x in range(len(line)/3))
        while indexes:
            index = indexes.popleft()
            yield line[index[0]:index[1]], 1

    def combiner(self, token, counts):
        yield token, sum(counts)

    def reducer(self, token, counts):
        yield token, sum(counts)


if __name__ == '__main__':
    TokenizeMRJob.run()