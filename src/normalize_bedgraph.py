#!/usr/bin/env python3

import sys
import gzip
import re
import collections
import argparse


class Bedgraph_element():

    def __init__(self, chrom, start, end, signal):
        assert(isinstance(chrom, str))
        assert(isinstance(start, int))
        assert(isinstance(end, int))
        assert(isinstance(signal, float) or isinstance(signal, int))
        self.chrom = chrom
        self.start = start
        self.end = end
        self.signal = float(signal)

    def __str__(self):
        return '\t'.join([str(x) for x in [self.chrom, self.start, self.end, self.signal]])

    def len(self):
        return self.end - self.start

    def get_chrom(self):
        return self.chrom

    def get_start(self):
        return self.start

    def get_end(self):
        return self.end

    def get_signal(self):
        return self.signal

    def set_signal(self, signal):
        self.signal = signal


def file_gzipped(f):
    # TODO: use the header, not the suffix.
    if re.search('.gz$', f):
        return True
    else:
        return False


def open_bedgraph(bedgraph_file):
    if file_gzipped(bedgraph_file):
        return gzip.open(bedgraph_file, 'rb')
    else:
        return open(bedgraph_file, 'r')


def parse_arguments():
    parser = argparse.ArgumentParser(prog='python normalize_bedgraph.py',
    description="""
    Normalize the signal in a bedgraph file.
    """
    )

    parser.add_argument('bedgraph', action='store', help='Path to the bedgraph file.')
    parser.add_argument('--to-mean-signal', action='store', metavar='TARGET_SIGNAL (float)', type=float, dest='to_mean_signal', help='Normalize the bedgraph such that the average signal per nucleotide is equal to this value.')
    parser.add_argument('--to-number-reads', action='store', metavar='NUMBER_READS (int)', type=int, dest='to_number_reads', help='Normalize the bedgraph to a given number of reads. This is useful if you are trying to keep the signal scale meaningful. Assumes that a signal value of 1.0 corresponds to a single read. To do this, the read length must be known. This can be passed manually (--read-length), or inferred from the bedgraph. Ignored if --to-mean-signal is passed.')
    parser.add_argument('--read-length', action='store', type=int, metavar='READ_LENGTH (int)', dest='read_length', help='If normalizing to a given library size, assume this read length (otherwise, inferred from the bedgraph file). Ignored if --to-mean-signal is passed.')

    return parser.parse_args()


def parse_line(line):
    line = line.rstrip().split('\t')
    chrom = line[0]
    start = int(line[1])
    end = int(line[2])
    signal = float(line[3])
    element = Bedgraph_element(chrom, start, end, signal)
    return element


def infer_read_length(bedgraph_file):
    # infer the read length
    # will do this by looking for stretches in the bedgraph file like this:
    #
    # chr13   20238340        20238589        0.00000
    # chr13   20238589        20238789        1.00000
    # chr13   20238789        20239638        0.00000
    #
    # i.e., a line with signal 0, then a line with non-zero signal,
    # followed by another line with signal 0. The end-start on the non-zero
    # signal line should be the read length

    potential_read_lengths = {}  # read_length --> number of above cases implying this read length
    read_length_queue = collections.deque()

    with open_bedgraph(bedgraph_file) as b:
        for line in b:
            read_length_queue.append(parse_line(line))
            if not len(read_length_queue) == 3:
                continue
            else:
                # look for 0, non-zero, 0
                if read_length_queue[0].get_signal() == 0.0 and read_length_queue[2].get_signal() == 0.0 and read_length_queue[1].get_signal() > 0:
                    potential_length = read_length_queue[1].len()
                    if potential_length not in potential_read_lengths:
                        potential_read_lengths[potential_length] = 0
                    potential_read_lengths[potential_length] += 1

                read_length_queue.popleft()

        if len(potential_read_lengths) == 0:
            sys.stderr.write('Cannot infer read length (try passing this information with --read-length); exiting.\n')
            sys.exit()

        # now find the read length with the most support
        max_support = max(potential_read_lengths.values())
        most_likely_read_lengths = []
        for read_length, support in potential_read_lengths.items():
            if support == max_support:
                most_likely_read_lengths.append(read_length)

        if len(most_likely_read_lengths) > 1:
            # 2+ read lengths are equally likely
            sys.stderr.write('Cannot infer read length (try passing this information with --read-length); exiting.\n')
            sys.exit()
        else:
            read_length = most_likely_read_lengths[0]
            sys.stderr.write('Read length inferred to be {}\n'.format(read_length))
            return read_length


if __name__ == '__main__':

    args = parse_arguments()

    if not (args.to_mean_signal or args.to_number_reads):
        sys.stderr.write('Error: either --to-mean-signal or --to-number-reads must be passed.\n')
        sys.exit()

    total_signal = 0
    bases = 0

    with open_bedgraph(args.bedgraph) as b:
        for line in b:
            element = parse_line(line)
            total_signal += (element.len() * element.get_signal())
            bases += element.len()

    sys.stderr.write('Current average signal per base: {}\n'.format(float(total_signal) / bases))

    normalization_factor = None

    if args.to_mean_signal:
        normalization_factor = args.to_mean_signal / (float(total_signal) / bases)
    elif args.to_number_reads:
        read_length = args.read_length or infer_read_length(args.bedgraph)
        number_reads = round(total_signal / read_length)
        sys.stderr.write('Number of reads detected: {}\n'.format(number_reads))
        normalization_factor = float(args.to_number_reads) / number_reads

    sys.stderr.write('Normalization factor: {}\n'.format(normalization_factor))

    with open_bedgraph(args.bedgraph) as b:
        for line in b:
            element = parse_line(line)
            element.set_signal(element.get_signal() * normalization_factor)
            print(element)
