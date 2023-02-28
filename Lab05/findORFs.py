#!/usr/bin/env python3
# Name: Your full name (CATS account username)
# Group Members: List full names (CATS usernames) or “None”

from sequenceAnalysis import Orf_Finder, FastAreader, ProteinParam, NucParams
########################################################################
# CommandLine
########################################################################
class CommandLine():
    '''
    Handle the command line, usage and help requests.
    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.
    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
    '''

    def __init__(self, inOpts=None):
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(
            description='Program prolog - a brief description of what this thing does',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
            )
        self.parser.add_argument('inFile', action='store', help='input file name')
        self.parser.add_argument('outFile', action='store', help='output file name')
        self.parser.add_argument('-lG', '--longestGene', action='store', nargs='?', const=True, default=True,
                                 help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices=range(0, 2000), action='store',
                                 help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action='append', nargs='?',
                                 help='start Codon')  # allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


########################################################################
# Main
# Here is the main program
# 
#
########################################################################
from sequenceAnalysis import NucParams, ProteinParam, FastAreader, Orf_Finder
import sys
def output(infile, outfile):

    fastaReader = FastAreader(infile)
    orfsFound = []

    f = open(outfile, 'w')
    sys.stdout = f

    for header, sequence in fastaReader.readFasta():
        print(header)

        '''find forward and reverse comp ORFs in sequence'''
        finder = Orf_Finder(sequence)
        orfList = finder.locate_ORFs()
        revOrfList = finder.locate_reverse_ORFs()

        '''convert orfList to clear tuples in single list for sorting'''
        for list in orfList:
            for tuple in list:
                frame = tuple[0] + 1
                startPos = tuple[1] + 1
                endPos = tuple[2]
                orfsFound.append((frame, startPos, endPos, tuple[3]))

        '''covert reverse complementary ORFs to top strand coordinates and add to single list'''
        for list in revOrfList:
            for tuple in list:
                frame = tuple[0] + 1
                startPos = tuple[1] + 1
                endPos = tuple[2]
                orfsFound.append((-frame, len(sequence)-endPos+1,
                                  len(sequence)-startPos+1, tuple[3]))

    '''sort by largest to smallest and by start position furthest to the left'''
    orfsFound.sort(key=lambda tup: (tup[3], tup[1]), reverse=True)

    '''print out list to file in proper format'''
    for orf in orfsFound:
        print('{:+d} {:>5d}..{:>5d} {:>5d}'.format(orf[0], orf[1], orf[2], orf[3]))

    f.close()

########################################################################
#
# Main
#
########################################################################

def main(myCommandLine=None):
    '''
    Implements the Usage exception handler that can be raised from anywhere in process.
    '''
    if myCommandLine is None:
        myCommandLine = CommandLine(['tass2.fa',
                                     'tass2ORFdata-ATG-100.txt',
                                     '--longestGene',
                                     '--start=ATG',
                                     '--minGene=100'])

        output(myCommandLine.args.inFile, myCommandLine.args.outFile)
    else:
        myCommandLine = CommandLine(myCommandLine)
        #output(myCommandLine.args.inFile, myCommandLine.args.outFile)
        ###### replace the code between comments.
        # myCommandLine.args.inFile has the input file name
        # myCommandLine.args.outFile has the output file name
        # myCommandLine.args.longestGene is True if only the longest Gene is desired
        # myCommandLine.args.start is a list of start codons
        # myCommandLine.args.minGene is the minimum Gene length to include
        #
        myReader = FastAreader(myCommandLine.args.inFile)

#######

if __name__ == "__main__":
    main()