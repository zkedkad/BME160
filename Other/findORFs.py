#!/usr/bin/env python3
# Name: Ziad Kedkad (zkedkad)
# Group Members: Qudsi Aljabiri, Justin Ngyuen, Aurko Mahesh, Emily Cheng

import Other.sequenceAnalysiss as sequenceAnalysiss
import sys

'''
    MAIN PSUEDOCODE
    orfFinder():
        for each frame (1,2,3):
            **remember to clear stop list from previous frame to search for new dangling stop cases**
            go through in positions in iterations of 3
                find codon at that position
                if start : add to startList
                if stop : add to stopList
                          if there is a start codon in start list:
                                find length from start to stop
                                if that meet minLength requirement
                                then add to orfList
                                **clear the startList if we are only looking for longest orf
                                we don't need the rest of the ones we found**
                          if we havent found any ORFs yet and this is our first stop
                          then we could potentially have a dangling orf:
                                check length from start of sequence to stop
                                if it meets the minLength requirement add to orfList
                if we reach the last position in sequence:
                    check if our startList still has starts
                        if it does then have a dangling orf:
                            check if length of that from start to end of sequence
                            meets the minLength requirement and if it does add to orfList
        return orfList
'''

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


class OrfFinder():

    orfsList = [[], [], []]
    startList = []
    stopList = []

    startCodonList = ["ATG"]
    stopCodonList = ["TAG", "TAA", "TGA"]

    def __init__(self, seq):
        self.seq = seq

    def findOrfs(self):

        for frame in range(3):
            ''' for each frame find the codon positions based on length of sequence'''
            self.stopList.clear()
            for pos in range(frame, len(self.seq), 3):
                ''' after finding codon append it to the right list if start or stop'''
                codon = self.seq[pos: pos+3]

                if codon in self.startCodonList:
                    self.startList.append(pos)

                if codon in self.stopCodonList:
                    self.stopList.append(pos)
                    '''if stop found after stop, check length and if greater than minLength add to orfList'''
                    if self.startList:
                        '''assume that first orf found is a dangling stop case'''
                        if not self.orfsList[frame] and len(self.stopList) == 1:
                            if pos+3 > 100:
                                self.saveOrf(0, pos+3, pos+3, frame)
                            self.startList.clear()
                        else:
                            length = (pos + 3) - self.startList[0]
                            if length > 100:
                                self.saveOrf(self.startList[0], (pos+3), length,frame)
                            '''clear the startList after if we are only looking for the greatest length ORF'''
                            self.startList.clear()

                    ''' if not ORFS were found and this is the first stop then its a dangling ORF'''
                    if not self.orfsList[frame] and len(self.stopList) == 1:
                        '''check length of that ORF add it to orfList if its meet the minLength requirement'''
                        if pos+3 > 100:
                            self.saveOrf(0, pos+3, pos+3, frame)

                '''we reached the end of the sequence'''
                if pos == len(self.seq) - 4:
                    '''if there is still start codons in the startList but no stops then we have a dangling ORF'''
                    if self.startList:
                        length = (len(self.seq) - 1) - self.startList[0]
                        ''' if it meets the minLength requirement add it to orfList '''
                        if length > 100:
                            self.saveOrf(self.startList[0], (len(self.seq) - 1), length, frame)
        return self.orfsList


    def findRevOrfs(self):

        '''clear ORFS we found in the forward reading frame '''
        self.orfsList = [[], [], []]
        self.startList.clear()
        self.stopList.clear()

        compDict = {"A": "T", "G": "C", "C": "G", "T": "A"}

        '''find the reverse complement of that sequence'''
        revSeq = list(self.seq)
        revSeq = reversed([compDict.get(base, base) for base in revSeq])
        revSeq = ''.join(revSeq)

        '''set it as our initial seq and just call findOrfs because we are lazy'''
        self.seq = revSeq
        return self.findOrfs()


    def saveOrf(self, start, stop, length, frame):
        ''' save the orf frame to orfList '''
        orf = (frame, start, stop, length)
        self.orfsList[frame].append(orf)

def output(infile, outfile):

    fastaReader = sequenceAnalysiss.FastAreader(infile)
    orfsFound = []

    f = open(outfile, 'w')
    sys.stdout = f

    for header, sequence in fastaReader.readFasta():
        print(header)

        '''find forward and reverse comp ORFs in sequence'''
        finder = OrfFinder(sequence)
        orfList = finder.findOrfs()
        revOrfList = finder.findRevOrfs()

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
        import Other.sequenceAnalysiss as sequenceAnalysiss
        myReader = sequenceAnalysiss.fastaReader(myCommandLine.args.inFile)

#######

if __name__ == "__main__":
    main()