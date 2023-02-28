#!/usr/bin/env python3
# Name: Your full name (CATS account username)
# Group Members: List full names (CATS usernames) or “None”

from sequenceAnalysis import Orf_Finder, FastAreader, ProteinParam, NucParams
########################################################################
# CommandLine
########################################################################
class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output'
                                             )
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices= (100,200,300,500,1000), default=100, action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', default = ['ATG'],nargs='?', 
                                 help='start Codon') #allows multiple list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = ['TAG','TGA','TAA'],nargs='?', help='stop Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

########################################################################
# Main
# Here is the main program
# 
#
########################################################################
   
def main(MP=None):
    '''
    Reads in a fasta file and outputs the ORFs frame, start, stop, 
    and length position on a output file.
    '''
    if MP is None:
        main_program = CommandLine()
        if main_program.args.longestGene:
            areader_file = FastAreader()
            for header, sequence in areader_file.readFasta():
                print(header)
                orfData = Orf_Finder(sequence)
                orfData.locate_ORFs()
                orfData.locate_reverse_ORFs()
                fil_list = filter(lambda orf: orf[3] > main_program.args.minGene, orfData.orfs)
                for reading_frame, start, stop, size, in sorted(fil_list, key=lambda orf:orf[3], reverse = True):
                    print('{:+d} {:>5d}..{:>5d} {:>5d}'.format(reading_frame, start, stop, size))
            else:
                main_program = CommandLine(MP)
            print(main_program.args)
    ###### replace the code between comments.
    # thisCommandLine.args.longestGene is True if only the longest Gene is desired
    # thisCommandLine.args.start is a list of start codons
    # thisCommandLine.args.stop is a list of stop codons
    # thisCommandLine.args.minGene is the minimum Gene length to include
    #
    #######
    
if __name__ == "__main__":
    main()