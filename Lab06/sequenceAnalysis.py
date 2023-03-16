import sys

class NucParams:

    rnaCodonTable = {
        # RNA codon table
        # U
        'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
        'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
        'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
        'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
        # C
        'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
        'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
        'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
        'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
        # A
        'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
        'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
        'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
        'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
        # G
        'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
        'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
        'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
        'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'   # GxG
    }
    dnaCodonTable = {key.replace('U', 'T'): value for key, value in rnaCodonTable.items()}

    '''Given accepted nucleotide set {ACGUTN}'''
    accepted = ['A', 'U', 'N', 'C', 'G', 'T']

    def __init__(self, sequence=""):
        '''initialize dictionaries for updates later in addSequence and calls from other fuctions'''
        self.nucCompDict = {'A': 0, 'T': 0, 'G': 0, "U": 0, 'C': 0, 'N': 0}
        self.aaCompDict = {aa: 0 for aa in NucParams.rnaCodonTable.values()}
        self.codonCompDict = {codon: 0 for codon in NucParams.rnaCodonTable.keys()}
        '''call addSequence if user overrides optional input for sequence'''
        self.addSequence(sequence)

    def addSequence(self, thisSequence):

        '''get rid of any white spaces and capitalize input'''
        l = ''.join(thisSequence).split()
        rawSeq = ''.join(l).upper()

        '''initialize codon array to collect codons found in sequence'''
        codonsFound = []

        #'''clean up any inputs that aren't accepted'''
        #for char in rawSeq:
        #    if char not in NucParams.accepted:
        #        rawSeq = rawSeq.replace(char, "")

        '''get nucleotide count for all valid nucleotides in {AGCTUN}'''
        for nuc in rawSeq:
            if nuc in NucParams.accepted:
                self.nucCompDict[nuc] += 1

        '''change input sequence to RNA'''
        rawSeq = rawSeq.replace("T", "U")

        '''get all codons in sequence'''
        for i in range(0, len(rawSeq), 3):
            codonsFound.append(rawSeq[i:i+3])

        '''remove any invalid codons'''
        for codon in codonsFound:
            if 'N' in codon:
                codonsFound.remove(codon)

        '''find aa count by translating codon'''
        for codon in codonsFound:
            if codon in NucParams.rnaCodonTable.keys():
                self.aaCompDict[NucParams.rnaCodonTable[codon]] += 1

        '''find codon count for composition'''
        for codon in codonsFound:
            if codon in NucParams.rnaCodonTable.keys():
                self.codonCompDict[codon] += 1

    def aaComposition(self):
        return self.aaCompDict

    def nucComposition(self):
        return self.nucCompDict

    def codonComposition(self):
        return self.codonCompDict

    def nucCount(self):
        return sum(self.nucCompDict.values())

class FastAreader:
    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen(self):
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):

        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence