import Other.sequenceAnalysiss as sequenceAnalysiss
import sys

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

class OrfFinder():

    def __init__(self, seq):
        self.seq = seq

    def findOrfs(self):
        '''
        find Orfs on top strand and return list of Orfs
        remember to handle the dangling start and stop cases
        '''
        orfsList = [[], [], []]
        startList = []
        stopList = []

        startCodonList = ["ATG"]
        stopCodonList = ["TAG", "TAA", "TGA"]

        for frame in range(3):
            stopList.clear()
            for pos in range(frame, len(self.seq), 3):
                codon = self.seq[pos: pos+3]

                if codon in startCodonList:
                    startList.append(pos)

                if codon in stopCodonList:
                    stopList.append(pos)

                    if startList:

                        if not orfsList[frame] and len(stopList) == 1:
                            if pos + 3 > 0:
                                orf = (frame, 0, pos + 3, pos + 3)
                                orfsList[frame].append(orf)
                            startList.clear()
                        else:
                            endPos = pos + 3
                            startPos = startList[0]
                            length = endPos - startPos
                            if length > 0:
                                orf = (frame, startPos, endPos, length)
                                orfsList[frame].append(orf)
                            startList.clear()

                    if not orfsList[frame] and len(stopList) == 1:
                        if pos+3 > 0:
                            orf = (frame, 0, pos+3, pos+3)
                            orfsList[frame].append(orf)

                if pos == len(self.seq) - 4:
                    if startList:
                        startPos = startList[0]
                        endPos = len(self.seq) - 1
                        length = endPos - startPos
                        if length > 0:
                            orf = (frame, startPos, endPos, length)
                            orfsList[frame].append(orf)
        return orfsList


    def findRevOrfs(self):
        '''
        find Orfs on the bottom strand and return that list of Orfs
        remember to fixup the orfList so that it r	efers to top strand coordinates and the rev frames
        '''
        compDict = {"A": "T", "G": "C", "C": "G", "T": "A"}

        revSeq = list(self.seq)
        revSeq = reversed([compDict.get(base, base) for base in revSeq])
        revSeq = ''.join(revSeq)

        self.seq = revSeq
        return self.findOrfs()


    def saveOrf(self, start, stop, length, frame):
        pass

def main():
    fastaReader = FastAreader("test.fa")
    orfsFound = []

    for header, sequence in fastaReader.readFasta():
        print(header)

        finder = OrfFinder(sequence)
        orfList = finder.findOrfs()
        revOrfList = finder.findRevOrfs()

        for list in orfList:
            for tuple in list:
                frame = tuple[0] + 1
                startPos = tuple[1] + 1
                endPos = tuple[2]
                print(frame, startPos, endPos, tuple[3])

        for list in revOrfList:
            for tuple in list:
                frame = tuple[0] + 1
                startPos = tuple[1] + 1
                endPos = tuple[2]
                print(-frame, len(sequence)-endPos+1,
                                  len(sequence)-startPos+1, tuple[3])

    #orfsFound.sort(key=lambda tup: (tup[3], tup[1]), reverse=True)

    #for orf in orfsFound:
    #    print('{:+d} {:>5d}..{:>5d} {:>5d}'.format(orf[0], orf[1], orf[2], orf[3]))

main()