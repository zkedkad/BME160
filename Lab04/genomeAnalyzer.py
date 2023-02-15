
from sequenceAnalysis import NucParams, ProteinParam, FastAreader
def main (fileName=None):
    myReader = FastAreader(fileName) 
    myNuc = NucParams()
    for head, seq in myReader.readFasta():
        myNuc.addSequence(seq)

    print('sequence length = {:0.2f} Mb \n'.format(myNuc.nucCount() / (1000000))) #print sequence length grabbed from function in sequence analysis and divide to make it into correct unit for Mb    
    GC_content = (myNuc.nucComp['G'] + myNuc.nucComp['C']) / (myNuc.nucCount())
    GC_percentage = GC_content * 100
    print("GC content = {:0.1f} \n".format(GC_percentage))

    emptylist = []

    for nuc in myNuc.codonComp:
        thisCodonComp = myNuc.rnaCodonTable[nuc]
        aa = myNuc.aaComp[myNuc.rnaCodonTable[nuc]]
        val = myNuc.codonComp[nuc] / aa
        emptylist.append([nuc, thisCodonComp, val * 100, myNuc.codonComp[nuc]]) #Add nucleotide, amino acid, percentage, and total count
        emptylist.sort(key = lambda x: (x[1], x[0]))
    
    for i in range(len(emptylist)):
        print('{:s} : {:s} {:5.1f} ({:6d})'.format(emptylist[i][0], emptylist[i][1], emptylist[i][2], emptylist[i][3])) #format the list of list by its place in the list
    # calculate relative codon usage for each codon and print
    # for nucI in nuc:
        # print ('{:s} : {:s} {:5.1f} ({:6d})'.format(nuc, aa, val*100, thisCodonComp[nuc]))

if __name__ == "__main__":
    main('testGenome.fa') # make sure to change this in order to use stdin