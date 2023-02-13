import sequenceAnalysis as seqA

def main (fileName=None):
    myReader = seqA.FastAreader(fileName) 
    myNuc = seqA.NucParams()
    
    for head, seq in myReader.readFasta() :
        codonTable, aaComp, nucComp, codonComp, nucCount = myNuc.addSequence(seq.replace('T','U'))
    
    #Calculates sequence length in Mb
    print(f'Sequence Length = {nucCount/10**6:.2f} Mb')
    print()
    
    #Gives percent GC content, rounded to nearest decimal
    gcContent = ((nucComp['G'] + nucComp['C'])/nucCount)*100
    print(f'GC content = {gcContent:.1f}%')
    print()
    
    #Calculates relative abundance of each codon
    for x in codonTable:
        pass
    for key, value in codonTable.items():
        print(f'{key} : {value} {codonComp[key]/aaComp[codonTable[key]]:.2f} ({codonComp[key]})')

if __name__ == "__main__":
    main('testGenome.fa') # make sure to change this in order to use stdin