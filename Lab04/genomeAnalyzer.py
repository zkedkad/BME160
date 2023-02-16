#!/usr/bin/env python3
# Name: Ziad Kedkad (zkedkad)
# Group Members: Qudsi Aljabiri, Justin Ngyuen, Aurko Mahesh, Emily Cheng
from sequenceAnalysis import NucParams, ProteinParam, FastAreader
def main (fileName=None):
    '''Function initializes all the data and functions from sequenceAnalysis to be used in this program. This also prints out
    and calculates the final data points which are displayed to the user. This program also functions to grab the file
    and fetch the data to be used within sequenceAnalysis.'''
    myReader = FastAreader(fileName)  #fetches FastAreader program and stores in variable
    myNuc = NucParams() #fetches NucParams and stores in variable myNuc
    for head, seq in myReader.readFasta(): 
        myNuc.addSequence(seq)

    print('sequence length = {:0.2f} Mb \n'.format(myNuc.nucCount() / (1000000))) #print sequence length grabbed from function in sequence analysis and divide to make it into correct unit for Mb    
    GC_content = (myNuc.nucComp['G'] + myNuc.nucComp['C']) / (myNuc.nucCount()) #calculates GC composition on all nucleotides within given sequence
    GC_percentage = GC_content * 100 #converts data point to percentage 
    print("GC content = {:0.1f} \n".format(GC_percentage)) #prints GC content

    emptylist = [] #creates empty list to be used to display codon composition data

    for nuc in myNuc.codonComp:
        thisCodonComp = myNuc.rnaCodonTable[nuc] #stores value in variable 
        aa = myNuc.aaComp[myNuc.rnaCodonTable[nuc]] #stores value in variable
        val = myNuc.codonComp[nuc] / aa #stores value in varaibale and divides over amino acids to produce value
        emptylist.append([nuc, thisCodonComp, val * 100, myNuc.codonComp[nuc]]) #Add nucleotide, amino acid, percentage, and total count
        emptylist.sort(key = lambda x: (x[1], x[0])) #lambda function used to sort data for further printing in next for statement
    
    for i in range(len(emptylist)):
        print('{:s} : {:s} {:5.1f} ({:6d})'.format(emptylist[i][0], emptylist[i][1], emptylist[i][2], emptylist[i][3])) #format the list of list by its place in the list

if __name__ == "__main__":
    main('testGenome.fa') # make sure to change this in order to use stdin