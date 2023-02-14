import sys
class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence



class NucParams:

    def __init__ (self, inString=''):
        
        '''
        Creates a sorted RNA Codon table, which is utilized as a baseline for intiailizing various
        default dictionaries. These are then modified in later methods.
        '''
        
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
        'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
        }
        
        
        
        sortedAAList = sorted(rnaCodonTable.values()) #Sorts the AAs from rnaCodonTable
        codonList = [] #List of codons from rnaCodonTable
        aaList = [] #List of AAs from rnaCodonTable
        for key, value in rnaCodonTable.items():
            codonList.append(key)
            aaList.append(value)
        
        '''
        Generates dictionary of codons sorted alphabetically by their AA values
        '''
        sortByAA = {}
        for x in sortedAAList: # x represents an AA from the sorted list (First is '-' )
            if x in aaList: #Checks to ensure that x is present in rnaCodonTable's values
                i = aaList.index(x) #Takes index of AA in the aaList
                aa = aaList[i] 
                codon = codonList[i]
                sortByAA[codon] = aa #Adds an entry into the dictionary with a codon key and AA value.
                del codonList[i], aaList[i] #Deletes the entered key, value from their respective lists to avoid double dipping (Multiple keys share the same value)
        
        '''
        Sorts the shared codons for each AA
        '''
        sortedAAList = sorted(set(sortedAAList)) #Removes duplicates from sortedAAList
        self.rnaCodonTable = {} #The final fully sorted RNA codon table
        for x in sortedAAList:
            temp = [] #Temporarily stores the all the codons associated with a specific AA
            for key in sortByAA.keys(): #If a key-value pair has the value X, the key is appended to temp
                if sortByAA[key] == x:
                    temp.append(key)
            temp = sorted(temp) #Alphabetically sorts all the codons associated with a specific AA
            for y in temp: #Assigns each codon in temp to their shared key, and adds them in order to the dictionary. 
                self.rnaCodonTable[y] = x

        self.rnaCodonComp = {} 
        for key in self.rnaCodonTable: #Sets default RNA codon count. Values = 0
            self.rnaCodonComp[key] = 0
        
        #Codon table for DNA. Replaces U for T in values from the RNA Codon Table
        self.dnaCodonTable = {key.replace('U','T'):value for key, value in self.rnaCodonTable.items()}
        
        self.dnaCodonComp = {} 
        for key in self.dnaCodonTable: #Sets default DNA codon count. Values = 0
            self.dnaCodonComp[key] = 0
        
        self.aaComp = {} #Default AA count. Values = 0
        for aa in self.dnaCodonTable.values():
            self.aaComp[aa] = 0
        
        self.dnaBases = 'ATCGN' #Possible DNA bases
        self.dnaNucDict = {} 
        for base in self.dnaBases: #Sets default DNA base counts. Values = 0
            self.dnaNucDict[base] = 0
        
        self.rnaBases = 'AUCGN' #Possible RNA bases
        self.rnaNucDict = {} 
        for base in self.rnaBases: #Sets default RNA base counts. Values = 0
            self.rnaNucDict[base] = 0
            
        self.baseCount = 0 #Utilized to count nucleotides
    
    def addSequence (self, inSeq):
        
        '''
        Takes codons from each inputted sequence and feeds them into various methods along with 
        relavent information based on whether sequence is DNA or RNA. 
        '''
        for x in range(0,len(inSeq),3): 
            codon = inSeq[x:x+3]
            
            if len(codon) != 3: #If a codon has a length < 3, prints error and exits program
                print('Error')
                break
        
            elif 'T' in inSeq: #In the event sequence is DNA:
                self.aaComposition(self.dnaCodonTable, codon) #aaComposition utilizes dnaCodonTable to determine AA of inputted codon
                self.nucComposition(self.dnaNucDict, codon) #nucComposition breaks the codon down into DNA nucleotides
                self.codonComposition(self.dnaCodonComp, codon) #codonComposition attempts to match codon with DNA codons
                self.nucCount(self.dnaNucDict) #Calls nucCount
                
            elif 'U' in inSeq: #In the event sequence is RNA:
                self.aaComposition(self.rnaCodonTable, codon) #aaComposition utilizes rnaCodonTable to determine AA of inputted codon
                self.nucComposition(self.rnaNucDict, codon) #nucComposition breaks the codon down into RNA nucleotides
                self.codonComposition(self.rnaCodonComp, codon) #codonComposition attempts to match codon with RNA codons
                self.nucCount(self.rnaNucDict)
        
        '''
        Returns proper sorted codon table, aa composition, nuc composition, codon composition, 
        and base count depending on type of sequence entered(DNA or RNA)
        
        '''
        if 'T' in inSeq: #When DNA sequences entered:
            return self.dnaCodonTable, self.aaComp, self.dnaNucDict, self.dnaCodonComp, self.baseCount
        elif 'U' in inSeq: #When RNA sequences entered:
            return self.rnaCodonTable, self.aaComp, self.rnaNucDict, self.rnaCodonComp, self.baseCount

    def aaComposition(self, codonTable, codon): #codonTable depends on whether seq is DNA or RNA. Refer to elif statements in addSequence()
        '''
        Takes codon from sequence, sees what AA it corresponds to using the relevant codonTable
        Adds 1 to that AA's count 
        Ignores sequences with N, they are not present in either rnaCodonTable or dnaCodonTable
        '''
        if codon in codonTable:
            aa = codonTable[codon]
            self.aaComp[aa] += 1
            
    def nucComposition(self, baseDict, codon): #baseDict depends on whether seq is DNA or RNA. Refer to elif statements in addSequence()
        '''
        Takes codon from sequence and iterates through it. 
        For each base that appears, 1 is added to that base's count
        '''
        for x in codon:
            if x in baseDict:
                baseDict[x] += 1
                
        
    def codonComposition(self, codonComp, codon):
        '''
        Takes codon and searches the relevant codonComp dictionary to see if it it present.
        If so, adds 1 to the codon's count
        Ignores sequences with N, since they are not present in either rnaCodonTable or dnaCodonTable
        '''
        if codon in codonComp:
            codonComp[codon] += 1
            
    def nucCount(self, bases):
        '''
        Takes the sum of all the bases calculated using nucComposition
        '''
        self.baseCount = sum(bases.values())
    

    
class ProteinParam:

    def __init__ (self, protein):
        #Dictionary of valid amino acids and their associated molecular weights (provided by professor)
        self.a2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

        '''Dictionary containing the valid AA's in the input string, and the number of times they appear (including zero)'''
        
        #Creates a "default" dictionary where every valid amino acid has a value of zero
        aaComposition = {}
        for x in self.a2mw:
            aaComposition[x] = 0       
        #Goes through the input string. If letter is a valid one, adds 1 to said letter's value. 
        protein = protein.upper()
        for x in protein:
            if x in aaComposition:
                aaComposition[x] += 1
        self.proteinComp = aaComposition 

    def aaCount (self):
        '''
        Returns the sum of the values in self.proteinComp(), which contains every AA along with the number of times
        they appear in the string
        
        '''
        return sum(self.proteinComp.values())

    def pI (self, min = 0, max = 14, pH = 7, precision = 2):
        '''
        Utilizes a binary search to calculate a pH between 0 and 14 that will result in a net charge closest to zero. 
        Calls upon the _charge_() function to calculate individual charges.
        Function has default parameters min = 0, max = 14, and pH = 7 to symbolize the range of calculations and the 
        original midpoint. This is done to avoid the need for an input array.
        
        (Note: Results calculated to the nearest 0.01 decimal)
        
        '''
        
        #Mid is set to the charge at the current pH (which is initially the pH midpoint of 7)
        #Upper is the charge at midpoint plus 10 to the negative precision (For default of 2, this would be 0.01)
        #Lower is the charge at midpoint minus 10 to the negative precision (For default of 2, this would be 0.01)
        
        mid = abs(self._charge_(pH))
        upper = abs(self._charge_(pH + 10**(-1*precision)))
        lower = abs(self._charge_(pH - 10**(-1*precision)))
        
        '''Compares the charge of Mid with Upper and Lower. '''
        #Returns current pH if its charge is lower than both the charge of the higher and lower pH (the isoelectric point)
        
        #If upper has a lower charge than current pH, we assume the desired point is in the upper half of data.
            #The minimum range is set to the current pH (initially 7)
            #A new pH is calculated (midpoint between the new min and the old max)
            #The function calls itself. Arguments: The new min, old max, and new pH (Initially 7, 14, and 10.5 respectively)
        
        #If lower has a lower charge than current pH, we assume desired point is in the lower half of data.
            #The maximum range is set to the current pH (initially 7)
            #A new pH is calculated (midpoint between the old min and new max)
            #The function calls itself. Arguments: The new min, old max, and new pH (Initially 0, 7, and 3.5 respectively)
       
        if mid < upper and mid < lower:
            return pH
        elif upper < mid:
            min = pH
            pH = pH + (max - pH)/2
            return self.pI(min, max, pH)
        elif lower < mid:
            max = pH
            pH = pH - ((pH - min)/2)
            return self.pI(min, max, pH)

    def aaComposition (self) :
        '''
        Returns self.proteinComp(), which contains every AA along with the number of times they appear in the string
        
        '''
        return self.proteinComp

    def _charge_ (self, pH):
        '''
        Calculates the net charge of the protein from given pKa and inputted pH values
        
        '''
        #pKa values provided by Professor
        aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
        aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
        aaNterm = 9.69
        aaCterm = 2.34

        posCharge = 0
        negCharge = 0
        for key, value in self.proteinComp.items():
            #Charge of a positive AA calculated, and multiplied by num of appearances. Total is added to posCharge
            if key in aa2chargePos and value > 0:
                charge = (10**aa2chargePos[key])/(10**aa2chargePos[key] + 10**pH)
                posCharge += (charge*value)
            #Charge of a negative AA calculated, and multiplied by num of appearances. Total is added to negCharge
            if key in aa2chargeNeg and value > 0:
                charge = (10**pH)/(10**aa2chargeNeg[key] + 10**pH)
                negCharge += (charge*value)
       
        #Charge of N-terminus and C-Terminus
        n_term = (10**aaNterm)/(10**aaNterm + 10**pH)
        c_term = (10**pH)/(10**aaCterm + 10**pH)

        #Total charge of the protein
        total_charge = abs((posCharge + n_term) - (negCharge + c_term))
        return total_charge

    def molarExtinction (self):
        '''
        Calculates the molar extinction coefficient of the protein using given absorbance values 
        
        '''
        #Absorbance values at 280nm, provided by Professor
        aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}
        
        #If value in string has absorbance, multiply by the number of times it appears and add it to the extinction coefficient 
        e_coeff = 0
        for key, value in self.proteinComp.items():
            if key in aa2abs280 and value > 0:
                e_coeff += (value*aa2abs280[key])
        return e_coeff

    def massExtinction (self):
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight (self):
        '''
        Takes the sum of all the molecular weights of present AAs. 
        Excludes water released during peptide bond formation. 
        
        '''
        total_mass = 0
        mwH2O = 18.015
        
        #Iterates through proteinComp() and pulls the molecular weight from a2mw utilizing the shared key. 
        #Multiplies the weight by the number of times the AA appears, and adds it to total_mass
        for key, value in self.proteinComp.items():
            total_mass += (value*self.a2mw[key])
        
        #Calculates the amount of water released (1 less than the number of AAs) and multiplies it by the weight of H2O
        #Subtracts the weight of H2O from the total mass, and returns it.
        return total_mass - ((self.aaCount() - 1)*mwH2O)

