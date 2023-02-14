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

class ProteinParam :
    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158, #stored molecular weight of amino acids
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }
    '''
    Initializes the objects used later on as well as comparing input to dictionary to ensure all charecters used is valid
    '''
    def __init__ (self, protein):
        self.mwH2O = 18.015
        self.protein = protein.upper() #uppercase inputed amino acid
        self.count = {}    #counts input and stores it
        for key,value in self.aa2mw.items(): #iterates every letter and stores key values
            self.count[key] = protein.count(key)            
    '''Counts every amino acid within the inputted value'''      
    def aaCount (self):
        return sum(self.count.values()) #returns the sum of the values within given input
    
    '''
        Function computes the optimal pH of the amino acid using binary search. Mid is the middle of the pH times 100
        due to the fact we need 2 decimal points. Add is the number that is subtracted from the middle value that creates
        a second mid point after the first search is done (quarters). 
        '''
    def pI_Binary(self, mid = 700, add = 350):
        if (abs(self._charge_(0)) < abs(self._charge_(0.01))): #if the pH is below 0.01 then it will be 0
            return 0
        elif (abs(self._charge_(14)) < abs(self._charge_(13.99))): #if the pH is bigger than 13.99 it is 14
            return 14
        elif (abs(self._charge_(mid / 100.0)) > abs(self._charge_((mid + 1) / 100.0))): #if charge is greater than greater, good, but if lower than 
            return self.pI((int)(mid + add),(int)(add / 2)) #it will trigger halfing the top half
        elif (abs(self._charge_(mid/100.0)) > abs(self._charge_((mid - 1) / 100.0))): #if charge is greater than lower, good, but if higher then
            return self.pI((int)(mid - add),(int)(add / 2)) #it will trigger the bottom half 
        else:
            return (mid / 100) #if nothing, then we settle with exactly 7 pH

        '''
        Function looks for optimal pH but instead iterating through every pH value linearly. It starts at initialized pH of 0.01. 
        It then scales each pH by 0.01. If pH of the charge is greater than the charge then the previous is correct, it subtracts by 0.01
        then stores the real value within charge in which the value is returned. 
        '''
    def pI (self):
        ph = 0.01 #initialize the variable to 0.01
        charge = self._charge_(0) #find the charge at pH 0 
        while ph <= 14:
            if abs(self._charge_(ph)) > abs(charge): #if charge of pH is greater than charge at 0 then the previous was optimum
                return ph - 0.01 #this subtracts to previous pH as it was optimum and sotres value
            charge = self._charge_(ph)  #stores the latest charge 
            ph += 0.01 #used to specify 2 decimal points
        return 14
    '''Calculates the composition of the inputted amino acid sequence '''
    def aaComposition (self) :
        return self.count #counts object inputed and returns aaComposition

    '''Using an equation to find the net charge at specfic pH value. Using the pH and current pH value, we can calculate the
    total charge based on previous methods using following equation.'''
    def _charge_ (self, ph):
        aaNterm = 9.69 #storing values within variable
        aaCterm = 2.34 #storing values within variable
        aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6} #storing values within variable
        aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10} #storing values within variable
        first = 0 #storing values within variable which will be used later within the equation to find charge 
        for key, value in aa2chargePos.items():
            first += self.count[key]*(10**aa2chargePos[key] / (10**aa2chargePos[key]+10**ph)) #calculates the first value
        first += (10**aaNterm / (10**aaNterm + 10**ph))
        second = 0 #initialize value of variable 
        for key,value in aa2chargeNeg.items():
            second += self.count[key]*(10**ph / (10**aa2chargeNeg[key] + 10**ph)) #calculates the second value
        second += (10**ph / (10**aaCterm + 10**ph))
        return first - second #difference between the two values provides the value for the charge 

    '''Based on amino acid composition we can calculate how much light it absorbs. We do this by using
    standard measures and using the count of the input and multiple it by the molarweight of the standards.'''
    def molarExtinction (self):
        aa2abs280= {'Y':1490, 'W': 5500, 'C': 125} #creating dictionary with specifc molarweight to calculate later on
        molar = 0 #initialized molar 
        for key, value in aa2abs280.items(): 
            molar += self.count[key]*aa2abs280[key] 
        return molar
    '''We find the mass provided based on the molecular weight and molar extinction devided amongst each other. This 
    is also based on molecular weight on light wavelength values calculated previous within molarExtinction.'''

    def massExtinction (self):
        molecular_weight =  self.molecularWeight() #stores molecular weight into variable 
        return self.molarExtinction() / molecular_weight if molecular_weight else 0.0 #return molar extension devided by molecular weight 

    '''molecular weight is calculated using an equation using the composition of amino acid and multiplying that
    by the standard weight provided by the dictionary and subtracting that by water. This is based on the equation 
    provided by the professor.'''
    def molecularWeight (self):
        total = self.mwH2O #store molecular weight of H2O in variable 
        for key,value in self.aa2mw.items():  #iterate through dictionary for every key and value of aa2mw
            total += self.aa2mw[key] * self.count[key] - self.mwH2O * self.count[key] #sum calculation of molecular weight minus weight of water
        return total

class NucParams:
    def __init__ (self, inString=''): 
        self.rnaCodonTable = {
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
        self.dnaCodonTable = {key.replace('U','T'):value for key, value in self.rnaCodonTable.items()}
        self.codonComposition = {codon:0 for codon in self.rnaCodonTable}
        self.aaComposition = {aa:0 for aa in self.rnaCodonTable.values()}
        self.nucComposition = {nuc:0 for nuc in 'ACGTNU'}

    def addSequence (self, inSeq): #inSeq for sequence
        for index in range (0, len(inSeq), 3): #start, seqnece length read(end), and final is step (increments of 3)
            codon = inSeq[index:index+3].upper()
            for i in codon:
                if i in self.nucComposition:
                    codon = {}
            if codon in self.rnaCodonTable:
                self.aaComp += 1
                self.codonComp += 1

        pass
    def aaComposition(self):
        return self.aaComp
    def nucComposition(self):
        return self.nucComp
    def codonComposition(self):
        return self.codonComp
    def nucCount(self):
        return sum(self.nucComp.values())
  