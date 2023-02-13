#!/usr/bin/env python3
# Name: Zarghaam Abbas
# Group Members: Andrew Matsiev, Aanad Verna, Samuel Hirsch, Rogelio Chavez

class ProteinParam :
    '''This class with its method is used to calculate
    
    Number of Amino Acids
    
    Molecular Weight
    
    molar Extinction coefficient
    
    mass Extinction coefficient
    
    Theoretical pI
    
    Binary Theoretical pI(For bonus marks)
    
    Amino acid composition 
    '''
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    def __init__ (self, protein):
        
        ''' This method is used for initializing the variables that will be used by other functions
        
        i)self.mwH2O This is the molecular weight of weight that is initialized here
        
        ii)self.aa2mw : this is a dictionary that contain the molecular weight of different amino acids
        
        iii)protein it will store the upper case converted equivalence of the input string
        
        iv)self.count it is a dictionary that will only store number of the valid amino acids from protein string
        
        '''
      #weight of water
        self.mwH2O = 18.015
      #The dictionary containing valid amino acid composition, Copied from above
        self.aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }
      #convert the string to upper case as it maybe lower case too
        self.protein=protein.upper()
      #This will create a dictionary,that will be used to store count of the amino acids
        self.count={}
      #WE WILL Iterate though the whole dictionary aa2mw
        for key,value in self.aa2mw.items():
        #count the times a specific amino acid comes and save it in the key of count dictionary
          self.count[key]=protein.count(key)

    def aaCount (self):
        ''' It will retrun the number of amino acids, 
        It basically retrun the sum of the values of the dictionary defined in self._INIT_'''
      #store the values of count in a variable
        value=self.count.values()
      #return the sum of the values
        return sum(value)
    def pIm(self,mid=700,add=350):
        '''
        THis method basically compute the optimium ph using binary search
        over the range of 1400(0.00-14.00)
        @Param1 mid: It is the mid*100 of ph(middle for binary search), as we need precision of 2 decimal points
        @Param2 add: It is the term that need to be added of subtracted from the mid, which result in the new mid
        It will work in a way that:
        it will see first for ph 7.00 but if 7.01 will have lower charge so, it will have new mid of 1050 
        else if 6.99 will have lower charge so, it will have new mid of 350 then correspondingly every search 
        will decrease the term by 1/2.
        '''
        #chechk if 0 is the optimium ph AS it will increase of decrease constantly
        if(abs(self._charge_(0))<abs(self._charge_(0.01))):
            return 0
        #check for 14 as optimium ph AS it will increase of decrease constantly
        elif (abs(self._charge_(14))<abs(self._charge_(13.99))):
            return 14
        #now if the charge decrease with increase of ph mean optimium point is not reached and we will see the second mid
        elif(abs(self._charge_(mid/100.0))>abs(self._charge_((mid+1)/100.0))):
            return self.pIm((int)(mid+add),(int)(add/2))
        #now if the charge decrease with decrease of ph mean optimium point is not reached and we will see the first mid
        elif (abs(self._charge_(mid/100.0))>abs(self._charge_((mid-1)/100.0))):
            return self.pIm((int)(mid-add),(int)(add/2))
        else:
            #neither increase with ph nor decrease with ph implies we have reached the optimium ph
            return mid/100
    def pI (self):
        '''
        This method also return the optimium ph where the charge gets close to 0, but it instead of 
        binary search iterate thoughout the range of 0.00-14.00
        '''
        #initialize ph to 0.01
        ph=0.01
        #find charge at ph 0
        charge=self._charge_(0)
        while ph<=14:
            #compare the charge with the charge at ph 0.01 less and 
            #if the charge increase then it mean the previous charge was optimium
            #as it work in a linear way
            if abs(self._charge_(ph))>abs(charge):
                return ph-0.01
            #store the latest charge
            charge=self._charge_(ph)
            #increase ph by 0.01 as we are only concerned to 2 decimal points
            ph+=0.01
        return 14

    def aaComposition (self) :
        '''This method return the number of time each valid amino acid has occured in the 
        string Protein
        '''
        #This will return the dictionary that is created in the _INIT_ 
        return self.count

    def _charge_ (self,ph):
        '''
        This method compute the charge on an amino acid at specific ph, 
        @Param1 ph: It is the ph at which this method return charge, it is used by method PI,and PIm
        it will use the formula 
        𝑛𝑒𝑡𝐶ℎ𝑎𝑟𝑔𝑒=[∑𝑎𝑎=(𝐴𝑟𝑔,𝐿𝑦𝑠,𝐻𝑖𝑠,𝑁𝑡𝑒𝑟𝑚𝑁𝑎𝑎 10^𝑝𝐾𝑎(𝑎𝑎)/10^𝑝𝐾𝑎(𝑎𝑎)+10^𝑝𝐻]−[∑𝑎𝑎=(𝐴𝑠𝑝,𝐺𝑙𝑢,𝐶𝑦𝑠,𝑇𝑦𝑟,𝐶𝑡𝑒𝑟𝑚𝑁𝑎𝑎 10𝑝𝐻/10^𝑝𝐾𝑎(𝑎𝑎)+10^𝑝𝐻]
        so, for that it will first calculate 
        [∑𝑎𝑎=(𝐴𝑟𝑔,𝐿𝑦𝑠,𝐻𝑖𝑠,𝑁𝑡𝑒𝑟𝑚𝑁𝑎𝑎 10^𝑝𝐾𝑎(𝑎𝑎)/10^𝑝𝐾𝑎(𝑎𝑎)+10^𝑝𝐻]
        and after that 
        [∑𝑎𝑎=(𝐴𝑠𝑝,𝐺𝑙𝑢,𝐶𝑦𝑠,𝑇𝑦𝑟,𝐶𝑡𝑒𝑟𝑚𝑁𝑎𝑎 10𝑝𝐻/10^𝑝𝐾𝑎(𝑎𝑎)+10^𝑝𝐻]
        then return their differnce
        '''
        #stroe the Nterminal pka
        aaNterm = 9.69
        #stroe the Cterminal pka
        aaCterm = 2.34  
        #Store the pka of the corresponding amino acids
        aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
        #Store the pka of the corresponding amino acids
        aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
        #It will strore the summition of the first term of the equation above
        first=0
        #Iterate though the first term pka
        for key,value in aa2chargePos.items():
            #Store correspoding sum in first
            first+=self.count[key]*(10**aa2chargePos[key]/(10**aa2chargePos[key]+10**ph))
        #also add Nterm pka term to first using the above equation
        first+=(10**aaNterm/(10**aaNterm+10**ph))
        #It will strore the summition of the second term of the equation above
        second=0
        #Iterate though the first term pka
        for key,value in aa2chargeNeg.items():
            #Store correspoding sum in first
            second+=self.count[key]*(10**ph/(10**aa2chargeNeg[key]+10**ph))
        #also add Cterm pka term to first using the above equation
        second+=(10**ph/(10**aaCterm+10**ph))
        #Return the difference according to the equation given(That is used for calculating charge at specific ph)
        return first-second

    def molarExtinction (self):
        '''This method return molar extinction of the amino acid sequence using the 
        extiniction cofficients
        we will use the equation
        𝐸=𝑁𝑌*𝐸𝑌+𝑁𝑊*𝐸𝑊+𝑁𝐶*𝐸𝐶
        AS we were given the repective extinction coffiecient and we already computed the
        cound of the amino acids in self.count
        '''
        
        #Extinction coffecients
        aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}
        #It will store the molarExtinction
        mol=0
        #Iterate though the dictionary contaning the molarExtinction cofficient
        for key,value in aa2abs280.items():
            #Multiply the count with its cofficent according to the equation above and sum in mol
            mol+=self.count[key]*aa2abs280[key]
        #return the mol
        return mol

    def massExtinction (self):
        '''This metho return massExtinction using molar extinction'''
        #Stroe molar weight in myMW
        myMW =  self.molecularWeight()
        #return molarExtinction devided by molar weight(strored in myMW)
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight (self):
        '''This method return the molar weight of the amino acid sequence given
        It just sum the count of each valid amino times it molar weight
        '''
      #store the total molecular weight
        sum=self.mwH2O
      #Iterate though the whole dictonary, for every key and value pair of aa2mw
        for key,value in self.aa2mw.items():
        #it will compute the summation of molecular weight minus molecular weight of water
          sum+=self.aa2mw[key]*self.count[key]-self.mwH2O*self.count[key]
        return sum

class NucParams:
    '''
    This class anylize any NUCLIC acid sequence passed to it
    It ignores all the symbols words and letters in the sequence
    It automatically handle upper case conversion
    It while calculating the AA ignore N and count remaining sequence
    New sequence can be added to it, 
    It caluculate the codon composition, Valid,
    It calculate the nucliotide composition
    It calculate amino-Acid composition

    '''
    def __init__ (self, inString=''):
        '''
        This is the most vital function where all the initializations are done
        Here self.countNuc is initiallized that is a dictionary with codon count that ignore N in codon
        self.Nuc='' is a string that contain the DNA sequence passed with valid
        Nucleotides(A,C,T,G,U,N)
        #Some of the dictionaries were pre defiened
        
        '''
        #Predefiened RNA table
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
        #Defining DNA Table
        self.dnaCodonTable = {key.replace('U','T'):value for key, value in self.rnaCodonTable.items()}
        #inString is convereted to upper case
        inString=inString.upper()
        #create a LIST with lenght of 1, basically convert the sequence to indivisual Amino Acid
        res = [inString[y - 1:y] for y in range(1, len(inString) + 1, 1)]
        #DEFIND possible nucleotides
        self.possibleNucleotides=['A','T','G','C','U','N']
        
        #String that contain all the possible nucleotide without N
        self.nuc=''
        #String that contain all the posible nucleotide with N
        self.Nuc=''
        for i in res:
            #CONDITIONING that the nucliotide is valid and ignoring N
            if i in self.possibleNucleotides and i!='N':
                self.nuc+=i
            #CONDITIONING that the nucliotide is valid
            if i in self.possibleNucleotides:
                self.Nuc+=i
            
        #This will create a list with out N and with length of 3 that will strore valid codons starting from 1
        codons = [self.nuc[y:y+3] for y in range(0, len(self.nuc), 3)]
        #Dictionary initialized to store the count of codons
        self.countNuc={}
        #Iterate though the rna dictionary and store number of times a codon occur ignoring N
        for key,value in self.rnaCodonTable.items():
            self.countNuc[key]=codons.count(key)
        #Iterate though the dna dictionary and store number of times a codon occur ignoring N
        for key,value in self.dnaCodonTable.items():
            self.countNuc[key]=codons.count(key)
    
    def addSequence (self, inSeq):
        '''
            It will add the new sequence to the vaibles declared previously in _init_
            New DNA sequence will be added
            Here self.countNuc is initiallized that is a dictionary with codon count that ignore N in codon
            self.Nuc='' is a string that contain the DNA sequence passed with valid
            Nucleotides(A,C,T,G,U,N)
        '''
        #inSeq is converted to Upper case
        inSeq=inSeq.upper()
        #A string is initialized that will store the amino acids with ignoring N
        ignoreN=''
        #iterate through the string passed and add the self.possibleNucleotides to self.Nuc
        for code in range(0, len(inSeq), 1):
            #Add the valid sequence
            if inSeq[code:code+1] in self.possibleNucleotides:
                self.Nuc+=inSeq[code:code+1]
            #ADd the sequence to a string with out N, THat will be used to calculate codon
            if inSeq[code:code+1] in self.possibleNucleotides and inSeq[code:code+1]!='N':
                ignoreN+=inSeq[code:code+1]
        #devide the whole sequence into 3 size chunks who's value will be added to self.countNuc
        for code in range(0, len(ignoreN), 3):
            if ignoreN[code:code+3] in self.countNuc.keys():
                try:
                    self.countNuc[ignoreN[code:code+3]]+=1
                except:
                    self.countNuc[ignoreN[code:code+3]]=1

    def aaComposition(self):
        '''
        This function return a dictionary that contain amino acid composition
        IT will look at the codon dictionary and find out what is the totall count of amino-acids
        '''
        #The dictianary that will store amino acid composition
        self.aaComp={}
        #Iterate though self.countNuc when there are stop codons put them else put the Amino Acid
        for key,value in self.countNuc.items():
            
            #If the codons are of RNA then the we will look at this condition
            if key in self.rnaCodonTable.keys():
                #If it is a stop codon it will have a key of the codon as it is not translated
                if self.rnaCodonTable[key]=='-':
                    #TRY else are here as for first time any new amino acid is added to dictionary
                    #it will give error and else it should be adde to previous amount of specific amino acid
                    try:
                        self.aaComp['-']=+value
                    except:
                        self.aaComp['-']=value


                    
                else:
                    #TRY else are here as for first time any new amino acid is added to dictionary
                    #it will give error and else it should be adde to previous amount of specific amino acid
                    try:
                        self.aaComp[self.rnaCodonTable[key]]+=value
                    except:
                        self.aaComp[self.rnaCodonTable[key]]=value
                        
             #If the codons are of RNA then the we will look at this condition
            else:
                #If it is a stop codon it will have a key of the codon as it is not translated
                if self.dnaCodonTable[key]=='-':
                    
                    #TRY else are here as for first time any new amino acid is added to dictionary
                    #it will give error and else it should be adde to previous amount of specific amino acid
                    try:
                        self.aaComp['-']+=value
                    except:
                        self.aaComp['-']=value
                else:
                    
                    #TRY else are here as for first time any new amino acid is added to dictionary
                    #it will give error and else it should be adde to previous amount of specific amino acid
                    
                    try:
                        self.aaComp[self.dnaCodonTable[key]]+=value
                    except:
                        self.aaComp[self.dnaCodonTable[key]]=value
        #Return the list of amino acids with the number of time they are translated from codon sequence
        return self.aaComp
    
    
    def nucComposition(self):
        '''This will basically count the number of time any specific Nucleotide has occured and 
        it will also sepratly tackle both DNA sequence without U and RNA without T, but the composition of 
        some sequence with both will also be returned
        '''
        #The self.nucComp store the number of valid nucletides in a DNA or RNA sequence
        self.nucComp={}
        #This dictinary is initialize as all the amino acids at start are 0 therefore, the values are 0
        nucComp={'A':0,'U':0,'C':0,'G':0,'T':0,'N':0}

        #THE case of DNA 
        if self.nuc.count('U')==0:
            #A dictionary that store the DNA sequence
            dna={}
            
            #Iterate though all the valid nucletides and store their count in the dictionary
            for key,value in nucComp.items():
                #THe U will not be return with the dictionary
                if key=='U':
                    continue
                dna[key]=self.Nuc.count(key)
                #THe nucComp will store the valid amino acids
                self.nucComp[key]=self.Nuc.count(key)
            #AS THE Sequence is of DNA the DNA will be returned
            return dna
        #CASE of RNA
        elif self.nuc.count('T')==0: 
            #A dictionary that store the RNA sequence
            rna={}
            #Iterate though all the valid nucletides and store their count in the dictionary
            for key,value in nucComp.items():
            #THe T will not be return with the dictionary
                if key=='T':
                    continue
                rna[key]=self.Nuc.count(key)
                #THe nucComp will store the valid amino acids
                self.nucComp[key]=self.Nuc.count(key)
            #AS THE Sequence is of DNA the DNA will be returned
            return rna
        
        #case where a sequence has both Thymine and uracil, 
        else:
            #Iterate though all the valid nucletides and store their count in the dictionary
            for key,value in nucComp.items():
                #THe nucComp will store the valid amino acids
                self.nucComp[key]=self.Nuc.count(key)
            #Return both uracil and thymine in the dictiory
            return self.nucComp
        
    def codonComposition(self):
        '''
        This will return the composition of the codons with ignoring N
        '''

        #return self.countNUc initialize in _init_ that was a dictioanry or codons for RNA or DNA
        delList=[]
        if self.Nuc.count('U')==0:
            for key,value in self.countNuc.items():
                if key.count('U')!=0:
                    delList.append(key)
        else:
            for key,value in self.countNuc.items():
                if key.count('T')!=0:
                    delList.append(key)
        for i in delList:
            self.countNuc.pop(i)
        return self.countNuc
    
    def nucCount(self):
        '''
        Return the Sum of the number of valid nucliotides
        '''
        #This will sum up the values of the nucleotides and return that
        
        return sum(self.nucComposition().values())