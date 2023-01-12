#!/usr/bin/env python3
# Name: Ziad Kedkad (zkedkad)
# Group Members: none
"""
This program inputs a DNA sequence from a user. It then counts
the total length of that sequence and outputs both the length
of the sequence and the number of each nucleotide present.
"""
# below I create the class dnaString which is defined as a string
# below the class I
class dnaString(str):

    def length(self):
        return(int(len(self)))

    def getA(self):
        num_A = self.count("A")
        return(num_A)

    def getT(self):
        num_T = self.count("T")
        return(num_T)

    def getG(self):
        num_G = self.count("G")
        return(num_G)

    def getC(self):
        num_C = self.count("C")
        return(num_C)

dna = input("Enter a dna sequence: ")
upperDNA = dnaString.upper(dna)
coolString = dnaString(upperDNA)
# Here the user is greeted with two statements with the calculated results
print("Your sequence is {0} ".format(coolString.length()) + " nucleotides long with the following break down of bases:")
print("number As = {0}".format(coolString.getA()) + " number Cs = {0}".format(coolString.getC()) + " number Gs = {0}".format(coolString.getG()) + " number Ts = {0}".format(coolString.getT()))