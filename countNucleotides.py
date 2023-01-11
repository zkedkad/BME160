#!/usr/bin/env python3
# Name: Ziad Kedkad (zkedkad)
# Group Members: none

class dnaString (str):
    def length (self):
        return (len(self))

    def getA (self):
        num_A = self.count("A")
        return (num_A)

    def getT (self):
        num_T = self.count("T")
        return (num_T)

    def getG (self):
        num_G = self.count("G")
        return (num_G)

    def getC (self):
        num_C = self.count("C")
        return (num_C)

dna = input("Enter a dna sequence: ")
upperDNA = dna.upper()
coolString = dnaString(upperDNA)
print("Your sequence is {0:0.1f}".format(coolString.length()) + " nucleotides long with the following break down of bases:")
print("number As = {0:0.1f}".format(coolString.getA()) + "number Cs = {0:0.1f}".format(coolString.getC()) + "number Gs = {0:0.1f}".format(coolString.getG()) + "number Ts = {0:0.1f}".format(coolString.getT()))