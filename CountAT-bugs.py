#!/usr/bin/env python3
# Name: Ziad Kedkad (zkedkad)
# Group Members: none

class dnaString (str):
    def length (self):
        return (len(self))

    def getAT (self):
        num_A = self.count("A")
        num_T = self.count("T")
        return ((num_A + num_T) / dnaString.length(self) )

dna = input("Enter a dna sequence: ")
upperDNA = dna.upper()
coolString = dnaString(upperDNA)
print("AT content = {0:0.1f}".format(coolString.getAT()))