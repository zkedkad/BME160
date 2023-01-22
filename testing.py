#!/usr/bin/env python3 
# Name: Your full name  
# Group Members: List full names or “None”

'''
Read a DNA string from user input and return a collapsed substring of embedded Ns to: {count}.

Example: 
 input: AaNNNNNNGTC
output: AA{6}GTC

Any lower case letters are converted to uppercase
'''


# Function to count length of DNA
class DNAstring:
    def __init__(self, dna):
        self.dna = dna

    def length(self):
        return (len(self))

        # Function to clean DNA

    def purify(self):
        upperdna = self.dna.upper()
        comma_data = upperdna.replace("", ",")
        list_data = comma_data.split(",")
        print(list_data)
        return

# Function that is called that runs data through purifier
def main():
    while (True):
        data = input('DNA data?')
        thisDNA = DNAstring(data)
        pureData = thisDNA.purify()
        print(pureData)


main()