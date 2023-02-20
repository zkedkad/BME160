class DNAstring(str):
    def length(self):
        return (length(self))

    ''' Return an upcased version of the string, collapsing a single run of Ns.'''
    def purify(self):
        '''Setting preset type for future use'''
        pureDNA = ''
        N_count = 0
        '''This takes dna and counts the number of N within statement'''
        for charecters in self:
            if charecters.upper() == 'N':
                '''This adds every N counted on each other'''
                N_count += 1
            else: #After all ns accounted for, move on the next if statement to compile the int into the string
                if N_count > 0:
                    pureDNA += '{' + str(N_count) + '}' #n integer is put within brackets
                    N_count = 0 #set n back to 0 or it will continue to run if statement
                pureDNA += charecters.upper() #int is compiled with full purified DNA statement
        if N_count > 0:
            pureDNA += '{' + str(N_count) + '}'
        return pureDNA


def main():
    ''' Get user DNA data and clean it up.'''
    while (True):
        data = input('DNA data?')
        thisDNA = DNAstring(data)
        pureData = thisDNA.purify()
        print(pureData)


main()
