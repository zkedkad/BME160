#aminoAcids = ['A', 'D', 'K', 'R']
#aminoAcids[1:1] = 'W'
#print(aminoAcids)
#del aminoAcids[3]
#print(aminoAcids)
counts = {'A':0, 'C':0, 'G':0, 'T':0}
dna = 'ACGTACGT'

for nuc in dna:
    counts[nuc] += 1
    counts[nuc] = counts[nuc] + 1



emptyDict = dict()       # same as {}
ThreeToOne = {'Ala':'A', 'Cys':'C'}
myAA = ThreeToOne['Ala'] # lookup of 1 letter code
ThreeToOne['Tyr'] = 'Y'  # adds a new AA to dict
print (ThreeToOne.get('Trp'))

answer = ThreeToOne.get(data) or nucleotide.get(data) or 'I dont know'