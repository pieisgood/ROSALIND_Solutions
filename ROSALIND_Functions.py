

#Rosalind Functions and Solutions

RNACodoneTable = [["UUU", "F"],["CUU", "L"],["AUU", "I"],["GUU", "V"],\
                  ["UUC", "F"],["CUC", "L"],["AUC", "I"],["GUC", "V"],\
                  ["UUA", "L"],["CUA", "L"],["AUA", "I"],["GUA", "V"],\
                  ["UUG", "L"],["CUG", "L"],["AUG", "M"],["GUG", "V"],\
                  ["UCU", "S"],["CCU", "P"],["ACU", "T"],["GCU", "A"],\
                  ["UCC", "S"],["CCC", "P"],["ACC", "T"],["GCC", "A"],\
                  ["UCA", "S"],["CCA", "P"],["ACA", "T"],["GCA", "A"],\
                  ["UCG", "S"],["CCG", "P"],["ACG", "T"],["GCG", "A"],\
                  ["UAU", "Y"],["CAU", "H"],["AAU", "N"],["GAU", "D"],\
                  ["UAC", "Y"],["CAC", "H"],["AAC", "N"],["GAC", "D"],\
                  ["UAA", "Stop"],["CAA", "Q"],["AAA", "K"],["GAA", "E"],\
                  ["UAG", "Stop"],["CAG", "Q"],["AAG", "K"],["GAG", "E"],\
                  ["UGU", "C"],["CGU", "R"],["AGU", "S"],["GGU", "G"],\
                  ["UGC", "C"],["CGC", "R"],["AGC", "S"],["GGC", "G"],\
                  ["UGA", "Stop"],["CGA", "R"],["AGA", "R"],["GGA", "G"],\
                  ["UGG", "W"],["CGG", "R"],["AGG", "R"],["GGG", "G"]]


def RabbitRecurrence( numMonths, numOffspring ):
    adultpairs = 0
    youngpairs = 1
    rabbits = [0, 1]
    for x in range(1, numMonths):
        rabbits = AtoY(rabbits, numOffspring)
        
    return rabbits[0] + rabbits[1]

def AtoY ( rab , mult):
    young = rab[0]*mult
    adults = rab[1] + rab[0]

    return [adults, young]

def RNAtoProteinS( RNAString ):
    ProteinS = ""
    stop = False
    total = 0
    while not stop:
        
        sRNA = StringGroup(RNAString, total, total + 2)
        sRNA = RNAProteinMap(sRNA)
        if sRNA != "Stop" and sRNA != True:
            ProteinS = ProteinS + sRNA
        else:
            stop = True
        total = total + 3
    return ProteinS

def RNAProteinMap( inRNA ):

    for x in RNACodoneTable:
        if inRNA == x[0]:
            return x[1]

    return True

                  
def SubStringLocations( inString , subString ):
    locations = []
    for x in range(0, len(inString) + 1):
        if subString == inString[x : len(subString) + x ]:
            locations.append(x + 1)
    return locations


def StringGroup( inString, start, stop ):
    outString = ""
    while start <= stop :
        outString = outString + inString[start]
        start = start + 1

    return outString

