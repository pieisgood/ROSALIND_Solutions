

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

def GCcontent( dna ):
    gcContent = 0.0
    
    for x in range(0, len(dna)):
        if dna[x] == 'G' or dna[x] == 'C':
            gcContent = gcContent + 1.0
    return  gcContent/len(dna) * 100

def PrintGCcontent( FASTA ):

    for x in FASTA:
        #print x
        content = GCcontent( x[1] )
        print x[0]
        print content

def CompareOverLapO3(tail, head):
    if tail[1][len(tail[1])-3:len(tail[1])] == head[1][0:3]:
        return True
    else:
        return False

def OverLapGraph( FASTA ):

    for x in FASTA:
        for y in FASTA:
            if x[0] != y[0]:
                edge = CompareOverLapO3(x, y)
                if edge:
                    print x[0] + ' ' + y[0]
                    

def FASTAImport( inFile ):
    dnaAssociation = []
    holder = open(inFile, 'r')
    notDone = True
    sum_name = ''
    rosalind = holder.readline()
    rosalind = rosalind[1:len(rosalind)-1]
    dna = ''
    
    while notDone:
        name = holder.readline()
        chain = name[0:len(name)-1]
        name = name[1:len(name)-1]
        if name != '' and name[0:8] != 'Rosalind':
            sum_name = sum_name + chain
        elif name != '' and name[0:8] == 'Rosalind':
            dna = [rosalind, sum_name]
            dnaAssociation.append(dna)
            rosalind = name
            sum_name = ''
            name = ''
        elif name == '':
            notDone = False
            dna = [rosalind, sum_name]
            dnaAssociation.append(dna)

    return dnaAssociation

OverLapGraph( FASTAImport('rosalind_grph.txt') )
            

