# This script goes through the prosite .dat file and create a dictionary with 
# accession pattern pairs
# also converts patterns to be compatible with python re module



# TODO : for the pattern of PS00012 we got by parsing
# [DEQGSTALMKRH][LIVMFYSTAC][GNQ][LIVMFYAG][DNEKHS]S[LIVMST]
# [ADEGHIKLMNQRSTVW][STAGCPQLIVMF][LIVMATN][DENQGTAKRHLM]
# [LIVMWSTA][LIVGSTACR][ACDEFGHKMNQRSTVW][ACDEFGHIKLMNPQRSTW][LIVMFA]
# but the actual is 
#[DEQGSTALMKRH]-[LIVMFYSTAC]-[GNQ]-[LIVMFYAG]-[DNEKHS]-S-[LIVMST]-{PCFY}-
# [STAGCPQLIVMF]-[LIVMATN]-[DENQGTAKRHLM]-[LIVMWSTA]-[LIVGSTACR]-{LPIY}-
# {VY}-[LIVMFA]
# it removes the {} brackets (for example {PCFY} and doesn't add them to the regex / changes it
# into the correct regex format, we need to fix that

import re

aaList =  ["A",
           "C",
           "D",
           "E",
           "F",
           "G",
           "H",
           "I",
           "K",
           "L",
           "M",
           "N",
           "P",
           "Q",
           "R",
           "S",
           "T",
           "V",
           "W",
           "Y"]

prositeFile = open('prosite.dat.txt','r')
outputFile = open('prosite_preprocessed.txt','a')

currentDict = {}

#currentPA holds the prosite pattern from dat file; if it spans more than one line, read multiple files and concatenate 
#pattern from consecutive lines with PA tag in the beginning of the line

#when a new line with AC tag is read, currentPA from previous accession is cleared
# AC is the row accession for the new protein (i.e. in that row accession of protein such as PS10203)


# PA is the row where the pattern is found 



for line in prositeFile:
    line = line.strip()
    if re.match("AC   ", line):
        currentAC = line.split()[1][:-1]
        if currentAC not in currentDict.keys():
            currentDict[currentAC] = []
            if 'currentPA' in globals():
                del currentPA
    elif re.match("PA   ", line): 
        currentPA = line.split()[1] # Gets the current Pattern
        if currentPA[-1] == ".":
            currentPA = currentPA[:-1]
        if lastLineTag == "PA   ": # check if in the last line thre was also a PA
            currentPA =  lastLinePA + currentPA # if yes, concatenate the patterns (both Rows PA)
            # since we iterate over all lines in the main loop, we catch all PA rows and concatenate everything with this code
    else:
        if not 'currentPA' in globals():
            lastLineTag = line[:5]
            continue # skip the current loop (i.e. go until we find again AC/PA)

        
        # TODO : for me this is just a bit confusing, shouldn't this code happen after the whole PA was found ?
        # TODO : because right now it is happening in the else case ? 
        currentPAList = currentPA.split("-")
        refinedCurrentPAList = []

        for n in range(len(currentPAList)):
        # Redefine the Pattern such that it works with Python Regex

            #change to range

            betweenCurlyList = []
            betweenSquaresList = []

            currentAAList = aaList.copy()
            if re.search("x", currentPAList[n]):
                currentPAList[n] = currentPAList[n].replace("x",".")
            if re.search("{" , currentPAList[n]):
                currentPAList[n] = currentPAList[n].replace("{","#")
            if re.search("}" , currentPAList[n]):
                currentPAList[n] = currentPAList[n].replace("}","%")
            if "(" in currentPAList[n]:
                currentPAList[n] = currentPAList[n].replace("(","{")
            if ")" in currentPAList[n]:
                currentPAList[n] = currentPAList[n].replace(")","}")
            if currentPAList[n][0] == "<":
                currentPAList[n] = currentPAList[n].replace("<","^")
            if currentPAList[n][-1] == ">":
                currentPAList[n] = currentPAList[n].replace(">","$")
            if re.search("#" , currentPAList[n]):
                element = currentPAList[n]
                betweenCurly =  element[element.find("#")+1 : element.find("%")]
                if len(betweenCurly)>1:
                    betweenCurlyList = list(betweenCurly)
                    for aa in betweenCurlyList:
                        if aa in currentAAList:
                            currentAAList.remove(aa)
                    betweenCurly = str("[") + "".join(currentAAList) + str("]")
                    #print(betweenCurly)
                    currentPAList[n] = re.sub(r'#.*%', betweenCurly,currentPAList[n])
                else:
                    if betweenCurly in currentAAList:
                            currentAAList.remove(betweenCurly)
                    betweenCurly = str("[") + "".join(currentAAList) + str("]")
                    #print(betweenCurly)
                    currentPAList[n] = re.sub(r'#.*%', betweenCurly,currentPAList[n])

            
        finalPA = "".join(currentPAList)
        if len(finalPA) > 1:
            if finalPA[-2] == ">":
                finalPA1 = finalPA[:-2] + "]"
                refinedCurrentPAList.append(finalPA1)
                finalPA2 = finalPA[:finalPA.rindex("[")] + "$"
                refinedCurrentPAList.append(finalPA2)
            else:
                refinedCurrentPAList.append(finalPA)

        currentDict[currentAC] = refinedCurrentPAList
        del currentPA
#save information from the present line for retrival while reading the next line
# useful while dealing with long Prosite patterns spanning multiple lines
    lastLineTag = line[:5]
    if re.match("PA   ", line):
        if 'currentPA' in globals():
            lastLinePA = currentPA
        else:
            lastLinePA = line.split()[1]

for key in currentDict.keys():
    if len(currentDict[key]) != 0:
        print(key,currentDict[key] )
        outputFile.write(key + "\t" + "\t&\t".join(currentDict[key])+ "\n")