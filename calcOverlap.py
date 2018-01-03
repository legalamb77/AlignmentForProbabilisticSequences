import sys

def main():
    filecount = int(sys.argv[1])
    queryLength = int(sys.argv[2])
    filenames = ["fullres"+str(num)+".txt" for num in range(filecount)]
    #Extract unmutated, start, and stop, and put in a list in order
    unmut = open("unmutated.txt", "r")
    allUnmut = []
    for line in unmut:
        word, start = line.split(",")
        start = int(start.strip())
        allUnmut.append((word.strip(),start))
    unmut.close()
    #Extract alignment start and stop
    count = 0
    scoringDict = dict()
    bestScoring = dict()
    for name in filenames:
        f = open(name, 'r')
        overlaps = []
        for line in f:
            if line[0] == "(":
                #Extract start and stop
                allvals = line.split(",")
                Astop = int(allvals[-1].rstrip(")\n"))
                Astart = int(allvals[-2].strip())
                Ostart = int(allUnmut[count][1])
                Ostop = Ostart+queryLength
                #Measure against corresponding start and stop of origin
                overlap = None
                sStart = None
                sStop = None
                if Astart > Ostop:
                    overlap = 0
                elif Ostart > Astop:
                    overlap = 0
                else:
                    # There is some overlap
                    sStart = max(Astart, Ostart)
                    sStop = min(Astop, Ostop)
                    overlap = float(sStop-sStart)/float(queryLength)
                overlaps.append(overlap)
                #Store in a dictionary, mapping from the unmutated word to the accuracy
            else:
                pass
        if sum(overlaps)==0:
            scoringDict[allUnmut[count][0]] = 0
        else:
            scoringDict[allUnmut[count][0]] = sum(overlaps)/float(len(overlaps))
        #Best scoring
        best = max(overlaps)
        bestScoring[allUnmut[count][0]] = best
        count += 1
        f.close()
    bestOut = open("bestOverlap.txt", "w")
    bestTotal = 0
    for val in bestScoring.items():
        bestOut.write(str(val[0])+","+str(val[1]))
        bestTotal+=val[1]
        bestOut.write("\n")
    bestCalc = bestTotal/(float(queryLength)*(count+1))
    bestOut.write("\n Total Percent Overlap: "+str(bestCalc))
    bestOut.close()
    out = open("overlap.txt", "w")
    totalOverlap = 0
    for val in scoringDict.items():
        out.write(str(val[0])+","+str(val[1]))
        totalOverlap+=val[1]
        out.write("\n")
    calc = totalOverlap/(float(queryLength)*(count+1))
    out.write("\n Total Percent Overlap: "+str(calc))
    out.close()

if __name__ == "__main__":
    main()
