#!/usr/bin/python
import sys

if __name__ == "__main__":
    # get inputs
    assert(len(sys.argv) == 2)
    fileName = sys.argv[1]

    # open file
    file = open(fileName,"r")

    # BGQ
    freq = (1./1.6)*10**(-9)

    
    # read file into rawData
    rawData = []
    for line in file:
        rawData.append(line)

    # add data into data
    data = []
    for i in range(1,len(rawData),2):
        #print "i = ",i
        #print rawData[i]
        
        # go through line of data and parse into a row
        row = []
        sbuffer = ""
        started = False
        for j in range(0,len(rawData[i])):
            # scanner hits a non space
            if ((rawData[i][j] != " ")):
                sbuffer = sbuffer+rawData[i][j]
                #print sbuffer
                started = True
            # scanner hits a space and sbuffer is building

            elif started:
                row.append(eval(sbuffer))
                sbuffer = ""
                started = False

        #print sbuffer
        row.append(sbuffer)
        data.append(row)

        pass

    # close file
    file.close()
    print '%5s %5s %5s %6s %6s %9s %10s %10s %10s %10s %10s %10s %15s' % (
        "Nx", "Ny", "Ngp", 
        "nIter", "nRanks", "nThreads", 
        "tInit", "tCalc", "tComm", "tbatch","tTotal","tC/tcal","(tC+tb)/tCal")
    for row in range(len(data)):
        print '%5d %5d %5d %6d %6d %9d %10.3e %10.3e %10.3e %10.3e %10e %10e %10e' % (
            data[row][0],data[row][1],data[row][2],
            data[row][3],data[row][4],data[row][5],
            data[row][6]/float((data[row][2]*data[row][3])),
            data[row][7]/float((data[row][2]*data[row][3])),
            data[row][8]/float((data[row][2]*data[row][3])),
            float(data[row][9])/float((data[row][2]*data[row][3])),
            (data[row][6]+data[row][7]+data[row][8]+float(data[row][9]))/float((data[row][2]*data[row][3])),
            data[row][8]/data[row][7],
            (data[row][8]+float(data[row][9]))/data[row][7],
        )
        
        pass

    
    pass
