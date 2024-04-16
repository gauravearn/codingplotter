def plotCDS(inputfile, outputfile):
    """
    Author: Gaurav Sablok
    Universitat Potsdam
    Date: 2024-4-16
    a coding region plotter for the protein annotations 
    given a reference genome with the proteome aligned,
    it will plot the coding regions.
    @inputfile = file aligned to the genome using the protein hints
    @outputfile = file to which the binary layout for the protein hints should be written. 
    """
    readfile = [i for i in open(inputgff, "r").readlines() if "#" not in i]
    with open(inputgff + ".coding.gff", "w") as writegff:
        writegff.write("col0 \t col1 \t col2 \t col3 \t col4 \t col5 \t col6 \t col7 \t col8 \t col9\n")
        for line in readfile:
            writegff.write(line)
        writegff.close()
    iterator = [i.strip().split() for i in open(inputgff + ".coding.gff").readlines() if i.strip().split()[2] == "CDS"]
    codinglength = []
    for i in range(len(readfile)):
        codinglength.append(int(readfile[i].strip().split()[4])-int(readfile[i].strip().split()[3]))
    binrange_500 = []
    binrange_100 = []
    binrange_1000 = []
    binrange_2000 = []
    for i in range(len(codinglength)):
        if codinglength[i] <= 100:
            binrange_100.append(codinglength[i])
        elif codinglength[i] <= 500 and codinglength[i] >= 100:
            binrange_500.append(codinglength[i])
        elif codinglength[i] <= 1000 and codinglength[i] >= 500:
            binrange_1000.append(codinglength[i])
        elif codinglength[i] <= 2000 and codinglength[i] >= 1000:
            binrange_2000.append(codinglength[i])
        else:
            pass
    lengthbins = [len(binrange_100), len(binrange_500), len(binrange_1000), len(binrange_2000)]
    figure = sns.histplot(pd.DataFrame(lengthbins, columns = ["lengthbins"])).get_figure()
    figure.savefig("save.png")
    with open(outputfile, "w") as writebins:
        writebins.write("The length distribution for the coding regions for the less than 100bp are:")
        for i in range(len(binrange_100)):
            writebins.write(f"{binrange_100[i]}\n")
        writebins.write("The length distribution for the coding regions for the less than 500bp are:")
        for i in range(len(binrange_500)):
            writebins.write(f"{binrange_500[i]}\n")
        writebins.write("The length distribution for the coding regions for the less than 1000bp are:")
        for i in range(len(binrange_1000)):
            writebins.write(f"{binrange_1000[i]}\n")
        writebins.write("The length distribution for the coding regions for the less than 2000bp are:")
        for i in range(len(binrange_2000)):
            writebins.write(f"{binrange_2000[i]}\n")
        writebins.close()
