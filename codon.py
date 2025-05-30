import os
import sys

from codon1.Bias import bias
from codon1.Correlation import correlation
from codon1.Correspondence import COA
from codon1.ENc import ENc
from codon1.Neutrality import neutrality
from codon1.Optimalcodon import optimalcodon
from codon1.Similar import similarcodon



def main():
    if len(sys.argv) < 2:
        print("")
    modulename = sys.argv[1]


    outpath = "codonOutputs/" + modulename + "/"
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    elif modulename == "bias":
        if len(sys.argv) != 4:
            print("Usage: python codon.py bias <cds.fasta> <figName>")
        figPath = outpath + sys.argv[3]
        bias.run(sys.argv[2], figPath, outpath)

    elif modulename == "coa":
        if len(sys.argv) != 4:
            print("Usage: python codon.py bias <cds.fasta> <figName>")
        figPath = outpath + sys.argv[3]
        COA.run(sys.argv[2], figPath, outpath)

    elif modulename == "enc":
        if len(sys.argv) != 4:
            print("Usage: python codon.py bias <cds.fasta> <figName>")
        figPath = outpath + sys.argv[3]
        ENc.run(sys.argv[2], figPath)

    elif modulename == "correlation":
        if len(sys.argv) != 5:
            print("Usage: python codon.py bias <cds.fasta> <figName>")
        figPath = outpath + sys.argv[4]
        correlation.run(sys.argv[2], sys.argv[3], figPath, outpath)

    elif modulename == "neutrality":
        if len(sys.argv) != 4:
            print("Usage: python codon.py bias <cds.fasta> <figName>")
        figPath = outpath + sys.argv[3]
        neutrality.run(sys.argv[2], figPath, outpath)

    elif modulename == "optimal":
        if len(sys.argv) != 5:
            print("Usage: python codon.py bias <cds.fasta> <figName>")
        figPath = outpath + sys.argv[4]
        optimalcodon.run(sys.argv[2], sys.argv[3], figPath, outpath)

    elif modulename == "similar":
        if len(sys.argv) != 4:
            print("Usage: python codon.py bias <cds.fasta> <figName>")

        similarcodon.run(sys.argv[2], sys.argv[3], outpath)

if __name__ == "__main__":
    main()