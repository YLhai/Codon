import os
import sys
from codon1.Bias import bias
from codon1.Correlation import correlation
from codon1.Correspondence import COA
from codon1.ENc import ENc
from codon1.Neutrality import neutrality
from codon1.Optimalcodon import optimalcodon
from codon1.Similar import similarcodon
from codon1.RSCU import rscu

def help(type):
    if type == 1:
        Commands = {
            "bias": "***",
            "coa": "***",
            "enc": "***",
            "optimal": "***",
            "neutrality": "***",
            "correlation": "***",
            "similar": "***",
            "rscu": "***",

        }
        print("Usage: \t python codon.py <command> <arguments>\n"
              "version: \t 1.1\n\n"
              "Description: Codon is a anilize tool of codon bias.\n\n"
              "Commands:"
              )
        for key, value in Commands.items():
            print(key + "\t" + value)
        print(
            "\nWriten by hai and peng in 2025."
        )


def main():

    Commands = {
        "bias":"***",
        "coa":"***",
        "enc":"***",
        "optimal": "***",
        "neutrality":"***",
        "correlation": "***",
        "similar": "***",
        "rscu": "***",

    }

    # 未输入命令
    if len(sys.argv) < 2:
        help(1)
        sys.exit(0)

    modulename = sys.argv[1]
    # 输入命令错误
    if modulename not in Commands.keys():
        help(1)
        sys.exit(0)

    outpath = "codonOutputs/" + modulename + "/"
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    elif modulename == "bias":
        if len(sys.argv) != 4:
            print("Usage: python codon.py bias <cds.fasta> <figName>\n\n"
                  "cds.fasta\t the fasta file of codons.\n"
                  "figName\t the name of the output figure.\n"
                  )
            sys.exit(0)
        figPath = outpath + sys.argv[3]
        bias.run(sys.argv[2], figPath, outpath)

    elif modulename == "coa":
        if len(sys.argv) != 4:
            print("Usage: python codon.py coa <cds.fasta> <figName>\n\n"
                  "cds.fasta\t the fasta file of codons.\n"
                  "figName\t the name of the output figure.\n"
                  )
            sys.exit(0)
        figPath = outpath + sys.argv[3]
        COA.run(sys.argv[2], figPath, outpath)

    elif modulename == "enc":
        if len(sys.argv) != 4:
            print("Usage: python codon.py enc <codonW.out> <figName>\n\n"
                  "codonW.out\t the outputfile of codonW using the codon.fasta.\n"
                  "figName\t the name of the output figure.\n"
                  )
            sys.exit(0)
        figPath = outpath + sys.argv[3]
        ENc.run(sys.argv[2], figPath)

    elif modulename == "correlation":
        if len(sys.argv) != 5:
            print("Usage: python codon.py correlation <cds.fasta> <codonW.out> <figName>\n\n"
                  "cds.fasta\t the fasta file of codons.\n"
                  "codonW.out\t the outputfile of codonW using the codon.fasta.\n"
                  "figName\t the name of the output figure.\n"
                  )
            sys.exit(0)
        figPath = outpath + sys.argv[4]
        correlation.run(sys.argv[2], sys.argv[3], figPath, outpath)

    elif modulename == "neutrality":
        if len(sys.argv) != 4:
            print("Usage: python codon.py neutrality <cds.fasta> <figName>\n\n"
                  "cds.fasta\t the fasta file of codons.\n"
                  "figName\t the name of the output figure.\n"
                  )
            sys.exit(0)
        figPath = outpath + sys.argv[3]
        neutrality.run(sys.argv[2], figPath, outpath)

    elif modulename == "optimal":
        if len(sys.argv) != 5:
            print("Usage: python codon.py optimal <cds.fasta> <codonW.out> <figName>\n\n"
                  "cds.fasta\t the fasta file of codons.\n"
                  "codonW.out\t the outputfile of codonW using the codon.fasta.\n"
                  "figName\t the name of the output figure.\n"
                  )
            sys.exit(0)
        figPath = outpath + sys.argv[4]
        optimalcodon.run(sys.argv[2], sys.argv[3], figPath, outpath)

    elif modulename == "similar":
        if len(sys.argv) != 4:
            print("Usage: python codon.py bias <cdsA.fasta> <cdsB.fasta>>\n\n"
                  "cdsA.fasta\t the fasta file of codons of speciesA\n"
                  "cdsB.fasta\t the fasta file of codons of speciesB\n"
                  )

        similarcodon.run(sys.argv[2], sys.argv[3], outpath)

    elif modulename == "rscu":
        if len(sys.argv) != 4:
            print("Usage: python codon.py rscu <cds.fasta> <figName>\n\n"
                  "cds.fasta\t the fasta file of codons.\n"
                  "figName\t the name of the output figure.\n"
                  )
            sys.exit(0)
        figPath = outpath + sys.argv[3]
        rscu.run(sys.argv[2], figPath, outpath)



if __name__ == "__main__":
    main()