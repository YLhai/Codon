#!/usr/bin/env python3

import os
import sys
from codonUsageBias.Bias import bias
from codonUsageBias.Correlation import correlation
from codonUsageBias.Correspondence import COA
from codonUsageBias.ENc import ENc
from codonUsageBias.Neutrality import neutrality
from codonUsageBias.Optimalcodon import optimalcodon
from codonUsageBias.Similar import similarcodon
from codonUsageBias.RSCU import rscu

Commands = {
    "bias": "PR2 bias analysis",
    "coa": " correspondence analysis(COA)",
    "enc": "ENC-GC3 plot analysis ",
    "optimal": "optimal codons analysis ",
    "neutrality": "neutrality plot analysis ",
    "correlation": "correlation analysis ",
    "similar": "similarity analysis of two species",
    "rscu": "RSCU values of synonymous codons analysis ",

}

def help(type):
    if type == 1:
        print("Usage: \t python CUB.py <command> <arguments>\n"
              "version: \t 1.1\n\n"
              "Description: CUB is a plot analysis tool of codon usage bias.\n\n"
              "Commands:"
              )
        for key, value in Commands.items():
            print(key + "\t" + value)
        print(
            "\nWriten by hai and peng in 2025.\n"
            "Github: https://github.com/YLhai/Codon"
        )


def main():

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

    if modulename == "bias":
        if len(sys.argv) < 4 or len(sys.argv) > 5:
            print("Usage: python CUB.py bias <cds.fasta> <figName>\n\n"
                  "cds.fasta\t the fasta file of codons.\n"
                  "figName\t the name of the output figure.\n"
                  )
            sys.exit(0)
        figPath = outpath + sys.argv[3]
        figType = sys.argv[4] if len(sys.argv) == 5 else "pdf"
        bias.run(sys.argv[2], figPath, outpath, figType)

    elif modulename == "coa":
        if len(sys.argv) <4 or len(sys.argv) > 5:
            print("Usage: python CUB.py coa <cds.fasta> <figName>\n\n"
                  "cds.fasta\t the fasta file of codons.\n"
                  "figName\t the name of the output figure.\n"
                  )
            sys.exit(0)
        figPath = outpath + sys.argv[3]
        figType = sys.argv[4] if len(sys.argv) == 5 else "pdf"
        COA.run(sys.argv[2], figPath, outpath,figType)

    elif modulename == "enc":
        if len(sys.argv) < 4 or len(sys.argv) > 5:
            print("Usage: python CUB.py enc <codonW.out> <figName>\n\n"
                  "codonW.out\t the outputfile of codonW using the codon.fasta.\n"
                  "figName\t the name of the output figure.\n"
                  )
            sys.exit(0)
        figPath = outpath + sys.argv[3]
        figType = sys.argv[4] if len(sys.argv) == 5 else "pdf"
        ENc.run(sys.argv[2], figPath,figType)

    elif modulename == "correlation":
        if len(sys.argv) <5 or len(sys.argv) > 6:
            print("Usage: python CUB.py correlation <cds.fasta> <codonW.out> <figName>\n\n"
                  "cds.fasta\t the fasta file of codons.\n"
                  "codonW.out\t the outputfile of codonW using the codon.fasta.\n"
                  "figName\t the name of the output figure.\n"
                  )
            sys.exit(0)
        figPath = outpath + sys.argv[4]
        figType = sys.argv[5] if len(sys.argv) == 6 else "pdf"
        correlation.run(sys.argv[2], sys.argv[3], figPath, outpath,figType)

    elif modulename == "neutrality":
        if len(sys.argv) <4 or len(sys.argv) > 5:
            print("Usage: python CUB.py neutrality <cds.fasta> <figName>\n\n"
                  "cds.fasta\t the fasta file of codons.\n"
                  "figName\t the name of the output figure.\n"
                  )
            sys.exit(0)
        figPath = outpath + sys.argv[3]
        figType = sys.argv[4] if len(sys.argv) == 5 else "pdf"
        neutrality.run(sys.argv[2], figPath, outpath,figType)

    elif modulename == "optimal":
        if len(sys.argv) <5 or len(sys.argv) > 6:
            print("Usage: python CUB.py optimal <cds.fasta> <codonW.out> <figName>\n\n"
                  "cds.fasta\t the fasta file of codons.\n"
                  "codonW.out\t the outputfile of codonW using the codon.fasta.\n"
                  "figName\t the name of the output figure.\n"
                  )
            sys.exit(0)
        figPath = outpath + sys.argv[4]
        figType = sys.argv[5] if len(sys.argv) == 5 else "pdf"
        optimalcodon.run(sys.argv[2], sys.argv[3], figPath, outpath,figType)

    elif modulename == "similar":
        if len(sys.argv) != 4:
            print("Usage: python CUB.py bias <cdsA.fasta> <cdsB.fasta>\n\n"
                  "cdsA.fasta\t the fasta file of codons of speciesA\n"
                  "cdsB.fasta\t the fasta file of codons of speciesB\n"
                  )

        similarcodon.run(sys.argv[2], sys.argv[3], outpath)

    elif modulename == "rscu":
        if len(sys.argv) <4 or len(sys.argv) > 5:
            print("Usage: python CUB.py rscu <cds.fasta> <figName>\n\n"
                  "cds.fasta\t the fasta file of codons.\n"
                  "figName\t the name of the output figure.\n"
                  )
            sys.exit(0)
        figPath = outpath + sys.argv[3]
        figType = sys.argv[4] if len(sys.argv) == 5 else "pdf"
        rscu.run(sys.argv[2], figPath, outpath,figType)


if __name__ == "__main__":
    main()