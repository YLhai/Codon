
CodonUsageBias(CUB) is a plot and analysis tool of codon bias. You can use this 
tool to analyse amino acid usage bias, COA, RSCU, ENc, correlation, neutrality and similar, 
and to generate plots to illustrate them. However, just two files(`CDS.fasta` and `codonW.out`) will be prepared.
# Install

## Quick Install CUB

You can install CodonUsageBias(CUB) using the following command:

```bash
git https://github.com/YLhai/Codon.git
cd CodonUsageBias
pip install -r request.txt
```

## Dependencies
### Python

version 3.11

### Python Packages

```commandline
bio==1.8.0
biopython==1.85
biothings_client==0.4.1
matplotlib==3.10.3
mygene==3.2.2
numpy==2.2.6
pandas==2.2.3
pillow==11.2.1
platformdirs==4.3.8
prince==0.7.1
pytz==2025.2
scipy==1.15.3
seaborn==0.13.2
```

# Commands

The parameters of all following commands are the same:
 
`<CDS.fasta>` is the CDS fasta file.

`<codonW.out>` is the output file of codonW by CDS fasta file.

`<figName>` is the name of output figures, you just need input the name without file extensions.

`[figFormat]` The default format is 'pdf' if you not appoint it. 

### bias

This command will help you analysis codon usage bias of ever single CDS and all CDS, 
and plot two figures to show them.

```bash

python CUB.py bias <CDS.fasta> <figName> [figFormat]

```

### coa

This command can generate a COA-plot.

```bash

python CUB.py coa <CDS.fasta> <figName> [figFormat]

```
### enc

This command can generate a ENC-plot.

```bash

python CUB.py enc <codonW.out> <figName> [figFormat]

```

### correlation

```bash

python CUB.py correlation <CDS.fasta> <codonW.out> <figName> [figFormat]

```
### neutrality

```bash

python CUB.py neutrality <CDS.fasta> <figName> [figFormat]

```
### optimal
This command can make a optimalcodon analysis and plot the result.

```bash

python CUB.py optimal <CDS.fasta> <codonW.out> <figName> [figFormat]

```
### similar

You can get the similarity of A and B species using their CDS.

```bash

python CUB.py similar <cdsA.fasta> <cdsB.fasta>

```
### rscu

This command can generate a stacked barplot  illustrating RSCU values of synonymous codons.
```bash

python CUB.py rscu <CDS.fasta> <figName> [figFormat]
```

