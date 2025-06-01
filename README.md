# Install

## Quick Install CUB

You can install CodonUsageBiad(CUB) using the following command:

```bash
git https://github.com/YLhai/Codon.git
cd CodonUsageBias
pip install -r request.txt
```

## Dependencies
### python

version 3.11

### modoule

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

### bias

```bash

python CUB.py bias <yourCDS.fasta> <figName>

```

### coa

```bash

python CUB.py coa <yourCDS.fasta> <figName>

```
### enc

```bash

python CUB.py enc <codonW.out> <figName>

```

### correlation

```bash

python CUB.py correlation <yourCDS.fasta> <codonW.out> <figName>

```
### neutrality

```bash

python CUB correneutrality <cds.fasta> <figName>

```
### optimal

```bash

python CUB.py optimal <cds.fasta> <codonW.out> <figName>

```
### similar

```bash

python CUB bias <cdsA.fasta> <cdsB.fasta>

```
### rscu

```bash

python CUB.py rscu <cds.fasta> <figName>

