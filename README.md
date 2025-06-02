# Install

## Quick Install CUB

You can install CodonUsageBiad(CUB) using the following command:

```bash
git https://github.com/YLhai/Codon.git
cd CodonUsageBias
pip install -r request.txt
```

## Dependencies
### Python

version 3.11

### Modules

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

python CUB.py bias <CDS.fasta> <figName>

```

### coa

This command can generate a COA-plot.

```bash

python CUB.py coa <CDS.fasta> <figName>

```
### enc

This command can generate a ENC-plot.

```bash

python CUB.py enc <codonW.out> <figName>

```

### correlation

```bash

python CUB.py correlation <CDS.fasta> <codonW.out> <figName>

```
### neutrality

```bash

python CUB.py neutrality <CDS.fasta> <figName>

```
### optimal

```bash

python CUB.py optimal <CDS.fasta> <codonW.out> <figName>

```
### similar

```bash

python CUB.py similar <cdsA.fasta> <cdsB.fasta>

```
### rscu

This command can generate a stacked barplot  illustrating RSCU values of synonymous codons.
```bash

python CUB.py rscu <CDS.fasta> <figName>
```

