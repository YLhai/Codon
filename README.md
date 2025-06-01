## 安装

python版本3.11

你可以使用一下命令安装CodonUsageBias:
```bash
git https://github.com/YLhai/Codon.git
cd CodonUsageBias
pip install -r request.txt
```

## 使用

### bias

```bash

codonUB bias <yourCDS.fasta> <figName>

```

### coa

```bash

codonUB coa <yourCDS.fasta> <figName>

```
### enc

```bash

codonUB enc <codonW.out> <figName>

```

### correlation

```bash

codonUB correlation <yourCDS.fasta> <codonW.out> <figName>

```
### neutrality

```bash

codonUB correneutrality <cds.fasta> <figName>

```
### optimal

```bash

codonUB optimal <cds.fasta> <codonW.out> <figName>

```
### similar

```bash

codonUB bias <cdsA.fasta> <cdsB.fasta>

```
### rscu

```bash

codonUB rscu <cds.fasta> <figName>

```
