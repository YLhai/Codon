from Bio import SeqIO
from Bio.SeqUtils import Seq


def GC123(in_fasta,out_results):

    def calculate_gc(seq):
        """手动计算GC含量百分比"""
        gc = seq.count('G') + seq.count('C')
        total = len(seq)
        return round((gc / total) * 100, 2) if total > 0 else 0.0



    def calculate_GC123(seq):
        GC1 = []
        GC2 = []
        GC3 = []
        for i in range(1, len(seq) + 1):
            if i % 3 == 1:
                GC1.append(str(seq[i - 1]))
            elif i % 3 == 2:
                GC2.append(str(seq[i - 1]))
            else:
                GC3.append(str(seq[i - 1]))
        GC1_seq = Seq(''.join(GC1))
        GC2_seq = Seq(''.join(GC2))
        GC3_seq = Seq(''.join(GC3))

        return [round(calculate_gc(GC1_seq), 2), round(calculate_gc(GC2_seq), 2), round(calculate_gc(GC3_seq), 2)]


    fw = open(out_results, 'w')
    fw.write("%s\t%s\t%s\t%s\t%s\n" % ("gene_name", "GC1", "GC2", "GC3", "GC_total"))

    for rec in SeqIO.parse(in_fasta, 'fasta'):
        GC_total = round(calculate_gc(rec.seq), 2)
        GC123 = calculate_GC123(rec.seq)
        GC_first = GC123[0]
        GC_second = GC123[1]
        GC_third = GC123[2]

        fw.write("%s\t%s\t%s\t%s\t%s\n" % (rec.id, str(GC_first), str(GC_second), str(GC_third), str(GC_total)))

    fw.close()