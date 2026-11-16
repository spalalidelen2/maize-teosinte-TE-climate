"""Genomic distribution of the teosinte SNP-temperature signal.

Counts SNPs significant at FDR < 0.05 (partial correlation with bio1
controlling for ancestry), their distribution over chromosomes and 5 Mb
windows, and the enrichment of significant SNPs in the chr9:120-150 Mb
candidate region and in the known inversion inv4m (hypergeometric test).
"""
from collections import Counter
import numpy as np
from scipy.stats import hypergeom
from common import SNPMatrix, climate, teosinte_ancestry, perlocus_p, bh_fdr

snp = SNPMatrix("teosinte")
pops = snp.pops
clim = climate()
bio1 = np.array([clim[p]["bio1"] for p in pops])
anc = teosinte_ancestry()
A = np.array([anc[p] for p in pops]).reshape(-1, 1)

p = perlocus_p(snp.Y, bio1, A)
q = bh_fdr(p)
sig = q < 0.05
print("SNPs significant at FDR < 0.05: %d of %d" % (sig.sum(), len(p)))
cc = Counter(snp.chrom[sig])
print("by chromosome:", dict(sorted(cc.items(), key=lambda x: int(x[0]))))

W = 5_000_000
cnt, sigcnt = Counter(), Counter()
for c, ps, s in zip(snp.chrom, snp.pos, sig):
    w = (c, ps // W)
    cnt[w] += 1
    if s:
        sigcnt[w] += 1
print("\n5 Mb windows with the most significant SNPs (chr : Mb : significant/total):")
for (c, b), sc in sigcnt.most_common(12):
    print("  chr%s : %d-%d Mb : %d/%d" % (c, b * 5, b * 5 + 5, sc, cnt[(c, b)]))

print("\nStrongest eight SNPs:")
for idx in np.argsort(p)[:8]:
    print("  chr%s:%s  p = %.2g" % (snp.chrom[idx], format(snp.pos[idx], ","), p[idx]))


def enrich(name, c, lo, hi):
    reg = (snp.chrom == c) & (snp.pos >= lo) & (snp.pos < hi)
    n_reg, s_reg = reg.sum(), (sig & reg).sum()
    exp = n_reg * sig.sum() / len(sig)
    pv = hypergeom.sf(s_reg - 1, len(sig), sig.sum(), n_reg)
    print("  %s (chr%s:%.0f-%.0f Mb): %d significant / %d SNPs (expected %.1f), enrichment %.1fx, p = %.2g"
          % (name, c, lo / 1e6, hi / 1e6, s_reg, n_reg, exp, s_reg / max(exp, 1e-9), pv))


print("\nEnrichment of significant SNPs in known regions:")
enrich("inv4m", "4", 171_771_502, 185_951_149)
enrich("chr9 candidate region", "9", 120e6, 150e6)
enrich("chr2 control window", "2", 120e6, 150e6)
n9 = (sig & (snp.chrom == "9")).sum()
print("\nchromosome 9 carries %d of %d significant SNPs (%.0f%%); chr9:120-150 Mb carries %d (%.0f%%)"
      % (n9, sig.sum(), 100 * n9 / sig.sum(), (sig & snp.chr9_block).sum(), 100 * (sig & snp.chr9_block).sum() / sig.sum()))
