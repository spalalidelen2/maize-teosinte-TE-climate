"""Genomic distribution of the maize SNP-temperature signal.

Per-SNP correlation with bio1 (unconditioned), FDR < 0.05 counts by
chromosome, and the fraction of the chromosome-4 hits that fall inside the
inv4m inversion (chr4:171.77-185.95 Mb, B73 RefGen_V4; Calfee et al. 2021).
"""
from collections import Counter
import numpy as np
from common import SNPMatrix, climate, perlocus_p, bh_fdr

snp = SNPMatrix("maize")
bio1 = np.array([climate()[p]["bio1"] for p in snp.pops])
p = perlocus_p(snp.Y, bio1, None)
q = bh_fdr(p)
sig = q < 0.05
print("maize SNPs at FDR < 0.05: %d of %d" % (sig.sum(), len(p)))
cc = Counter(snp.chrom[sig])
print("by chromosome:", dict(sorted(cc.items(), key=lambda x: int(x[0]))))
c4 = (snp.chrom == "4") & sig
inv4m = c4 & (snp.pos >= 171_771_502) & (snp.pos < 185_951_149)
print("chromosome 4: %d significant SNPs, of which %d (%.0f%%) lie inside inv4m (chr4:171.8-186.0 Mb)"
      % (c4.sum(), inv4m.sum(), 100 * inv4m.sum() / max(c4.sum(), 1)))
print("chromosome-4 significant SNPs span %.1f-%.1f Mb" % (snp.pos[c4].min() / 1e6, snp.pos[c4].max() / 1e6))
blk = snp.chr9_block
print("chr9:120-150 Mb: %d significant SNPs" % (sig & blk).sum())
