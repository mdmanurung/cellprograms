
## DEV-014 — MR-2 substrate rebuild: Y scale correction (2026-10-23)

The original MR-2 semi_base was lost with the dead sandbox. The rebuild from
E-MTAB-10026 (ArrayExpress, h5ad md5 35428b921da91b1205d8025f31740c15) first
stored raw counts as Y; the sim_cmp.R contract is Y = log2(CPM + 1) (the
ambient scenario inverts via 2^Y - 1; planting amplitude is amp x msd on the
log scale). The raw-count substrate made planted signals vanishingly small
and crashed the ambient transform (2^counts = Inf; 12/96 datasets ERR).
Fixed (commit da5a932 + Ylog patch); all fits from the raw-count base were
discarded. Additional rebuild deviations, all disclosed: (a) corrupt gzip
block in X/indices [550M, 560M) zeroed (data = 0); (b) 192 AB_* CITE-seq ADT
columns per type dropped (centered negative protein values, not RNA counts);
(c) HSC (15 donors) and RBC (17 donors) dropped: MIN_DONORS = 32 required by
the C4 residual-df guard (added 2026-10-16, after the v1.2 run) at the
registered max_factors = 30; (d) HSC fine subtypes mapped to the HSC major
label before the donor filter. Session ended before the corrected-base cache
completed; partial results only (see MR-2 section of the final report).
