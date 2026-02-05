### for super tier cell typing

# =============================================================================
# Broad Cell Type Markers (Non-Neuronal)
# =============================================================================

# Astrocytes
mAstro <- unique(c("AGT", "ALDOC", "AQP4", "GFAP", "ALDH1L1", "SLC1A2"))

# Oligodendrocytes
mOligo <- unique(c("PLP1", "OPALIN", "MOG", "MOBP", "MBP", "CNP"))

# Oligodendrocyte Precursor Cells (OPC)
mOPC <- unique(c("PDGFRA", "OLIG2", "CSPG4"))

# Microglia
mMicro <- unique(c("PTPRC", "C1QC", "HEXB", "CTSS", "P2RY12", "CX3CR1"))

# Ependymal Cells
mEpendy <- unique(c("FOXJ1", "DTHD1", "RASSF9", "ADGB", "CFAP126",
                    "TMEM212", "DNAH5")) 

# Vascular Cells (including Endothelial)
mVas <- unique(c("RGS5", "CLDN5", "FLT1", "PECAM1"))

# Fibroblasts
mFibro <- unique(c("LUM", "DCN", "COL1A2", "COL1A1", "VIM"))

# =============================================================================
# Neuronal Markers (with refined subtypes)
# =============================================================================

# Pan-neuronal
mNeurons <- unique(c("SNAP25", "SYT1"))

# General Excitatory Neurons (Glutamatergic)
mExN <- unique(c("SLC17A7", "GRIA1", "GRIN2B", "SLC17A6", "TCF7L2"))

# General Inhibitory Neurons (GABAergic)
mInN <- unique(c("PVALB", "LHX6", "DLX1", "CCK", "ERBB4"))

# Thalamic Excitatory Neurons
mTEx <- unique(c("RORB", "FOXP2", "TBR1", "BCL11B", "SATB2", "POU4F1", "CALB1"))

# Excitatory Neuron subtype 1 (EN1)
mEN1 <- unique(c("NTS", "S100A6", "PROX1", "CALB2", "FOXP2",
                 "ROBO2", "KALRN", "RBFOX1"))

# Excitatory Neuron subtype 2 (EN2)
mEN2 <- unique(c("PENK", "CRTAC1", "RNF220", "NECAB1"))

# Midbrain-derived Neurons
mMidbrainN <- unique(c("OTX2", "LMO1"))

# TNR Neurons 
mTNR_N <- unique(c("ISL1", "GAD2", "PAX6", "SIX3", "SST", "EBF1", "EBF2", "EBF3"))

# Pretectum-derived Neurons
mPretectumN <- unique(c("ADGRL2", "ISLR2"))

# Telencephalic-derived Neurons
mTelencephalonN <- unique(c("FOXG1", "VIP", "PDZRN3", "NPY", "CRABP1", "MAF", "RARB", "NR2F2"))
