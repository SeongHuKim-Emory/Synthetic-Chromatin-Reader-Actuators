## README.md

### Repository title
Synthetic-Chromatin-Reader-Actuators

### Associated manuscript
Decoding Gene Responsiveness to Synthetic Chromatin Reader-Actuators with Multi-Modal Epigenomic Profiling

### Journal
NPJ Systems Biology and Applications

### Manuscript status
Submitted

### Authors
Seong Hu Kim

Isioma Enwerem-Lackland

Natecia L. Williams

Rachel Fisher

Christopher Plaisier

Karmella A. Haynes

### Corresponding contact
Karmella A. Haynes
kahayne@emory.edu

### Keywords
epigenome engineering, polycomb, chromatin, gene regulation, topologically associating domains, machine learning, breast cancer

### Repository URL
github.com/SeongHuKim-Emory/Synthetic-Chromatin-Reader-Actuators

### Last updated
March 8th, 2026

### 1. Overview

This repository contains code used for machine learning, classification, and downstream analyses associated with the manuscript "Decoding Gene Responsiveness to Synthetic Chromatin Reader-Actuators with Multi-Modal Epigenomic Profiling."

The repository includes scripts required to reproduce the computational analyses and associated figures presented in the study.

The study combines ChIP-seq, time-course RNA-seq, public bioinformatics datasets, and machine learning analyses to investigate how chromatin context shapes responsiveness to synthetic reader-actuators in transgenic MCF7 breast cancer cells.


### 2. Repository contents

Directory layout is shown below.
```text
├───public
│   ├───ChIPseq
│   │       SHK_2024_06_11_MCF7_DBN021_0_0_ug_mL_Rep1.broadPeak
│   │       SHK_2024_06_11_MCF7_DBN021_0_0_ug_mL_Rep2.broadPeak
│   │       SHK_2024_06_11_MCF7_DBN021_0_0_ug_mL_Rep2_ppois.bigWig
│   │       SHK_2024_06_11_MCF7_DBN021_0_1_ug_mL_Rep1.broadPeak
│   │       SHK_2024_06_11_MCF7_DBN021_0_1_ug_mL_Rep2.broadPeak
│   │       SHK_2024_06_11_MCF7_DBN021_0_1_ug_mL_Rep2_ppois.bigWig
│   │       SHK_2024_06_11_MCF7_DBN021_0_5_ug_mL_Rep1.broadPeak
│   │       SHK_2024_06_11_MCF7_DBN021_0_5_ug_mL_Rep2.broadPeak
│   │       SHK_2024_06_11_MCF7_DBN021_0_5_ug_mL_Rep2_ppois.bigWig
│   │       SHK_2024_06_11_MCF7_DBN021_1_0_ug_mL_Rep1.broadPeak
│   │       SHK_2024_06_11_MCF7_DBN021_1_0_ug_mL_Rep2.broadPeak
│   │       SHK_2024_06_11_MCF7_DBN021_1_0_ug_mL_Rep2_ppois.bigWig
│   │       SHK_2024_06_11_MCF7_DBN025_0_5_ug_mL_Rep1.broadPeak
│   │       SHK_2024_06_11_MCF7_DBN025_0_5_ug_mL_Rep2.broadPeak
│   │       SHK_2024_06_11_MCF7_DBN025_1_0_ug_mL_Rep1.broadPeak
│   │       SHK_2024_06_11_MCF7_DBN025_1_0_ug_mL_Rep2.broadPeak
│   │       SHK_2024_06_11_MCF7_DBN025_1_0_ug_mL_Rep2_ppois.bigWig
│   │
│   ├───Multiomics
│   │   ├───Gene Enhancers
│   │   ├───Step18_ML_Prediction
│   │   ├───Step18_ML_Prediction_H1
│   │   ├───Step18_ML_Prediction_K562
│   │   ├───Step18_ML_Prediction_SRA
│   │   ├───Step19_ML_UpDEG_Prediction_SRA
│   │   └───Step26_TAD_Insulation_Score
│   │
│   ├───Public Dataset
│   │   ├───ATACseq
│   │   │       MCF7 ATACseq GSE169929 ENCFF976UNK hg38.bigWig
│   │   │
│   │   ├───ChIPseq
│   │   │   │   ENCSR636HFF GRCh38_unified_blacklist.bed
│   │   │   │   GSE174945 ENCFF195FSD hg38 MCF7 H3K36me3 ChIPseq.bed
│   │   │   │   GSE86714 ENCFF991HJA hg38 MCF7 H3K4me1 ChIPseq.bed
│   │   │   │   GSE95898 ENCFF348DEB hg38 MCF7 H3K9ac ChIPseq.bed
│   │   │   │   GSE96283 ENCFF714DEQ hg38 MCF7 H4K20me1 ChIPseq.bed
│   │   │   │   GSE96352 ENCFF491LQY hg38 MCF7 H3K27ac ChIPseq.bed
│   │   │   │   GSE96363 ENCFF669NUD hg38 MCF7 H3K27me3 ChIPseq.bed
│   │   │   │   GSE96439 ENCFF188VRU hg38 MCF7 H3K4me2 ChIPseq.bed
│   │   │   │   GSE96506 ENCFF268RXB hg38 MCF7 H3K4me3 ChIPseq.bed
│   │   │   │   GSE96517 ENCFF501UHK hg38 MCF7 H3K9me3 ChIPseq.bed
│   │   │   │   MCF7 ARID3A ChIPseq GSE91595 ENCFF823HUJ hg38.bed
│   │   │   │   MCF7 ARID3A ChIPseq GSE91595 ENCFF907DWT hg38.bigWig
│   │   │   │   MCF7 ATF7 ChIPseq GSE106025 ENCFF261DZR hg38.bed
│   │   │   │   MCF7 ATF7 ChIPseq GSE106025 ENCFF703ZRP hg38.bigWig
│   │   │   │   MCF7 BMI1 ChIPseq GSE105933 ENCFF846CPP hg38.bed
│   │   │   │   MCF7 BMI1 ChIPseq GSE105933 ENCFF977TJK hg38.bigWig
│   │   │   │   MCF7 CEBPB ChIPseq GSM1010889 ENCFF240ZHW hg38.bed
│   │   │   │   MCF7 CEBPB ChIPseq GSM1010889 ENCFF464JOC hg38.bigWig
│   │   │   │   MCF7 CHD1 ChIPseq GSE91698 ENCFF587BZB hg38.bed
│   │   │   │   MCF7 CHD1 ChIPseq GSE91698 ENCFF923MDX hg38.bigWig
│   │   │   │   MCF7 CLOCK ChIPseq GSE127438 ENCFF045RXD hg38.bigWig
│   │   │   │   MCF7 CLOCK ChIPseq GSE127438 ENCFF594YBN hg38.bed
│   │   │   │   MCF7 CLOCK ChIPseq GSE127640 ENCFF065FJX hg38.bigWig
│   │   │   │   MCF7 CLOCK ChIPseq GSE127640 ENCFF622TWI hg38.bed
│   │   │   │   MCF7 COPS2 ChIPseq GSE105356 ENCFF102ABZ hg38.bigWig
│   │   │   │   MCF7 COPS2 ChIPseq GSE105356 ENCFF395PJC hg38.bed
│   │   │   │   MCF7 CREB1 ChIPseq GSE105525 ENCFF599FYT hg38.bed
│   │   │   │   MCF7 CREB1 ChIPseq GSE105525 ENCFF954WHX hg38.bigWig
│   │   │   │   MCF7 CREB1 ChIPseq GSE106068 ENCFF030DLI hg38.bed
│   │   │   │   MCF7 CREB1 ChIPseq GSE106068 ENCFF958DMB hg38.bigWig
│   │   │   │   MCF7 CSDE1 ChIPseq GSE105291 ENCFF576UBB hg38.bed
│   │   │   │   MCF7 CSDE1 ChIPseq GSE105291 ENCFF722GUE hg38.bigWig
│   │   │   │   MCF7 CTBP1 ChIPseq GSE91938 ENCFF882ZZZ hg38.bed
│   │   │   │   MCF7 CTBP1 ChIPseq GSE91938 ENCFF922FKY hg38.bigWig
│   │   │   │   MCF7 CTCF ChIPseq GSE123219 ENCFF278FNP hg38.bed
│   │   │   │   MCF7 CTCF ChIPseq GSE123219 ENCFF507CRU hg38.bigWig
│   │   │   │   MCF7 CTCF ChIPseq GSM1010734 ENCFF596MGE hg38.bed
│   │   │   │   MCF7 CTCF ChIPseq GSM1022663 ENCFF203RMD hg38.bigWig
│   │   │   │   MCF7 CTCF ChIPseq GSM1022663 ENCFF948CYD hg38.bed
│   │   │   │   MCF7 CTCF ChIPseq GSM822305 ENCFF603QKR hg38.bigWig
│   │   │   │   MCF7 CUX1 ChIPseq GSE91415 ENCFF184PHK hg38.bed
│   │   │   │   MCF7 CUX1 ChIPseq GSE91415 ENCFF663YXL hg38.bigWig
│   │   │   │   MCF7 DDX20 ChIPseq GSE105517 ENCFF565XMS hg38.bed
│   │   │   │   MCF7 DDX20 ChIPseq GSE105517 ENCFF895IXP hg38.bigWig
│   │   │   │   MCF7 DPF2 ChIPseq GSE91598 ENCFF142KHL hg38.bigWig
│   │   │   │   MCF7 DPF2 ChIPseq GSE91598 ENCFF870ORM hg38.bed
│   │   │   │   MCF7 E2F8 ChIPseq GSE127615 ENCFF096JXF hg38.bed
│   │   │   │   MCF7 E2F8 ChIPseq GSE127615 ENCFF919YCR hg38.bigWig
│   │   │   │   MCF7 E4F1 ChIPseq GSE127582 ENCFF434IDL hg38.bed
│   │   │   │   MCF7 E4F1 ChIPseq GSE127582 ENCFF555HAI hg38.bigWig
│   │   │   │   MCF7 EGR1 ChIPseq GSM1010844 ENCFF427RSA hg38.bed
│   │   │   │   MCF7 EGR1 ChIPseq GSM1010844 ENCFF487DKA hg38.bigWig
│   │   │   │   MCF7 ELF1 ChIPseq GSE105503 ENCFF751TAA hg38.bigWig
│   │   │   │   MCF7 ELF1 ChIPseq GSE105503 ENCFF989DLF hg38.bed
│   │   │   │   MCF7 ELF1 ChIPseq GSM1010764 ENCFF091NQE hg38.bed
│   │   │   │   MCF7 ELF1 ChIPseq GSM1010764 ENCFF289KXH hg38.bigWig
│   │   │   │   MCF7 ELK1 ChIPseq GSE91713 ENCFF045XJM hg38.bigWig
│   │   │   │   MCF7 ELK1 ChIPseq GSE91713 ENCFF427XNV hg38.bed
│   │   │   │   MCF7 EP300 ChIPseq GSM1010800 ENCFF290DJX hg38.bed
│   │   │   │   MCF7 EP300 ChIPseq GSM1010800 ENCFF708NMR hg38.bigWig
│   │   │   │   MCF7 ESRRA ChIPseq GSE127591 ENCFF560XEP hg38.bed
│   │   │   │   MCF7 ESRRA ChIPseq GSE127591 ENCFF650XLA hg38.bigWig
│   │   │   │   MCF7 ESRRA ChIPseq GSE92187 ENCFF059ZSC hg38.bigWig
│   │   │   │   MCF7 ESRRA ChIPseq GSE92187 ENCFF214ICQ hg38.bed
│   │   │   │   MCF7 FOS ChIPseq GSE105734 ENCFF327LTA hg38.bigWig
│   │   │   │   MCF7 FOS ChIPseq GSE105734 ENCFF519VZH hg38.bed
│   │   │   │   MCF7 FOSL2 ChIPseq GSM1010768 ENCFF704HAX hg38.bigWig
│   │   │   │   MCF7 FOXA1 ChIPseq GSE105305 ENCFF112JVK hg38.bed
│   │   │   │   MCF7 FOXA1 ChIPseq GSE105305 ENCFF512UGW hg38.bigWig
│   │   │   │   MCF7 FOXK2 ChIPseq GSE91785 ENCFF222VDW hg38.bed
│   │   │   │   MCF7 FOXK2 ChIPseq GSE91785 ENCFF477VZD hg38.bigWig
│   │   │   │   MCF7 FOXM1 ChIPseq GSM1010769 ENCFF405RBL hg38.bigWig
│   │   │   │   MCF7 FOXM1 ChIPseq GSM1010769 ENCFF563MHZ hg38.bed
│   │   │   │   MCF7 GABPA ChIPseq GSM1010864 ENCFF494GAR hg38.bigWig
│   │   │   │   MCF7 GABPA ChIPseq GSM1010864 ENCFF678DJM hg38.bed
│   │   │   │   MCF7 GATA3 ChIPseq GSE127656 ENCFF217KJR hg38.bed
│   │   │   │   MCF7 GATA3 ChIPseq GSE127656 ENCFF971AZB hg38.bigWig
│   │   │   │   MCF7 GATA3 ChIPseq GSM1010783 ENCFF120BNO hg38.bigWig
│   │   │   │   MCF7 GATA3 ChIPseq GSM825711 ENCFF336XFC hg38.bed
│   │   │   │   MCF7 GATA3 ChIPseq GSM825711 ENCFF668UPC hg38.bigWig
│   │   │   │   MCF7 GATA3 ChIPseq GSM935445 ENCFF381ULU hg38.bed
│   │   │   │   MCF7 GATA3 ChIPseq GSM935445 ENCFF951CUN hg38.bigWig
│   │   │   │   MCF7 GATAD2B ChIPseq GSE105980 ENCFF304HIM hg38.bigWig
│   │   │   │   MCF7 GATAD2B ChIPseq GSE105980 ENCFF327XIO hg38.bed
│   │   │   │   MCF7 GATAD2B ChIPseq GSE91723 ENCFF401YUZ hg38.bed
│   │   │   │   MCF7 GATAD2B ChIPseq GSE91723 ENCFF560WGB hg38.bigWig
│   │   │   │   MCF7 GTF2F1 ChIPseq GSE91869 ENCFF736KUY hg38.bigWig
│   │   │   │   MCF7 GTF2F1 ChIPseq GSE91869 ENCFF819AEG hg38.bed
│   │   │   │   MCF7 H3K27ac ChIPseq GSE96352 ENCFF138YNG hg38.bigWig
│   │   │   │   MCF7 H3K27me3 ChIPseq GSE96363 ENCFF163QKN hg38.bigWig
│   │   │   │   MCF7 H3K36me3 ChIPseq GSE174945 ENCFF910BRP hg38.bigWig
│   │   │   │   MCF7 H3K4me1 ChIPseq GSE86714 ENCFF763NCP hg38.bigWig
│   │   │   │   MCF7 H3K4me2 ChIPseq GSE96439 ENCFF442RRY hg38.bigWig
│   │   │   │   MCF7 H3K4me3 ChIPseq GSE96506 ENCFF163MXP hg38.bigWig
│   │   │   │   MCF7 H3K9ac ChIPseq GSE95898 ENCFF327XJC hg38.bigWig
│   │   │   │   MCF7 H3K9me3 ChIPseq GSE96517 ENCFF481DZL hg38.bigWig
│   │   │   │   MCF7 H4K20me1 ChIPseq GSE96283 ENCFF366GLZ hg38.bigWig
│   │   │   │   MCF7 HCFC1 ChIPseq GSE91992 ENCFF103ABY hg38.bigWig
│   │   │   │   MCF7 HCFC1 ChIPseq GSE91992 ENCFF796EMB hg38.bed
│   │   │   │   MCF7 HDAC2 ChIPseq GSM1010825 ENCFF178WAQ hg38.bigWig
│   │   │   │   MCF7 HDAC2 ChIPseq GSM1010825 ENCFF661VXV hg38.bed
│   │   │   │   MCF7 HDGF ChIPseq GSE91567 ENCFF178LPH hg38.bigWig
│   │   │   │   MCF7 HDGF ChIPseq GSE91567 ENCFF480XAU hg38.bed
│   │   │   │   MCF7 HES1 ChIPseq GSE105287 ENCFF132MJX hg38.bed
│   │   │   │   MCF7 HES1 ChIPseq GSE105287 ENCFF509HTJ hg38.bigWig
│   │   │   │   MCF7 HSF1 ChIPseq GSE91444 ENCFF505ZBM hg38.bigWig
│   │   │   │   MCF7 HSF1 ChIPseq GSE91444 ENCFF766UTT hg38.bed
│   │   │   │   MCF7 JUN ChIPseq GSE91550 ENCFF730TVS hg38.bigWig
│   │   │   │   MCF7 JUN ChIPseq GSE91550 ENCFF769TUY hg38.bed
│   │   │   │   MCF7 JUND ChIPseq GSM1010892 ENCFF398OOD hg38.bed
│   │   │   │   MCF7 JUND ChIPseq GSM1010892 ENCFF990FGN hg38.bigWig
│   │   │   │   MCF7 LARP7 ChIPseq GSE105864 ENCFF148SFF hg38.bigWig
│   │   │   │   MCF7 LARP7 ChIPseq GSE105864 ENCFF734FBR hg38.bed
│   │   │   │   MCF7 MAFK ChIPseq GSE127366 ENCFF174HAP hg38.bigWig
│   │   │   │   MCF7 MAFK ChIPseq GSE127366 ENCFF425GWD hg38.bed
│   │   │   │   MCF7 MAX ChIPseq GSM1010863 ENCFF315FIO hg38.bigWig
│   │   │   │   MCF7 MAX ChIPseq GSM1010863 ENCFF624CRN hg38.bed
│   │   │   │   MCF7 MAZ ChIPseq GSE91633 ENCFF290CNC hg38.bed
│   │   │   │   MCF7 MAZ ChIPseq GSE91633 ENCFF297WGG hg38.bigWig
│   │   │   │   MCF7 MBD2 ChIPseq GSE96480 ENCFF517UQD hg38.bigWig
│   │   │   │   MCF7 MBD2 ChIPseq GSE96480 ENCFF853ULZ hg38.bed
│   │   │   │   MCF7 MLLT1 ChIPseq GSE91754 ENCFF654ISN hg38.bigWig
│   │   │   │   MCF7 MLLT1 ChIPseq GSE91754 ENCFF972FNO hg38.bed
│   │   │   │   MCF7 MNT ChIPseq GSE105652 ENCFF278RMF hg38.bigWig
│   │   │   │   MCF7 MNT ChIPseq GSE105652 ENCFF810GPR hg38.bed
│   │   │   │   MCF7 MNT ChIPseq GSE91968 ENCFF623IRB hg38.bigWig
│   │   │   │   MCF7 MNT ChIPseq GSE91968 ENCFF655RVZ hg38.bed
│   │   │   │   MCF7 MTA1 ChIPseq GSE91687 ENCFF070JYN hg38.bed
│   │   │   │   MCF7 MTA1 ChIPseq GSE91687 ENCFF883HBF hg38.bigWig
│   │   │   │   MCF7 MTA2 ChIPseq GSE91864 ENCFF556NWZ hg38.bed
│   │   │   │   MCF7 MTA2 ChIPseq GSE91864 ENCFF557XGU hg38.bigWig
│   │   │   │   MCF7 MTA3 ChIPseq GSE91727 ENCFF486PVE hg38.bigWig
│   │   │   │   MCF7 MTA3 ChIPseq GSE91727 ENCFF975UOJ hg38.bed
│   │   │   │   MCF7 NBN ChIPseq GSE91899 ENCFF065ASC hg38.bed
│   │   │   │   MCF7 NBN ChIPseq GSE91899 ENCFF242QKO hg38.bigWig
│   │   │   │   MCF7 NCOA3 ChIPseq GSE105681 ENCFF770MLI hg38.bigWig
│   │   │   │   MCF7 NCOA3 ChIPseq GSE105681 ENCFF778ZBT hg38.bed
│   │   │   │   MCF7 NCOA3 ChIPseq GSE105740 ENCFF342NSW hg38.bed
│   │   │   │   MCF7 NCOA3 ChIPseq GSE105740 ENCFF546FXD hg38.bigWig
│   │   │   │   MCF7 NEUROD1 ChIPseq GSE127346 ENCFF501HJF hg38.bigWig
│   │   │   │   MCF7 NEUROD1 ChIPseq GSE127346 ENCFF796KBR hg38.bed
│   │   │   │   MCF7 NFIB ChIPseq GSE105751 ENCFF079JBB hg38.bed
│   │   │   │   MCF7 NFIB ChIPseq GSE105751 ENCFF631HNZ hg38.bigWig
│   │   │   │   MCF7 NFIB ChIPseq GSE105805 ENCFF440SEE hg38.bigWig
│   │   │   │   MCF7 NFIB ChIPseq GSE105805 ENCFF489IWB hg38.bed
│   │   │   │   MCF7 NFRKB ChIPseq GSE105388 ENCFF411WFM hg38.bigWig
│   │   │   │   MCF7 NFRKB ChIPseq GSE105388 ENCFF532IPJ hg38.bed
│   │   │   │   MCF7 NFXL1 ChIPseq GSE105498 ENCFF142NAV hg38.bed
│   │   │   │   MCF7 NFXL1 ChIPseq GSE105498 ENCFF537TBX hg38.bigWig
│   │   │   │   MCF7 NONO ChIPseq GSE92159 ENCFF115WXJ hg38.bed
│   │   │   │   MCF7 NONO ChIPseq GSE92159 ENCFF875CIS hg38.bigWig
│   │   │   │   MCF7 NR2F2 ChIPseq GSM1010837 ENCFF386FQQ hg38.bed
│   │   │   │   MCF7 NR2F2 ChIPseq GSM1010837 ENCFF678MPN hg38.bigWig
│   │   │   │   MCF7 NRF1 ChIPseq GSE91522 ENCFF113MFU hg38.bigWig
│   │   │   │   MCF7 NRF1 ChIPseq GSE91522 ENCFF970MPY hg38.bed
│   │   │   │   MCF7 PAX8 ChIPseq GSE105903 ENCFF490MXJ hg38.bed
│   │   │   │   MCF7 PAX8 ChIPseq GSE105903 ENCFF775POG hg38.bigWig
│   │   │   │   MCF7 PKNOX1 ChIPseq GSE92210 ENCFF043NUC hg38.bed
│   │   │   │   MCF7 PKNOX1 ChIPseq GSE92210 ENCFF737NUU hg38.bigWig
│   │   │   │   MCF7 PML ChIPseq GSM1010838 ENCFF569JWM hg38.bed
│   │   │   │   MCF7 PML ChIPseq GSM1010838 ENCFF911FVQ hg38.bigWig
│   │   │   │   MCF7 POLR2A ChIPseq GSM822295 ENCFF235UTX hg38.bed
│   │   │   │   MCF7 POLR2A ChIPseq GSM822295 ENCFF827YIP hg38.bigWig
│   │   │   │   MCF7 PPP1R10 ChIPseq GSE105558 ENCFF167JMR hg38.bigWig
│   │   │   │   MCF7 PPP1R10 ChIPseq GSE105558 ENCFF428GSG hg38.bed
│   │   │   │   MCF7 RAD21 ChIPseq GSM1010791 ENCFF728MPA hg38.bed
│   │   │   │   MCF7 RAD21 ChIPseq GSM1010791 ENCFF775EKJ hg38.bigWig
│   │   │   │   MCF7 RAD51 ChIPseq GSE105597 ENCFF062KLW hg38.bed
│   │   │   │   MCF7 RAD51 ChIPseq GSE105597 ENCFF522WFI hg38.bigWig
│   │   │   │   MCF7 RCOR1 ChIPseq GSE91726 ENCFF160KBR hg38.bed
│   │   │   │   MCF7 RCOR1 ChIPseq GSE91726 ENCFF536KZN hg38.bigWig
│   │   │   │   MCF7 REST ChIPseq GSM1010891 ENCFF068TEY hg38.bigWig
│   │   │   │   MCF7 REST ChIPseq GSM1010891 ENCFF680JMZ hg38.bed
│   │   │   │   MCF7 RFX1 ChIPseq GSE91448 ENCFF310ZPD hg38.bigWig
│   │   │   │   MCF7 RFX1 ChIPseq GSE91448 ENCFF782YIG hg38.bed
│   │   │   │   MCF7 RFX1 ChIPseq GSE92059 ENCFF160KCP hg38.bed
│   │   │   │   MCF7 RFX1 ChIPseq GSE92059 ENCFF319WVD hg38.bigWig
│   │   │   │   MCF7 RFX5 ChIPseq GSE105874 ENCFF085XIV hg38.bigWig
│   │   │   │   MCF7 RFX5 ChIPseq GSE105874 ENCFF753LKY hg38.bed
│   │   │   │   MCF7 SIN3A ChIPseq GSE91789 ENCFF402DKZ hg38.bed
│   │   │   │   MCF7 SIN3A ChIPseq GSE91789 ENCFF807QKK hg38.bigWig
│   │   │   │   MCF7 SIN3A ChIPseq GSM1010862 ENCFF792DRZ hg38.bigWig
│   │   │   │   MCF7 SIN3A ChIPseq GSM1010862 ENCFF803WYI hg38.bed
│   │   │   │   MCF7 SIX4 ChIPseq GSE91630 ENCFF032MSQ hg38.bigWig
│   │   │   │   MCF7 SIX4 ChIPseq GSE91630 ENCFF728MNA hg38.bed
│   │   │   │   MCF7 SMARCA5 ChIPseq GSE105722 ENCFF297JNO hg38.bigWig
│   │   │   │   MCF7 SMARCA5 ChIPseq GSE105722 ENCFF526KAO hg38.bed
│   │   │   │   MCF7 SMARCE1 ChIPseq GSE105546 ENCFF375LRP hg38.bigWig
│   │   │   │   MCF7 SMARCE1 ChIPseq GSE105546 ENCFF889QAN hg38.bed
│   │   │   │   MCF7 SNIP1 ChIPseq GSE105188 ENCFF231LYV hg38.bigWig
│   │   │   │   MCF7 SNIP1 ChIPseq GSE105188 ENCFF285KBA hg38.bed
│   │   │   │   MCF7 SP1 ChIPseq GSE92014 ENCFF367HVY hg38.bigWig
│   │   │   │   MCF7 SP1 ChIPseq GSE92014 ENCFF932UJX hg38.bed
│   │   │   │   MCF7 SREBF1 ChIPseq GSE91561 ENCFF571VYR hg38.bed
│   │   │   │   MCF7 SREBF1 ChIPseq GSE91561 ENCFF791DFW hg38.bigWig
│   │   │   │   MCF7 SRF ChIPseq GSM1010839 ENCFF131TYZ hg38.bed
│   │   │   │   MCF7 SRF ChIPseq GSM1010839 ENCFF534AJF hg38.bigWig
│   │   │   │   MCF7 SUZ12 ChIPseq GSE105981 ENCFF488VNT hg38.bed
│   │   │   │   MCF7 SUZ12 ChIPseq GSE105981 ENCFF598JRC hg38.bigWig
│   │   │   │   MCF7 TAF1 ChIPseq GSM1010811 ENCFF075TIU hg38.bed
│   │   │   │   MCF7 TAF1 ChIPseq GSM1010811 ENCFF418WWV hg38.bigWig
│   │   │   │   MCF7 TARDBP ChIPseq GSE105812 ENCFF065RQM hg38.bigWig
│   │   │   │   MCF7 TARDBP ChIPseq GSE105812 ENCFF666QVW hg38.bed
│   │   │   │   MCF7 TCF12 ChIPseq GSM1010861 ENCFF031XVC hg38.bigWig
│   │   │   │   MCF7 TCF12 ChIPseq GSM1010861 ENCFF987QLI hg38.bed
│   │   │   │   MCF7 TCF7L2 ChIPseq GSM816438 ENCFF332AZX hg38.bed
│   │   │   │   MCF7 TCF7L2 ChIPseq GSM816438 ENCFF711CKB hg38.bigWig
│   │   │   │   MCF7 TEAD4 ChIPseq GSM1010860 ENCFF751VAZ hg38.bed
│   │   │   │   MCF7 TEAD4 ChIPseq GSM1010860 ENCFF977ZGZ hg38.bigWig
│   │   │   │   MCF7 TOE1 ChIPseq GSE105831 ENCFF318YIK hg38.bigWig
│   │   │   │   MCF7 TRIM22 ChIPseq GSE127607 ENCFF100OCC hg38.bigWig
│   │   │   │   MCF7 TRIM22 ChIPseq GSE127607 ENCFF238QAG hg38.bed
│   │   │   │   MCF7 YBX1 ChIPseq GSE92114 ENCFF033VON hg38.bigWig
│   │   │   │   MCF7 YBX1 ChIPseq GSE92114 ENCFF745UTS hg38.bed
│   │   │   │   MCF7 ZBTB1 ChIPseq GSE105482 ENCFF062ZHZ hg38.bigWig
│   │   │   │   MCF7 ZBTB1 ChIPseq GSE105482 ENCFF192BPQ hg38.bed
│   │   │   │   MCF7 ZBTB11 ChIPseq GSE95952 ENCFF328ZJV hg38.bed
│   │   │   │   MCF7 ZBTB11 ChIPseq GSE95952 ENCFF510SQD hg38.bigWig
│   │   │   │   MCF7 ZBTB33 ChIPseq GSE91596 ENCFF506CTQ hg38.bigWig
│   │   │   │   MCF7 ZBTB40 ChIPseq GSE91661 ENCFF160TJG hg38.bed
│   │   │   │   MCF7 ZBTB40 ChIPseq GSE91661 ENCFF539ZHK hg38.bigWig
│   │   │   │   MCF7 ZBTB7B ChIPseq GSE105418 ENCFF073COQ hg38.bigWig
│   │   │   │   MCF7 ZBTB7B ChIPseq GSE105418 ENCFF339BDC hg38.bed
│   │   │   │   MCF7 ZFX ChIPseq GSE105562 ENCFF482RYD hg38.bigWig
│   │   │   │   MCF7 ZFX ChIPseq GSE105562 ENCFF861DOL hg38.bed
│   │   │   │   MCF7 ZHX2 ChIPseq GSE96441 ENCFF443WEJ hg38.bed
│   │   │   │   MCF7 ZHX2 ChIPseq GSE96441 ENCFF724GPW hg38.bigWig
│   │   │   │   MCF7 ZKSCAN1 ChIPseq GSE91769 ENCFF239BWH hg38.bigWig
│   │   │   │   MCF7 ZKSCAN1 ChIPseq GSE91769 ENCFF856VYW hg38.bed
│   │   │   │   MCF7 ZNF207 ChIPseq GSE91475 ENCFF668FJE hg38.bigWig
│   │   │   │   MCF7 ZNF207 ChIPseq GSE91475 ENCFF761KBK hg38.bed
│   │   │   │   MCF7 ZNF217 ChIPseq GSE105662 ENCFF197CVM hg38.bed
│   │   │   │   MCF7 ZNF217 ChIPseq GSE105662 ENCFF491CXP hg38.bigWig
│   │   │   │   MCF7 ZNF217 ChIPseq GSE127625 ENCFF343WRL hg38.bigWig
│   │   │   │   MCF7 ZNF217 ChIPseq GSE127625 ENCFF607VIS hg38.bed
│   │   │   │   MCF7 ZNF217 ChIPseq GSM935563 ENCFF083LCM hg38.bed
│   │   │   │   MCF7 ZNF217 ChIPseq GSM935563 ENCFF092KFI hg38.bigWig
│   │   │   │   MCF7 ZNF24 ChIPseq GSE105531 ENCFF258VTL hg38.bed
│   │   │   │   MCF7 ZNF24 ChIPseq GSE105531 ENCFF505ORH hg38.bigWig
│   │   │   │   MCF7 ZNF444 ChIPseq GSE127367 ENCFF101ION hg38.bigWig
│   │   │   │   MCF7 ZNF444 ChIPseq GSE127367 ENCFF896QHA hg38.bed
│   │   │   │   MCF7 ZNF507 ChIPseq GSE127652 ENCFF659LTF hg38.bigWig
│   │   │   │   MCF7 ZNF507 ChIPseq GSE127652 ENCFF987CPY hg38.bed
│   │   │   │   MCF7 ZNF512B ChIPseq GSE105987 ENCFF320KEW hg38.bed
│   │   │   │   MCF7 ZNF512B ChIPseq GSE105987 ENCFF413BLF hg38.bigWig
│   │   │   │   MCF7 ZNF512B ChIPseq GSE127363 ENCFF692RAY hg38.bigWig
│   │   │   │   MCF7 ZNF512B ChIPseq GSE127363 ENCFF760TRJ hg38.bed
│   │   │   │   MCF7 ZNF574 ChIPseq GSE127638 ENCFF281CNU hg38.bed
│   │   │   │   MCF7 ZNF574 ChIPseqGSE127638 ENCFF914VHT hg38.bigWig
│   │   │   │   MCF7 ZNF579 ChIPseq GSE105224 ENCFF223BRJ hg38.bed
│   │   │   │   MCF7 ZNF579 ChIPseq GSE105224 ENCFF913XOW hg38.bigWig
│   │   │   │   MCF7 ZNF592 ChIPseq GSE91425 ENCFF566KDV hg38.bed
│   │   │   │   MCF7 ZNF592 ChIPseq GSE91425 ENCFF950SDY hg38.bigWig
│   │   │   │   MCF7 ZNF592 ChIPseq GSE91993 ENCFF271UJG hg38.bigWig
│   │   │   │   MCF7 ZNF592 ChIPseq GSE91993 ENCFF645NLC hg38.bed
│   │   │   │   MCF7 ZNF687 ChIPseq GSE92150 ENCFF609HVM hg38.bed
│   │   │   │   MCF7 ZNF687 ChIPseq GSE92150 ENCFF910XUA hg38.bigWig
│   │   │   │   MCF7 ZNF8 ChIPseq GSE127545 ENCFF435UDQ hg38.bed
│   │   │   │   MCF7 ZNF8 ChIPseq GSE127545 ENCFF550FAH hg38.bigWig
│   │   │   │
│   │   │   └───Other Cell Lines
│   │   │       ├───H1
│   │   │       │       H1 CBX8 ChIPseq GSE123227 ENCFF483UZG hg38.bed
│   │   │       │       H1 H3K27ac ChIPseq GSM733718 ENCFF771GNB hg38.bigWig
│   │   │       │       H1 H3K27me3 ChIPseq GSM433167 ENCFF927FVH hg38.bigWig
│   │   │       │       H1 H3K36me3 ChIPseq GSM733725 ENCFF985CVI hg38.bigWig
│   │   │       │       H1 H3K4me1 ChIPseq GSM733782 ENCFF706CHK hg38.bigWig
│   │   │       │       H1 H3K4me2 ChIPseq GSM733670 ENCFF556VEB hg38.bigWig
│   │   │       │       H1 H3K4me3 ChIPseq GSE96392 ENCFF698DKQ hg38.bigWig
│   │   │       │       H1 H3K9ac ChIPseq GSM537685 ENCFF225ZYN hg38.bigWig
│   │   │       │       H1 H3K9me3 ChIPseq GSM1003585 ENCFF441SLB hg38.bigWig
│   │   │       │       H1 H4K20me1 ChIPseq GSM733687 ENCFF674EST hg38.bigWig
│   │   │       │
│   │   │       └───K562
│   │   │               K562 CBX8 ChIPseq GSM1003569 ENCFF522HZT hg38.bed
│   │   │               K562 H3K27ac ChIPseq GSM733656 ENCFF469JMR hg38.bigWig
│   │   │               K562 H3K27me3 ChIPseq GSM733658 ENCFF665RDD hg38.bigWig
│   │   │               K562 H3K36me3 ChIPseq GSM733714 ENCFF296TSL hg38.bigWig
│   │   │               K562 H3K4me1 ChIPseq GSM733692 ENCFF100FDI hg38.bigWig
│   │   │               K562 H3K4me2 ChIPseq GSM733651 ENCFF054RSU hg38.bigWig
│   │   │               K562 H3K4me3 ChIPseq GSE96303 ENCFF405ZDL hg38.bigWig
│   │   │               K562 H3K9ac ChIPseq GSM733778 ENCFF239EBH hg38.bigWig
│   │   │               K562 H3K9me3 ChIPseq GSM733776 ENCFF632NQA hg38.bigWig
│   │   │               K562 H4K20me1 ChIPseq GSM733675 ENCFF694ODT hg38.bigWig
│   │   │
│   │   ├───DNase_seq
│   │   │       MCF7 DNase_seq GSM1008565 ENCFF137FTO hg38.bigWig
│   │   │
│   │   ├───GeneHancer
│   │   │       GeneHancer_AnnotSV_elements_v5.25.txt
│   │   │       GeneHancer_AnnotSV_gene_association_scores_v5.25.txt
│   │   │       GeneHancer_AnnotSV_ReadMe.txt
│   │   │       GeneHancer_AnnotSV_tissues_v5.25.txt
│   │   │       GeneHancer_Genes_Elements.csv
│   │   │       GeneHancer_Protein_Coding_Genes_Elements.csv
│   │   │       GeneHancer_TFBSs_v5.25.txt
│   │   │       GeneHancer_TFBSs_v5.25_MCF7.csv
│   │   │       GeneHancer_Tissues_v5.25.txt
│   │   │       GeneHancer_v5.25.gff
│   │   │
│   │   ├───HiC
│   │   │   │   MCF7 HiC Contact Domains GSE237722 ENCFF164AGX hg38 Original.bedpe
│   │   │   │   MCF7 HiC Contact Domains GSE237722 ENCFF164AGX hg38.bedpe
│   │   │   │   MCF7 HiC Loops GSE237722 ENCFF797NNQ hg38 Original.bedpe
│   │   │   │   MCF7 HiC Loops GSE237722 ENCFF797NNQ hg38.bedpe
│   │   │   │
│   │   │   └───GSE66733_MCF7_HiC
│   │   │           hg19ToHg38.over.chain
│   │   │           HiCStein-MCF7-WT__hg19__chr10__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           HiCStein-MCF7-WT__hg19__chr11__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           HiCStein-MCF7-WT__hg19__chr12__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           HiCStein-MCF7-WT__hg19__chr13__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           HiCStein-MCF7-WT__hg19__chr14__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           HiCStein-MCF7-WT__hg19__chr15__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           HiCStein-MCF7-WT__hg19__chr16__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           HiCStein-MCF7-WT__hg19__chr17__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           HiCStein-MCF7-WT__hg19__chr18__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           HiCStein-MCF7-WT__hg19__chr19__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           HiCStein-MCF7-WT__hg19__chr1__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           HiCStein-MCF7-WT__hg19__chr20__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           HiCStein-MCF7-WT__hg19__chr21__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           HiCStein-MCF7-WT__hg19__chr22__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           HiCStein-MCF7-WT__hg19__chr2__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           HiCStein-MCF7-WT__hg19__chr3__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           HiCStein-MCF7-WT__hg19__chr4__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           HiCStein-MCF7-WT__hg19__chr5__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           HiCStein-MCF7-WT__hg19__chr6__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           HiCStein-MCF7-WT__hg19__chr7__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           HiCStein-MCF7-WT__hg19__chr8__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           HiCStein-MCF7-WT__hg19__chr9__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           HiCStein-MCF7-WT__hg19__chrX__C-40000-iced.is1000000.ids240000.insulation.boundaries
│   │   │           MCF7_hg19_hg38_TAD_Boundaries.csv
│   │   │           MCF7_hg19_TAD_Boundaries.csv
│   │   │           MCF7_hg38_TADs.csv
│   │   │           MCF7_hg38_TAD_Boundaries.bed
│   │   │
│   │   └───RNAseq
│   │           GSE175204 ENCFF721BRA hg38 MCF7 total RNAseq.tsv
│   │            
│   └───RNAseq
│           gexp_counts.csv
│
└───src
    │   hg19ToHg38.over.chain
    │   Step11_ChIPseq_Heatmaps.R
    │   Step12_Visualize_Genomic_Annotation.R
    │   Step14_RNAseq_Plots.R
    │   Step16_Linux_ChIPseq_Overlap.R
    │   Step1_GeneCountMatrix_AddGeneSymbol.R
    │   Step26_1_TAD_Insulation_Score.R
    │   Step26_2_TAD_Insulation_Score_Summary.R
    │   Step26_3_TAD_UpDEG_Category.R
    │   Step2_RNAseq_Sample_Analysis.R
    │   Step3_DEG_Categorization.R
    │   Step4_Public_HiC_Preprocessing.R
    │   Step5_ChIPseq_Preprocessing.R
    │   Step6_Multiomics.py
    │
    ├───Step18_ML_Prediction
    │       hg38.chrom.sizes
    │       Step18_1_ML_Prediction_hg38_bin.py
    │       Step18_2_ML_Prediction_pos_neg_sets.py
    │       Step18_3_ML_Prediction_histone_bigwig.py
    │       Step18_4_0_ML_Prediction_Model_Training.py
    │       Step18_4_1_ML_Prediction_Model_Training_H3K27me3_Only.py
    │       Step18_5_ML_Prediction_Test_Result_bed.py
    │       Step18_6_Custome_SHAP_Beeswarm.py
    │
    ├───Step18_ML_Prediction_H1
    │       hg38.chrom.sizes
    │       Step18_1_ML_Prediction_hg38_bin.py
    │       Step18_2_ML_Prediction_pos_neg_sets.py
    │       Step18_3_ML_Prediction_histone_bigwig.py
    │       Step18_4_0_ML_Prediction_Model_Training.py
    │       Step18_4_1_ML_Prediction_Model_Training_H3K27me3_Only.py
    │       Step18_5_ML_Prediction_Test_Result_bed.py
    │       Step18_6_Custome_SHAP_Beeswarm.py
    │
    ├───Step18_ML_Prediction_K562
    │       hg38.chrom.sizes
    │       Step18_1_ML_Prediction_hg38_bin.py
    │       Step18_2_ML_Prediction_pos_neg_sets.py
    │       Step18_3_ML_Prediction_histone_bigwig.py
    │       Step18_4_0_ML_Prediction_Model_Training.py
    │       Step18_4_1_ML_Prediction_Model_Training_H3K27me3_Only.py
    │       Step18_5_ML_Prediction_Test_Result_bed.py
    │       Step18_6_Custome_SHAP_Beeswarm.py
    │
    ├───Step18_ML_Prediction_SRA
    │       hg38.chrom.sizes
    │       Step18_1_ML_Prediction_hg38_bin.py
    │       Step18_2_ML_Prediction_pos_neg_sets.py
    │       Step18_3_ML_Prediction_histone_bigwig.py
    │       Step18_4_0_ML_Prediction_Model_Training.py
    │       Step18_4_1_ML_Prediction_Model_Training_H3K27me3_Only.py
    │       Step18_5_ML_Prediction_Test_Result_bed.py
    │       Step18_6_Custome_SHAP_Beeswarm.py
    │   
    └───Step19_ML_UpDEG_Prediction_SRA
            Step19_10_0_ML_UpDEG_Prediction_Feature_Matrix_TSS_SRA_Max_Pooling_Transcript_GO.py
            Step19_11_ML_UpDEG_Prediction_Feature_Matrix_TSS_SRA_Max_Pooling_Histone_Sets.py
            Step19_12_ML_UpDEG_Prediction_Feature_Filtering.py
            Step19_13_ML_UpDEG_Prediction_Final_Training.py
            Step19_1_ML_UpDEG_Prediction_Merge_Enhancer_RNAseq.py
            Step19_2_ML_UpDEG_Prediction_Enhancer_DEG_Category.py
            Step19_3_ML_UpDEG_Prediction_Feature_Matrix_TSS.py
            Step19_4_ML_UpDEG_Prediction_Feature_Matrix_TSS_CleanUp.py
            Step19_5_ML_UpDEG_Prediction_Feature_Matrix_TSS_SRA.py
            Step19_6_ML_UpDEG_Prediction_Feature_Matrix_TSS_SRA_Max_Pooling.py
            Step19_7_ML_UpDEG_Prediction_Feature_Matrix_TSS_SRA_Max_Pooling_Promoter.py
            Step19_8_ML_UpDEG_Prediction_Feature_Matrix_TSS_SRA_Max_Pooling_Histone.py
            Step19_9_ML_UpDEG_Prediction_Feature_Matrix_TSS_SRA_Max_Pooling_Transcript.py
```

Directory description:

- public/: all non-code files
- public/ChIPseq/: ChIP-seq experiment files
- public/Multiomics/: Output files from computational analyses
- public/Public Dataset/: Experiment files obtained from public datasets
- public/RNAseq/: RNA-seq experiment file
- src/: all source code files


### 3. What this repository reproduces

This repository is intended to reproduce the computational analyses and associated figures presented in the manuscript.

Main figures reported in the manuscript:
- Figure 2B, 2C
- Figure 3A
- Figure 4
- Figure 5C
- Figure 6A
- Figure 7

Supplementary items mentioned in the manuscript:
- Supplemental Figure S3
- Supplemental Figure S4


### 4. Software environment

Analyses reported in the manuscript used the following software and package versions where stated.

Software package versions:
- Cutadapt v4.1
- STAR 2.7.1
- htseq-count v0.11.1
- DESeq2 v1.44.0
- ggplot2 v4.0.0
- GenomicRanges v1.56.2
- regioneR v1.36.0
- ChIPseeker v1.40.0
- rtracklayer v1.64.0
- pyBigWig v0.3.24
- XGBoost v3.0.1
- LightGBM v4.6.0
- CatBoost v1.2.8
- scikit-learn v1.7.2
- GraphPad v10.6.1

Operating system:
- Windows 11 Version 25H2 (OS Build 26200.7840)
- Linux 24.04.3 LTS

Python version:
Python 3.12.0

R version:
R 4.4.0


### 5. Data availability

Gene-level raw count matrices and differential expression results generated in this study are being deposited in the NCBI Gene Expression Omnibus and will be publicly available upon publication; accession is pending.

RNA-seq raw sequencing files are unavailable.

The processed data are sufficient to reproduce all RNA-seq analyses presented in the manuscript.

ChIP-seq datasets are being deposited in GEO and will be made publicly available upon publication.

Publicly available datasets analyzed in the study are listed in the Methods and Supplemental sections of the manuscript including:
- ENCODE histone ChIP-seq datasets for MCF7 cells
- ENCODE transcription factor ChIP-seq datasets for MCF7 cells
- ENCODE DNase-seq and ATAC-seq datasets for MCF7 cells
- ENCODE MCF7 total RNA-seq dataset ENCFF721BRA
- MCF7 Hi-C data GSE66733
- GeneHancer v5.25 enhancer-gene interaction database


Access notes:
- GEO accession for newly generated processed RNA-seq data: pending publication
- GEO accession for newly generated ChIP-seq data: pending publication
- RNA-seq raw sequencing files: unavailable


### 6. Installation.

Document the actual repository setup here, for example:
1. Clone the repository.
2. Replicate the directory indicated in Repository contents
3. Download RNA-seq, ChIP-seq, and public dataset files then place the files as shown in the Repository contents
4. Run the analysis scripts in the order listed below.


### 7. Execution order

All the source code filenames include step numberings. 
Follow the step numbers. 

Missing steps, e.g. step 7, was not used in the final manuscript


### 8. Figure reproduction map

Figure 2: Genomic distributions and chromatin feature overlaps for SRA and PCD-RFP ChIP-seq peaks.

Figure 2B
- Step1_GeneCountMatrix_AddGeneSymbol.R
- Step2_RNAseq_Sample_Analysis.R
- Step3_DEG_Categorization.R
- Step5_ChIPseq_Preprocessing.R
- Step11_ChIPseq_Heatmaps.R
- Using .csv outputs, plotted heatmap using GraphPad Psism (v10.6.1)

Figure 2C
- Step5_ChIPseq_Preprocessing.R
- Step16_Linux_ChIPseq_Overlap.R
- Use .jpg outputs directly or use .csv files and plot using GraphPad Psism (v10.6.1)

Figure 3: Annotations of regions with PCD-fusion ChIP-seq peaks.

Figure 3A
- Step12_Visualize_Genomic_Annotation.R
- Output .jpg files 

Figure 4: XGBoost prediction of PCD-fusion genome-wide binding.

Figure 4A, 4B
- Step5_ChIPseq_Preprocessing.R
- Step18_ML_Prediction/Step18_1_ML_Prediction_hg38_bin.py
- Step18_ML_Prediction/Step18_2_ML_Prediction_pos_neg_sets.py
- Step18_ML_Prediction/Step18_3_ML_Prediction_histone_bigwig.py
- Step18_ML_Prediction/Step18_4_0_ML_Prediction_Model_Training.py
- Use .jpg outputs for PCD-RFP

- Step5_ChIPseq_Preprocessing.R
- Step18_ML_Prediction_SRA/Step18_1_ML_Prediction_hg38_bin.py
- Step18_ML_Prediction_SRA/Step18_2_ML_Prediction_pos_neg_sets.py
- Step18_ML_Prediction_SRA/Step18_3_ML_Prediction_histone_bigwig.py
- Step18_ML_Prediction_SRA/Step18_4_0_ML_Prediction_Model_Training.py
- Use .jpg outputs for SRA

Figure 4C
- Step5_ChIPseq_Preprocessing.R
- Step18_ML_Prediction_SRA/Step18_1_ML_Prediction_hg38_bin.py
- Step18_ML_Prediction_SRA/Step18_2_ML_Prediction_pos_neg_sets.py
- Step18_ML_Prediction_SRA/Step18_3_ML_Prediction_histone_bigwig.py
- Step18_ML_Prediction_SRA/Step18_4_0_ML_Prediction_Model_Training.py
- Use .jpg outputs

Figure 4D
- Step5_ChIPseq_Preprocessing.R
- Step18_ML_Prediction_SRA/Step18_1_ML_Prediction_hg38_bin.py
- Step18_ML_Prediction_SRA/Step18_2_ML_Prediction_pos_neg_sets.py
- Step18_ML_Prediction_SRA/Step18_3_ML_Prediction_histone_bigwig.py
- Step18_ML_Prediction_SRA/Step18_4_0_ML_Prediction_Model_Training.py
- Step18_ML_Prediction_SRA/Step18_6_Custome_SHAP_Beeswarm.py
- Use .jpg outputs

Figure 4E
- Step5_ChIPseq_Preprocessing.R
- Step18_ML_Prediction_SRA/Step18_1_ML_Prediction_hg38_bin.py
- Step18_ML_Prediction_SRA/Step18_2_ML_Prediction_pos_neg_sets.py
- Step18_ML_Prediction_SRA/Step18_3_ML_Prediction_histone_bigwig.py
- Step18_ML_Prediction_SRA/Step18_4_0_ML_Prediction_Model_Training.py
- Step18_ML_Prediction_SRA/Step18_4_1_ML_Prediction_Model_Training_H3K27me3_Only.py
- Change H3K27me3 to other histone marks for the complete table in figure 4E
- Use Google Spreadsheet to tabulate the results 

Figure 5: Time-course transcriptome profiling of MCF7 cells with PCD-fusion proteins.

Figure 5C
- Step1_GeneCountMatrix_AddGeneSymbol.R
- Step2_RNAseq_Sample_Analysis.R
- Step3_DEG_Categorization.R
- Step14_RNAseq_Plots.R
- Output .jpg files

Figure 6: UpDEG and ChIP-seq peak distribution within topologically insulated regions.

Figure 6A
- Step1_GeneCountMatrix_AddGeneSymbol.R
- Step2_RNAseq_Sample_Analysis.R
- Step3_DEG_Categorization.R
- Step4_Public_HiC_Preprocessing.R
- Step5_ChIPseq_Preprocessing.R
- Step6_Multiomics.py
- Step26_1_TAD_Insulation_Score.R
- Step26_2_TAD_Insulation_Score_Summary.R
- Step26_3_TAD_UpDEG_Category.R
- Use output .csv file to plot the Venn Diagram

Figure 7: Predictive features at UpDEGs vs. non-DEGs with enhancers that are bound by PCD-RFP.

Figure 7C, 7D, 7E, 7F
- Step1_GeneCountMatrix_AddGeneSymbol.R
- Step5_ChIPseq_Preprocessing.R
- Step19_ML_UpDEG_Prediction_SRA/Step19_1_ML_UpDEG_Prediction_Merge_Enhancer_RNAseq.py
- Step19_ML_UpDEG_Prediction_SRA/Step19_2_ML_UpDEG_Prediction_Enhancer_DEG_Category.py
- Step19_ML_UpDEG_Prediction_SRA/Step19_3_ML_UpDEG_Prediction_Feature_Matrix_TSS.py
- Step19_ML_UpDEG_Prediction_SRA/Step19_4_ML_UpDEG_Prediction_Feature_Matrix_TSS_CleanUp.py
- Step19_ML_UpDEG_Prediction_SRA/Step19_5_ML_UpDEG_Prediction_Feature_Matrix_TSS_SRA.py
- Step19_ML_UpDEG_Prediction_SRA/Step19_6_ML_UpDEG_Prediction_Feature_Matrix_TSS_SRA_Max_Pooling.py
- Step19_ML_UpDEG_Prediction_SRA/Step19_7_ML_UpDEG_Prediction_Feature_Matrix_TSS_SRA_Max_Pooling_Promoter.py
- Step19_ML_UpDEG_Prediction_SRA/Step19_8_ML_UpDEG_Prediction_Feature_Matrix_TSS_SRA_Max_Pooling_Histone.py
- Step19_ML_UpDEG_Prediction_SRA/Step19_9_ML_UpDEG_Prediction_Feature_Matrix_TSS_SRA_Max_Pooling_Transcript.py
- Step19_ML_UpDEG_Prediction_SRA/Step19_10_0_ML_UpDEG_Prediction_Feature_Matrix_TSS_SRA_Max_Pooling_Transcript_GO.py
- Step19_ML_UpDEG_Prediction_SRA/Step19_11_ML_UpDEG_Prediction_Feature_Matrix_TSS_SRA_Max_Pooling_Histone_Sets.py
- Step19_ML_UpDEG_Prediction_SRA/Step19_12_ML_UpDEG_Prediction_Feature_Filtering.py
- Step19_ML_UpDEG_Prediction_SRA/Step19_13_ML_UpDEG_Prediction_Final_Training.py
- Use output .jpg files

Supplementary Figure S3: PCD-RFP SHAP-related supplemental analysis referenced in the Figure 4 section.

Figure S3A
- Step5_ChIPseq_Preprocessing.R
- Step18_ML_Prediction/Step18_1_ML_Prediction_hg38_bin.py
- Step18_ML_Prediction/Step18_2_ML_Prediction_pos_neg_sets.py
- Step18_ML_Prediction/Step18_3_ML_Prediction_histone_bigwig.py
- Step18_ML_Prediction/Step18_4_0_ML_Prediction_Model_Training.py
- Use .jpg outputs

Figure S3B
- Step5_ChIPseq_Preprocessing.R
- Step18_ML_Prediction/Step18_1_ML_Prediction_hg38_bin.py
- Step18_ML_Prediction/Step18_2_ML_Prediction_pos_neg_sets.py
- Step18_ML_Prediction/Step18_3_ML_Prediction_histone_bigwig.py
- Step18_ML_Prediction/Step18_4_0_ML_Prediction_Model_Training.py
- Step18_ML_Prediction/Step18_6_Custome_SHAP_Beeswarm.py
- Use .jpg outputs

Figure S3C
- Step5_ChIPseq_Preprocessing.R
- Step18_ML_Prediction/Step18_1_ML_Prediction_hg38_bin.py
- Step18_ML_Prediction/Step18_2_ML_Prediction_pos_neg_sets.py
- Step18_ML_Prediction/Step18_3_ML_Prediction_histone_bigwig.py
- Step18_ML_Prediction/Step18_4_0_ML_Prediction_Model_Training.py
- Step18_ML_Prediction/Step18_4_1_ML_Prediction_Model_Training_H3K27me3_Only.py
- Change H3K27me3 to other histone marks for the complete table in figure 4E
- Use Google Spreadsheet to tabulate the results 

Supplementary Figure S4: CBX8 ENCODE-based XGBoost analysis referenced in the Figure 4 section.

Figure S4A
- Step5_ChIPseq_Preprocessing.R
- Step18_ML_Prediction_H1/Step18_1_ML_Prediction_hg38_bin.py
- Step18_ML_Prediction_H1/Step18_2_ML_Prediction_pos_neg_sets.py
- Step18_ML_Prediction_H1/Step18_3_ML_Prediction_histone_bigwig.py
- Step18_ML_Prediction_H1/Step18_4_0_ML_Prediction_Model_Training.py
- Use .jpg outputs for H1

- Step5_ChIPseq_Preprocessing.R
- Step18_ML_Prediction_K562/Step18_1_ML_Prediction_hg38_bin.py
- Step18_ML_Prediction_K562/Step18_2_ML_Prediction_pos_neg_sets.py
- Step18_ML_Prediction_K562/Step18_3_ML_Prediction_histone_bigwig.py
- Step18_ML_Prediction_K562/Step18_4_0_ML_Prediction_Model_Training.py
- Use .jpg outputs for K562

Figure S4B
- Step5_ChIPseq_Preprocessing.R
- Step18_ML_Prediction_H1/Step18_1_ML_Prediction_hg38_bin.py
- Step18_ML_Prediction_H1/Step18_2_ML_Prediction_pos_neg_sets.py
- Step18_ML_Prediction_H1/Step18_3_ML_Prediction_histone_bigwig.py
- Step18_ML_Prediction_H1/Step18_4_0_ML_Prediction_Model_Training.py
- Step18_ML_Prediction_H1/Step18_6_Custome_SHAP_Beeswarm.py
- Use .jpg outputs for H1

- Step5_ChIPseq_Preprocessing.R
- Step18_ML_Prediction_K562/Step18_1_ML_Prediction_hg38_bin.py
- Step18_ML_Prediction_K562/Step18_2_ML_Prediction_pos_neg_sets.py
- Step18_ML_Prediction_K562/Step18_3_ML_Prediction_histone_bigwig.py
- Step18_ML_Prediction_K562/Step18_4_0_ML_Prediction_Model_Training.py
- Use .jpg outputs for K562

### 9. Hardware specification

Windows PC (mainly used for R)
- CPU: Intel(R) Core(TM) i5-10505 CPU @ 3.20GHz (3.20 GHz)
- RAM: 16.0 GB
- GPU: Intel(R) UHD Graphics 630

Linux PC (mainly used for Python)
- CPU: Intel(R) Xeon(R) Gold 6126 CPU @ 2.60GHz
- RAM: 754 Gb
- GPU: NVIDIA Corporation TU102 [GeForce RTX 2080 Ti]


### 10. Known limitations

RNA-seq raw sequencing files are unavailable.

Newly generated GEO accessions are pending publication.

Accurate enhancer-gene pairing remains a major limitation and GeneHancer provides a working approximation rather than a definitive cell-type-specific map.

Some UpDEGs may reflect secondary induction rather than direct SRA targeting.
