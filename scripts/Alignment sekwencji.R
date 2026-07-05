#Dependancies
library("Biostrings")
library("msa")
library("DECIPHER")

#Wczytanie plików FASTA do R
#Sekwencje referencyjne
Ref_file <- "data/Seq. ref Ex2_Intr2.txt"
Seq_ref <- readDNAStringSet(Ref_file, format = "fasta")
#Zmiana nazwy do wizualizacji
names(Seq_ref) <- c("Ref_alt_Boxer_Ex2In2","Ref_Cfam_Ex2In2")

#Wczytanie sekwencji badanych (po oczyszczaniu)
Rodzina1_F_file <- "data/Chromogramy_59-63/59-63 F clean.fasta"
Rodzina1_R_file <- "data/Chromogramy_59-63/59-63 R clean.fasta"

Rodzina2_F_file <- "data/Chromogramy_64-67/64-67 F clean.fasta"
Rodzina2_R_file <- "data/Chromogramy_64-67/64-67 R clean.fasta"

Rodzina2_F_0.01_file <- "data/Chromogramy_64-67/64-67 F clean 0.01.fasta"

Kohorta_consensus_file <- "data/Consensus SangerContig EX2IN2_56-67.fasta"

Seq_Rodzina1_F <- readDNAStringSet(Rodzina1_F_file, format = "fasta")
Seq_Rodzina1_R <- readDNAStringSet(Rodzina1_R_file, format = "fasta")
Reverse_Seq_Rodzina1_R <- reverseComplement(Seq_Rodzina1_R)

Seq_Rodzina2_F <- readDNAStringSet(Rodzina2_F_file, format = "fasta")
Seq_Rodzina2_R <- readDNAStringSet(Rodzina2_R_file, format = "fasta")
Reverse_Seq_Rodzina2_R <- reverseComplement(Seq_Rodzina2_R)

Seq_Rodzina2_F_0.01 <- readDNAStringSet(Rodzina2_F_0.01_file, format = "fasta")

Seq_kohorta <- readDNAStringSet(Kohorta_consensus_file, format = "fasta")
# Zmiana nazw do wizualizacji - ewentualnie
# names(Seq_Rodzina1_F) <- c()

#RODZINA 1 FORWARD - MUSCLE
#Rodzina 1 F wobec Ref_Cfam MUSCLE
seqs_r1F_muscle_cfam <- c(Seq_ref[2],Seq_Rodzina1_F)
alignment_r1F_muscle_cfam <- msa(seqs_r1F_muscle_cfam,
        method = "Muscle",
        order = "input",
        verbose = T)
aln_r1F_muscle_cfam <- as(alignment_r1F_muscle_cfam,
        "DNAStringSet")
BrowseSeqs(aln_r1F_muscle_cfam)

#Rodzina 1 F wobec Alt Boxer MUSCLE
seqs_r1F_muscle_alt <- c(Seq_ref[1],Seq_Rodzina1_F)
alignment_r1F_muscle_alt <- msa(seqs_r1F_muscle_alt,
        method = "Muscle",
        order = "input",
        verbose = T)
aln_r1F_muscle_alt <- as(alignment_r1F_muscle_alt,
        "DNAStringSet")
BrowseSeqs(aln_r1F_muscle_alt)

#RODZINA 1 FORWARD - CLUSTAL OMEGA
#Rodzina 1 F wobec Ref_Cfam CLUSTAL OMEGA
seqs_r1F_clustalO_cfam <- c(Seq_ref[2],Seq_Rodzina1_F)
alignment_r1F_clustalO_cfam <- msa(seqs_r1F_clustalO_cfam,
        method = "ClustalOmega",
        order = "input",
        verbose = T)
aln_r1F_clustalO_cfam <- as(alignment_r1F_clustalO_cfam,
        "DNAStringSet")
BrowseSeqs(aln_r1F_clustalO_cfam)

#Rodzina 1 F wobec Alt Boxer CLUSTAL OMEGA
seqs_r1F_clustalO_alt <- c(Seq_ref[1],Seq_Rodzina1_F)
alignment_r1F_clustalO_alt <- msa(seqs_r1F_clustalO_alt,
        method = "ClustalOmega",
        order = "input",
        verbose = T)
aln_r1F_clustalO_alt <- as(alignment_r1F_clustalO_alt,
        "DNAStringSet")
BrowseSeqs(aln_r1F_clustalO_alt)

#RODZINA 1 REVERSE - MUSCLE
Reverse_Seq_Rodzina1_R <- reverseComplement(Seq_Rodzina1_R)

#Rodzina 1 REVERSE wobec Ref_Cfam MUSCLE
seqs_Rev_r1R_muscle_cfam <- c(Seq_ref[2],Reverse_Seq_Rodzina1_R)
alignment_Rev_r1R_muscle_cfam <- msa(seqs_Rev_r1R_muscle_cfam,
        method = "Muscle",
        order = "input",
        verbose = T)
aln_Rev_r1R_muscle_cfam <- as(alignment_Rev_r1R_muscle_cfam,
        "DNAStringSet")
BrowseSeqs(aln_Rev_r1R_muscle_cfam)

#Rodzina 1 REVERSE wobec Alt Boxer MUSCLE
seqs_Rev_r1R_muscle_alt <- c(Seq_ref[1],Reverse_Seq_Rodzina1_R)
alignment_Rev_r1R_muscle_alt <- msa(seqs_Rev_r1R_muscle_alt,
        method = "Muscle",
        order = "input",
        verbose = T)
aln_Rev_r1R_muscle_alt <- as(alignment_Rev_r1R_muscle_alt,
        "DNAStringSet")
BrowseSeqs(aln_Rev_r1R_muscle_alt)

#RODZINA 1 REVERSE - CLUSTAL OMEGA
#Rodzina 1 REVERSE wobec Ref_Cfam CLUSTAL OMEGA
seqs_Rev_r1R_clustalO_cfam <- c(Seq_ref[2],Reverse_Seq_Rodzina1_R)
alignment_Rev_r1R_clustalO_cfam <- msa(seqs_Rev_r1R_clustalO_cfam,
        method = "ClustalOmega",
        order = "input",
        verbose = T)
aln_Rev_r1R_clustalO_cfam <- as(alignment_Rev_r1R_clustalO_cfam,
        "DNAStringSet")
BrowseSeqs(aln_Rev_r1R_clustalO_cfam)

#Rodzina 1 REVERSE wobec Alt Boxer CLUSTAL OMEGA
seqs_Rev_r1R_clustalO_alt <- c(Seq_ref[1],Reverse_Seq_Rodzina1_R)
alignment_Rev_r1R_clustalO_alt <- msa(seqs_Rev_r1R_clustalO_alt,
        method = "ClustalOmega",
        order = "input",
        verbose = T)
aln_Rev_r1R_clustalO_alt <- as(alignment_Rev_r1R_clustalO_alt,
        "DNAStringSet")
BrowseSeqs(aln_Rev_r1R_clustalO_alt)

#RODZINA 2 FORWARD - MUSCLE
#Rodzina 2 F wobec Ref_Cfam MUSCLE
seqs_r2F_muscle_cfam <- c(Seq_ref[2],Seq_Rodzina2_F)
alignment_r2F_muscle_cfam <- msa(seqs_r2F_muscle_cfam,
        method = "Muscle",
        order = "input",
        verbose = T)
aln_r2F_muscle_cfam <- as(alignment_r2F_muscle_cfam,
        "DNAStringSet")
BrowseSeqs(aln_r2F_muscle_cfam)

#Rodzina 2 F wobec Alt Boxer MUSCLE
seqs_r2F_muscle_alt <- c(Seq_ref[1],Seq_Rodzina2_F)
alignment_r2F_muscle_alt <- msa(seqs_r2F_muscle_alt,
        method = "Muscle",
        order = "input",
        verbose = T)
aln_r2F_muscle_alt <- as(alignment_r2F_muscle_alt,
        "DNAStringSet")
BrowseSeqs(aln_r2F_muscle_alt)

#RODZINA 2 FORWARD - CLUSTAL OMEGA
#Rodzina 2 F wobec Ref_Cfam CLUSTAL OMEGA
seqs_r2F_clustalO_cfam <- c(Seq_ref[2],Seq_Rodzina2_F)
alignment_r2F_clustalO_cfam <- msa(seqs_r2F_clustalO_cfam,
        method = "ClustalOmega",
        order = "input",
        verbose = T)
aln_r2F_clustalO_cfam <- as(alignment_r2F_clustalO_cfam,
        "DNAStringSet")
BrowseSeqs(aln_r2F_clustalO_cfam)

#Rodzina 2 F wobec Alt Boxer CLUSTAL OMEGA
seqs_r2F_clustalO_alt <- c(Seq_ref[1],Seq_Rodzina2_F)
alignment_r2F_clustalO_alt <- msa(seqs_r2F_clustalO_alt,
        method = "ClustalOmega",
        order = "input",
        verbose = T)
aln_r2F_clustalO_alt <- as(alignment_r2F_clustalO_alt,
        "DNAStringSet")
BrowseSeqs(aln_r2F_clustalO_alt)

#RODZINA 2 REVERSE - MUSCLE
Reverse_Seq_Rodzina2_R <- reverseComplement(Seq_Rodzina2_R)

#Rodzina 2 REVERSE wobec Ref_Cfam MUSCLE
seqs_Rev_r2R_muscle_cfam <- c(Seq_ref[2],Reverse_Seq_Rodzina2_R)
alignment_Rev_r2R_muscle_cfam <- msa(seqs_Rev_r2R_muscle_cfam,
        method = "Muscle",
        order = "input",
        verbose = T)
aln_Rev_r2R_muscle_cfam <- as(alignment_Rev_r2R_muscle_cfam,
        "DNAStringSet")
BrowseSeqs(aln_Rev_r2R_muscle_cfam)

#Rodzina 2 REVERSE wobec Alt Boxer MUSCLE
seqs_Rev_r2R_muscle_alt <- c(Seq_ref[1],Reverse_Seq_Rodzina2_R)
alignment_Rev_r2R_muscle_alt <- msa(seqs_Rev_r2R_muscle_alt,
        method = "Muscle",
        order = "input",
        verbose = T)
aln_Rev_r2R_muscle_alt <- as(alignment_Rev_r2R_muscle_alt,
        "DNAStringSet")
BrowseSeqs(aln_Rev_r2R_muscle_alt)

#RODZINA 2 REVERSE - CLUSTAL OMEGA
#Rodzina 2 REVERSE wobec Ref_Cfam CLUSTAL OMEGA
seqs_Rev_r2R_clustalO_cfam <- c(Seq_ref[2],Reverse_Seq_Rodzina2_R)
alignment_Rev_r2R_clustalO_cfam <- msa(seqs_Rev_r2R_clustalO_cfam,
        method = "ClustalOmega",
        order = "input",
        verbose = T)
aln_Rev_r2R_clustalO_cfam <- as(alignment_Rev_r2R_clustalO_cfam,
        "DNAStringSet")
BrowseSeqs(aln_Rev_r2R_clustalO_cfam)

#Rodzina 1 REVERSE wobec Alt Boxer CLUSTAL OMEGA
seqs_Rev_r2R_clustalO_alt <- c(Seq_ref[1],Reverse_Seq_Rodzina2_R)
alignment_Rev_r2R_clustalO_alt <- msa(seqs_Rev_r2R_clustalO_alt,
        method = "ClustalOmega",
        order = "input",
        verbose = T)
aln_Rev_r2R_clustalO_alt <- as(alignment_Rev_r2R_clustalO_alt,
        "DNAStringSet")
BrowseSeqs(aln_Rev_r2R_clustalO_alt)


###OBYDWIE RODZINY
#OBYDWIE - CLUSTAL OMEGA
#Obydwie F wobec Ref_Cfam CLUSTAL OMEGA
seqs_F_clustalO_cfam <- c(Seq_ref[2],Seq_Rodzina1_F,Seq_Rodzina2_F)
alignment_F_clustalO_cfam <- msa(seqs_F_clustalO_cfam,
        method = "ClustalOmega",
        order = "input",
        verbose = T)
aln_F_clustalO_cfam <- as(alignment_F_clustalO_cfam,
        "DNAStringSet")
BrowseSeqs(aln_F_clustalO_cfam, title = "ClustalOmega")

#Obydwie R wobec Ref_Cfam CLUSTAL OMEGA
reverse_seq_both <- c(Reverse_Seq_Rodzina1_R,Reverse_Seq_Rodzina2_R)

seqs_R_clustalO_cfam <- c(Seq_ref[2],reverse_seq_both)
alignment_R_clustalO_cfam <- msa(seqs_R_clustalO_cfam,
        method = "ClustalOmega",
        order = "input",
        verbose = T)
aln_R_clustalO_cfam <- as(alignment_R_clustalO_cfam,
        "DNAStringSet")
BrowseSeqs(aln_R_clustalO_cfam)

#Obydwie F wobec Ref_Cfam MUSCLE
seqs_F_muscle_cfam <- c(Seq_ref[2],Seq_Rodzina1_F,Seq_Rodzina2_F)
alignment_F_muscle_cfam <- msa(seqs_F_muscle_cfam,
        method = "Muscle",
        order = "input",
        verbose = T)
aln_F_muscle_cfam <- as(alignment_F_muscle_cfam,
        "DNAStringSet")
BrowseSeqs(aln_F_muscle_cfam, title = "MUSCLE")

#Obydwie R wobec Ref_Cfam MUSCLE
seqs_R_clustalO_muscle <- c(Seq_ref[2],reverse_seq_both)
alignment_R_clustalO_muscle <- msa(seqs_R_clustalO_muscle,
        method = "Muscle",
        order = "input",
        verbose = T)
aln_R_clustalO_muscle <- as(alignment_R_clustalO_muscle,
        "DNAStringSet")
BrowseSeqs(aln_R_clustalO_muscle)


#26.03.2026
#Obydwie F (R1 0.05, R2 0.01) wobec Ref_Cfam MUSCLE
seqs_F_muscle_cfam <- c(Seq_ref[2],Seq_Rodzina1_F,Seq_Rodzina2_F_0.01)
alignment_F_muscle_cfam <- msa(seqs_F_muscle_cfam,
        method = "Muscle",
        order = "input",
        verbose = T)
aln_F_muscle_cfam <- as(alignment_F_muscle_cfam,
        "DNAStringSet")
BrowseSeqs(aln_F_muscle_cfam)

#Rodzina 2 F 0.01  wobec Ref_Cfam MUSCLE
seqs_F_muscle_cfam <- c(Seq_ref[2],Seq_Rodzina2_F_0.01)
alignment_F_muscle_cfam <- msa(seqs_F_muscle_cfam,
        method = "Muscle",
        order = "input",
        verbose = T)
aln_F_muscle_cfam <- as(alignment_F_muscle_cfam,
        "DNAStringSet")
BrowseSeqs(aln_F_muscle_cfam)

##KOHORTA
seqs_kohorta_cfam <- c(Seq_ref[2], Seq_kohorta)
alignment_muscle_cfam <- msa(seqs_kohorta_cfam,
        method = "Muscle",
        order = "input",
        verbose = T)
aln_muscle_cfam <- as(alignment_muscle_cfam, 
        "DNAStringSet")
BrowseSeqs(aln_muscle_cfam)

##KOHORTA
Kohorta_con_popr_file <- "data/Consensus SangerContig EX2IN2_56-67 Popr.fasta"
Seq_kohorta_popr <- readDNAStringSet(Kohorta_con_popr_file, format="fasta")

seqs_kohorta_cfam <- c(Seq_ref[2], Seq_kohorta_popr)
alignment_muscle_cfam <- msa(seqs_kohorta_cfam,
        method = "Muscle",
        order = "input",
        verbose = T)
aln_muscle_cfam <- as(alignment_muscle_cfam, 
        "DNAStringSet")
BrowseSeqs(aln_muscle_cfam, title = "Muscle")

intron2_file <- "data/Seq. ref. Int2.txt"
seq_ref_in2 <- readDNAStringSet(intron2_file, format = "fasta")

seqs_kohorta_cfam <- c(seq_ref_in2, Seq_kohorta_popr)
alignment_muscle_cfam <- msa(seqs_kohorta_cfam,
        method = "Muscle",
        order = "input",
        verbose = T)
aln_muscle_cfam <- as(alignment_muscle_cfam, 
        "DNAStringSet")
BrowseSeqs(aln_muscle_cfam)