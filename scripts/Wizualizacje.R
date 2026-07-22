#Dependancies
library("Biostrings")
library("msa")
library("ggmsa")
library("ggplot2")
library("DECIPHER")

#Wczytanie plików FASTA do R
#Sekwencje referencyjne
Ref_file <- "data/Seq. ref Ex2_Intr2.txt"
Seq_ref <- readDNAStringSet(Ref_file, format = "fasta")
#Zmiana nazwy do wizualizacji
names(Seq_ref) <- c("Ref_alt_Boxer_Ex2In2","Ref_Cfam_Ex2In2")

RefIn_file <- "data/Seq. ref. Int2.txt"
Seq_intron2 <- readDNAStringSet(RefIn_file, format = "fasta")
names(Seq_intron2) <- c("Ref_alt_Boxer_Intron2","Ref_Cfam_Intron2")

Kohorta_consensus_file <- "data/Consensus SangerContig EX2IN2_56-67 Popr.fasta"
Seq_kohorta <- readDNAStringSet(Kohorta_consensus_file, format = "fasta")
Seq_Rodzina1 <- Seq_kohorta[1:5]
names(Seq_Rodzina1) <- c("59 szczenię R1", "60 szczenię R1", "61 szczenię R1",
                                "62 szczenię R1", "63 matka R1")
Seq_Rodzina1_ref <- c(Seq_ref[2], Seq_Rodzina1)

Seq_Rodzina2 <- Seq_kohorta[6:9]
names(Seq_Rodzina2) <- c("64 szczenię R2", "65 szczenię R2", "66 matka R2",
                                "67 ojciec R2")
Seq_Rodzina2_ref <- c(Seq_ref[2],Seq_Rodzina2)

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

ggmsa(aln_r1F_muscle_cfam, 
      start = 20, 
      end = 40,
      char_width = 0.5,
      seq_name = TRUE,
      color = "Chemistry_NT")+
  geom_msaBar()


#Obydwie F wobec Ref_Cfam MUSCLE
seqs_F_cfam <- c(Seq_ref[2],Seq_Rodzina1_F,Seq_Rodzina2_F)
sample_ids_F <- paste("EX2IN2", 59:67, sep="-")
sample_ids_F <- paste(sample_ids_F, "F", sep="_")
names(seqs_F_cfam) <-c("Ref_Cfam_EX2IN2", 
                                sample_ids_F)

alignment_F_muscle_cfam <- msa(seqs_F_cfam,
        method = "Muscle",
        order = "input",
        verbose = T)
aln_F_muscle_cfam <- as(alignment_F_muscle_cfam,
        "DNAStringSet")
BrowseSeqs(aln_F_muscle_cfam, title = "MUSCLE")

names(aln_F_muscle_cfam) <- c("Ref_Cfam_EX2IN2", 
                                sample_ids_F)


custom_palette <- data.frame(
        names = c("A","T","C","G"),
        color = c("mediumseagreen","tomato2","royalblue3","gray65"))


plot_F_Muscle <- ggmsa(aln_F_muscle_cfam, 
      start = 0, 
      end = 35,
      char_width = 0.2,
      seq_name = TRUE,
      custom_color = custom_palette)+
        theme_minimal()+
        labs(title = "MSA algorytm MUSCLE",
                subtitle = "Sekwencje Forward")

#jeżeli chcesz użyć geom_msaBar() + title to musi być po labs(title=...),
#inaczej nie działa

path <- "data/Wizualizacje"
ggsave("MSA_F_Muscle.png",
        plot_F_Muscle,
        path = path, 
        dpi = 300,
        width = 15, height = 7)

#Obydwie F wobec Ref_Cfam CLUSTAL OMEGA
alignment_F_clustalO_cfam <- msa(seqs_F_cfam,
        method = "ClustalOmega",
        order = "input",
        verbose = T)
aln_F_clustalO_cfam <- as(alignment_F_clustalO_cfam,
        "DNAStringSet")
BrowseSeqs(aln_F_clustalO_cfam, title = "ClustalOmega")

plot_F_ClustalO <- ggmsa(aln_F_clustalO_cfam, 
      start = 0, 
      end = 35,
      char_width = 0.2,
      seq_name = TRUE,
      custom_color = custom_palette)+
        theme_minimal()+
        labs(title = "MSA algorytm ClustalOmega",
                subtitle = "Sekwencje Forward")

ggsave("MSA_F_ClustalOmega.png",
        plot_F_ClustalO,
        path = path, 
        dpi = 300,
        width = 15, height = 7)

#Obydwie R wobec Ref_Cfam MUSCLE
seqs_R_cfam <- c(Seq_ref[2],reverse_seq_both)
sample_ids_R <- paste("EX2IN2", 59:67, sep="-")
sample_ids_R <- paste(sample_ids_R, "R", sep="_")
names(seqs_R_cfam) <-c("Ref_Cfam_EX2IN2", 
                                sample_ids_R)

alignment_R_muscle_cfam <- msa(seqs_R_cfam,
        method = "Muscle",
        order = "input",
        verbose = T)
aln_R_muscle_cfam <- as(alignment_R_muscle_cfam,
        "DNAStringSet")
BrowseSeqs(aln_R_muscle_cfam, title = "Muscle")

plot_R_Muscle <- ggmsa(aln_R_muscle_cfam, 
      start = 123, 
      end = 157,
      char_width = 0.2,
      seq_name = TRUE,
      custom_color = custom_palette)+
        theme_minimal()+
        labs(title = "MSA algorytm MUSCLE",
                subtitle = "Sekwencje Reverse po reverse complement")
ath <- "data/Wizualizacje"
ggsave("MSA_R_Muscle.png",
        plot_R_Muscle,
        path = path, 
        dpi = 300,
        width = 15, height = 7)


#Rodzina 1 consensus wobec Ref_Cfam MUSCLE
alignment_R1 <- msa(Seq_Rodzina1_ref,
                                method = "Muscle",
                                order = "input",
                                verbose = T)
aln_R1 <- as(alignment_R1, "DNAStringSet")
BrowseSeqs(aln_R1, title = "Rodzina 1")

Seq_Rodzina1_in2 <- c(Seq_intron2[2],Seq_Rodzina1)
alignment_R1_In2 <- msa(Seq_Rodzina1_in2,
                                method = "Muscle",
                                order = "input",
                                verbose = T)
aln_R1_In2 <- as(alignment_R1_In2, "DNAStringSet")
BrowseSeqs(aln_R1_In2, title = "Rodzina 1 ref intron 2")

#Rodzina 2 consensus wobec Ref_Cfam MUSCLE
alignment_R2 <- msa(Seq_Rodzina2_ref,
                                method = "Muscle",
                                order = "input",
                                verbose = T)
aln_R2 <- as(alignment_R2, "DNAStringSet")
BrowseSeqs(aln_R2, title = "Rodzina 2")

Seq_Rodzina2_in2 <- c(Seq_intron2[2],Seq_Rodzina2)
alignment_R2_In2 <- msa(Seq_Rodzina2_in2,
                                method = "Muscle",
                                order = "input",
                                verbose = T)
aln_R2_In2 <- as(alignment_R2_In2, "DNAStringSet")
BrowseSeqs(aln_R2_In2, title = "Rodzina 2 ref intron 2")

#PLOTS
custom_palette <- data.frame(
        names = c("A","T","C","G"),
        color = c("mediumseagreen","tomato2","royalblue3","gray65"))
#Rodzina 1
plot_R1 <- ggmsa(aln_R1, 
      start = 240, 
      end = 264,
      char_width = 0.2,
      seq_name = TRUE,
      custom_color = custom_palette)+
        theme_minimal()+
        labs(title = "Rodzina 1",
                subtitle = "MSA algorytm MUSCLE, sekwencje konsensusowe, intron 2 fragment")
path <- "data/Wizualizacje"
ggsave("Rodzina1_Muscle_Intron2.png",
        plot_R1,
        path = path, 
        dpi = 300,
        width = 15, height = 7)
#Rodzina 2
plot_R2 <- ggmsa(aln_R2, 
      start = 242, 
      end = 266,
      char_width = 0.2,
      seq_name = TRUE,
      custom_color = custom_palette)+
        theme_minimal()+
        labs(title = "Rodzina 2",
                subtitle = "MSA algorytm MUSCLE, sekwencje konsensusowe, intron 2 fragment")
ggsave("Rodzina2_Muscle_Intron2.png",
        plot_R2,
        path = path, 
        dpi = 300,
        width = 15, height = 7)

