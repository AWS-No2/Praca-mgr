#Dependancies
library("Biostrings")
library("pwalign")
library("msa")
library("DECIPHER")

file_seq_E2_ref <- "data/Seq. ref. Ex2.txt"
file_seq_I2_ref <- "data/Seq. ref. Int2.txt"

file_seq_kohorta <- "data/Consensus SangerContig EX2IN2_56-67 Popr.fasta"

file_seq_E2_species <- "data/Species Ex2.txt"
file_seq_I2_species <- "data/Species In2.txt"

seq_E2_ref <- readDNAStringSet(file_seq_E2_ref, format="fasta")
seq_I2_ref <- readDNAStringSet(file_seq_I2_ref, format="fasta")
seq_E2_species <- readDNAStringSet(file_seq_E2_species,format="fasta")
seq_I2_species <- readDNAStringSet(file_seq_I2_species, format="fasta")

seq_kohorta <- readDNAStringSet(file_seq_kohorta, format="fasta")

raport <- data.frame(
      Nazwa_Ref = character(),
      Nazwa_Sekwencji = character(),
      Procent_conserved= c(),
      Metoda = c()
)

seqs <- seq_kohorta
method <- "overlap"
for (i in seq_along(seqs)){
  #Alignment globalny wobec Cfam
  aln_global <- pwalign::pairwiseAlignment(seq_I2_ref[2],
        seqs[i],
        substitutionMatrix = "BLOSUM62",
        type = method
      )
  #Obliczenie % tożsamości
  pid <- pwalign::pid(aln_global,
        type = "PID1"
      )
  #Tworzenie raportu
  wiersz_raport <- data.frame(
    Nazwa_Ref = names(seq_I2_ref[2]),
    Nazwa_Sekwencji = names(seqs[i]),
    Procent_conserved = round(pid,2),
    Metoda = method
  )
  raport <- rbind(wiersz_raport, raport)
}

raport_dir <- "results/Konserwatywność sekwencji w kohorcie.csv"
write.csv2(raport, raport_dir, row.names = FALSE)


seqs_kohorta_cfam <- c(seq_I2_ref[2],seq_kohorta)
alignment_muscle_cfam <- msa(seqs_kohorta_cfam,
        method = "Muscle",
        order = "input",
        verbose = T)
aln_muscle_cfam <- as(alignment_muscle_cfam, 
        "DNAStringSet")
BrowseSeqs(aln_muscle_cfam)


raport_species <- data.frame(
      Nazwa_Ref = character(),
      Nazwa_Sekwencji = character(),
      Procent_conserved = c(),
      Metoda = method
)

seqs <- seq_E2_species
method <- "overlap"

for (i in seq_along(seqs)){
  #Alignment globalny wobec Cfam
  aln_global <- pwalign::pairwiseAlignment(seq_E2_ref[2],
        seqs[i],
        substitutionMatrix = "BLOSUM62",
        type = method
      )
  #Obliczenie % tożsamości
  pid <- pwalign::pid(aln_global,
        type = "PID1"
      )
  #Tworzenie raportu
  wiersz_raport <- data.frame(
    Nazwa_Ref = names(seq_E2_ref[2]),
    Nazwa_Sekwencji = names(seqs[i]),
    Procent_conserved = round(pid,2),
    Metoda = method
  )
  raport_species <- rbind(wiersz_raport, raport_species)
}

raport_species_dir <- "results/Konserwatywność sekwencji u różnych gat.csv"
write.csv2(raport_species, raport_species_dir, row.names = FALSE)


seqs_kohorta_cfam <- c(seq_I2_ref[2],seq_I2_species)
alignment_muscle_cfam <- msa(seqs_kohorta_cfam,
        method = "Muscle",
        order = "input",
        verbose = T)
aln_muscle_cfam <- as(alignment_muscle_cfam, 
        "DNAStringSet")
BrowseSeqs(aln_muscle_cfam)