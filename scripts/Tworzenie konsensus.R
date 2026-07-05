#Dependancies
library(sangeranalyseR)
library(Biostrings)

parent_dir <- "data/Chromogramy razem"
sample_ids <- paste("EX2IN2",59:67, sep="-")

final_consensus_list <- DNAStringSet()

for (id in sample_ids){
  contig <- SangerContig(
      ABIF_Directory = parent_dir,
      REGEX_SuffixForward = "_F.ab1",
      REGEX_SuffixReverse = "_R.ab1",
      contigName = id,
      TrimmingMethod = "M1",
      M2CutoffQualityScore = 13
    )
  consensus_raw <- contig@contigSeq #format DNAString
  #Konwersja na DNAStringSet (aby dodać nazwę)
  consensus_named <- DNAStringSet(consensus_raw)
  names(consensus_named) <- id
  final_consensus_list <- c(final_consensus_list, consensus_named)
}
fasta_file <- "data/Consensus SangerContig EX2IN2_56-67.fasta"
writeXStringSet(final_consensus_list,fasta_file) #Zapis powstałych sekwencji konsensusowych w pliku

alphabetFrequency(final_consensus_list) #Sprawdzanie, czy są nukleotydy niejednoznaczne
# launchApp(contig) #można obejrzeć chromatogramy za pomocą interfejsu shiny

#Generowanie raportu (omijamy shiny)
dir.create("results/Raporty konsensus")

for (id in sample_ids){
    contig <- SangerContig(
      ABIF_Directory = parent_dir,
      REGEX_SuffixForward = "_F.ab1",
      REGEX_SuffixReverse = "_R.ab1",
      contigName = id,
      TrimmingMethod = "M1",
      M2CutoffQualityScore = 13
    )
  generateReport(contig, outputDir = "results/Raporty konsensus")
}

#Ręczne poprawki
#Poszukiwanie pozycji niejednoznacznych nukleotydów
iupac <- c("R", "Y", "K", "M", "S", "W", "B", "D", "H", "V", "N")
find_iupac_positions <- function(DNA_set) {
  results_list <- list()
  for (i in 1:length(DNA_set)){
    seq_name <- names(DNA_set)[i]
    seq_char <- as.character(DNA_set[[i]])
    seq_vec <- strsplit(seq_char,"")[[1]] #rozbijamy sekwencję na pojedyncze znaki
    #szukamy pozycji, gdzie znak należący do listy IUPAC się znajduje
    pos <- which(seq_vec %in% iupac)

    if (length(pos)>0){
      results_list[[seq_name]] <- data.frame(
        Osobnik = seq_name,
        Pozycja = pos,
        Kod = seq_vec[pos]
      )
    }
  }
  do.call(rbind, results_list)
}
mapa_pozycji_iupac <- find_iupac_positions(final_consensus_list)

# # Przykład: w 3. sekwencji (szczenie 59) na pozycji 142 zmieniamy 'R' na 'A'
# subseq(final_consensus_list[3], start=142, end=142) <- "A"
# final_consensus_list[1] <- chartr("Y", "T", final_consensus_list[1])

#Ex2In2-59 [1]
final_consensus_list[1] <- chartr("Y", "T", final_consensus_list[1])
final_consensus_list[1] <- chartr("R", "A", final_consensus_list[1])

#EX2IN2-60 [2]
subseq(final_consensus_list[2], start=51, end=51) <- "A"
subseq(final_consensus_list[2], start=504, end=504) <- "A"
subseq(final_consensus_list[2], start=492, end=492) <- "T"
#41R to najprawdopodobniej ta 4 adenina - zmieniam
subseq(final_consensus_list[2], start=41, end=41) <- "A"

# EX2IN2-61 [3]
#38K, 39R, 40W, 42K, 43R, 494Y
#przesunięcie spowodowane poślizgiem przy adeninach
subseq(final_consensus_list[3], start=38, end=40) <- "GAT"
subseq(final_consensus_list[3], start=42, end=43) <- "GA"
subseq(final_consensus_list[3], start=494, end=494) <- "T"

# EX2IN2-62 [4]
#40W, 42K, 43R, 494Y
#przesunięcie spowodowane poślizgiem przy adeninach
subseq(final_consensus_list[4], start=40, end=40) <- "T"
subseq(final_consensus_list[4], start=42, end=43) <- "GA"
#tutaj powtarzalny błąd, początek sekwencji R, ale na F wyraźnie widoczne 2 piki T
subseq(final_consensus_list[4], start=494, end=494) <- "T"

# EX2IN2-63 [5]
#38K, 39R, 40W, 42K, 43R, 494Y
#przesunięcie spowodowane poślizgiem przy adeninach
subseq(final_consensus_list[5], start=38, end=40) <- "GAT"
subseq(final_consensus_list[5], start=42, end=43) <- "GA"
#tutaj powtarzalny błąd, początek sekwencji R, ale na F wyraźnie widoczne 2 piki T
subseq(final_consensus_list[5], start=494, end=494) <- "T"

# EX2IN2-64 [6]
#44W
subseq(final_consensus_list[6], start=44, end=44) <- "A"

# EX2IN2-65 [7]
#40W, 42K, 43R, 53R, 494Y
subseq(final_consensus_list[7], start=40, end=43) <- "TTGA"
subseq(final_consensus_list[7], start=53, end=53) <- "A"
subseq(final_consensus_list[7], start=494, end=494) <- "T"

# EX2IN2-66 [8]
#46R, 461K, 497W
subseq(final_consensus_list[8], start=46, end=46) <- "A"
subseq(final_consensus_list[8], start=461, end=461) <- "T"
subseq(final_consensus_list[8], start=497, end=497) <- "T"

# EX2IN2-67 [9]
#39K, 40R, 41W, 43K, 44R
subseq(final_consensus_list[9], start=39, end=44) <- "GATTGA"


fasta_file_popr <- "results/Consensus SangerContig EX2IN2_56-67 Popr.fasta"
writeXStringSet(final_consensus_list,fasta_file_popr)
