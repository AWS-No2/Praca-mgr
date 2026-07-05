#Dependancies
# BiocManager::install(c("sangerseqR","Biostrings"))
library("sangerseqR")
library("Biostrings")

#1. Ustawienie ścieżki do folderu z plikami w formacie .ab1
input_dir_r1_f <- "data/Chromogramy_59-63/Forward"
input_dir_r1_r <- "data/Chromogramy_59-63/Reverse"

input_dir_r2_f <- "data/Chromogramy_64-67/Forward"
input_dir_r2_r <- "data/Chromogramy_64-67/Reverse"

#Ustawienie pliku dla output file z sekwencjami po obróbce
output_file_r1_f <- "results/Chromogramy_59-63/59-63 F clean.fasta"
output_file_r1_r <- "results/Chromogramy_59-63/59-63 R clean.fasta"

output_file_r2_f <- "results/Chromogramy_64-67/64-67 F clean.fasta"
output_file_r2_r <- "results/Chromogramy_64-67/64-67 R clean.fasta"

output_file_r2_f_0.01 <- "results/Chromogramy_64-67/64-67 F clean 0.01.fasta"
#Utworzenie obiektu zawierającego sekwencje .ab1
files_r1_f <- list.files(path = input_dir_r1_f, pattern = "\\.ab1$", full.names = TRUE)
files_r1_r <- list.files(path = input_dir_r1_r, pattern = "\\.ab1$", full.names = TRUE)

files_r2_f <- list.files(path = input_dir_r2_f, pattern = "\\.ab1$", full.names = TRUE)
files_r2_r <- list.files(path = input_dir_r2_r, pattern = "\\.ab1$", full.names = TRUE)


#Aby nie powtarzać pętli 4 razy przypisz nowy obiekt 
files <- files_r1_f
#Utwórz pustą listę do przechowywania oczyszczonych sekwencji
cleaned_list <- DNAStringSet()

#Utwórz pusty data frame do przechowywania danych z processingu
#UWAGA - czyszczenie raportu nie jest konieczne między runami
raport <- data.frame(
     Nazwa_Sekwencji = character(),
     Dlugosc_Raw = c(),
     Srednia_Jakosc_Q_przed = c(),
     Dlugosc_Trim = c(),
     Srednia_Jakosc_Q_po = c(),
     Ilosc_przycietych = c(),
     Procent_straty = c(),
     stringsAsFactors = FALSE
)

#UWAGA - przed puszczeniem pętli sprawdź
#a. Czy odpowiednie pliki są podpisane pod obiekt files? linijka 28
#b. Czy wyczyściłeś listę? linijka 30
#c. Czy zmieniłeś output file na odpowiedni? linijka 114

## Ustawienie błędu maksymalnego
error_threshold <- 0.01 #0.01 (Q20), 0.05 (Q13), 0.10 (Q10)

#Pętla do automatycznego trimmingu
for (f in files){
  #Wczytujemy plik .ab1
  ab1 <- read.abif(f)
  sanger <- sangerseq(ab1)
  #Wyciągamy surową sekwencję i jakości przypisane do odczytów
  seq_raw <- primarySeq(sanger)
  quals <- ab1@data$PCON.1

  len_raw <- length(seq_raw) #długość surowej sekwencji przed przycięciem
  mean_q <- round(mean(quals), 1) #średnia jakość odczytów przed przycięciem
  ## Automatyczny trimming - metoda Motta
  diff <- error_threshold - (10^(-quals/10))
  cumsum_diff <- cumsum(diff)

  #Ustalenie punktów startu i końca o najwyższej sumarycznej jakości (te co zostają)
  start_pos <- which.max(cumsum_diff <= min(cumsum_diff))
  end_pos <- which.max(cumsum_diff)

  #Wycięcie fragmentu o najwyższej jakości
  trimmed_seq <- seq_raw[start_pos:end_pos]

  #Średnia jakość odczytów po przycięciu
  mean_q_processed <- round(mean(quals[start_pos:end_pos]),1)
  
  #Długość przetworzonej sekwencji
  len_trim <- length(trimmed_seq)

  #Zmiana typu pliku na DNAString
  trimmed_string <- DNAStringSet(trimmed_seq)

  #Usunięcie rozszerzenia z nagłówka
  names(trimmed_string) <- gsub(".ab1","",basename(f))

  #Dodaj przetworzoną sekwencję do listy zbiorczej
  cleaned_list <- c(cleaned_list, trimmed_string)
  
  strata <- len_raw-len_trim
  procent_straty <- round((((len_raw-len_trim)/len_raw)*100),1)

  message(paste(
     "Przetworzono:", basename(f), 
     "|Długość oryginalna", len_raw,
     "|Średnia jakość odczytu przed przycięciem", mean_q,
     "|Długość po przycięciu:",len_trim,
     "|Średnia jakość odczytu po przycięciu", mean_q_processed,
     "|Przycięto:", strata,
     "|Procent straty:", procent_straty
     )
  )
  #Tworzenie raportu z informacjami nt. dokonanych zmian
  wiersz_raport <- data.frame(
     Nazwa_Sekwencji = basename(f),
     Dlugosc_Raw = len_raw,
     Srednia_Jakosc_Q_przed = mean_q,
     Dlugosc_Trim = len_trim,
     Srednia_Jakosc_Q_po = mean_q_processed,
     Ilosc_przycietych = strata,
     Procent_straty = procent_straty
  )
  raport <- rbind(wiersz_raport, raport)
}

#2. Zapis oczyszczonych sekwencji do pliku .fasta
output_file <- output_file_r2_f_0.01 #Aby nie powtarzać kodu przypisz nowy obiekt
writeXStringSet(cleaned_list, filepath = output_file)

#Zapis raportu w formacie csv
raport_dir_r1 <- "results/Chromogramy_59-63/Raport R1 59-63.csv"
raport_dir_r1R <- "results/Chromogramy_59-63/Raport R1 59-63 R.csv"
raport_dir_r2 <- "results/Chromogramy_64-67/Raport R2 64-67.csv"

raport_dir_r2_0.01 <- "results/Chromogramy_64-67/Raport R2 64-67 0.01.csv"

raport_dir <- raport_dir_r2_0.01
write.csv2(raport, raport_dir, row.names = FALSE)

### Wykresy jakości
pdf_dir <- "results/Raporty jakości chromogramów.pdf"
pdf(pdf_dir, width = 11, height = 8)

files_all <- c(files_r1_f,files_r1_r,files_r2_f,files_r2_r)
for (f in files_all){
  #Wczytanie pliku i quals do rysowania
  ab1_temp <- read.abif(f)
  quals_temp <- ab1_temp@data$PCON.1
  #Powtórzenie logiki cięcia
  #ważne - error_treshold musi być taki sam (wyciągnięty przed pętlę)
  diff_temp <- error_threshold - (10^(-quals_temp/10))
  cumsum_diff_temp <- cumsum(diff_temp)

  #Ustalenie punktów startu i końca o najwyższej sumarycznej jakości (te co zostają)
  s_pos <- which.max(cumsum_diff_temp <= min(cumsum_diff_temp))
  e_pos <- which.max(cumsum_diff_temp)

  #Rysowanie wykresu
  plot(quals_temp, type="l", col="dodgerblue",
            main=paste("Analiza jakości:", basename(f)),
            xlab="Pozycja (bp)", ylab="Phred Score",ylim=c(0,50)
         )
   abline(h=13, col="red",lty=3) #h - próg Q13 dla error_threshold 0.05
   abline(v=c(s_pos, e_pos), col="darkorange", lwd=2, lty=2) #linie cięcia
  
  mtext(paste("Trimmed range:", s_pos,"-",e_pos, "bp"), side=3, cex=0.8)
}
dev.off()
