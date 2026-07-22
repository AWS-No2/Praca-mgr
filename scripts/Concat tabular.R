# Dependancies
# install.packages(c("dplyr", "purrr", "readr", "writexl", "paletteer", "ggplot2", "tidyr"))

library(dplyr)
library(purrr)
library(readr)
library(writexl)
library(paletteer)
library(ggplot2)
library(tidyr)

# 1. Wczytanie pliku głównego (Samdept z już nazwanymi kolumnami)
setwd("data/Samdepth tabular")
samdept <- read_tsv("Samdept SRR35635694,SRR35635695,SRR35635698,SRR35635699,SRR37655155.tabular")

# 2. Utworzenie wektora ze ścieżkami do pozostałych plików
# Wpisz tutaj nazwy swoich plików z kolumną "0"
file_list <- c("Samdepth SRR38085904.tabular", "Samdepth SRR38085905.tabular", 
        "Samdepth SRR38085906.tabular", "Samdepth SRR38085907.tabular", 
        "Samdepth SRR38085908.tabular")

# 3. Zdefiniowanie właściwych numerów prób dla powyższych plików
# Skrypt podmieni nagłówek "0" na te właśnie nazwy
sample_names <- c("SRR38085904", "SRR38085905", "SRR38085906", "SRR38085907", "SRR38085908")

# 4. Wczytanie mniejszych plików i nadanie im właściwych nazw kolumn
other_files <- map2(file_list, sample_names, function(file_path, new_name) {
  df <- read_tsv(file_path)
  
  # Zmiana nazwy trzeciej kolumny (pokrycia) na nazwę z wektora 'sample_names'
  colnames(df)[3] <- new_name
  return(df)
})

# 5. Zgrupowanie pliku głównego oraz pozostałych plików w jedną listę
all_datasets <- c(list(samdept), other_files)

# 6. Połączenie wszystkich ramek danych po kolumnach `#CHROM` i `POS`
# Używamy full_join, co gwarantuje, że żadna pozycja POS nie zniknie, 
# nawet jeśli występuje tylko w jednym pliku.
final_data <- reduce(all_datasets, full_join, by = c("#CHROM", "POS"))

# Opcjonalnie: Jeśli pozycja występuje w pliku A, ale nie w pliku B, 
# R wstawi w to miejsce 'NA'. Poniższa funkcja zamienia 'NA' na '0' (brak pokrycia).
final_data <- final_data %>% mutate(across(everything(), ~replace_na(.x, 0)))

# 7. Zapisanie połączonego pliku wynikowego jako .csv
write.csv2(final_data, "merged_coverage.csv")

sum(is.na(final_data)==TRUE)

samdept_all_full <- read.csv2("merged_coverage.csv")

#Pozycje zaczerpnięte z wykonanego MSA w MEGA
#start ekson 2 - poz 52
samdept_all[52,]
#POS: 21884490
#koniec ekson 2 - 245
samdept_all[245,]
#POS: 21884683

#start intron 2 - 246
samdept_all[246,]
#POS: 21884684
#koniec intron 2 - 572
samdept_all[572,]
#POS: 21885010

samdept_all <- samdept_all_full[52:572,] #subset tylko z ekson 2 i intron 2
samdept_ekson2 <- samdept_all_full[52:245,] #subset tylko ekson 2
samdept_intron2 <- samdept_all_full[246:572,] #subset tylko intron 2


#Raport całość
raport <- data.frame(
  Sample = c(),
  Srednia = c(),
  Odch.st. = c(),
  Mediana = c(),
  Min = c(),
  Max = c()
)

for (i in 4:13){
  srednia = mean(samdept_all[,i])
  odch.st = sd(samdept_all[,i])
  mediana = median(samdept_all[,i])
  min_val = min(samdept_all[,i])
  max_val = max(samdept_all[,i])
  wiersz_raport <- data.frame(
    Sample = colnames(samdept_all)[i],
    Srednia = srednia,
    Odch.st. = odch.st,
    Mediana = mediana,
    Min = min_val,
    Max = max_val
  )
  raport <- rbind(raport, wiersz_raport)
}
#Raport ekson 2
samdept_all <- samdept_ekson2

raport_ekson <- data.frame(
  Sample = c(),
  Srednia = c(),
  Odch.st. = c(),
  Mediana = c(),
  Min = c(),
  Max = c()
)

for (i in 4:13){
  srednia = mean(samdept_all[,i])
  odch.st = sd(samdept_all[,i])
  mediana = median(samdept_all[,i])
  min_val = min(samdept_all[,i])
  max_val = max(samdept_all[,i])
  wiersz_raport <- data.frame(
    Sample = colnames(samdept_all)[i],
    Srednia = srednia,
    Odch.st. = odch.st,
    Mediana = mediana,
    Min = min_val,
    Max = max_val
  )
  raport_ekson <- rbind(raport_ekson, wiersz_raport)
}
#Raport intron 2
samdept_all <- samdept_intron2

raport_intron <- data.frame(
  Sample = c(),
  Srednia = c(),
  Odch.st. = c(),
  Mediana = c(),
  Min = c(),
  Max = c()
)

for (i in 4:13){
  srednia = mean(samdept_all[,i])
  odch.st = sd(samdept_all[,i])
  mediana = median(samdept_all[,i])
  min_val = min(samdept_all[,i])
  max_val = max(samdept_all[,i])
  wiersz_raport <- data.frame(
    Sample = colnames(samdept_all)[i],
    Srednia = srednia,
    Odch.st. = odch.st,
    Mediana = mediana,
    Min = min_val,
    Max = max_val
  )
  raport_intron <- rbind(raport_intron, wiersz_raport)
}

raporty_do_excela <- list(
  "Calosc" = raport,
  "Ekson 2" = raport_ekson,
  "Intron 2" = raport_intron
)
# write.csv2(raport, "Raport_WGS.csv")
write_xlsx(raporty_do_excela, "Raport_WGS.xlsx")


#WYKRESY

# 1. Upewniamy się, że dane są w formacie długim (niezbędne dla ggplot2)
long_data <- samdept_all %>%
  pivot_longer(cols = starts_with("SRR"), 
               names_to = "Sample", 
               values_to = "Depth")

# 2. Definiujemy granice eksonu i intronu
exon_start <- 21884490
exon_end <- 21884683
intron_start <- 21884684
intron_end <- 21885010

# 3. Generowanie wykresu
plot <- ggplot(long_data, aes(x = POS, y = Depth, color = Sample)) +
  
  # Rysowanie krzywych pokrycia
  geom_line(linewidth = 0.8, alpha = 0.8) +
  
  # Jedna pionowa linia oddzielająca
  geom_vline(xintercept = exon_end, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  # Napis "Ekson" po lewej stronie linii
  # Zmień wartość '2' na większą (np. 10 lub 50), jeśli napis będzie za blisko linii 
  # w zależności od tego, jak szeroki jest Twój wycinek osi X
  annotate("text", x = exon_end - 10, y = Inf, label = "Ekson", 
           hjust = 1, vjust = 1.5, fontface = "bold", size = 5, color = "black") +
  
  # Napis "Intron" po prawej stronie linii
  annotate("text", x = exon_end + 10, y = Inf, label = "Intron", 
           hjust = 0, vjust = 1.5, fontface = "bold", size = 5, color = "black") +
  
  # Wymuszenie startu osi Y DOKŁADNIE od 0
  # expand = expansion(mult = c(0, 0.05)) sprawia, że na dole nie ma marginesu (0),
  # a na górze zostaje 5% luzu, żeby wykres i napisy nie dotykały samej krawędzi okna.
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) +
  
  #Zmiana koloru linii - jeżeli potrzebne
  # scale_color_paletteer_d("ggthemes::gdoc") +
  # Etykiety i kosmetyka
  labs(
    title = "Głębokość pokrycia WGS",
    x = "Pozycja na chromosomie (POS)",
    y = "Głębokość pokrycia (X)",
    color = "Próba"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    
    # Dodanie ciągłej, czarnej linii dla osi X (żeby ładnie podkreślić wartość 0)
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    
    #Odsunięcie podpisu osi X w dół
    axis.title.x = element_text(margin = margin(t = 20)),
    #Odsunięcie podpisu osi Y w lewo
    axis.title.y = element_text(margin = margin(r = 20)),
    panel.grid.minor = element_blank()
  )
path <- "data/Wizualizacje"
ggsave("WGS_depth.png",
        plot,
        path = path, 
        dpi = 300,
        width = 15, height = 7)