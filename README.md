[![DOI](https://zenodo.org/badge/1290169035.svg)](https://doi.org/10.5281/zenodo.21208697)
Autor: Alicja Sójka
Instytucja: Szkoła Główna Gospodarstwa Wiejskiego w Warszawie, Wydział Hodowli, Bioinżynierii i Ochrony Zwierząt
Rok: 2026

O PROJEKCIE:
Niniejsze repozytorium zawiera zbiór autorskich skryptów napisanych w języku R, które stanowią integralną część
mojej pracy magisterskiej pt. "Analiza sekwencji genu DLX6 w celu identyfikacji mutacji mogących 
wywoływać rozszczep wargi i/lub podniebienia u psów rasy buldog francuski".

Kod obejmuje pełen workflow dla danych pochodzących z sekwencjonowania metodą Sangera:
od kontroli jakości i oczyszczania surowych chromatogramów, poprzez budowę sekwencji konsensusowych,
aż po wielokrotne dopasowanie sekwencji (MSA) i analizę konserwatywności (obliczanie procentu tożsamości).

STRUKTURA REPOZYTORIUM I OPIS SKRYPTÓW:
Repozytorium składa się z czterech głównych skryptów analitycznych:
1. Oczyszczanie chromogramów.R
   Skrypt służący do wczytywania plików .ab1, ekstrakcji odczytów oraz automatycznego przycinania sekwencji
   (trimming), osobno Forward i Reverse, na podstawie Phred scores (Quality score) przy wykorzystaniu
   algorytmu Motta. Generuje również raport podsumowujący parametry jakościowe przed i po przycięciu,
   a także wykresy.
2. Tworzenie konsensus.R
   Narzędzie do składania sekwencji (odczyty Forward i Reverse) w pojedyncze sekwencje konsensusowe dla
   badanych prób. Skrypt zawiera również pętlę umożliwiającą odnalezienie pozycji nukleotydów niejednoznacznych
   oraz ręczne ich poprawienie. Plik ten zawiera również informację dot. zamienionych nukleotydów niejednoznacznych
   w badanych sekwencjach konsensusowych.
3. Analiza sekwencji.R
   Przeprowadza dopasowanie sekwencji (Multiple Sequence Alignment) z wykorzystaniem algorytmów MUSCLE oraz
   ClustalOmega, porównując uzyskane sekwencje Forward, Reverse i konsensusowe z sekwencjami referencyjnymi.
4. Procent tożsamości.R
   Oblicza procent tożsamości badanych sekwencji (w tym sekwencji referencyjnych dla różnych gatunków)
   w stosunku do wybranej sekwencji referencyjnej, eksportując wyniki do pliku .csv w formie raportu.
   Pozwala na zastosowanie zarówno metody globalnej, jak i overlap algorytmu Needleman-Wunscha.

DOSTĘPNOŚĆ DANYCH (DATA AVAILABILITY):
**Surowe chromatogramy (pliki .ab1) nie zostały udostępnione w tym repozytorium** ze względu na zasady
dot. nieopublikowanych danych z sekwencjonowania. Skrypty "Oczyszczanie chromogramów.R" oraz "Tworzenie konsensus.R"
zostały upublicznione wyłącznie w celach metodologicznych, aby umożliwić weryfikację logiki analizy.

W folderze data/ udostępniono natomiast publicznie dostępne sekwencje referencyjne i gatunkowe
(Seq. ref Ex2_Intr2.txt, Seq. ref. Ex2.txt, Seq. ref. Int2.txt, Species Ex2.txt, Species In2.txt),
które pozwalają na przetestowanie działania pozostałych części kodu (dopasowania sekwencji i obliczania tożsamości).
Użytkownicy chcący przetestować pełen workflow mogą podstawić własne pliki .ab1, zachowując odpowiednią strukturę
folderów opisaną w kodzie.

WYMAGANIA TECHNICZNE I PAKIETY:
Do uruchomienia skryptów wymagana jest instalacja środowiska R oraz następujących pakietów 
(głównie z repozytorium Bioconductor):
Biostrings
DECIPHER
msa
pwalign
sangerseqR
sangeranalyseR

LICENCJA:
Ten projekt jest udostępniany na warunkach licencji MIT. Szczegóły znajdują się w pliku LICENSE.





