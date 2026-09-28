# Input data

This directory is used as input for downstream analysis. The structure of this repository is:

```         
|-- README.md
|-- metadata.tsv
|-- usher_barcodes.feather # Freyja SARS-CoV-2 barcodes of 6 Feb 2026, available from the Freyja-data GitHub repository: https://github.com/andersen-lab/Freyja-data
|-- coverage # .tsv files containing genome coverage (in fractions) at a depths 1, 10, 100, etc. (columns: index    depth_threshold coverage)
|   |-- <sample_hash>_genome.tsv
|-- freyja_demix
|   |-- <sample_hash>_lineages.tsv
|-- min_depths_10 # Consensus sequences at minimum depth 10, generated using iVar Consensus
|   |-- <sample_hash>.fasta
|-- mpileup_depths # Coverage depth files generated using samtools mpileup.
|   |-- <sample_hash>.tsv
|-- snv # SNV files generated using iVar variants
    |-- <sample_hash>_ivar.tsv
```

The sample hashes correspond to those listed in the metadata file.

Karthikeyan S, Levy JI, De Hoff P, Humphrey G, Birmingham A, Jepsen K, Farmer S, Tubb HM, Valles T, Tribelhorn CE, Tsai R, Aigner S, Sathe S, Moshiri N, Henson B, Mark AM, Hakim A, Baer NA, Barber T, Belda-Ferre P, Chacón M, Cheung W, Cresini ES, Eisner ER, Lastrella AL, Lawrence ES, Marotz CA, Ngo TT, Ostrander T, Plascencia A, Salido RA, Seaver P, Smoot EW, McDonald D, Neuhard RM, Scioscia AL, Satterlund AM, Simmons EH, Abelman DB, Brenner D, Bruner JC, Buckley A, Ellison M, Gattas J, Gonias SL, Hale M, Hawkins F, Ikeda L, Jhaveri H, Johnson T, Kellen V, Kremer B, Matthews G, McLawhon RW, Ouillet P, Park D, Pradenas A, Reed S, Riggs L, Sanders A, Sollenberger B, Song A, White B, Winbush T, Aceves CM, Anderson C, Gangavarapu K, Hufbauer E, Kurzban E, Lee J, Matteson NL, Parker E, Perkins SA, Ramesh KS, Robles-Sikisaka R, Schwab MA, Spencer E, Wohl S, Nicholson L, McHardy IH, Dimmock DP, Hobbs CA, Bakhtar O, Harding A, Mendoza A, Bolze A, Becker D, Cirulli ET, Isaksson M, Schiabor Barrett KM, Washington NL, Malone JD, Schafer AM, Gurfield N, Stous S, Fielding-Miller R, Garfein RS, Gaines T, Anderson C, Martin NK, Schooley R, Austin B, MacCannell DR, Kingsmore SF, Lee W, Shah S, McDonald E, Yu AT, Zeller M, Fisch KM, Longhurst C, Maysent P, Pride D, Khosla PK, Laurent LC, Yeo GW, Andersen KG, Knight R. Wastewater sequencing reveals early cryptic SARS-CoV-2 variant transmission. Nature. 2022 Sep;609(7925):101-108. doi: 10.1038/s41586-022-05049-6. Epub 2022 Jul 7. PMID: 35798029; PMCID: PMC9433318.

Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352. Epub 2009 Jun 8. PMID: 19505943; PMCID: PMC2723002.

Grubaugh ND, Gangavarapu K, Quick J, Matteson NL, De Jesus JG, Main BJ, Tan AL, Paul LM, Brackney DE, Grewal S, Gurfield N, Van Rompay KKA, Isern S, Michael SF, Coffey LL, Loman NJ, Andersen KG. An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar. Genome Biol. 2019 Jan 8;20(1):8. doi: 10.1186/s13059-018-1618-7. PMID: 30621750; PMCID: PMC6325816.
