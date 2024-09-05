# Working directory: code

setwd("your/path/code")  

# example: setwd("~/lipidomics/raton/000chronics/code")

rm (list = ls ())


# Class annotation was conducted using the RefMet database and compared with the
# LIPID MAPS database. The classification is hierarchical. As an initial step in
# this division, lipids were divided into several principal categories
# ("super classes") containing distinct main classes and sub classes of
# molecules, devising a standard manner of representing the chemical structures
# of individual lipids and their derivatives.


# Load package
library(dplyr)

# Previously we have downloaded a csv from the RefMet website
# https://www.metabolomicsworkbench.org/databases/refmet/browse.php

# Load database
refmet <-
  read.csv("03ClassAnnotation/refmet_22062022.csv")

# Mice negative lipids ----
## Load data
df_MN_c <- read.csv(
  "00Data/all_entities_raton_cronico_negativo_normalizados_217100017.csv",
  comment.char = "#",
  row.names = 1
)

# We only want the lipid names, abundance values in this case is not necessary
lip_MN <- as.vector(rownames(df_MN_c))

# We create a data.frame extracting our lipids from refmet

lip_MN_annot <-
  data.frame(refmet[refmet$refmet_name == lip_MN[1], ])

for (lipido in lip_MN[2:length(lip_MN)]) {
  lip_MN_annot <-
    rbind(lip_MN_annot, refmet[refmet$refmet_name == lipido, ])
}

lip_MN_df <- data.frame(lip_MN)

colnames(lip_MN_df) <- "refmet_name"

lip_MN_annot_df <-
  dplyr::full_join(lip_MN_df, lip_MN_annot, by = "refmet_name")

rm(lip_MN_annot, lip_MN_df, lip_MN, df_MN_c)

lip_MN_annot_df$super_class <-
  as.vector(lip_MN_annot_df$super_class)
lip_MN_annot_df$main_class <- as.vector(lip_MN_annot_df$main_class)
lip_MN_annot_df$sub_class <- as.vector(lip_MN_annot_df$sub_class)

# Since not all lipids have been annotated, we adding manual annotations based
# on lipids which are annotated and compared with lipid maps data base.

indices_fa <-
  which((sub("\\s.*", "", lip_MN_annot_df$refmet_name) == "FA") == TRUE)

lip_MN_annot_df[indices_fa,]$super_class <- "Fatty Acyls"

lip_MN_annot_df[indices_fa,]$main_class <- "Fatty acids"

lip_MN_annot_df[indices_fa,]$sub_class <- "Unsaturated FA"


indices_fahfa <-
  which((sub("\\s.*", "", lip_MN_annot_df$refmet_name) == "FAHFA") == TRUE)

lip_MN_annot_df[indices_fahfa,]$super_class <- "Fatty Acyls"

lip_MN_annot_df[indices_fahfa,]$main_class <- "Fatty esters"

lip_MN_annot_df[indices_fahfa,]$sub_class <- "FAHFA"


indices_Cer_NS <-
  which((sub("\\s.*", "", lip_MN_annot_df$refmet_name) == "Cer_NS") == TRUE)

lip_MN_annot_df[indices_Cer_NS,]$super_class <- "Sphingolipids"

lip_MN_annot_df[indices_Cer_NS,]$main_class <- "Ceramides"

lip_MN_annot_df[indices_Cer_NS,]$sub_class <- "Cer_NS"



indices_Cer_NDS <-
  which((sub("\\s.*", "", lip_MN_annot_df$refmet_name) == "Cer_NDS") == TRUE)

lip_MN_annot_df[indices_Cer_NDS,]$super_class <- "Sphingolipids"

lip_MN_annot_df[indices_Cer_NDS,]$main_class <- "Ceramides"

lip_MN_annot_df[indices_Cer_NDS,]$sub_class <- "Cer_NDS"


indices_Cer_AP <-
  which((sub("\\s.*", "", lip_MN_annot_df$refmet_name) == "Cer_AP") == TRUE)

lip_MN_annot_df[indices_Cer_AP,]$super_class <- "Sphingolipids"

lip_MN_annot_df[indices_Cer_AP,]$main_class <- "Ceramides"

lip_MN_annot_df[indices_Cer_AP,]$sub_class <- "Cer_AP"


indices_Cer_AS <-
  which((sub("\\s.*", "", lip_MN_annot_df$refmet_name) == "Cer_AS") == TRUE)

lip_MN_annot_df[indices_Cer_AS,]$super_class <- "Sphingolipids"

lip_MN_annot_df[indices_Cer_AS,]$main_class <- "Ceramides"

lip_MN_annot_df[indices_Cer_AS,]$sub_class <- "Cer_AS"


indices_Cer_NP <-
  which((sub("\\s.*", "", lip_MN_annot_df$refmet_name) == "Cer_NP") == TRUE)

lip_MN_annot_df[indices_Cer_NP,]$super_class <- "Sphingolipids"

lip_MN_annot_df[indices_Cer_NP,]$main_class <- "Ceramides"

lip_MN_annot_df[indices_Cer_NP,]$sub_class <- "Cer_NP"


indices_Cer_ADS <-
  which((sub("\\s.*", "", lip_MN_annot_df$refmet_name) == "Cer_ADS") == TRUE)

lip_MN_annot_df[indices_Cer_ADS,]$super_class <- "Sphingolipids"

lip_MN_annot_df[indices_Cer_ADS,]$main_class <- "Ceramides"

lip_MN_annot_df[indices_Cer_ADS,]$sub_class <- "Cer_ADS"


indices_HexCer_NS <-
  which((sub("\\s.*", "", lip_MN_annot_df$refmet_name) == "HexCer_NS") ==
          TRUE)

lip_MN_annot_df[indices_HexCer_NS,]$super_class <- "Sphingolipids"

lip_MN_annot_df[indices_HexCer_NS,]$main_class <-
  "Glycosphingolipids"

lip_MN_annot_df[indices_HexCer_NS,]$sub_class <-
  "HexCer_NS" # hexosylceramide non-hydroxyfatty acid-sphingosines




indices_HexCer_AP <-
  which((sub("\\s.*", "", lip_MN_annot_df$refmet_name) == "HexCer_AP") ==
          TRUE)

lip_MN_annot_df[indices_HexCer_AP, ]$super_class <- "Sphingolipids"

lip_MN_annot_df[indices_HexCer_AP, ]$main_class <-
  "Glycosphingolipids"

lip_MN_annot_df[indices_HexCer_AP, ]$sub_class <-
  "HexCer_AP" # hexosylceramide non-hydroxyfatty acid-sphingosines



indices_HexCer_NDS <-
  which((sub("\\s.*", "", lip_MN_annot_df$refmet_name) == "HexCer_NDS") ==
          TRUE)

lip_MN_annot_df[indices_HexCer_NDS,]$super_class <- "Sphingolipids"

lip_MN_annot_df[indices_HexCer_NDS,]$main_class <-
  "Glycosphingolipids"

lip_MN_annot_df[indices_HexCer_NDS,]$sub_class <-
  "HexCer_NDS" # hexosylceramide non-hydroxyfatty acid-sphingosines


indices_sm <-
  which((sub("\\s.*", "", lip_MN_annot_df$refmet_name) == "SM") == TRUE)

lip_MN_annot_df[indices_sm,]$super_class <- "Sphingolipids"

lip_MN_annot_df[indices_sm,]$main_class <- "Sphingomyelins"

lip_MN_annot_df[indices_sm,]$sub_class <- "SM"


indices_lpc <-
  which((sub("\\s.*", "", lip_MN_annot_df$refmet_name) == "LPC") == TRUE)

lip_MN_annot_df[indices_lpc,]$super_class <- "Glycerophospholipids"

lip_MN_annot_df[indices_lpc,]$main_class <- "Glycerophosphocholines"

lip_MN_annot_df[indices_lpc,]$sub_class <- "LPC"


indices_PA <-
  which((sub("\\s.*", "", lip_MN_annot_df$refmet_name) == "PA") == TRUE)

lip_MN_annot_df[indices_PA,]$super_class <- "Glycerophospholipids"

lip_MN_annot_df[indices_PA,]$main_class <- "Glycerophosphates"

lip_MN_annot_df[indices_PA,]$sub_class <- "PA"


indices_PC <-
  which((sub("\\s.*", "", lip_MN_annot_df$refmet_name) == "PC") == TRUE)

lip_MN_annot_df[indices_PC,]$super_class <- "Glycerophospholipids"

lip_MN_annot_df[indices_PC,]$main_class <- "Glycerophosphocholines"

lip_MN_annot_df[indices_PC,]$sub_class <- "PC"


indices_PE <-
  which((sub("\\s.*", "", lip_MN_annot_df$refmet_name) == "PE") == TRUE)

lip_MN_annot_df[indices_PE,]$super_class <- "Glycerophospholipids"

lip_MN_annot_df[indices_PE,]$main_class <-
  "Glycerophosphoethanolamines"

lip_MN_annot_df[indices_PE,]$sub_class <- "PE"


indices_EherPE <-
  which((sub("\\s.*", "", lip_MN_annot_df$refmet_name) == "EtherPE") == TRUE)

lip_MN_annot_df[indices_EherPE,]$super_class <-
  "Glycerophospholipids"

lip_MN_annot_df[indices_EherPE,]$main_class <-
  "Glycerophosphoethanolamines"

lip_MN_annot_df[indices_EherPE,]$sub_class <- "PE-O"


indices_EherPC <-
  which((sub("\\s.*", "", lip_MN_annot_df$refmet_name) == "EtherPC") == TRUE)

lip_MN_annot_df[indices_EherPC,]$super_class <-
  "Glycerophospholipids"

lip_MN_annot_df[indices_EherPC,]$main_class <-
  "Glycerophosphocholines"

lip_MN_annot_df[indices_EherPC,]$sub_class <- "PC-O"


indices_PI <-
  which((sub("\\s.*", "", lip_MN_annot_df$refmet_name) == "PI") == TRUE)

lip_MN_annot_df[indices_PI, ]$super_class <- "Glycerophospholipids"

lip_MN_annot_df[indices_PI, ]$main_class <- "Glycerophosphoinositols"

lip_MN_annot_df[indices_PI, ]$sub_class <- "PI"





# Summary
table(lip_MN_annot_df$super_class)
table(lip_MN_annot_df$main_class)
table(lip_MN_annot_df$sub_class)

# Save
save(lip_MN_annot_df, file = "03ClassAnnotation/lip_MN_annot.RData")


# Mice positive lipids ----
## Load data
df_MP_c <-
  read.csv(
    "00Data/all_entities_raton_cronico_positivo_normalizados_217100017.csv",
    comment.char = "#",
    row.names = 1
  )


# We only want the lipid names, abundance values in this case is not necessary
lip_MP <- as.vector(rownames(df_MP_c))

# We create a data.frame extracting our lipids from refmet

lip_MP_annot <-
  data.frame(refmet[refmet$refmet_name == lip_MP[1], ])

for (lipido in lip_MP[2:length(lip_MP)]) {
  lip_MP_annot <-
    rbind(lip_MP_annot, refmet[refmet$refmet_name == lipido, ])
}

lip_MP_df <- data.frame(lip_MP)

colnames(lip_MP_df) <- "refmet_name"

lip_MP_annot_df <-
  dplyr::full_join(lip_MP_df, lip_MP_annot, by = "refmet_name")

rm(lip_MP_annot, lip_MP_df, lip_MP, df_MP_c)

lip_MP_annot_df$super_class <-
  as.vector(lip_MP_annot_df$super_class)
lip_MP_annot_df$main_class <- as.vector(lip_MP_annot_df$main_class)
lip_MP_annot_df$sub_class <- as.vector(lip_MP_annot_df$sub_class)

# Since not all lipids have been annotated, we adding manual annotations based
# on lipids which are annotated and compared with lipid maps data base.


indices_sm <-
  which((sub("\\s.*", "", lip_MP_annot_df$refmet_name) == "SM") == TRUE)

lip_MP_annot_df[indices_sm,]$super_class <- "Sphingolipids"

lip_MP_annot_df[indices_sm,]$main_class <- "Sphingomyelins"

lip_MP_annot_df[indices_sm,]$sub_class <- "SM"



indices_lpc <-
  which((sub("\\s.*", "", lip_MP_annot_df$refmet_name) == "LPC") == TRUE)

lip_MP_annot_df[indices_lpc,]$super_class <- "Glycerophospholipids"

lip_MP_annot_df[indices_lpc,]$main_class <- "Glycerophosphocholines"

lip_MP_annot_df[indices_lpc,]$sub_class <- "LPC"



indices_Cer_NS <-
  which((sub("\\s.*", "", lip_MP_annot_df$refmet_name) == "Cer_NS") == TRUE)

lip_MP_annot_df[indices_Cer_NS,]$super_class <- "Sphingolipids"

lip_MP_annot_df[indices_Cer_NS,]$main_class <- "Ceramides"

lip_MP_annot_df[indices_Cer_NS,]$sub_class <- "Cer_NS"



indices_TG <-
  which((sub("\\s.*", "", lip_MP_annot_df$refmet_name) == "TG") == TRUE)

lip_MP_annot_df[indices_TG,]$super_class <- "Glycerolipids"

lip_MP_annot_df[indices_TG,]$main_class <- "Triradylglycerols"

lip_MP_annot_df[indices_TG,]$sub_class <- "TAG"



indices_ACar <-
  which((sub("\\s.*", "", lip_MP_annot_df$refmet_name) == "ACar") == TRUE)

lip_MP_annot_df[indices_ACar,]$super_class <- "Fatty Acyls"

lip_MP_annot_df[indices_ACar,]$main_class <- "Fatty esters"

lip_MP_annot_df[indices_ACar,]$sub_class <- "CAR"



indices_PC <-
  which((sub("\\s.*", "", lip_MP_annot_df$refmet_name) == "PC") == TRUE)

lip_MP_annot_df[indices_PC,]$super_class <- "Glycerophospholipids"

lip_MP_annot_df[indices_PC,]$main_class <- "Glycerophosphocholines"

lip_MP_annot_df[indices_PC,]$sub_class <- "PC"


indices_PE <-
  which((sub("\\s.*", "", lip_MP_annot_df$refmet_name) == "PE") == TRUE)

lip_MP_annot_df[indices_PE, ]$super_class <- "Glycerophospholipids"

lip_MP_annot_df[indices_PE, ]$main_class <-
  "Glycerophosphoethanolamines"

lip_MP_annot_df[indices_PE, ]$sub_class <- "PE"



# Summary
table(lip_MP_annot_df$super_class)
table(lip_MP_annot_df$main_class)
table(lip_MP_annot_df$sub_class)

# Save
save(lip_MP_annot_df, file = "03ClassAnnotation/lip_MP_annot.RData")

