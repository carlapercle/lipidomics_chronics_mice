# Working directory: code
setwd("~/lipidomics/raton/000chronics")  # example: "~/lipidomics/chronics/code"

# This script might appear complex, but its purpose is to reformat our data 
# correctly for uploading it to the LINEX2 web tool with various nomenclatures. 
# Trying the different nomenclatures,  comparing the failed information and 
# combining in the same file the nomebclature, we get the maximum number of 
# lipids that were recognised by the tool.

rm (list = ls ())


# Load data
load("00Data/datos_cronicos_mice_wo.RData")

## File obtained using different sources and doing a manual revision
nom <- read.csv("00Data/Nomenclatura_lipidos_Descargar.csv") 
colnames(nom)
colnames(nom) <- c("My_data","Input_name","MetaboAnalyst",
                   "LipidlynxX_pasando.Input_name",
                   "LipidlynxX_pasando.MetaboAnalyst",
                   "Standarized_name",
                   "Annotations")

df_M_c_wo$lip <- sapply(strsplit(rownames(df_M_c_wo),"_neg"),`[`, 1)
df_M_c_wo$lip <- sapply(strsplit(df_M_c_wo$lip,"_pos"),`[`, 1)
df_M_c_wo$Input_name <- sapply(strsplit(df_M_c_wo$lip,"  RT"),`[`, 1)

df_M_c_wo$Input_name[!(df_M_c_wo$Input_name %in% nom$Input_name)]

summary(unique(df_M_c_wo$Input_name))
summary(unique(nom$Input_name))

intersect(df_M_c_wo$Input_name,nom$Input_name) # 310 coincidencias
summary(intersect(df_M_c_wo$Input_name,nom$Input_name)) # 310 coincidencias

class(df_M_c_wo$Input_name)
class(nom$Input_name)

df_M_c_wo <- df_M_c_wo[order(df_M_c_wo$Input_name),]
nom <- nom[order(nom$Input_name),]  

summary(df_M_c_wo$Input_name == nom$Input_name)

# Otra opcion
# library(dplyr)
# a <- left_join(df_M_c_wo, nom, by = "Input_name")
# a <- data.frame(df_M_c_wo[,"Input_name"],nom[,"Input_name"])

# Otra
match(df_M_c_wo$Input_name,nom$Input_name)

df_M_c_wo$Input_name <- nom$Input_name

# Two columns are equal to merge
abun_nom <- cbind(df_M_c_wo, nom)


# Now, I have the lipid name in different formats and the aundance values

summary(unique(abun_nom$Input_name)) #310 check


# Standarized name

abun_SN <- abun_nom[,c(56,1:48)]

summary(unique(abun_SN$Standarized_name)) # 310

write.csv(abun_SN, file = "05LINEX/abun_SN.csv",
          row.names = F)


# RefMet
# 
# abun_RefMet <- abun_nom[,c(56,1:49)]
# 
# summary(unique(abun_RefMet$Refmet)) # 452:

# # I lose two lipids in this nomenclature (452/454).
# # Cer 42:2;O3 --> Cer_ADS d42:2 Cer_AS d42:2
# # Cer 44:1;O2 --> Cer_NDS d44:1 Cer_NS d44:1
# 
# # Median
# abun_RefMet_median <- stats::aggregate(abun_RefMet, by = list(Refmet = abun_RefMet$Refmet), FUN = median)
# 
# abun_RefMet_median <- abun_RefMet_median[,c(1,3:25)]
# 
# write.csv(abun_RefMet_median, file = "06LINEX/data_steps/abun_RefMet_median.csv",
#           row.names = F)
# 
# # Lynxx
# 
# abun_Lynx <- abun_nom[,c(27,2:24)]
# 
# summary(unique(abun_Lynx$LipidlynxX)) # 451:
# 
# # I lose three lipids in this nomenclature (451/454).
# # Cer 42:2;O3 --> Cer_ADS d42:2 Cer_AS d42:2
# # Cer 44:1;O2 --> Cer_NDS d44:1 Cer_NS d44:1
# # OxPC 18:1_18:3+1O
# # OxPC 18:1_20:3+1O
# 
# write.csv(abun_RefMet_median, file = "06LINEX/data_steps/abun_RefMet_median.csv",
#           row.names = F)
#         
# # After processing the data in LINEX web tool, I obtained the files containing 
# # the species that were not recognized using LipidLynxX nomencature.

LINEX_molspec_failed <-
  read.csv(
    "05LINEX/26_09/LINEX_molspec_failed_246_202309191351.txt",
    header = FALSE,
    sep = "\t"
  )

LINEX_lynx_failed <-
  read.csv(
    "05LINEX/26_09/LINEX_lynx_failed_246_202309191351.txt",
    header = FALSE,
    sep = "\t"
  )

failed_refmet <- rbind(LINEX_lynx_failed,LINEX_molspec_failed)
failed_refmet <- failed_refmet$V1
summary(unique(failed_refmet))


quiero_Met <- nom[nom$MetaboAnalyst %in% failed_refmet,] 

quiero_Lynx_Input <- nom[nom$LipidlynxX_pasando.Input_name %in% failed_refmet,]

quiero_Lynx_Met <- nom[nom$LipidlynxX_pasando.MetaboAnalyst %in% failed_refmet,]

quiero <- rbind(quiero_Met, quiero_Lynx_Input, quiero_Lynx_Met)

quiero <- quiero[!duplicated(quiero),]



# Trying with different versions
abun_nom_failed <- abun_nom[abun_nom$My_data %in% quiero$My_data,]

abun_Lynx_IN <- abun_nom_failed[,c(54,1:48)]

abun_Lynx_IN_median <- stats::aggregate(abun_Lynx_IN, by = list(LipidlynxX = abun_Lynx_IN$LipidlynxX_pasando.Input_name), FUN = median)
abun_Lynx_IN_median <- abun_Lynx_IN_median[c(2:25),c(1,3:50)]

write.csv(abun_Lynx_IN_median, file = "05LINEX/abun_Lynx_median_IN.csv",
          row.names = F)
# Fallan todos jeje

abun_Met <-  abun_nom_failed[,c(53,1:48)]

abun_Met_median <- stats::aggregate(abun_Met, by = list(MetaboAnalyst = abun_Met$MetaboAnalyst), FUN = median)
abun_Met_median <- abun_Met_median[,c(1,3:50)]

write.csv(abun_Met_median, file = "05LINEX/abun_Met_median.csv",
          row.names = F)

# Como solo me reconoce 1, añado dos que me reconoce
A <- abun_SN[abun_SN$Standarized_name == c("Acar(16:0)","Acar(18:0)"),]

colnames(A) <- colnames(abun_Met_median)
# rownames(A) <- NULL
abun_Met_median_plus2 <- rbind(abun_Met_median,A)



write.csv(abun_Met_median_plus2 , file = "05LINEX/abun_Met_plus2_median.csv",
          row.names = F)


# Lynx_con met
abun_Lynx_met <- abun_nom_failed[,c(55,1:48)]

abun_Lynx_met_median <- stats::aggregate(abun_Lynx_met, by = list(LipidlynxX = abun_Lynx_met$LipidlynxX_pasando.MetaboAnalyst), FUN = median)
abun_Lynx_met_median <- abun_Lynx_met_median[c(2:5),c(1,3:50)]

write.csv(abun_Lynx_IN_median, file = "05LINEX/abun_Lynx_median_met.csv",
          row.names = F)



# From abun_nom, I take all lipids with Refmet except these 19 with lipidLynx


abun_SN_median <- stats::aggregate(abun_SN, by = list(Lipid = abun_SN$Standarized_name), FUN = median)
abun_SN_median <- abun_SN_median[,c(1,3:50)]

colnames(abun_SN_median)

# colnames(abun_SN_median) <- c("Lipid","C_M1","C_M2","C_M3","C_M4","C_M5","C_M6",
#                              "C_F1","C_F2","C_F3","C_F4","C_F5","C_F6",
#                              "AUD_M1","AUD_M2","AUD_M3","AUD_M4","AUD_M5","AUD_M6",
#                              "AUD_F1","AUD_F2","AUD_F3","AUD_F4","AUD_F5")

abun_SN_f_median <- abun_SN_median[,c("Lipid",
                                       "MF_KO_ET1","MF_KO_ET2", "MF_KO_ET3" ,"MF_KO_ET4" ,"MF_KO_ET5","MF_KO_ET6",
                                       "MM_KO_ET1","MM_KO_ET2", "MM_KO_ET3" ,"MM_KO_ET4" ,"MM_KO_ET5","MM_KO_ET6",
                                       
                                       "MF_WT_ET1","MF_WT_ET2", "MF_WT_ET3" ,"MF_WT_ET4" ,"MF_WT_ET5","MF_WT_ET6",
                                       "MM_WT_ET1","MM_WT_ET2", "MM_WT_ET3" ,"MM_WT_ET4" ,"MM_WT_ET5","MM_WT_ET6",
                                       
                                       "MF_KO_C1","MF_KO_C2", "MF_KO_C3" ,"MF_KO_C4" ,"MF_KO_C5","MF_KO_C6",
                                       "MM_KO_C1","MM_KO_C2", "MM_KO_C3" ,"MM_KO_C4" ,"MM_KO_C5","MM_KO_C6",
                                       
                                       "MF_WT_C1","MF_WT_C2", "MF_WT_C3" ,"MF_WT_C4" ,"MF_WT_C5","MF_WT_C6",
                                       "MM_WT_C1","MM_WT_C2", "MM_WT_C3" ,"MM_WT_C4" ,"MM_WT_C5","MM_WT_C6")]

# This is the file that we upload to LINEX2 web tool
write.csv(abun_SN_f_median, file = "05LINEX/abun_SN_f_median.csv",
          row.names = F)


