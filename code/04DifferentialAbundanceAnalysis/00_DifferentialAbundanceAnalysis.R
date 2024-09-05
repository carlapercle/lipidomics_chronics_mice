# Working directory: code

setwd("your/path/code")  

# example: setwd("~/lipidomics/raton/000chronics/code")

rm (list = ls ())

# Lipid abundance levels between groups were compared using the limma R package.
# P-values were adjusted using the Benjamini & Hochberg (BH) procedure, and
# significant lipids were considered when the BH-adjusted p-value ≤ 0.05.

# To perform the analysis we need an object (matrix or data frame) with the
# lipid abundance values and a data frame object with the sample characteristics

# It is important to check that the samples of the two variables follow the same
# order.

# Contrasts:

# A comparison between ethanol-exposed and untreated mice was performed in both 
# sexes (F: female and M: male) for the two genotypes of mice, resulting in the 
# following comparisons, where IF and IM indicate the impact of the ethanol 
# exposure in females and males, respectively:

# · Ethanol exposure vs. control in female WT mice (IF WT)
# · Ethanol exposure vs. control in male WT mice (IM WT)
# · Ethanol exposure vs. control in female TLR4-KO mice (IF KO)
# · Ethanol exposure vs. control in male TLR4-KO mice (IM KO)

# NOTE: 
# The statistics used to measure the differential patterns were the logarithm of
# fold change (LFC) to quantify the effect of differential lipid abundance 
# analysis. A positive statistical sign indicates a higher mean for the variable
# in the first element of the comparison, whereas a negative statistical sign 
# indicates a higher mean value for the second element. In the study, the first 
# element corresponds to ethanol exposure, and the second element is used to
# control.

#_______________________________________________________________________________

# Load packages
library(knitr)
library(stringr)
library(limma)
library(dplyr)
library(DT)

# Function

dif_exp <- function(GSEmat,
                    # df of lipid abundance values
                    contraste,
                    # comparison (contrast)
                    grupos = NULL,
                    # variable with the contrast
                    lotes = NULL,
                    # bacth effect variable
                    trend = FALSE) {
  grupos <- as.factor(grupos) # groups have to be a factor
  
  if (is.null(lotes)) {
    design <- model.matrix( ~ 0 + grupos)
    colnames(design) <- levels(grupos)
  } else {
    lotes <- as.factor(lotes)
    design <- model.matrix( ~ 0 + grupos + lotes)
    colnames(design) <- c(levels(grupos), levels(lotes)[-2])
    #colnames(design) <- c(levels(grupos), levels(lotes))
  }
  
  cont.matrix <-
    makeContrasts(contrasts = contraste, levels = design)
  
  fit <- lmFit(GSEmat, design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  
  genes_limma.res <- eBayes(fit2, trend = trend)
  
  return(genes_limma.res)
}

## Load data ####

### Lipid abundance values  ####

load("00Data/chronic_mice_wo.RData")

### Metadata  ####
samplename <- colnames(df_M_c_wo)

interaccion <- c(substr(colnames(df_M_c_wo),1,nchar(colnames(df_M_c_wo))-1))


datos <- as.data.frame(cbind(samplename, interaccion), as.is = T)

### Check   ####
all(colnames(df_M_c_wo) == datos$samplename)

## IF WT  ####
AD_IF_WT <-
  dif_exp(df_M_c_wo, "(MF_WT_ET - MF_WT_C)", datos$interaccion)

save(AD_IF_WT,
     file = "04DifferentialAbundanceAnalysis/results/AD_IF_WT.Rdata")

AD_IF_WT_sig <-
  topTable(AD_IF_WT , number = Inf, p.value = 0.05)

save(AD_IF_WT_sig,
     file = "04DifferentialAbundanceAnalysis/results/sig/AD_IF_WT_sig.Rdata")

## IM WT  ####
AD_IM_WT <-
  dif_exp(df_M_c_wo, "(MM_WT_ET - MM_WT_C)", datos$interaccion)

save(AD_IM_WT,
     file = "04DifferentialAbundanceAnalysis/results/AD_IM_WT.Rdata")

AD_IM_WT_sig <-
  topTable(AD_IM_WT, number = Inf, p.value = 0.05)

save(AD_IM_WT_sig,
     file = "04DifferentialAbundanceAnalysis/results/sig/AD_IM_WT_sig.Rdata")


## IF KO  ####
AD_IF_KO <-
  dif_exp(df_M_c_wo, "(MF_KO_ET - MF_KO_C)", datos$interaccion)

save(AD_IF_KO,
     file = "04DifferentialAbundanceAnalysis/results/AD_IF_KO.Rdata")

AD_IF_KO_sig <-
  topTable(AD_IF_KO, number = Inf, p.value = 0.05)

save(AD_IF_KO_sig,
     file = "04DifferentialAbundanceAnalysis/results/sig/AD_IF_KO_sig.Rdata")

## IM KO  ####
AD_IM_KO <-
  dif_exp(df_M_c_wo, "(MM_KO_ET - MM_KO_C)", datos$interaccion)

save(AD_IM_KO,
     file = "04DifferentialAbundanceAnalysis/results/AD_IM_KO.Rdata")

AD_IM_KO_sig <-
  topTable(AD_IM_KO, number = Inf, p.value = 0.05)

save(AD_IM_KO_sig,
     file = "04DifferentialAbundanceAnalysis/results/sig/AD_IM_KO_sig.Rdata")

