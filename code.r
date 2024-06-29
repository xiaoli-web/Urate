# Code for the analysis and plot for the study: Sex-specific epigenetic signatures of circulating urate and its increase after BCG vaccination

# Quality control of the DNA methylation data was described in paper: https://github.com/CiiM-Bioinformatics-group/BCG_methylation_project

# Below is the code for perfomeing epigenome-wide association analysis and plotting
# Load required libraries
library(minfi)
library(ggplot2)
library(ggfortify)
library(RColorBrewer)
library(data.table)
library(MASS)
library(sandwich)
library(lmtest)
library(parallel)
library(R.utils)
library(openxlsx)
library(qqman)
library(dplyr)
library(tidyverse)
library(ggtext)
library(normentR)
library(readxl)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Load M-values and perform outlier removal
M.val <- readRDS("mVals_filtered.rds")

removeOutliers <- function(probes) {
  require(matrixStats)
  if (nrow(probes) < ncol(probes)) warning("expecting probes as rows (long dataset)")
  rowIQR <- rowIQRs(probes, na.rm = TRUE)
  row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = TRUE)
  maskL <- probes < row2575[, 1] - 3 * rowIQR 
  maskU <- probes > row2575[, 2] + 3 * rowIQR 
  initial_NAs <- rowSums(is.na(probes))
  probes[maskL] <- NA
  removed_lower <- rowSums(is.na(probes)) - initial_NAs
  probes[maskU] <- NA
  removed_upper <- rowSums(is.na(probes)) - removed_lower - initial_NAs
  N_for_probe <- rowSums(!is.na(probes))
  Log <- data.frame(initial_NAs, removed_lower, removed_upper, N_for_probe)
  return(list(probes, Log))
}

system.time(OutlierResults <- removeOutliers(M.val))
M.val <- OutlierResults[[1]]
Log <- OutlierResults[[2]]

save(M.val, file = "input/mVals_filtered_trimmed.Rdata")
save(Log, file = "input/mVals_filtered_trimmed_log.Rdata")

# Load required data for EWAS
load("input/m1.trim.rdata")
load("input/uv1.rdata")
load("input/pdv1.rdata")

# Prepare phenotype data
pdv1$age <- as.numeric(pdv1$age)
mv1.t <- t(m1.trim)

# Define the robust linear model (RLM) test function
RLMtest <- function(meth_matrix, methcol, urate, age, gender, plate, CD8, CD4, NK, B, Mono, Neu) {
  mod <- try(rlm(meth_matrix[, methcol] ~ urate + age + gender + plate + CD8 + CD4 + NK + B + Mono + Neu, maxit = 200))
  cf <- try(coeftest(mod, vcov = vcovHC(mod, type = "HC0")))
  
  if (class(cf) == "try-error") {
    bad <- as.numeric(rep(NA, 3))
    names(bad) <- c("Estimate", "Std. Error", "Pr(>|z|)")
    bad
  } else {
    cf[2, c("Estimate", "Std. Error", "Pr(>|z|)")]
  }
}

# Perform the EWAS
res <- mclapply(setNames(seq_len(ncol(mv1.t)), dimnames(mv1.t)[[2]]), RLMtest, 
                meth_matrix = mv1.t, urate = uv1$con, age = pdv1$age, 
                gender = pdv1$gender, plate = pdv1$Sample_Plate, 
                CD8 = pdv1$CD8T, CD4 = pdv1$CD4T, NK = pdv1$NK, 
                B = pdv1$Bcell, Mono = pdv1$Mono, Neu = pdv1$Neu)

setattr(res, 'class', 'data.frame')
setattr(res, "row.names", c(NA_integer_, 4))
setattr(res, "names", make.names(names(res), unique = TRUE))
probelistnamesB <- names(res)
result <- t(data.table(res))
result <- data.table(result)
result[, probeID := probelistnamesB]
setnames(result, c("BETA", "SE", "P_VAL", "probeID"))
setcolorder(result, c("probeID", "BETA", "SE", "P_VAL"))
result$padj <- p.adjust(result$P_VAL, method = "BH")

# Calculate lambda
lambda <- median(qchisq(as.numeric(as.character(result$P_VAL)), df = 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1)
lambda 

# Save EWAS results
write.xlsx(result, file = "output/v1_mval_trimmed_model4_remove_plate11.rlm.xlsx")

# Interaction analysis
RLMtest_interaction <- function(meth_matrix, methcol, urate, age, gender, plate, CD8, CD4, NK, B, Mono, Neu) {
  mod <- try(rlm(meth_matrix[, methcol] ~ urate + age + gender + plate + CD8 + CD4 + NK + B + Mono + Neu + urate * gender, maxit = 200))
  cf <- try(coeftest(mod, vcov = vcovHC(mod, type = "HC0")))
  
  if (class(cf) == "try-error") {
    bad <- as.numeric(rep(NA, 3))
    names(bad) <- c("Estimate", "Std. Error", "Pr(>|z|)")
    bad
  } else {
    cf[20, c("Estimate", "Std. Error", "Pr(>|z|)")]
  }
}

res <- lapply(setNames(seq_len(ncol(mv1.t)), dimnames(mv1.t)[[2]]), RLMtest_interaction, 
              meth_matrix = mv1.t, urate = uv1$con, age = pdv1$age, 
              gender = pdv1$gender, plate = pdv1$Sample_Plate, 
              CD8 = pdv1$CD8T, CD4 = pdv1$CD4T, NK = pdv1$NK, 
              B = pdv1$Bcell, Mono = pdv1$Mono, Neu = pdv1$Neu)

setattr(res, 'class', 'data.frame')
setattr(res, "row.names", c(NA_integer_, 4))
setattr(res, "names", make.names(names(res), unique = TRUE))
probelistnamesB <- names(res)
result <- t(data.table(res))
result <- data.table(result)
result[, probeID := probelistnamesB]
setnames(result, c("BETA", "SE", "P_VAL", "probeID"))
setcolorder(result, c("probeID", "BETA", "SE", "P_VAL"))
result$padj <- p.adjust(result$P_VAL, method = "BH")

lambda <- median(qchisq(as.numeric(as.character(result$P_VAL)), df = 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1)
lambda 

write.xlsx(result, file = "output/v1_m4_interaction.rlm.xlsx")

# Sex-stratified analysis (Female)
female <- pdv1 %>% filter(sex == "female") %>% pull(id)
uv1.female <- uv1[which(uv1$id %in% female), ]
pdv1 <- pdv1 %>% filter(sex == "female")
mv1.t <- mv1.t[which(paste0("X", rownames(mv1.t)) %in% uv1.female$id2), ]

pdv1$Sample_Plate <- droplevels(pdv1$Sample_Plate)
str(pdv1)

sum(is.na(uv1))
sum(is.na(mv1.t))
sum(paste0("X", rownames(mv1.t)) == pdv1$id2)
sum(paste0("X", rownames(mv1.t)) == uv1.female$id2)

RLMtest_sex_stratified <- function(meth_matrix, methcol, urate, age, plate, CD8, CD4, NK, B, Mono, Neu, smoker) {
  mod <- try(rlm(meth_matrix[, methcol] ~ urate + age + plate + CD8 + CD4 + NK + B + Mono + Neu + smoker, maxit = 200))
  cf <- try(coeftest(mod, vcov = vcovHC(mod, type = "HC0")))
  
  if (class(cf) == "try-error") {
    bad <- as.numeric(rep(NA, 3))
    names(bad) <- c("Estimate", "Std. Error", "Pr(>|z|)")
    bad
  } else {
    cf[2, c("Estimate", "Std. Error", "Pr(>|z|)")]
  }
}

res <- lapply(setNames(seq_len(ncol(mv1.t)), dimnames(mv1.t)[[2]]), RLMtest_sex_stratified, 
              meth_matrix = mv1.t, urate = uv1.female$con, age = pdv1$age, 
              plate = pdv1$Sample_Plate, CD8 = pdv1$CD8T, CD4 = pdv1$CD4T, 
              NK = pdv1$NK, B = pdv1$Bcell, Mono = pdv1$Mono, Neu = pdv1$Neu, 
              smoker = pdv1$smoker)

setattr(res, 'class', 'data.frame')
setattr(res, "row.names", c(NA_integer_, 4))
setattr(res, "names", make.names(names(res), unique = TRUE))
probelistnamesB <- names(res)
result <- t(data.table(res))
result <- data.table(result)
result[, probeID := probelistnamesB]
setnames(result, c("BETA", "SE", "P_VAL", "probeID"))
setcolorder(result, c("probeID", "BETA", "SE", "P_VAL"))
result$padj <- p.adjust(result$P_VAL, method = "BH")

write.xlsx(result, file = "output/rlm_v1_female_mval_trimmed_methy_dependent_smoker.xlsx")

# Manhattan plot
gwas_data_load <- read_excel("output/EWAS_result.xlsx")
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annotdata <- data.frame(ann850k[c('chr', 'pos', 'UCSC_RefGene_Name', 'Relation_to_Island')])
annotdata$id <- rownames(annotdata)

anno <- merge(gwas_data_load, annotdata, by.x = "probeID", by.y = "id")
anno$chr <- as.numeric(substr(anno$chr, 4, 5))
anno$pos <- as.numeric(anno$pos)

sig_data <- anno %>% filter(padj < 0.05)
notsig_data <- anno %>% filter(padj >= 0.05) %>% group_by(chr) %>% sample_frac(0.1)
gwas_data <- bind_rows(sig_data, notsig_data)

data_cum <- gwas_data %>% group_by(chr) %>% summarise(max_bp = max(pos)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% select(chr, bp_add)

gwas_data <- gwas_data %>% inner_join(data_cum, by = "chr") %>% 
  mutate(bp_cum = pos + bp_add)

axis_set <- gwas_data %>% group_by(chr) %>% summarize(center = mean(bp_cum))
ylim <- gwas_data %>% filter(P_VAL == min(P_VAL)) %>% 
  mutate(ylim = abs(floor(log10(P_VAL))) + 2) %>% pull(ylim)

sig <- 5e-08 # Modified based on the p-value that FDR = 0.05

manhplot <- ggplot() +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") +
  geom_hline(yintercept = log10(sig), color = "grey40", linetype = "dashed") +
  geom_point(data = gwas_data %>% filter(padj >= 0.05), 
             aes(x = bp_cum, y = -log10(P_VAL) * sign(BETA), color = as_factor(chr), size = -log10(P_VAL)), alpha = 0.75) +
  geom_point(data = gwas_data %>% filter(padj < 0.05), 
             aes(x = bp_cum, y = -log10(P_VAL) * sign(BETA), size = -log10(P_VAL)), color = "#276FBF", alpha = 0.75) +
  geom_hline(yintercept = 0, color = "grey40") +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_color_manual(values = rep(c("grey", "lightgrey"), unique(length(axis_set$chr)))) +
  scale_size_continuous(range = c(0.5, 3)) +
  labs(x = NULL, y = "-log<sub>10</sub>(p) * sign(coef)") + 
  theme_minimal() +
  theme(legend.position = "none", panel.border = element_blank(), 
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
        axis.title.y = element_markdown(), axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5))

pdf("plot/mahantan_plot.pdf", width = 7, height = 3.5)
print(manhplot)
dev.off()





