library(TwoSampleMR)
library(MVMR)
library(RMVMR)
library(MendelianRandomization)

perform_mvmr <- function(exposure_datas_fine, exposure_datas, outcome_data) {
  SNP_strong <- c()
  for (i in 1:length(exposure_datas_fine)){
    SNP_strong <- c(SNP_strong, exposure_datas_fine[[i]]$SNP)
  }
  SNP_comb <- intersect(outcome_data$SNP, SNP_strong)
  for (i in 1:length(exposure_datas)){
    SNP_comb <- intersect(SNP_comb, exposure_datas[[i]]$SNP)
  }
  SNP_clump <- ieugwasr::ld_clump(dplyr::tibble(rsid=SNP_comb),
                                          plink_bin = "F:/GWAS/plink_win64_20231211/plink.exe",
                                          bfile = "F:/GWAS/1kg.v3/EUR",
                                          clump_r2 = 0.001,
                                          clump_kb = 10000)$rsid
  SNP_df <- data.frame("SNP" = SNP_clump)
  outcome_beta_colname <- sprintf("outcome.beta.%s", outcome_data$id.outcome[1])
  outcome_se_colname <- sprintf("outcome.se.%s", outcome_data$id.outcome[1])
  SNP_df[,outcome_beta_colname] <- merge(data.frame(SNP=SNP_clump), outcome_data, by="SNP")$beta.outcome
  SNP_df[,outcome_se_colname] <- merge(data.frame(SNP=SNP_clump), outcome_data, by="SNP")$se.outcome
  exposure_beta_colnames <- c()
  exposure_se_colnames <- c()
  for (i in 1:length(exposure_datas)){
    exposure_beta_colnames <- c(exposure_beta_colnames, sprintf("exposure.beta.%s",exposure_datas[[i]]$id.exposure[1]))
    exposure_se_colnames <- c(exposure_se_colnames, sprintf("exposure.se.%s", exposure_datas[[i]]$id.exposure[1]))
    SNP_df[,sprintf("exposure.beta.%s",exposure_datas[[i]]$id.exposure[1])] <- merge(data.frame(SNP=SNP_clump), exposure_datas[[i]], by="SNP")$beta.exposure
    SNP_df[,sprintf("exposure.se.%s", exposure_datas[[i]]$id.exposure[1])] <- merge(data.frame(SNP=SNP_clump), exposure_datas[[i]], by="SNP")$se.exposure
  }
  mv_data_fine <- format_mvmr(BXGs = SNP_df[,exposure_beta_colnames],
                              BYG = SNP_df[,outcome_beta_colname],
                              seBXGs = SNP_df[,exposure_se_colnames],
                              seBYG = SNP_df[,outcome_se_colname],
                              RSID = SNP_df[,"SNP"])
  mvmr_res <- mvmr(r_input= mv_data_fine)
  sres <- strength_mvmr(r_input = mv_data_fine, gencov = 0)
  pres <- pleiotropy_mvmr(r_input = mv_data_fine, gencov = 0)
  res <- ivw_mvmr(r_input = mv_data_fine)
  # mvmrcovmatrix<-matrix(c(1,-0.1,-0.05,-0.1,1,0.2,-0.05,0.2,1), nrow = 3, ncol = 3)
  # res1 <- qhet_mvmr(mv_data_fine, mvmrcovmatrix, CI = F, iterations = 100)
  
  res <- as.data.frame(res)
  result_df <- data.frame()
  for (i in 1:length(exposure_datas_fine)){
    result_df[i,"Exposure"] <- exposure_datas[[i]]$exposure[1]
    result_df[i,"Outcome"] <- outcome_data$outcome[1]
    result_df[i,"F_Statistic"] <- sres[,sprintf("exposure%s", i)]
    result_df[i,"Q_Strength"] <- mvmr_res$Q_strength[,sprintf("exposure%s", i)]
    result_df[i,"Q_Statistic"] <- pres$Qstat
    result_df[i,"Q_Pval"] <- pres$Qpval
    result_df[i,"Estimate"] <- res[i,"Estimate"]
    result_df[i,"Std. Error"] <- res[i,"Std. Error"]
    result_df[i,"t Value"] <- res[i,"t value"]
    result_df[i,"Pr(>|t|)"] <- res[i,"Pr(>|t|)"]
  }
  write.table(result, file="./result/table.Multi_result.csv", row.names=F, col.names = T, sep=',', append = T)
}
