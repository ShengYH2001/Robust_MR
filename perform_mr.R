library(TwoSampleMR)

# Version: 240928.1
perform_mr <- function(exposure_dat = NULL,
                       outcome_dat = NULL,
                       harmonise_dat = NULL,
                       plot.svg = TRUE,
                       round = 0){
  ifelse(!dir.exists("./result"),dir.create("./result"),FALSE)
  
  # Initialize dataframe
  if (is.null(harmonise_dat)) {
    tryCatch({exposure <- exposure_dat$exposure[1]}, error = function(e) {exposure <<- "NA"})
    tryCatch({id.exposure <- exposure_dat$id.exposure[1]}, error = function(e) {id.exposure <<- "NA"})
    tryCatch({outcome <- outcome_dat$outcome[1]}, error = function(e) {outcome <<- "NA"})
    tryCatch({id.outcome <- outcome_dat$id.outcome[1]}, error = function(e) {id.outcome <<- "NA"})
  } else {
    tryCatch({exposure <- harmonise_dat$exposure[1]}, error = function(e) {exposure <<- "NA"})
    tryCatch({id.exposure <- harmonise_dat$id.exposure[1]}, error = function(e) {id.exposure <<- "NA"})
    tryCatch({outcome <- harmonise_dat$outcome[1]}, error = function(e) {outcome <<- "NA"})
    tryCatch({id.outcome <- harmonise_dat$id.outcome[1]}, error = function(e) {id.outcome <<- "NA"})
  }
  res_str <- "-"
  res_df <- data.frame(
    Exposure = exposure,
    Exposure_ID = id.exposure,
    Outcome = outcome,
    Outcome_ID = id.outcome,
    nSNP = "NA",
    Main.Method = "NA",
    Main.Beta = "NA",
    Main.SE = "NA",
    Main.OR = "NA",
    Main.OR_CI95 = "NA",
    Main.Pval = "NA",
    Leave_one_out.Pval = "NA",
    Pleiotropy.Pval = "NA",
    Heterogeneity.IVW.Q_Pval = "NA",
    Heterogeneity.MR_egger.Q_Pval = "NA",
    MR_Presso.Global_Test.Pval = "NA",
    MR_Result = 0,
    Main.OR_CI95_Low = "NA",
    Main.OR_CI95_Up = "NA",
    IVW.Beta = "NA",
    IVW.SE = "NA",
    IVW.OR = "NA",
    IVW.OR_CI95 = "NA",
    IVW.OR_CI95_Low = "NA",
    IVW.OR_CI95_Up = "NA",
    IVW.Pval = "NA",
    IVW_MRE.Beta = "NA",
    IVW_MRE.SE = "NA",
    IVW_MRE.OR = "NA",
    IVW_MRE.OR_CI95 = "NA",
    IVW_MRE.OR_CI95_Low = "NA",
    IVW_MRE.OR_CI95_Up = "NA",
    IVW_MRE.Pval = "NA",
    MR_Egger.Beta = "NA",
    MR_Egger.SE = "NA",
    MR_Egger.OR = "NA",
    MR_Egger.OR_CI95 = "NA",
    MR_Egger.OR_CI95_Low = "NA",
    MR_Egger.OR_CI95_Up = "NA",
    MR_Egger.Pval = "NA",
    Wald_Ratio.Beta = "NA",
    Wald_Ratio.SE = "NA",
    Wald_Ratio.OR = "NA",
    Wald_Ratio.OR_CI95 = "NA",
    Wald_Ratio.OR_CI95_Low = "NA",
    Wald_Ratio.OR_CI95_Up = "NA",
    Wald_Ratio.Pval = "NA",
    Pleiotropy.Egger_intercept = "NA",
    Pleiotropy.SE = "NA",
    Heterogeneity.IVW.Q = "NA",
    Heterogeneity.IVW.Q_DF = "NA",
    Heterogeneity.MR_egger.Q = "NA",
    Heterogeneity.MR_egger.Q_DF = "NA",
    MR_Presso.Global_Test.RSSobs = "NA",
    Round = round
  )
  run_df <- data.frame(
    Exposure = exposure,
    Exposure_ID = id.exposure,
    Outcome = outcome,
    Outcome_ID = id.outcome,
    Status = "NA",
    Round = round
  )
  tryCatch({
    # Save information of SNP used for MR analysis
    if (is.null(harmonise_dat)) {
      exposure_dat_output <- exposure_dat[,c("exposure","id.exposure","SNP","chr.exposure","pos.exposure","effect_allele.exposure","other_allele.exposure","beta.exposure","eaf.exposure","se.exposure","pval.exposure","samplesize.exposure","rsq.exposure","F")]
      write.table(exposure_dat_output, file="./result/table.SNP.csv", col.names =!file.exists("./result/table.SNP.csv"), row.names=F, sep=',', append = T)
      # Harmonise data
      harmonise_dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
      harmonise_dat <- harmonise_dat[which(harmonise_dat$mr_keep==TRUE),]
    }
    # Perform MR analysis and estimate 95% CI
    res <- mr(harmonise_dat, method_list = c('mr_wald_ratio','mr_egger_regression','mr_weighted_median','mr_ivw','mr_ivw_mre','mr_simple_mode','mr_weighted_mode'))
    odds <- generate_odds_ratios(res)
    odds <- transform(odds, or_ci95 = sprintf("%.2f (%.2f~%.2f)", or, or_lci95, or_uci95))
    
    # Save results of MR analysis into matrix form
    res_df$nSNP <- odds[1,"nsnp"]
    
    if (length(odds[odds$method == "Inverse variance weighted", ]$or) > 0) {
      res_df$IVW.Beta <- odds[odds$method == "Inverse variance weighted", ]$b
      res_df$IVW.SE <- odds[odds$method == "Inverse variance weighted", ]$se
      res_df$IVW.OR <- odds[odds$method == "Inverse variance weighted", ]$or
      res_df$IVW.OR_CI95 <- odds[odds$method == "Inverse variance weighted", ]$or_ci95
      res_df$IVW.OR_CI95_Low <- odds[odds$method == "Inverse variance weighted", ]$or_lci95
      res_df$IVW.OR_CI95_Up <- odds[odds$method == "Inverse variance weighted", ]$or_uci95
      res_df$IVW.Pval <- odds[odds$method == "Inverse variance weighted", ]$pval
    }
    
    if (length(odds[odds$method == "Inverse variance weighted (multiplicative random effects)", ]$or) > 0) {
      res_df$IVW_MRE.Beta <- odds[odds$method == "Inverse variance weighted (multiplicative random effects)", ]$b
      res_df$IVW_MRE.SE <- odds[odds$method == "Inverse variance weighted (multiplicative random effects)", ]$se
      res_df$IVW_MRE.OR <- odds[odds$method == "Inverse variance weighted (multiplicative random effects)", ]$or
      res_df$IVW_MRE.OR_CI95 <- odds[odds$method == "Inverse variance weighted (multiplicative random effects)", ]$or_ci95
      res_df$IVW_MRE.OR_CI95_Low <- odds[odds$method == "Inverse variance weighted (multiplicative random effects)", ]$or_lci95
      res_df$IVW_MRE.OR_CI95_Up <- odds[odds$method == "Inverse variance weighted (multiplicative random effects)", ]$or_uci95
      res_df$IVW_MRE.Pval <- odds[odds$method == "Inverse variance weighted (multiplicative random effects)", ]$pval
    }
    
    if (length(odds[odds$method == "MR Egger", ]$or) > 0) {
      res_df$MR_Egger.Beta <- odds[odds$method == "MR Egger", ]$b
      res_df$MR_Egger.SE <- odds[odds$method == "MR Egger", ]$se
      res_df$MR_Egger.OR <- odds[odds$method == "MR Egger", ]$or
      res_df$MR_Egger.OR_CI95 <- odds[odds$method == "MR Egger", ]$or_ci95
      res_df$MR_Egger.OR_CI95_Low <- odds[odds$method == "MR Egger", ]$or_lci95
      res_df$MR_Egger.OR_CI95_Up <- odds[odds$method == "MR Egger", ]$or_uci95
      res_df$MR_Egger.Pval <- odds[odds$method == "MR Egger", ]$pval
    }
    
    if (length(odds[odds$method == "Wald ratio", ]$or) > 0) {
      res_df$Wald_Ratio.Beta <- odds[odds$method == "Wald ratio", ]$b
      res_df$Wald_Ratio.SE <- odds[odds$method == "Wald ratio", ]$se
      res_df$Wald_Ratio.OR <- odds[odds$method == "Wald ratio", ]$or
      res_df$Wald_Ratio.OR_CI95 <- odds[odds$method == "Wald ratio", ]$or_ci95
      res_df$Wald_Ratio.OR_CI95_Low <- odds[odds$method == "Wald ratio", ]$or_lci95
      res_df$Wald_Ratio.OR_CI95_Up <- odds[odds$method == "Wald ratio", ]$or_uci95
      res_df$Wald_Ratio.Pval <- odds[odds$method == "Wald ratio", ]$pval
    }
    
    # Pleiotropy test
    res_ple <- mr_pleiotropy_test(harmonise_dat)
    if (length(res_ple$egger_intercept) > 0) {
      tryCatch({
        res_df$Pleiotropy.Egger_intercept = res_ple$egger_intercept
        res_df$Pleiotropy.SE = res_ple$se
        res_df$Pleiotropy.Pval = res_ple$pval
      }, error=function(e){e})
    }
    
    # Heterogeneity test
    res_het <- mr_heterogeneity(harmonise_dat)
    if (length(res_het$method) > 0) {
      tryCatch({
        res_df$Heterogeneity.IVW.Q = res_het[res_het$method == "Inverse variance weighted", ]$Q
        res_df$Heterogeneity.IVW.Q_DF = res_het[res_het$method == "Inverse variance weighted", ]$Q_df
        res_df$Heterogeneity.IVW.Q_Pval = res_het[res_het$method == "Inverse variance weighted", ]$Q_pval
        res_df$Heterogeneity.MR_egger.Q = res_het[res_het$method == "MR Egger", ]$Q
        res_df$Heterogeneity.MR_egger.Q_DF = res_het[res_het$method == "MR Egger", ]$Q_df
        res_df$Heterogeneity.MR_egger.Q_Pval = res_het[res_het$method == "MR Egger", ]$Q_pval
      }, error=function(e){e})
    }
    
    # Choose main method
    res_df$Main.Method = ifelse(res_df$nSNP == 1, "Wald_Ratio", ifelse(res_df$Pleiotropy.Pval > 0.05, ifelse(res_df$Heterogeneity.IVW.Q_Pval > 0.05, "IVW", "IVW_MRE"), "MR_Egger"))
    if (is.na(res_df$Main.Method)) {res_df$Main.Method="IVW"}
    res_df$Main.Beta = eval(parse(text = sprintf("res_df$%s.Beta", res_df$Main.Method)))
    res_df$Main.SE = eval(parse(text = sprintf("res_df$%s.SE", res_df$Main.Method)))
    res_df$Main.OR = eval(parse(text = sprintf("res_df$%s.OR", res_df$Main.Method)))
    res_df$Main.OR_CI95 = eval(parse(text = sprintf("res_df$%s.OR_CI95", res_df$Main.Method)))
    res_df$Main.OR_CI95_Low = eval(parse(text = sprintf("res_df$%s.OR_CI95_Low", res_df$Main.Method)))
    res_df$Main.OR_CI95_Up = eval(parse(text = sprintf("res_df$%s.OR_CI95_Up", res_df$Main.Method)))
    res_df$Main.Pval = eval(parse(text = sprintf("res_df$%s.Pval", res_df$Main.Method)))
    
    if (res_df$Main.Pval < 0.05) {res_df$MR_Result <- 1; res_str <- "+"}
    
    # Presso
    res_presso <- run_mr_presso(harmonise_dat)
    write.table(res_presso, file="./result/table.Presso_result.csv", col.names = !file.exists("./result/table.Presso_result.csv"), row.names=F, sep=',', append = T)
    res_df$MR_Presso.Global_Test.RSSobs = res_presso[[1]]$`MR-PRESSO results`$`Global Test`$RSSobs
    res_df$MR_Presso.Global_Test.Pval = res_presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
    run_df$Status <- "Finish"
  }, error=function(e){
    run_df$Status <<- "Error"
    if (exists("res_presso")) {
      try({res_df$MR_Presso.Global_Test.RSSobs <<- res_presso[[1]]$`MR-PRESSO results`$`Global Test`$RSSobs}, silent = TRUE)
      try({res_df$MR_Presso.Global_Test.Pval <<- res_presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue}, silent = TRUE)
      run_df$Status <<- "Finish"
    }
    
  })
  
  tryCatch({
    # Save forest plot & funnel plot
    res_single <- mr_singlesnp(harmonise_dat)
    # Save PNG
    png(sprintf(".\\result\\figure.%s→%s%s.forest.png", make.names(exposure), make.names(outcome), res_str), width = 1080, height = 1080, res=300)
    fig <- mr_forest_plot(res_single); print(fig); dev.off()
    
    # Save leave-one-out plot
    res_loo <- mr_leaveoneout(harmonise_dat)
    res_df$Leave_one_out.Pval = res_loo[res_loo$SNP == "All", ]$p
    if (plot.svg) {
      svg(sprintf(".\\result\\figure.%s→%s%s.scatter.svg", make.names(exposure), make.names(outcome), res_str), width = 6, height = 6)
      fig <- mr_scatter_plot(res, harmonise_dat); print(fig); dev.off()
      svg(sprintf(".\\result\\figure.%s→%s%s.forest.svg", make.names(exposure), make.names(outcome), res_str), width = 6, height = max(6, 0.1*res_df$nSNP))
      fig <- mr_forest_plot(res_single); print(fig); dev.off()
      svg(sprintf(".\\result\\figure.%s→%s%s.funnel.svg", make.names(exposure), make.names(outcome), res_str), width = 6, height = 6)
      fig <- mr_funnel_plot(res_single); print(fig); dev.off()
      svg(sprintf(".\\result\\figure.%s→%s%s.leaveoneout.svg", make.names(exposure), make.names(outcome), res_str), width = 6, height = 6)
      fig <- mr_leaveoneout_plot(res_loo); print(fig); dev.off()
    }
  }, error=function(e) {})
  
  write.table(run_df, "./result/table.MR_run.csv", col.names =!file.exists("./result/table.MR_run.csv"), row.names=F, sep=',', append = T)
  write.table(res_df, file="./result/table.MR_result_df_all.csv", col.names = !file.exists("./result/table.MR_result_df_all.csv"), row.names=F, sep=',', append = T)
  
  if (res_df$Pleiotropy.Pval != "NA" & res_df$Pleiotropy.Pval < 0.05 & exists("res_presso")) {
    if (any(res_presso[[1]]$`MR-PRESSO results`$`Outlier Test`$Pvalue < 0.05)) {
      harmonise_dat <- harmonise_dat[res_presso[[1]]$`MR-PRESSO results`$`Outlier Test`$Pvalue > 0.05,]
      perform_mr(harmonise_dat = harmonise_dat, plot.svg = plot.svg, round = round+1)
      return()
    }
    if (any(res_presso[[1]]$`MR-PRESSO results`$`Outlier Test`$Pvalue < 1)) {
      harmonise_dat <- harmonise_dat[res_presso[[1]]$`MR-PRESSO results`$`Outlier Test`$Pvalue == 1,]
      perform_mr(harmonise_dat = harmonise_dat, plot.svg = plot.svg, round = round+1)
      return()
    }
  }
  
  tryCatch({
    write.table(res_df, file="./result/table.MR_result_df.csv", col.names = !file.exists("./result/table.MR_result_df.csv"), row.names=F, sep=',', append = T)
    # Save information of SNP harmonized
    harmonise_dat_output <- harmonise_dat[,c("exposure","id.exposure","outcome","id.outcome","SNP","chr.exposure","pos.exposure","effect_allele.exposure","other_allele.exposure","beta.exposure","eaf.exposure","se.exposure","pval.exposure","samplesize.exposure","rsq.exposure","F","effect_allele.outcome","other_allele.outcome","beta.outcome","eaf.outcome","se.outcome","pval.outcome","samplesize.outcome")]
    write.table(harmonise_dat_output, file="./result/table.SNP_Harmonise.csv", col.names = !file.exists("./result/table.SNP_Harmonise.csv"), row.names=F, sep=',', append = T)
    # Save all OR results of each method of MR analysis
    odds_df <- odds[,c("exposure","id.exposure","outcome","id.outcome","nsnp","method","or","or_lci95","or_uci95","or_ci95","pval","b","se","lo_ci","up_ci")]
    colnames(odds_df) <- c('Exposure','Exposure_ID','Outcome','Outcome_ID','nSNP','Method','OR','OR.Low_CI95','OR.Up_CI95','OR_CI95','Pval','Beta','SE','Low_CI','Up_CI')
    odds_df$MR_Result <- res_df$MR_Result
    write.table(odds_df, file="./result/table.MR_result.csv", col.names = !file.exists("./result/table.MR_result.csv"), row.names=F, sep=',', append = T)
  }, error=function(e) {})
}
