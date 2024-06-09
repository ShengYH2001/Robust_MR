perform_mr <- function(exposure_data, outcome_data){
  # Extract the name and ID of exposure and outcome data
  exposure <- exposure_data$exposure[1]
  id.exposure <- exposure_data$id.exposure[1]
  outcome <- outcome_data$outcome[1]
  id.outcome <- outcome_data$id.outcome[1]
  # Set up strings and dataframes for storing the result information
  res_str <- "-" # A string used to indicate whether there is a positive result, with the content being '+' or '-'
  res_df <- data.frame( # A dataframe used to store the result
    Exposure = exposure,
    Exposure_ID = id.exposure,
    Outcome = outcome,
    Outcome_ID = id.outcome,
    nSNP = NA,
    IVW.OR_CI95 = NA,
    IVW.Pval = NA,
    IVW_MRE.OR_CI95 = NA,
    IVW_MRE.Pval = NA,
    MR_Egger.OR_CI95 = NA,
    MR_Egger.Pval = NA,
    Wald_Ratio.OR_CI95 = NA,
    Wald_Ratio.Pval = NA,
    MR_Result = 0,
    Leave_one_out.Pval = NA,
    Pleiotropy.Egger_intercept = NA,
    Pleiotropy.SE = NA,
    Pleiotropy.Pval = NA,
    Heterogeneity.IVW.Q = NA,
    Heterogeneity.IVW.Q_DF = NA,
    Heterogeneity.IVW.Q_Pval = NA,
    Heterogeneity.MR_egger.Q = NA,
    Heterogeneity.MR_egger.Q_DF = NA,
    Heterogeneity.MR_egger.Q_Pval = NA,
    MR_Presso.Global_Test.RSSobs = NA,
    MR_Presso.Global_Test.Pval = NA
  )
  run_df <- data.frame( # A dataframe used to store the situation of analysis
    Exposure = exposure,
    Exposure_ID = id.exposure,
    Outcome = outcome,
    Outcome_ID = id.outcome,
    Status = NA
  )
  tryCatch({
    # Save information of SNP used for MR analysis
    snp_output <- exposure_data[,c("exposure","id.exposure","SNP","chr.exposure","pos.exposure","effect_allele.exposure","other_allele.exposure","beta.exposure","eaf.exposure","se.exposure","pval.exposure","samplesize.exposure","rsq.exposure","F")]
    write.table(snp_output, file="./result/table.SNP.csv", col.names =!file.exists("./result/table.SNP.csv"), row.names=F, sep=',', append = T)
    # Harmonise data
    har_data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)
    # Save information of SNP harmonised
    har_data_output <- har_data[,c("exposure","id.exposure","SNP","chr.exposure","pos.exposure","effect_allele.exposure","other_allele.exposure","beta.exposure","eaf.exposure","se.exposure","pval.exposure","samplesize.exposure","rsq.exposure","F")]
    write.table(har_data_output, file="./result/table.SNP_Harmonise.csv", col.names = !file.exists("./result/table.SNP_Harmonise.csv"), row.names=F, sep=',', append = T)
    # Perform MR analysis and estimate 95% CI
    res <- mr(har_data, method_list = c('mr_wald_ratio','mr_egger_regression','mr_weighted_median','mr_ivw','mr_ivw_mre','mr_simple_mode','mr_weighted_mode'))
    odds <- generate_odds_ratios(res)
    odds <- transform(odds, or_ci95 = sprintf("%.3f~%.3f", or_lci95, or_uci95))
    # Save results of MR analysis into matrix form
    res_df$nSNP <- odds[1,"nsnp"]
    res_df$IVW.OR_CI95 = ifelse(length(odds[odds$method == "Inverse variance weighted", ]$or_ci95) == 0, "NA", odds[odds$method == "Inverse variance weighted", ]$or_ci95)
    res_df$IVW.Pval =ifelse(length(odds[odds$method == "Inverse variance weighted", ]$pval) == 0, "NA", odds[odds$method == "Inverse variance weighted", ]$pval)
    res_df$IVW_MRE.OR_CI95 = ifelse(length(odds[odds$method == "Inverse variance weighted (multiplicative random effects)", ]$or_ci95) == 0, "NA", odds[odds$method == "Inverse variance weighted (multiplicative random effects)", ]$or_ci95)
    res_df$IVW_MRE.Pval =ifelse(length(odds[odds$method == "Inverse variance weighted (multiplicative random effects)", ]$pval) == 0, "NA", odds[odds$method == "Inverse variance weighted (multiplicative random effects)", ]$pval)
    res_df$MR_Egger.OR_CI95 = ifelse(length(odds[odds$method == "MR Egger", ]$or_ci95) == 0, "NA", odds[odds$method == "MR Egger", ]$or_ci95)
    res_df$MR_Egger.Pval = ifelse(length(odds[odds$method == "MR Egger", ]$pval) == 0, "NA", odds[odds$method == "MR Egger", ]$pval)
    res_df$Wald_Ratio.OR_CI95 = ifelse(length(odds[odds$method == "Wald ratio", ]$or_ci95) == 0, "NA", odds[odds$method == "Wald ratio", ]$or_ci95)
    res_df$Wald_Ratio.Pval = ifelse(length(odds[odds$method == "Wald ratio", ]$pval) == 0, "NA", odds[odds$method == "Wald ratio", ]$pval)
    if (res_df$IVW.Pval != "NA" & res_df$IVW.Pval < 0.05) {res_df$MR_Result <- 1; res_str <- "+"}
    if (res_df$Wald_Ratio.Pval != "NA" & res_df$Wald_Ratio.Pval < 0.05) {res_df$MR_Result <- 1; res_str <- "+"}
    # Save all OR results of each method of MR analysis
    odds_df <- odds[,c("exposure","id.exposure","outcome","id.outcome","nsnp","method","or","or_lci95","or_uci95","or_ci95","pval","b","se","lo_ci","up_ci")]
    colnames(odds_df) <- c('Exposure','Exposure_ID','Outcome','Outcome_ID','nSNP','Method','OR','OR.Low_CI95','OR.Up_CI95','OR_CI95','Pval','Beta','SE','Low_CI','Up_CI')
    odds_df$MR_Result <- res_df$MR_Result
    write.table(odds_df, file="./result/table.MR_result.csv", col.names = !file.exists("./result/table.MR_result.csv"), row.names=F, sep=',', append = T)
    # Save scatter plot
    svg(sprintf(".\\result\\figure.%s→%s%s.scatter.svg", exposure, outcome, res_str), width = 4, height = 4)
    fig <- mr_scatter_plot(res, har_data)
    print(fig)
    dev.off()
    # Save forest plot & funnel plot
    res_single <- mr_singlesnp(har_data)
    svg(sprintf(".\\result\\figure.%s→%s%s.forest.svg", exposure, outcome, res_str), width = 4, height = 4)
    fig <- mr_forest_plot(res_single)
    print(fig)
    dev.off()
    svg(sprintf(".\\result\\figure.%s→%s%s.funnel.svg", exposure, outcome, res_str), width = 4, height = 4)
    fig <- mr_funnel_plot(res_single)
    print(fig)
    dev.off()
    png(sprintf(".\\result\\figure.%s→%s%s.forest.png", exposure, outcome, res_str), width = 1080, height = 1080, res=300)
    fig <- mr_forest_plot(res_single)
    print(fig)
    dev.off()
    # Save leave-one-out plot
    res_loo <- mr_leaveoneout(har_data)
    svg(sprintf(".\\result\\figure.%s→%s%s.leaveoneout.svg", exposure, outcome, res_str), width = 4, height = 4)
    fig <- mr_leaveoneout_plot(res_loo)
    print(fig)
    dev.off()
    res_df$Leave_one_out.Pval = res_loo[res_loo$SNP == "All", ]$p
    # Pleiotropy test
    res_ple <- mr_pleiotropy_test(har_data)
    res_df$Pleiotropy.Egger_intercept = res_ple$egger_intercept
    res_df$Pleiotropy.SE = res_ple$se
    res_df$Pleiotropy.Pval = res_ple$pval
    # Heterogeneity test
    res_het <- mr_heterogeneity(har_data)
    res_df$Heterogeneity.IVW.Q = res_het[res_het$method == "Inverse variance weighted", ]$Q
    res_df$Heterogeneity.IVW.Q_DF = res_het[res_het$method == "Inverse variance weighted", ]$Q_df
    res_df$Heterogeneity.IVW.Q_Pval = res_het[res_het$method == "Inverse variance weighted", ]$Q_pval
    res_df$Heterogeneity.MR_egger.Q = res_het[res_het$method == "MR Egger", ]$Q
    res_df$Heterogeneity.MR_egger.Q_DF = res_het[res_het$method == "MR Egger", ]$Q_df
    res_df$Heterogeneity.MR_egger.Q_Pval = res_het[res_het$method == "MR Egger", ]$Q_pval
    # Presso
    res_presso <- run_mr_presso(har_data)
    write.table(res_presso, file="./result/table.Presso_result.csv", col.names = !file.exists("./result/table.Presso_result.csv"), row.names=F, sep=',', append = T)
    res_df$MR_Presso.Global_Test.RSSobs = res_presso[[1]]$`MR-PRESSO results`$`Global Test`$RSSobs
    res_df$MR_Presso.Global_Test.Pval = res_presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
    run_df$Status <- "Finish"
  }, error=function(e){
    run_df$Status <- "Error"
  })
  write.table(run_df, "./result/table.MR_run.csv", col.names =!file.exists("./result/table.MR_run.csv"), row.names=F, sep=',', append = T)
  write.table(res_df, file="./result/table.MR_result_df.csv", col.names = !file.exists("./result/table.MR_result_df.csv"), row.names=F, sep=',', append = T)
}
