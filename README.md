# perform_mr()
This R function, perform_mr, is designed to conduct Mendelian Randomization (MR) analysis, a statistical method used to assess causal relationships by utilizing genetic variants as instrumental variables. Key features of this function include:
1. Data Integration: It extracts necessary information from provided exposure and outcome datasets and harmonizes the data for analysis.
2. Multiple MR Methods: Performs MR analysis using various methods such as the Wald ratio, Egger regression, weighted median, inverse variance weighted (IVW), multiplicative random effects model for IVW, and both simple and weighted modes.
3. Result Output: Saves the results of the MR analysis into a data frame, detailing odds ratios, confidence intervals, and p-values for different methods.
4. Graphical Output: Generates scatter plots, forest plots, funnel plots, and leave-one-out plots, saving them in both SVG and PNG formats.
5. Error Handling: Employs tryCatch to capture and handle potential errors that may occur during the MR analysis process.
6. Heterogeneity and Pleiotropy Tests: Conducts tests for heterogeneity and pleiotropy to assess the robustness of the MR findings.
7. MR-PRESSO Sensitivity Analysis: Applies the MR-PRESSO method for sensitivity analysis to detect potential outliers or data manipulation.
8. Data Saving: Stores information about SNPs, harmonized data, MR analysis results, graphical outputs, and the status of the run in CSV files.
This function
