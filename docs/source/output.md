# Output

Pydamage generates both a tabular and a visual output.

The tabular output is a comma-separated file (`.csv`) with the following information for all analysed references:

  * `reference`: name of the reference genome/contig
  * `null_model_p0`: parameter `p0` of the null model
  * `null_model_p0_stdev`: standard error of the null model paramater `p0`
  * `damage_model_p`: parameter `p` of the damage model
  * `damage_model_p_stdev`: standard error of the parameter `p` of the damage model
  * `damage_model_pmin`: paramater `p_min` of the damage model
  * `damage_model_pmin_stdev`: standard error of the paramater `p_min` of the damage model
  * `damage_model_pmax`: paramater `p_max` of the damage model
  * `damage_model_pmax_stdev`: standard error of the paramater `p_max` of the damage model
  * `pvalue`: p-value calculated from the likelihood-ratio test-statistic using a chi-squared distribution
  * `qvalue`: p-value corrected for multiple testing using Benjamini-Hochberg procedure
  * `RMSE`: residual mean standard error of the model fit of the damage model
  * `nb_reads_aligned`: number of aligned reads
  * `coverage`: average coverage along the reference genome

Finally, the remaining columns indicate the frequency of C to T and G to A transitions at every read position observed in the sequencing data.

The visual output are PNG files, one per reference contig. They show the frequency of observed C to T, and G to A transitions at the 5' end of the sequencing data and overlay it with the fitted models for both the null and the damage model, including 95% confidence intervals. Furthermore, it provides a "residuals versus fitted" plot to help evaluate the fit of the pydamage damage model. Finally, the plot contains informtion on the average coverage along the reference and the p-value calculated from the likelihood-ratio test-statistic using a chi-squared distribution.

> The visual output is only produced when using the `--plot` flag 