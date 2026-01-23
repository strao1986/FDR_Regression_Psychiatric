# FDR_Regression_Psychiatry
Here we presented an analytic framework based on the FDR regression (FDRreg) approach to identify additional association signals, by integrating association evidence from GWAS summary statistics of related traits/disorders, and a wealth of biological information, which may be derived from real experimental data. 

For instance, the regression model could incorporate evidence from drug targets or animal models. This statistical method not only effectively controlled the false discovery rate, but also extracted association evidence from multiple correlated traits related to the target phenotype. This significantly increased the likelihood of discovering novel risk loci or associated genes/pathways. 

We validate the FDRreg approach through extensive simulations that confirm proper FDR control, and demonstrate that discoveries from a smaller sample GWAS are replicated in a larger, higher powered study.

Our simulations reveal that pre analysis variable selection, whether by marginal association test or regularized regression, not only maintains FDR control but also greatly improves the computational efficiency. The variable selection approach also outperforms standard FDRreg in terms of FDR control when the number of hypotheses tested is moderate. This opens a new approach for the incorporation of high dimensional covariates in FDR regression.

We first applied the strategy to investigate several common types of psychiatric disorders/traits, while this proposed strategy holds potential for discovering risk loci/genes in other complex diseases. 
