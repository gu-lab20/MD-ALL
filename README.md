# ShinyALL
In the quantification and classification of B-cell Acute Lymphoblastic Leukemia (B-ALL), a total of approximately 20,000 gene characteristics have been recorded along with the categorization of 27 different subtypes. In order for future predictions to be made based on patient samples and to improve visualization of gene expression profiles that suit each subtype, a t-distributed stochastic neighbor embedding (t-SNE) plot is required. This model will serve to compress and transform high-dimensional data such as that from RNA-seq into a two or three-dimensional map, where each data point represents one patient sample in the form of an x-y plot. To do this, we have ensured that the median absolute deviation is around 800 for most variable genes so as to achieve distinct clusters of data yet reduce noise simultaneously. The generated plot will be based on RNA-seq performed on 1409 samples, each from a unique B-ALL patient, where 19 majorly different subtypes are noted. An adjustable parameter is included in the form of "Perplexity" in order to visualize the potential patient's data ranging from local to global relativity. This Shinyapp is designed to be coding-free and still easily accessible for proper patient diagnosis (hence, appropriate treatment) as well as for general research use.
###data preprocessing
###


![adf]
