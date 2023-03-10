#################################
Analysis of Metagenomic data in R
#################################

Now that we have formatted out abundance table we can import it in R.

******************
Install R packages
******************

.. highlight:: r


We will use the R package ``Phyloseq`` to explore the microbiome data in the normalized abundance table. Once imported in R, the table can be tweaked to focus on a subset of samples. 
Install Phyloseq from `Bioconductor`_.

  .. _Bioconductor: https://bioconductor.org/packages/release/bioc/html/phyloseq.html


.. _installed before:

Let's first set up the work-directory containing the table that we generated: 
::

  setwd("path/to/folder")

If needed, install the following R packages: 
::
  
  install.packages("phyloseq")
  install.packages("ggplot2")
  install.packages("vegan")
  
The packages **microbiome** and **DESeq2** must be installed via **BiocManager**
::

  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
  BiocManager::install("microbiome")
  BiocManager::install("DESeq2")


***************************************
Importing abundance tables in Phyloseq
***************************************

Phyloseq works with ``biom`` objects that store all the information required for the analyses: abundance of taxa, sample metadata, taxonomy
First you can import the table with the species abundances and the full taxonomic ranks. This table has the taxa in the rows and the samples in the columns.
Note that we will set the column #20, corresponding to the species names, as row names. 
::

  table = as.matrix(read.delim("abundance_table_IDs.merged.norm.taxonomy.final.Archaea_Bacteria", header=T, fill=T, row.names=20, sep="\t"))

From this table, we will extract as separate objects the taxonomy and the species abundances.  
::

  # extract the taxonomy by selecting the corresponding columns.
  taxonomy = (table[,c(3:19)])
  # extract the otu by removing the columns corresponding to the taxonomy and other unwanted columns. 
  otu = (table[,-c(1:23)])
  # make the matrix values as numeric
  class(otu) <- "numeric"

Then you can generate the ``biom`` object with phyloseq. 
::

  library(phyloseq)
  library(ggplot2)
  library(vegan)
  # assing otu and tax tables
  OTU = otu_table(otu, taxa_are_rows = TRUE)
  TAX = tax_table(taxonomy)
  my_biom = phyloseq(OTU, TAX)

As additional step, you can add to the biom file the metadata associated with the samples. These will be useful to describe the samples and select them in the downstream analyses (e.g. to work on subset of samples).
The metadata are stored in a mapping file, a tab-delimited text file where you have in the first column the names of your samples (matching the names in the columns of your OTU abundance table), and in the other columns as many metadata as you want (e.g. sample source, literature ID etc.). 
You can download the metadata of the dataset that we are analysing: `metadata_chimps_2Gb.txt`_

  .. _metadata_chimps_2Gb.txt: https://drive.google.com/file/d/1oedAI5qA_0RsVEL20oe8Yf80WmGnMzYR/view?usp=sharing
  
The metadata can be incorporated in your biom file as follows:
::

  # import mapping file with metadata
  biom_metadata <- import_qiime_sample_data("metadata_chimps_2Gb.txt")
  # merge data
  my_biom <- merge_phyloseq(my_biom, biom_metadata)

Once the biom object has been generated, you can try some commands to explore it: 
::

  # check sample names:
  sample_data(my_biom)
  # check taxa:
  tax_table(my_biom)
  # check rank names:
  colnames(tax_table(my_biom))

.. note::
  If needed, it is convenient to rename the taxa columns more appropriately (some tools to retrieve the full taxonomy may return "Rank1" etc.)
  ::
  
    # check the labels used for rank names
    head(tax_table(my_biom))
    # change them if needed
    colnames(tax_table(my_biom)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

In the next step we'll remove low-abundance taxa below a desired threshold (most commonly 0.02% is chosen).
This is done because the very rare taxa (below the threshold) may be false positive resulting from taxonomic misassignments.
To do that we will create various objects and use a function to remove the taxa below the desired threshold. 
::

  minTotRelAbun = 0.0002
  x = taxa_sums(my_biom)
  keepTaxa = taxa_names(my_biom)[which((x / sum(x)) > minTotRelAbun)]
  my_biom_flt = prune_taxa(keepTaxa, my_biom)

Finally, we convert the absolute abundance values in relative abundance. This works as normalization for the different sequencing depths of the datasets that we are analysing. 
To do that we use the function ``transform_sample_counts`` in Phyloseq
::

  my_biom_rel = transform_sample_counts(my_biom_flt, function(x) x / sum(x))


**************************
Barplot of taxa abundances
**************************
Now that our table is formatted as biom object, we can start to analyse it by generating barplots to display the relative abundance of each sample at different taxonomic ranks. 
To do that we will first merge the abundance of the taxa of our table (e.g. species) to a higher taxonomic rank (e.g. Phylum). We will do this with the function `tax_glom` of Phyloseq. The command will be run on the filtered biom file, and then we will convert phyla abundances to relative frequencies.
::
  
  # merge taxa at the Phylum rank with tax_glom (in Phyloseq)
  my_biom_phylum = tax_glom(my_biom_flt, "phylum", NArm=FALSE)
  # convert to relative freq
  my_biom_phylum_rel = transform_sample_counts(my_biom_phylum, function(x) x / sum(x))
  # plot the relative abundances (you can use ggplot2 themes to improve the chart)
  plot_bar(my_biom_phylum_rel, fill = "phylum") + theme(legend.position="bottom", axis.text.x = element_text(size = 4, angle = 90, hjust = 1))

.. warning::
  In the ``tax_glom`` command the label of the rank must match the one that is defined in your original biom file, as reported in the column names of the tax_table command (see above). 

.. warning::
  You may see that in one sample Kraken2 did not return any taxonomic classification (the bar is empty). This is because the database is too small (2Gb!), hence, there may not be enough resolution to make the taxonomic classification in some instances. 
  Another reason of this result is that in our Bracken analysis, we used a threshold of 50 reads. This means that taxa below that threshold are removed. You can try to use a lower threshold, for example 10. 
  In our dataset, the sample that did not return any taxonomic classification is **Met0199**, which corresponds to the column 85 of our imported table. Let's remove it from the original table in R and create again the biom file. To do that we will use the function `subset_samples` of phyloseq.
  ::
  
    # remove the sample identified as "Met0199" in the "Sample_short" column of the metadata
    my_biom = subset_samples(my_biom, Sample_short != "Met0199")
    # filter again 
    minTotRelAbun = 0.0002
    x = taxa_sums(my_biom)
    keepTaxa = taxa_names(my_biom)[which((x / sum(x)) > minTotRelAbun)]
    my_biom_flt = prune_taxa(keepTaxa, my_biom)
    # convert to relative frequencies again
    my_biom_rel = transform_sample_counts(my_biom_flt, function(x) x / sum(x))
 
   
  You can see the difference in the results by repeating the analyses on the same dataset (including more chimps) generated by using a Kraken2 database of 8Gb (see section 6). 
  

Try to generate barplots at different taxonomic ranks (e.g. family). 
To save the plot in the pdf format you can use this command: 
::

  dev.print(pdf, 'filename.pdf)')

You can also focus your analysis on a subset of samples and use the metadata contained in your original mapping file. It would be convenient to add metadata the way it more convenient to select and compare your samples. Here we subsample for some groups of individuals as defined in the metadata column **Group1**.
::

  # subsample single group phyla
  subsample = subset_samples(my_biom_phylum_rel, Group1=="Dental calculus chimp")
  # subsample multiple groups phyla
  subgroup = c("Dental_calculus_chimp","Dental_calculus","Oral_plaque")
  subsample = subset_samples(my_biom_phylum_rel, Group1 %in% subgroup)
  # barplots in groups with facet_grid
  plot_bar(subsample, fill="phylum") +
  facet_grid(~Group1, scales= "free_x", space = "free_x") +
  theme(legend.position="bottom", axis.text.x = element_text(size = 6, angle = 90, hjust = 1))

Try now to generate barplots at the family level. You will end-up with more labels in the legend, to tweak it play with the options ``legend.position="right"``, ``legend.text = element_text(size = 4)``, ``legend.key.size = unit(0.5,"line")``, inside ``theme`` in ggplot2.  
::
  
  plot_bar(my_biom_family_rel, fill = "family") + 
  theme(legend.position="right", legend.text = element_text(size = 4), legend.key.size = unit(0.5,"line"), axis.text.x = element_text(size = 4, angle = 90, hjust = 1))



***************
Alpha diversity
***************
Alpha diversity (??-diversity) is defined as the mean diversity of species in different sites or habitats within a local scale, it is used to describe the "within-sample" diversity. 
A number of metrics are commonly-used to measure alpha diversity, such as the Shannon index and the Simpson index (read more about `Alpha diversity`_).
It is advisable to use untrimmed, non-normalized count data for meaningful results, as many of alpha diversity estimates are highly dependent on the number of singletons (taxa with count=1). You can always trim the data later on if needed, just not before using this function.

  .. _Alpha diversity: https://docs.onecodex.com/en/articles/4136553-alpha-diversity
  
We will estimate the Shannon index and the Simpson index on our dataset and plot them with the function ``plot_richness`` of Phyloseq. 
For alpha diversity estimation, Phyloseq accepts only integer values, so first we'll round out table to integers and then plot the alpha diversity values for each sample and for group of samples (as identified in the Group1 metadata).
::

  # make the data in the otu table integer
  otu_table(my_biom) = round(otu_table(my_biom))
  # plot two diversity indices as dots for each individual and customize legend and axis-x
  plot_richness(my_biom, measures=c("Simpson", "Shannon"), color="Group1") +
  	theme(legend.position="right", axis.text.x = element_text(size = 4, angle = 90, hjust = 1))
  # plot richness based on group variable
  plot_richness(my_biom, x="Group1", measures=c("Simpson", "Shannon"), color="Group1")
  # plot richness as violin plot
  plot_richness(my_biom, x="Group1", measures=c("Simpson", "Shannon"), color="Group1") +
  	geom_violin()
  # plot richness by including fots in the violin plots
  plot_richness(my_biom, x="Group1", measures=c("Simpson", "Shannon"), color="Group1") +
  	geom_violin() +
  	geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
  # plot richness as boxplots - you can also combine violin and boxplot
  plot_richness(my_biom, x="Group1", measures=c("Simpson", "Shannon"), color="Group1") +
      geom_boxplot(width=0.3)


**************************
Note on compositional data
**************************

An important concept to keep in mind for microbiome analysis from NGS data is that **microbiome datasets are compositional**, which means that they have an arbitrary total imposed by the sequencing instrument. This is due to the fact that high-throughput sequencing (HTS) instruments can deliver reads only up to the capacity of the instrument.
Thus, the total read count observed in an NGS run is a fixed-size, random sample of the relative abundance of the molecules in the underlying ecosystem (e.g. oral microbiota). 
Several methods applied include count-based strategies (normalized to a constant area or volume, e.g. the sequencing lane output) such as Bray-Curtis dissimilarity, zero-inflated Gaussian models and negative binomial models, but these do not account for the limitations imposed by the instrument and the so-called principle of true indpendence of species in ecology. 
Read more about compositonal data in the paper by `Gloor et al. (2017)`_ and `Quinn et al. (2019)`_.

  .. _Gloor et al. (2017): https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full
  .. _Quinn et al. (2019): https://academic.oup.com/gigascience/article/8/9/giz107/5572529
  
Due to the compositional nature of microbiome datasets, the centered log-ratio (clr) transformation introduced by Aitchison (1986) is often used. The clr-transformed values are scale-invariant, which means that the same ratio is expected to be obtained in a sample with few read counts or an identical sample with many read counts, only the precision of the clr estimate is affected. 
The clr transformation uses the geometric mean g(x) of the sample vector as reference. 
Given an observation vector of D ???counted??? features (taxa, operational taxonomic units or OTUs, genes, etc.) in a sample, x = [x1, x2, ???xD], the geometric mean g(x) is: 

.. math::

  g(x) = \sqrt[D]{x_{i} * ... *  x_{D}}

The log-ratio transformation is applied to each subject vector *x* with *D* features (OTUs). Basically, the analysis of clr data will reveal how OTUs behave relative to the per-sample average.

.. math::

  clr(x) = [\ln{\frac {x_{i}}{g(x)}},...,\ln{\frac {x_{D}}{g(x)}}]

.. note::
  Compositional methods depend on logarithms that **do not compute for zeros**. Therefore, zeros in the abundance table must be addressed prior to, or during, the pipeline. 
  Two different methods of zero-handling can be used: zero-replacement using the pseudo-counts method from the R package `zCompositions`_ version 1.3.3 and the use of an arbitrary "offset", which means replacing the zeroes with ones. 
  It is good-praxis to always demonstrate that the removal or modification of zero-laden features does not change the overall interpretation of the results (`Quinn et al. (2019)`_).
  
    .. _zCompositions: https://www.sciencedirect.com/user/identity/landing?code=B0sdK2bMVyKF904thdCXKg1gPy2uFRWS5qoj4UPg&state=retryCounter%3D0%26csrfToken%3Deb9aa274-f4e8-4f8c-9adc-04d3770e7d74%26idpPolicy%3Durn%253Acom%253Aelsevier%253Aidp%253Apolicy%253Aproduct%253Ainst_assoc%26returnUrl%3D%252Fscience%252Farticle%252Fpii%252FS0169743915000490%26prompt%3Dnone%26cid%3Darp-25746d25-eda8-434f-bad4-f41bcdbb30ca

To account for the compositional nature of microbiome datasets, in the following analyses we will use Aitchinson-based (namely, clr-transformed) abundance values. 


********************************************************
Principal Coordinate Analysis (Multidimensional Scaling)
********************************************************
The Principal Coordinate Analysis (PCoA), also referred to as metric Multidimensional Scaling (MDS) is a multivariate reduction method performed on distance (or dissimilarity) indexes. 
In fact, PCoA can handle (dis)similarity matrices calculated from quantitative, semi-quantitative, qualitative, and mixed variables. 
Here we will use the Aitchinson distance, which corresponds to the Euclidean distance of the CLR-transformed dataset of relative abundances. To make of the clr-transformation of the abundance data we will use the library ``microbiome``. 
Read more here about the MDS and so-called `ordination methods`_

  .. _ordination methods: https://ourcodingclub.github.io/tutorials/ordination/

.. note::
  The CLR transformation of the package microbiome automatically deals with zeroes by appliying a pseudocount of min(relative abundance)/2 to exact zero relative abundance entries in the OTU table before taking logs.

::

  # calculate Aitchinson distance (the Aitchison distance is a sub-compositionally coherent distance defined as the Euclidean distance after the clr-transformation of the compositions). 
  library(microbiome)
  my_biom_rel_clr = transform(my_biom_rel, "clr")
  distAIT = distance(my_biom_rel_clr, method = "euclidean")


To make the PCoA we will use the ``ordinate`` function in Phyloseq. 
::

  # make PCoA ordination (MDS) of Aitchison distances
  ordAIT = ordinate(my_biom_rel_clr, method = "PCoA", distance = distAIT)
  # plot PCoA (MDS)
  plot_ordination(my_biom_rel_clr, ordAIT, color = "Group1")

We can also investigate how much of the total distance structure we will capture in the first few axes. We can do this graphically with a **scree plot**, an ordered barplot of the relative fraction of the total eigenvalues associated with each axis.
::
  
  plot_scree(ordAIT, "Scree plot for PCoA of Aitchinson distances")

Like we already did above, you can make the MDS on a subset of samples: 
::
  
  # subsample multiple groups phyla
  subgroup = c("Dental_calculus_chimp","Dental_calculus","Oral_plaque")
  subsample = subset_samples(my_biom_rel, Group1 %in% subgroup)
  # calculate Aitchinson distance 
  library(microbiome)
  subsample_clr = transform(subsample, "clr")
  distAIT = distance(subsample_clr, method = "euclidean")
  # make the PCoA ordination (MDS) setting the colors of the points automatically as defined in the metadata "Group1".
  ordAIT = ordinate(subsample, method = "PCoA", distance = distAIT)
  # two-dimension plot PCoA (MDS) setting the colors of the points automatically as defined in the metadata "Group2".
  plot_ordination(subsample, ordAIT, color = "Group1")
  

***********************************
Non-Metric Multidimensional Scaling
***********************************
Unlike other ordination techniques that rely on (primarily Euclidean) distances, such as Principal Coordinates Analysis, nMDS uses rank orders (so not the abundances), for this reason it is non-metric.
You can read more on the nMDS in the blogs curated by ecologists: `link1`_

To begin, NMDS requires a matrix of dissimilarities. Raw Euclidean distances are not ideal for this purpose: they???re sensitive to total abundances. 
Previously, NMDS plots based on Bray-Curtis distances were most commonly used with rarefied and relative abundance transformations, but Euclidean distances with clr-transformed data are becoming more common now. See for example `Sisk-Hackworth & Kelley (2020)`_.
NMDS arranges points to maximize rank-order correlation between real-world distance and ordination space distance. 

  .. _link1: https://archetypalecology.wordpress.com/2018/02/18/non-metric-multidimensional-scaling-nmds-what-how/
  .. _Sisk-Hackworth & Kelley (2020): https://academic.oup.com/nargab/article/2/4/lqaa079/5917299?login=true


The technique uses a *trial and error* to find the best positioning in dimensional space. Goodness-of-fit is measured by **stress** ??? a measure of rank-order disagreement between observed and fitted distances. 
A stress values below 0.2 ic considered acceptable. 
You can follow the same steps as above, just change the ordinatiom method supported in the command. 
::

  # make PCoA ordination (MDS) of Aitchison distances
  ordAIT = ordinate(my_biom_rel_clr, method = "NMDS", distance = distAIT)
  # plot the NMDS setting the colors of the points automatically as defined in the metadata "Group1".
  plot_ordination(my_biom_rel_clr, ordAIT, color = "Group1")


.. _DESeq2:

*********************************************
Differential taxonomic abundances with DESeq2
*********************************************
A common goal of microbiome studies is to identify differentially abundant taxa (species, genera etc.) between different groups of samples.
One of the tools used to do that is DESeq2, which was originally developed to identify differentially expressed genes in RNAseq data, but it is commonly adopted also in microbiome studies. You can read more about DESeq2 in the `original publication`_ and in the `dedicated page`_ of Bioconductor.
You can install `DESeq2 from Bioconductor`_. 

  .. _original publication: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8
  .. _dedicated page: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
  .. _DESeq2 from Bioconductor: https://bioconductor.org/packages/release/bioc/html/phyloseq.html

Run DESeq2 on the raw abundace data (`here is why`_), filtered for low-abundance taxa. You can choose to work on a subset of sample, as done above. For example, we can select only the oral samples.

  .. _here is why: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts

::

  library(DESeq2)
  # subsample oral microbiomes 
  subgroup = c("Dental_calculus_chimp","Dental_calculus","Oral_plaque")
  subsample = subset_samples(my_biom, Group1 %in% subgroup)

Then, we filter out low-abundance taxa with a mean < 0.02%. 
::

  minTotRelAbun = 0.0002
  x = taxa_sums(subsample)
  keepTaxa = taxa_names(subsample)[which((x / sum(x)) > minTotRelAbun)]
  subsample_flt = prune_taxa(keepTaxa, subsample)

.. warning::
  It is better to remove the low-abundant taxa only **after** subsampling the abudance table. Doing this before will keep unwanted taxa that are more common in other microbiota (e.g. soil) and irrelevant for describing oral envrionments (or the specific environment you are analysing with the PCA).  


DESeq2 does not tolerate zero-counts, for this reason an offset (+1) is commonly applied to remove all zeroes from your table.
After that, run DESeq2 on your table and contrat two sets of samples to find differential taxa abundances. 
::

  # make a +1 offset to run Deseq2 (which does not tolerate zero counts)
  otu_table(subsample_flt) = otu_table(subsample_flt)+1  
  # make deseq object from phyloseq object
  ds = phyloseq_to_deseq2(subsample_flt, ~ Group1)
  # Run DESeq2
  dds.data = DESeq(ds)
  # With the 'contrast' function you screen two different set of samples (based on your metadata) for differential taxa. 
  res = results(dds.data, contrast=c("Group1","Dental_calculus","Oral_plaque"))
  
Then you can explore your results, filter and sort the differential taxa detected based on a False Dicovery Rate (FDR) threshold (e.g. set to 0.01), the log2-fold-change and the base mean. Read more about the `FDR threshold here`_.

  .. _FDR threshold here: https://www.biostars.org/p/209118/
  
::

  # sort based on p-value adjusted:
  resOrdered = res[order(res$padj, na.last=NA), ]
  # set a threshold value for the False Discovery Rate (FDR):
  alpha = 0.01
  # get only significant taxa based on p-value adjusted (the FDR):
  resSig <- subset(resOrdered, padj < alpha)
  # sort significant values based on the log2-fold-change:
  resfc = resSig[order(resSig$log2FoldChange),]
  # sort significant values based on abundance (the base mean):
  resbm = resSig[order(resSig$baseMean),]
  # save the the tables 
  write.table(as.data.frame(resbm), file="deseq2.tsv", sep="\t", row.names=TRUE, col.names=TRUE, quote=F)
  
.. note::
  A positive log2-fold-change for a comparison of A vs B means that the feature (OTU) in A is more abundant than in B (and viceversa).
  For example, a log2-fold-change of ???1 means that in A the OTU is of 2^???1 = 0.5 less abundant than in B.


You can plot the log2-fold-changes to visualize for example the differential distribution of species based on their phyla (or genera) in the two groups contrasted. You could also color each differential taxon based on the phylum it belongs to.  
::

  # convert the data in a dataframe
  resSigMod = cbind(as(resSig, "data.frame"), as(tax_table(subsample_flt)[rownames(resSig), ], "matrix"))
  # Plot log-fold-changes of the OTUs based on Phylum
  ggplot(resSigMod, aes(x=phylum, y=log2FoldChange, color=phylum)) +
      geom_jitter(size=3, width = 0.2) +
      theme(axis.text.x = element_text(size = 10, angle = -90, hjust = 0, vjust=0.5)) +
      ggtitle("Dental calculus vs Oral plaque")
  # Plot log-fold-changes of the OTUs based on genera
  ggplot(resSigMod, aes(x=genus, y=log2FoldChange, color=phylum)) +
      geom_jitter(size=3, width = 0.2) +
      theme(axis.text.x = element_text(size = 5, angle = -90, hjust = 0, vjust=0.5)) +
      ggtitle("Dental calculus vs Oral plaque")

You can also generate a heatmap of the transformed abundance data. For visualization or clustering purposes it might be useful to work with transformed versions of the count data.
Various `transformations`_ may be applied on the data: 
1) ntd: Log transformation of the the data
2) vsd: Variance Stabilizing Transformation 
3) rld: Regularized log transformation 

  .. _transformations: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-transformations-and-visualization

When dealing with many samples it is useful to use the vst. To make the various transformatinos run the folowing commands: 
::

  # variance stabilizing transformation
  vsd = varianceStabilizingTransformation(dds.data, blind=FALSE)
  # log transformation
  ntd <- normTransform(dds.data)
  # regularized log transformation 
  rld <- rlog(dds.data, blind=FALSE)

We will plot the heatmap of the transformed data only for a reduced number of taxa (with decreasing counts).
::

  library("pheatmap")
  #adjust max number of taxa to display.
  select <- order(rowMeans(counts(dds.data,normalized=TRUE)),
                  decreasing=TRUE)[1:50]
  # create a dataframe importing using the metadata (e.g. Group1) and incorporating the taxa names as row names.
  df <- as.data.frame(colData(dds.data)[,"Group1"])
  rownames(df) <- colnames(subsample_flt@otu_table) 
  # plot heatmap (on ntd transformation)
  pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
           cluster_cols=FALSE, annotation_col=df, fontsize=4)


**********************************
Principal Component Analysis (PCA)
**********************************
The Principal Component Analysis (PCA) is another multivariate reduction method used for data visualization. 
It does not use a distance matrix and the benefit of that is that you can plot loadings onto your PCA axes. 
You can read more about the general principles of PCA here: `link3`_, `link4`_

  .. _link3: https://builtin.com/data-science/step-step-explanation-principal-component-analysis
  .. _link4: https://archetypalecology.wordpress.com/2018/02/17/principal-component-analysis-in-r/
  
  
We will make the PCA only on the oral microbiomes (dental calculus and plaque). So we will use the subsampled table that we already used for :ref:`DESeq2<DESeq2>`. 
Make sure you have the subsampled table, or do that again as follows: 
::

  # subsample oral microbiomes 
  subgroup = c("Dental_calculus_chimp","Dental_calculus","Oral_plaque")
  subsample = subset_samples(my_biom, Group1 %in% subgroup)
  #filter out low-abundance taxa with a mean < 0.02%. 
  minTotRelAbun = 0.0002
  x = taxa_sums(subsample)
  keepTaxa = taxa_names(subsample)[which((x / sum(x)) > minTotRelAbun)]
  subsample_flt = prune_taxa(keepTaxa, subsample)
  # Convert to relative frequencies. 
  subsample_rel = transform_sample_counts(subsample_flt, function(x) x / sum(x))

.. warning::
  It is better to remove the low-abundant taxa only **after** subsampling the abudance table. Doing this before may keep unwanted taxa that are more common in other microbiota (e.g. soil) and irrelevant for describing oral envrionments (or the specific environment you are analysing with the PCA).  


We will do the clr-transformation of the otu table in the biom object with the package `microbiome`, that we already :ref:`installed before<installed before>`
Activate the library and transform the data with the following command: 
::
  
  library(microbiome)
  subsample_rel_clr = transform(subsample_rel,"clr")

Finally, make the PCA with the command `ordinate`, adopting the ordination method `RDA` (redundancy analysis, a constrained ordination method), which without constraints correspond to the PCA in phyloseq.
::

  ord_clr = ordinate(subsample_rel_clr, "RDA")

You can plot the variance associated with each PC: 
::
  
  plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
  
And finally plot the PCA in a 2D chart. 
::

  plot_ordination(subsample_rel_clr, ord_clr, color="Group1") + 
    geom_point(size = 2)
 
To display the labels of each symbol: 
::

  plot_ordination(subsample_rel_clr, ord_clr, color="Group1", label="Sample_short") + 
    geom_point(size = 2)


.. note::
  To refine the analysis you can remove unwanted samples (e.g. those that appeared to be contaminated or outliers). Let's remove these four samples in the oral plaque group and repeat the PCA. 
  ::
  
    remove = c("SRR1564193","SRR1564223","SRR1564208","SRR1564211")
    subsample_pruned = prune_samples(!(subsample_rel_clr@sam_data$Sample_short %in% remove), subsample_rel_clr)
    ord_clr_pruned = ordinate(subsample_pruned, "RDA")
    plot_ordination(subsample_pruned, ord_clr_pruned, color="Group1") + 
    geom_point(size = 2)


