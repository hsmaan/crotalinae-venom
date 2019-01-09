Phylogenetic Analysis of Venom Proteome Composition and Toxicity for the Crotalinae Subfamily of Snakes
================
Author: Hassaan Maan
University of Guelph

----------------------------------------------------------------------

#### 1. INTRODUCTION: THE UTILITY OF SNAKE VENOM

Venomous snakes use an integrated combination of proteins, most commonly referred to as venom, to subdue their prey and ward off predators (Tasoulis & Isbister, 2017). Different snakes have venom of varying protein composition, which in turn leads to varying toxicological properties (Tasoulis & Isbister, 2017). However, many of the proteins found in snake venom have pharmacological and medicinal applications as well (Ojeda et al., 2018). Toxins extracted from snake venom have been used to develop various pharmaceutical compounds, such as captopril for the treatment of hypertension (Ojeda et al., 2018). Therefore, an understanding of the properties of snake venom could have great implications for drug development. A thorough knowledge of the proteomics, toxicology, and evolution of snake venom are all important considerations for the study of snake venom toxins. Venomous snakes that have implications for pharmacological development are typically in the following families: Atractaspidae, Elapidae, and Viperidae (Tasoulis & Isbister, 2017). Vipers, snakes from the Viperidae family, are further classified into two subfamilies based on the presence of a heat-sensing and temperature regulating pit gland (Bakken et al., 2018). Vipers that have this gland are classified into the Crotalinae family (Pit vipers), and those that do not are classified into the Viperinae family (True Vipers) (Tasoulis & Isbister, 2017).

Previous evolutionary studies of the snake venom proteome have allowed for many significant insights, such as findings related to which genes were recruited for the creation of the various venom protein toxins (Fry, 2005). Technological advancement for proteome analysis has led to a large influx of snake venom proteomes being published, and Tasoulis and Isbister conducted a review aggregating snake venom proteomes published over the last ten years (Tasoulis & Isbister, 2017).

In this exploratory analysis, the composition data from the aforementioned study for the Crotalinae subfamily will be used to investigate the correlation of genome evolution with venom composition and venom toxicity. Only Crotalinae species will be used, as traditional marker gene delineation may lead to unbalanced trees at higher taxonomic levels, and therefore a subfamily level will be considered to alleviate these effects. This analysis will be done in the following steps:

1.  Determination of the optimal marker gene for phylogenetic tree reconstruction using maximum likelihood estimation and bootstrapping.

2.  Phylogenetic tree comparison using various tree distance measures, between the marker-based evolution tree, and trees determined based on venom composition and toxicity level.

3.  Phylogenetic signalling analysis for venom composition and toxicity using Abouheif Moran's I spatial auto-correlation, as well as a phylogenetic trait mapping for venom toxicity level.

4.  Considering phylogenetic signalling effects, phylogenetic generalized least squares will be used to determine correlation between venom composition factors and toxicity level. If signalling effects are not significant, other correlation methods may be considered.

This study aims to provide insight into the correlation between Crotalinae evolution, venom composition, and toxicity. Insights into how venom composition changes with evolution could be useful, as further analysis could closely examine these changes to find novel genes associated with venom proteome composition. Correlation between venom composition and toxicity could provide toxicological insights into how protein composition and protein-protein interactions contribute to overall toxicity. Therefore, this analysis could serve as a starting point for investigations of a greater scope.

Let's load some necessary packages. A few key packages will be discussed further when explaining the choice of software tools used for the analysis. If you do not have any of these packages in your library, they can be installed using 'install.packages()'.

``` r
library(tidyverse)
library(rvest)
library(stringr)
library(picante)
library(ape)
library(adephylo)
library(ade4)
library(phylobase)
library(geiger)
library(phytools)
library(rentrez)
library(seqinr)
library(muscle)
library(Biostrings)
library(DECIPHER)
library(phangorn)
library(stringi)
library(ggplot2)
library(rlist)
library(distory)
library(geiger)
library(nlme)
library(corrr)
library(knitr)
library(data.table)
library(markdown)
```

------------------------------------------------------------------------

#### 2. DATASETS: CROTALINAE VENOM COMPOSITION AND VENOM TOXICITY

As outlined in the introduction, one of the main data sets involves protein composition data taken from a review study done by Tasoulis and Isbister (Tasoulis & Isbister, 2017). The proteome compositions of the venom are based on a certain number of abundant protein families, such as SVSP - snake venom serine protease. For the purposes of this analysis, only the abbreviations of these protein families will be used. The complete list of protein families considered with their full names is included in the Supplementary Information - Appendix A. This data set contains the protein families as an abundance percentage of the total venom proteome, based on the species of Pit vipers (Crotalinae). This data set was not able to web-scraped using the rvest package, as the XPATH for the HTML table was not clearly defined on PubMed, and therefore it was copied and pasted into an Excel csv file on December 3, 2018 from the following url: "<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5618223/table/toxins-09-00290-t003/?report=objectonly>".

The data set for LD50 toxicity of Crotalinae species was taken from "snakedatabase.org", specifically "<http://snakedatabase.org/pages/ld50.php>". The values for LD50, or lethal dose to kill 50% of a mice population, are taken from various studies on venomous snakes, each of which are cited within the data table. This analysis uses intravenous (IV) LD50 (mg/kg), as it is the most comprehensive type of toxicity data available in the database. The data set is sorted based on scientific name/common name of the snake species, followed by some physical data such as fang length, and various data for toxicity in mice, such as different LD50s based on method of venom injection. For the purposes of this analysis, only the species name and IV50 LD50 data will be used. The data set was web-scraped using the rvest package on December 3, 2018.

Genomic marker data for phylogenetic tree reconstruction will be taken from the NCBI Nucleotide database, and the following genes will be analyzed to determine the most appropriate for Crotalinae tree reconstruction: COI (Cytochrome c oxidase subunit I), 16S rRNA (ribosomal RNA), and 18s rRNA. Further explanation for this data will done in the 'additional' analysis section after the toxicity and venom composition data quality control.

------------------------------------------------------------------------

#### 3. DATA EXPLORATION AND QUALITY CONTROL

The venom composition data set had to be cleaned and imputed in various ways beforehand, and the steps taken are outlined in Supplementary Materials- Appendix B. We can import the cleaned data set into our environment, as long as the file is in our working directory. Please ensure this is the case before this step:

``` r
comp_data <- read.csv("comp_data.csv")

View(comp_data)
```

The first ten results from the first five columns of the composition data table are shown as follows:

``` r
kable(comp_data[1:10,1:5], caption = "Proteome Composition Dataframe")
```

| species\_name                   |   PLA2|   SVSP|   SVMP|   LAAO|
|:--------------------------------|------:|------:|------:|------:|
| Calloselasma rhodostoma         |   4.40|  14.90|  41.20|   7.00|
| Cryptelytrops purpureomaculatus |   8.00|  12.00|  35.00|  10.00|
| Gloydius brevicaudus            |  25.00|   3.70|  64.40|   0.90|
| Gloydius intermedius            |   9.90|  36.20|   2.60|  13.10|
| Ovophis okinavensis             |   0.65|  93.10|   4.20|   0.62|
| Protobothrops elegans           |  77.10|  10.40|   8.00|   0.50|
| Protobothrops flavoviridis      |  55.50|  11.80|  17.30|   3.10|
| Protobothrops mucrosquamatus    |  22.50|  10.40|  43.00|   2.00|
| Viridovipera stejnegeri         |  24.50|  11.00|  43.10|   3.30|
| Agkistrodon bilineatus          |  38.15|  12.25|  27.65|   3.75|

As we can see, the composition data is a table that has the Crotalinae species as row-names, and percentage abundance for protein families in the venom as columns. The last column, 'X.WV', summarizes what percentage of whole venom the listed components add up to. There are a significant number of zero values in this data, lets determine how many there are per given protein family.

``` r
comp_zero_count <- as.data.frame(sapply(comp_data, 
                   function(y) sum(length(which(y == 0))))) %>%
  setDT(keep.rownames = TRUE) %>%
  `colnames<-` (c("comp_column","zero_count"))
```

Let's plot this data using ggplot to get a visualization.

``` r
ggplot(data = comp_zero_count, aes(x = comp_column, y = zero_count)) + 
  geom_point(color = 'orange', size = 4) +
  geom_text(aes(label = comp_column), position = position_nudge(x = 0, y = -2.5), size = 3.5) +
  labs(title = "Zero Counts for Venom Protein Composition Data of Crotalinae Subfamily",
       y = "Number of Zero Values", x = "Composition Data Column") +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=9)) +
  theme(axis.text.x = element_text(size = 6))
```

<img src="Software_Final_Proj-Hassaan_Maan_files/figure-markdown_github/unnamed-chunk-5-1.png" width="100%" style="display: block; margin: auto;" />

As we can see from this plot, there are two columns in the composition data that have a very high number of zeros - the 'DEF' (defensins) and 'MPi' (snake venom metalloprotease inhibitor) columns. Since we are doing correlation analyses, and correlations can be severely impacted by high zero-value distributions, we will remove these two protein families from our data set. We can also go ahead and remove the 'X.WV' column, as we do not need a percentage of whole venom for our correlations. This can be done using the subset() function.

``` r
comp_data_sub <- subset(comp_data, select = -c(DEF,MPi,X.WV))
```

Now that we have our working data set for venom proteome composition, we need to get the venom toxicity data from "snakedatabase.org". We can use the R package rvest for this, as the website data has a clear XPATH for the table we need. Let's begin by reading the table into a variable.

``` r
tox_url <- "http://snakedatabase.org/pages/ld50.php"

ld50 <- tox_url %>%
  read_html() %>%
  html_nodes(xpath = '//*[@id="sortHeader"]') %>%
  html_table() 
  
ld50 <- ld50[[1]]

View(ld50)
```

As we can see, ld50 is a data frame containing the various parameters from the toxicity database. We are only concerned with the 'Scientific NameCommon Name' and 'IV(mg/kg)' columns. We will also have to clean up these columns a bit, as they both contain some unnecessary information and residual text. We will be using the word() function from the stringr package, which lets us manipulate strings in various types of data. Let's go ahead and do that.

``` r
iv_toxo <- str_remove(ld50$`IV(mg/kg)`,"ref")

species <- ld50$`Scientific NameCommon Name` %>%
  as.character() %>%
  word(1, sep = "\n\t\t\t\t\t\t\t\n\t\t\t\t\t\t\t\n\t\t\t\t\t\t\t\n\t\t\t\t\t\t\t") %>%
  word(1) %>%
  as.character()
```

Now we can put this information back into a data frame and correct for one last issue - non breaking spaces present in the species name. We can use another function from stringr - str\_replace, which looks for a string based on a pattern and replaces it with the second given string.

``` r
iv_50 <- data.frame(species, iv_toxo) %>%
  `colnames<-` (c("species_name", "iv50")) 

iv_50$species_name <- as.character(iv_50$species_name) %>%
  str_replace("\u00A0", " ") # '\u00A0' was a non-breaking space..
```

We also need to convert the toxicity data from factor to numeric, as leaving it as factor will complicate our downstream analyses because of inability to manipulate it as numerical data.

``` r
iv_50$iv50 <- as.numeric(as.character(iv_50$iv50))

View(iv_50)
```

Now we can get of the 'NA' rows. We can do this using the na.omit() function from the data.table library.

``` r
iv50 <- na.omit(iv_50)
```

We've eliminated all of the 'NA' species values for toxicity level. However, let's plot this data to determine if there are any outliers. A special function being used here, and one that will be used later, is the ifelse() function, which is essentially a one-line if/else statement. It takes a test, and then returns values based on the provided if/else settings.

``` r
ggplot(data = iv50, aes(x = species_name, y = iv50)) + 
  geom_point(color = 'red', size = 2) +
  geom_text(aes(label = ifelse(iv50>10, as.character(species_name),'')), 
            position = position_nudge(x = 0, y = 1)) +
  labs(title = "Intravenous LD50 Toxicity Levels for Venomous Snakes 
       from the Crotalinae Subfamily",
       y = "Intravenous LD50 (mg/kg)", x = "Crotalinae Species") +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=9), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank())
```

<img src="Software_Final_Proj-Hassaan_Maan_files/figure-markdown_github/unnamed-chunk-12-1.png" width="90%" style="display: block; margin: auto;" />

From this plot, we can see we have quite a few species that don't follow the trend of having a &lt;5 mg/kg IV50 toxicity level. Of these, there are some severe outliers that have a &gt;10 IV50 toxicity level, the most severe being *Natrix tessellata*, which has an IV50 of about 25 mg/kg. Going back to the database itself, there is a primary reference for this toxicity level as well as secondary. Examining the primary reference, this value is in fact supported in the literature (Mackessy, 2002). Therefore, this severe outlier will be kept in the data set, and if matched with the composition and genomic data, might even provide some insight as to why the toxicity level for this snake species is particularly low.

Great, we have our cleaned and processed toxicity data set. After selection of the marker gene, some more quality control will be done for the marker gene data itself. The additional section on marker gene selection will follow the software tool description section.

------------------------------------------------------------------------

#### 4. MAIN SOFTWARE TOOLS: METHODOLOGY

There are many different R packages used in this analysis (R Core Team, 2018). Specialty packages, such as rvest, are used for unique tasks throughout the script. However, for the main analysis, there are 3 packages that are used extensively throughout and/or have potential replacements: phangorn, phytools, and muscle. A short summary of the rationale behind choosing each package is outlined:

-   **muscle** - This package is used for sequence alignment when determining the best marker for phylogenetic reconstruction of the Crotalinae subfamily. This package contains a very simple, lightweight, and fast alignment tool in muscle(), and these were the primary reasons it was chosen over an alternative such as DECIPHER's alignment. Speed is an important consideration, as we will see later that a function call to alignment is made multiple times, and therefore it is imperative to have an alignment algorithm that is fast and accurate. (Edgar, 2004)

-   **phangorn** - Phangorn is a comprehensive package for phylogenetic reconstruction of trees using optimizations such as maximum likelihood. It also contains methods for manipulation of trees and distance analysis, the latter point being the main reason this package was chosen. Since an analysis of tree distance between the marker tree, venom composition tree, and toxicity level tree will be done, it is best to use a comprehensive library that will not only be able to create the optimized trees, but also do the distance measurements necessary for the analysis. An alternative could have been to use the BEAST (Bayesian Evolutionary Analysis Sampling Trees) software and import the trees, but since this analysis is being done in R, this was not considered. (Schliep, 2011)

-   **phytools** - The phytools package is used for quantitative trait mapping. There are a many alternatives to trait mapping in R, and the primary one that was considered was ggtree, because of the high customizability. However, there is currently no stable release of ggtree for R version 3.5.1, which is the version being used for this analysis. An attempt was made install the package using the Bioconductor software repository, but this destabilized the dependencies for many other packages, so a decision was made to use phytools. (Revell, 2012)

It is worth mentioning that adephylo also has an alternative for determining phylogenetic signalling in the ape package, however since the package imports from ape this is a trivial issue.

------------------------------------------------------------------------

#### 5. MARKER SELECTION AND PHYLOGENETIC RECONSTRUCTION

We want to determine the best marker to use for our phylogenetic analysis. In this analysis, we will consider the cytochrome c oxidase subunit I mitochondrial gene (COI), the 16S mitochondrial ribosomal RNA (16S), and the 12S mitochondrial ribosomal RNA (12S). Availability of these marker sequences for the Crotalinae subfamily was assessed beforehand, and appropriate search terms for the NCBI Nucleotide Database were developed. Unfortunately, other potential marker sequences are not readily available for Crotalinae species.

To make this analysis quick and simple, we will be using a function that conducts the search for the given term using the rentrez package function entrez\_search(), fetches the fasta files, matches the species with those in the venom composition and venom toxicity data sets, performs multiple sequence alignment using the muscle package, and creates a bootstrapped maximum likelihood tree using the phangorn package. We will be using visual inspection of the trees to assess the best marker to use for this analysis, based on the phylogeny itself and the confidence values from bootstrapping. However, once we have chosen our marker, we will do some more quality control before proceeding.

The search terms for the three marker regions are the following:

``` r
marker_coi <- "Crotalinae[ORGN] AND 600:700[Sequence Length] AND COI 
               NOT WHOLE NOT COMPLETE NOT UNVERIFIED"

marker_16s <- "Crotalinae[ORGN] AND 400:500[Sequence Length] AND 16S 
               NOT WHOLE NOT COMPLETE NOT UNVERIFIED"

marker_12s <- "Crotalinae[ORGN] AND 400:500[Sequence Length] AND 12S 
               NOT WHOLE NOT COMPLETE NOT UNVERIFIED"
```

The following function, 'tree\_generator', will perform the computation of the maximum likelihood bootstrapped trees for our three markers. The function itself is thoroughly commented. It will return a list of important data for use later on, such as the composite dataset, and it will output the bootstrapped and optimized tree.

``` r
tree_generator <- function(marker_name, search_term) {
  
  # The variable 'datalist' will be used to store informative variables created
  # throughout the function, such as the distance matrix, and will return a list
  # containing those variables after the execution of the function.
  
  data_list <- list()
  
  # We can use the rentrez package, and the function entrez_search to retrieve
  # the data given the search term, and use a web_history token (in case the
  # search returns a very large number of hits) to fetch the data.
  
  search <- entrez_search(db = "nuccore", term = search_term, 
                          retmax = 10000, use_history = TRUE)
  fetch <- entrez_fetch(db = "nuccore", web_history = search$web_history,
                        rettype = "fasta")
  
  # Now that we have the data, we can write the results into a fasta file and
  # create a dataframe.
  
  write(fetch, "marker.fasta", sep = "\n") 
  dna <- readDNAStringSet("marker.fasta") 
  marker <- data.frame(data = names(dna), seq = paste(dna))
  
  # Using the word function, we can clean our data, and get the columns we need
  # - species, id, and sequence.
  
  marker$species_name <- word(marker$data, 2L, 3L)
  marker$uniq_id <- word(marker$data, 1L)
  marker <- marker[, c("species_name", "uniq_id", "seq")]
  
  # We want only unique species, and not multiple sequences for each. Therefore,
  # we can shuffle and pick a unique sequence for each species.
  
  marker_shuff <- marker[sample(nrow(marker)),]
  marker_uniq <- distinct(marker_shuff, species_name, .keep_all = TRUE)
  
  # Let's join our data with the venom composition and toxicity dataframes,
  # since we only want to consider species that are in the union of the three
  # datasets. We can use inner_join() from the tidyverse package, which will
  # keep only species that are in both the joined datasets.
  
  marker_match <- marker_uniq %>%
    inner_join(comp_data_sub, by = "species_name") %>%
    inner_join(iv_50, by = "species_name")
  
  # Set the dataframe with matching composition, toxicity and genomic data as a
  # variable in data_list. This will be done throughout the function for
  # different variables.
  
  data_list[[1]] = marker_match
  
  
  # We will be using the muscle library for alignment with no optimizations for
  # the best accuracy. Rationale for using muscle alignment over other alignment
  # algorithms is outlined in Section 4: Main Software Tools.

  marker_align <- DNAStringSet(muscle::muscle(DNAStringSet(marker_match$seq)))
  names(marker_align) <- marker_match$species_name
  marker_dist <- dist.dna((as.DNAbin(marker_align)), model = "K80", as.matrix = TRUE, 
                             pairwise.deletion = TRUE)
  marker_phy <- as.phyDat(as.DNAbin(marker_align))
  
  # Adding some more important variables to the data_list. These variables will
  # be used for downstream analyses, an example being using the phyDat for
  # phylogenetic signalling analysis.
  
  data_list[[2]] = marker_align
  data_list[[3]] = marker_dist
  data_list[[4]] = marker_phy
  
  # Using the phangorn library, we can create a rudimentary neighbor joining
  # tree, and optimize it using the best model from a model test, and then use
  # maximum likelihood estimation followed by bootstrapping.
  
  tree=NJ(marker_dist)
  
  
  # This model test will give us the best fitting model for our data, using the
  # AIcc goodness of fit measure.
  
  mt <- modelTest(object = marker_phy, tree = tree)
  mt[order(mt$AICc), ]
  bestmodel <- mt$Model[which.min(mt$AICc)] 
  
  # Now we can use the optim.pml() function to optimize our tree using maximum
  # likelihood methods.
  
  env = attr(mt, "env")
  fitStart = eval(get(bestmodel, env), env)
  fit = optim.pml(fitStart, rearrangement = "stochastic",
                  optGamma = TRUE, optInv = TRUE, model = "GTR")
  
  # Bootstrapping will allow us to have confidence values for our trees. We will
  # perform 100 bootstrap samplings.
  
  bs = bootstrap.pml(fit, bs = 100, optNni = TRUE)
  
  # Now we can plot the optimized and bootstrapped tree with the appropriate
  # titles. This plot will be an output of the function.
  
  plotBS(fit$tree, bs, p = 50, type = "phylogram", cex = 0.70, 
         main = paste(marker_name,"ML Optimized Tree for Crotaline Subfamily 
         with Bootstrapped Confidence Values"))
  
  data_list[[5]] = fit$tree
  
  # The assign() function with 'envir = .GlobalEnv' will return our data_list
  # for the current marker to our working environment. This allows us to use
  # this data when the function is complete.
  
  assign(paste("marker", marker_name, "datalist", sep = "_"), data_list, envir = .GlobalEnv)
  
}
```

Now using our three search parameters, we can create bootstrapped maximum likelihood trees using the three marker genes:

``` r
tree_generator("COI", marker_coi)
```

<img src="Software_Final_Proj-Hassaan_Maan_files/figure-markdown_github/unnamed-chunk-15-1.png" width="100%" style="display: block; margin: auto;" />

``` r
tree_generator("16S_rRNA", marker_16s)
```

<img src="Software_Final_Proj-Hassaan_Maan_files/figure-markdown_github/unnamed-chunk-15-2.png" width="100%" style="display: block; margin: auto;" />

``` r
tree_generator("12S_rRNA", marker_12s)
```

<img src="Software_Final_Proj-Hassaan_Maan_files/figure-markdown_github/unnamed-chunk-15-3.png" width="100%" style="display: block; margin: auto;" />

From the bootstrapped COI based phylogenetic tree, we can see that there are not that many species overall as compared to 16S and 12S. Since we are doing an analysis on a subfamily, we want as much coverage as possible to make the phylogenetic tree as stable as possible, and therefore we will not consider using the COI marker for this analysis. The number of species available from the 12S and 16S markers are similar, although 12S has more species, but they are fairly different in terms of the species they actually have. Therefore, a combination of the two markers is out of question, as we will lose a lot of data in that case, and the statistical power of our downstream tests will decrease. Because the bootstrapping confidence values between 12S and 16S are fairly similar, and because 12S contains more species (24 vs. 26), we will proceed with using 12S as the marker for our phylogenetic analysis.

However, from the actual 12S tree, we can see that there are some peculiar branchings with low confidence values. Overall, the genuses are grouped into their own clades very well, but there is one exception - *Atropoides*. *A. nummifer* and *A. picadoi* are in separate clades, which is not in line with the rest of the tree. There are also a number of species from a unique genus in the tree, however, since there is no conflict with the non-unique genus species, these will be left as is. For the tree we will be using, we will drop the two *Atropoides* genus members. Instead of pruning, we can do a complete reconstruction by adding an exception to the search term.

``` r
marker_12s_sub <- "Crotalinae[ORGN] AND 400:500[Sequence Length] AND 12S 
NOT WHOLE NOT COMPLETE NOT UNVERIFIED NOT ATROPOIDES" 

tree_generator("12S_sub", marker_12s_sub)
```

<img src="Software_Final_Proj-Hassaan_Maan_files/figure-markdown_github/unnamed-chunk-16-1.png" width="100%" style="display: block; margin: auto;" />

We have our optimized and bootstrapped tree. We have the data necessary for proceeding with our main analysis in the 'marker\_12S\_sub\_datalist' list object.

------------------------------------------------------------------------

### 6. MAIN ANALYSES

To begin, we want to create trees for the toxicity and venom proteome data sets. This can be done by creating distance matrices for these data sets and then creating rudimentary neighbor joining trees using the phangorn package. First let's extract our full working data set from the 'marker\_12s\_sub\_datalist' list of objects.

``` r
cro_comp_tox <- marker_12S_sub_datalist[[1]]

View(cro_comp_tox)
```

As we can see, our data frame contains all of our necessary data, including sequence data for 12S rRNA, venom proteome protein family composition, and intravenous LD50 toxicity, all sorted by the Crotalinae species in question. We can go ahead and extract our tree from the list of lists as well.

``` r
marker_tree <- marker_12S_sub_datalist[[5]]
```

Now we want to create distance matrices for the IV50 and proteome composition data. For both toxicity and proteome composition, this can be done using the dist() function from the stats library in R. We will use euclidean distance for both trees, as we have no data-related reason to use any other specialized distance, such as Manhattan. Let's go ahead and do this for both data sets. We can pipe the output to create a matrix from these distance objects with the appropriate species row names and column names.

``` r
tox_dist <- as.matrix(dist(cro_comp_tox$iv50, method = "euclidean")) %>%
  `colnames<-` (cro_comp_tox$species_name) %>%
  `rownames<-` (cro_comp_tox$species_name)

comp_dist <- as.matrix(dist(cro_comp_tox[4:11], method = "euclidean")) %>%
  `colnames<-` (cro_comp_tox$species_name) %>%
  `rownames<-` (cro_comp_tox$species_name)
```

Now we can go ahead and create the neighbor joining trees using phangorn for both distance matrices.

``` r
tox_tree = NJ(tox_dist)

comp_tree = NJ(comp_dist)
```

We'll be testing for similarities between the three trees in a quantitative fashion. The phangorn package will be used again, as it contains a comprehensive list of tree distance determination algorithms. These algorithms are encompassed in the treedist() function and the ones that will be used include: the Robinson-Foulds distance (difference in topology of trees), the weighted Robinson-Foulds distance (factors in edge length), and the path distance (only considers pathing) (Kuhner & Yamato, 2015). The interpretation of these different distance computation methods will be expanded upon in the discussion section.

To begin, we can put our trees in a multiPhylo object so that we can perform the distance measurements for all three using one function call.

``` r
tree_list <- as.multiPhylo(c(comp_tree, tox_tree, marker_tree))
```

To compute the distances, we will use the specific distance parameter (in the first case, Robinson-Foulds as RF.dist), and the multiPhylo tree list. We can pipe this as a matrix to label our columns and rows, to have the appropriate names corresponding to their respective trees.

``` r
rf_dist <- RF.dist(tree_list, check.labels = TRUE, rooted = FALSE) %>%
  as.matrix %>% 
  `rownames<-` (c("comp_tree", "tox_tree", "marker_tree")) %>%
  `colnames<-` (c("comp_tree", "tox_tree", "marker_tree"))
```

And similarly for weighted Robinson-Foulds (wRF.dist) and path distance (path.dist).

``` r
wrf_dist <- wRF.dist(tree_list, check.labels = TRUE, rooted = FALSE) %>%
  as.matrix %>% 
  `rownames<-` (c("comp_tree", "tox_tree", "marker_tree")) %>%
  `colnames<-` (c("comp_tree", "tox_tree", "marker_tree"))

path_dist <- path.dist(tree_list, check.labels = TRUE) %>%
  as.matrix %>% 
  `rownames<-` (c("comp_tree", "tox_tree", "marker_tree")) %>%
  `colnames<-` (c("comp_tree", "tox_tree", "marker_tree"))
```

Now let's examine the results for each distance metric.

``` r
View(rf_dist)
```

``` r
kable(rf_dist, caption = "Robinson-Foulds Distance")
```

|              |  comp\_tree|  tox\_tree|  marker\_tree|
|--------------|-----------:|----------:|-------------:|
| comp\_tree   |           0|         42|            40|
| tox\_tree    |          42|          0|            42|
| marker\_tree |          40|         42|             0|

As we can see, for the Robinson-Foulds distance, there is an approximately equal distance between each of the three trees, in the range of 40-42 arbitrary units. Therefore, using this distance metric, the topology of the three trees is approximately equally distant.

``` r
View(wrf_dist)
```

``` r
kable(wrf_dist, caption = "Weighted Robinson-Foulds Distance")
```

|              |  comp\_tree|  tox\_tree|  marker\_tree|
|--------------|-----------:|----------:|-------------:|
| comp\_tree   |      0.0000|   327.3624|      324.2220|
| tox\_tree    |    327.3624|     0.0000|       12.5611|
| marker\_tree |    324.2220|    12.5611|        0.0000|

For the weighted Robinson-Foulds distance, we can see that there are large distances between the composition tree and toxicity tree, and the composition tree and marker tree. However, the distance is VERY small for the marker tree and toxicity tree as compared to the other two combinations. Therefore, when considering edge lengths, the Robinson-Foulds distance between the toxicity tree and marker based phylogenetic tree indicates that they are quite similar. This result will be further analyzed in the discussion.

``` r
View(path_dist)
```

``` r
kable(path_dist, caption = "Path Based Distance")
```

|              |  comp\_tree|  tox\_tree|  marker\_tree|
|--------------|-----------:|----------:|-------------:|
| comp\_tree   |     0.00000|   110.4129|      77.87169|
| tox\_tree    |   110.41286|     0.0000|     102.11268|
| marker\_tree |    77.87169|   102.1127|       0.00000|

For path distance, we see similar results as when we computed the Robinson-Foulds distance. There are large distances between all three trees when considering the pathing. An interesting note is that the path distance for the marker tree and toxicity tree is greater than the two others.

Now before we move on to considering whether or not to use phylogenetic generalized least squares to compute correlations between venom proteome components and toxicity level, we will do an analysis of the phylogenetic signalling itself. We will be creating a phylogenetic trait mapping of toxicity level over-top our marker tree phylogeny, and then computing the Abouheif Moran's I coefficient for phylogenetic signalling using the adephylo package. At 0, the Moran's I coefficient indicates that species that are highly related are as similar for a given trait as predicted by random chance (usually correlated with Brownian motion) (Pavoine et al., 2008). If the coefficient is less than 0, it indicates that species that are highly related have greater dissimilarity for the trait as compared to chance (Pavoine et al., 2008). If greater than 0, then the similarity is greater than as predicted by random chance (Pavoine et al., 2008).

Let's begin by creating a phylogenetic trait mapping for toxicity to get a visualization of toxicity level signalling. We will use the phytools package for this task.

First we need to create a variable that has our toxicity values as a matrix and set the species names for the matrix using the names() function.

``` r
tox_map <- as.matrix(cro_comp_tox$iv50)
names(tox_map) <- cro_comp_tox$species_name
```

Now we can go ahead and use the contMap() function from phytools to create a trait mapping plot for toxicity level.

``` r
par(mar = c(1,1,1,1), oma = c(1,1,3,0))
wings_map <- contMap(midpoint(marker_tree), tox_map)
par(mar = c(1,1,2,1), oma = c(1,1,1,0))
title(main = "Continuous Trait Mapping for Toxicity Level of 
      Crotalinae Subfamily Species")
```

<img src="Software_Final_Proj-Hassaan_Maan_files/figure-markdown_github/unnamed-chunk-31-1.png" width="100%" style="display: block; margin: auto;" />

As we can see, there is some degree of evolutionary conservation between similar species for toxicity level. For example, the *Bothrops* genus clade members seem to all share relatively high levels of toxicity (low IV50). However, there is also evidence to the contrary, such as the *Crotalus* genus clade, where although for most members the toxicity is relatively high, there is an outlier in *Crotalus basiliscus* that has a very low toxicity level (high IV50). Perhaps this species is highly diverged in terms of toxicity due to geographical isolation or some other evolutionary phenomenon.

Now we can go ahead and perform quantitative phylogenetic signalling analysis for all the proteome composition traits and toxicity level. We will be using the adephylo package for this task. First let's go ahead and extract the necessary columns names from our main data set.

``` r
comp_cols <- colnames(cro_comp_tox[4:12])
```

We want to use a for loop to compute the Moran's I coefficient for each trait. Before we do that, let's create a data frame to store the values from the for loop. We can do that by first creating a matrix of length equal to the traits we have (9), and then transposing this into a data frame with the appropriate column names.

``` r
abo_moran <- data.frame(t(matrix(vector(), 1, 9,
                               dimnames=list(c(), c(comp_cols))))) %>%
  `colnames<-` ("Moran's I Coefficient")

View(abo_moran)
```

We're all set to loop through the traits and get the Moran's I coefficient for each trait. The first three lines of the loop are simply to perform a manipulation similar to the one we did for the toxicity level trait mapping. We need to extract the values for each trait, and append it to the numerical vector or matrix with the names of the object being the species names. Fortunately, since we have the column names extracted, we can loop over this to produce the necessary numerical vectors. Then we can use a phylo4d object from the phylobase package and the phylo4d() function, to combine the trait values and marker tree data. This allows the abouheif.moran() function from the adephylo package to perform the multivariate analysis necessary to test for phylogenetic signalling. There is also a counter being used to append the correct values to the appropriate place in the abo\_moran data frame.

``` r
counter <- 1  
for (x in comp_cols) {
  
  # Get numerical vectors for trait data 
  a=select(cro_comp_tox, x)
  a <- as.numeric(unlist(a))
  names(a) <- cro_comp_tox$species_name
  
  # Create phylo4d object and perform multivariate analysis to determine Moran's
  # I coefficient
  four_d <- phylo4d(marker_tree, a)
  moran <- abouheif.moran(four_d, method = "Abouheif")
  
  # Append coefficient to abo_moran dataframe
  abo_moran[counter,] <- moran$obs
  counter <- counter + 1
}

View(abo_moran)
```

``` r
kable(abo_moran, caption = "Moran's I Coefficient for Crotalinae Traits")
```

|             |  Moran's I Coefficient|
|-------------|----------------------:|
| PLA2        |              0.0912147|
| SVSP        |             -0.0174976|
| SVMP        |             -0.1334267|
| LAAO        |              0.2590137|
| CRiSP       |              0.0206043|
| CTL.SNACLEC |             -0.0246197|
| DIS         |              0.3677768|
| NP          |              0.2117609|
| iv50        |              0.0185482|

We can go ahead and plot this data to visualize the quantitative phylogenetic signalling for each trait. Let's go ahead and convert the data frame to a data table to get the two columns we need for ggplot.

``` r
abo_moran_dt <- abo_moran %>%
  setDT(keep.rownames = TRUE) %>%
  `colnames<-` (c("protein_family", "moran_coeff")) 
```

Now we can go ahead and create the plot using ggplot.

``` r
ggplot(data = abo_moran_dt, aes(x = protein_family, y = moran_coeff)) + 
  geom_segment(aes(x = protein_family, xend = protein_family, y= 0, 
                   yend = moran_coeff)) +
  geom_point( color="green", size=5) +
  labs(title = "Moran's I Phylogenetic Signalling Coefficient for 
       Venom Traits of Crotalinae Subfamily",
       y = "Abouheif Moran's I Coefficient", x = "Crotalinae Trait") +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title = element_text(size=9))
```

<img src="Software_Final_Proj-Hassaan_Maan_files/figure-markdown_github/unnamed-chunk-37-1.png" width="100%" style="display: block; margin: auto;" />

From this plot, we can clearly see that no given trait (protein family or toxicity level) has a strong phylogenetic signalling in either direction. There are certain protein family traits that have a weak positive signal, such as DIS, however this may well not be significant, as the statistical power for this phylogenetic analysis is lacking, which is a topic that will be elaborated in the discussion. Now we can go ahead and proceed with the correlation of venom proteome components and toxicity level.

A preliminary analysis of the combined data set and the optimized tree revealed certain issues when trying to perform phylogenetic generalized least squares regression. The model had a tendency to over-fit the data, even after a lot of fine tuning. This is likely because of the low statistical power of the model itself (24 observations). Therefore, we must continue with a regular correlation test, and this is adequate because we learned that there were no strong phylogenetic signalling effects for any of the traits in question from the analysis of Abouheif Moran's I coefficient. Thus, a regular correlation that does not consider phylogenetic signalling will be apt.

We will be using Pearson's correlation, as our data is interval in nature and not ordinal, and therefore it is best to avoid the rank-order correlations. We can use the corrr package to correlate the series of venom composition variables to the singular toxicity variable. Let's start by sub-setting our data frame to only our required variables - venom composition protein families and toxicity.

``` r
comp_tox <- cro_comp_tox[4:12]

View(comp_tox)
```

Now we can create a correlation table for the entire data frame. We will use the correlate() function from corrr with the default method (Pearson) and then pipe this output to the focus() function, to get our variable of interest (IV50). This will return to us the correlations for the proteome data with the toxicity. Let's go ahead and do that.

``` r
tox_corr <- comp_tox %>%
  correlate(quiet = TRUE) %>% 
  focus(iv50) %>%
  `colnames<-` (c("prot_family", "iv50_corr")) 

tox_corr
```

``` r
kable(tox_corr, caption = "Venom Component Correlation with IV50 Toxicity (Pearson's)")
```

| prot\_family |  iv50\_corr|
|:-------------|-----------:|
| PLA2         |  -0.1409751|
| SVSP         |  -0.0536340|
| SVMP         |   0.2808405|
| LAAO         |  -0.3196048|
| CRiSP        |  -0.1260457|
| CTL.SNACLEC  |  -0.0910557|
| DIS          |   0.3034187|
| NP           |   0.1889330|

We have our correlations. To get a visualization of this data, we can use ggplot once again. Before we do this, would also like to create a column to indicate if the correlation is positive or negative. We can do this using the ifelse() function mentioned before, with the combination of the with() function.

``` r
tox_corr$corr_sign <- with(tox_corr, ifelse(iv50_corr>0, "positive", "negative"))
```

Now let's go ahead and plot the correlations.

``` r
ggplot(data = tox_corr, aes(x = prot_family, y = iv50_corr)) + 
  geom_bar(stat = "identity", aes(fill = corr_sign)) +
  labs(title = "Pearson's Coefficient for Correlation of Crotalinae Venom Protein
       Families and Toxicity Level",
       y = "Pearson's Correlation Coefficient", x = "Venom Proteome Component") +
  scale_fill_manual(name="Coefficient Sign", 
                    labels = c("Positive", "Negative"), 
                    values = c("positive"="#00ba38", "negative"="#f8766d")) +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=9)) +
  coord_flip(ylim = c(-1, 1))
```

<img src="Software_Final_Proj-Hassaan_Maan_files/figure-markdown_github/unnamed-chunk-42-1.png" width="100%" style="display: block; margin: auto;" />

The x-axis limits of this plot are set to the boundary of possible values of Pearson's correlation \[-1,1\], and thus we can see that none of the protein families have a particularly strong correlation one way or another to the toxicity level. The three strongest correlations seem to be LAAO (negative), DIS (positive), and SVMP (positive). These results will be elaborated further in the discussion.

------------------------------------------------------------------------

### 7. INTERPRETATION AND DISCUSSION

The overarching theme of this study was to provide exploratory insights into relationships between Crotalinae subfamily venom proteome composition, venom toxicity, and genetic evolution. Through an analysis of three markers (12S rRNA, 16S rRNA, and COI), it was determined that taxonomic delineation using 12S rRNA was most appropriate for this analysis due to high confidence values for phylogenetic reconstruction and sample size. Three metrics were considered when analyzing the relationship between the venom proteome composition, venom toxicity level, and marker-based evolution: a quantitative comparison of trees based on distances from the three data sets, phylogenetic signalling analysis for venom proteome protein families and toxicity level using the 12S rRNA based phylogenetic reconstruction, and a correlation analysis between the venom proteome composition components and the toxicity level.

In terms of distance between trees, both path distance and strictly topology-based Robinson-Foulds distance did not show any similarities between the three trees, as both metrics returned high distance values. However, weighted Robinson-Foulds distance, which considers edge length as a factor, showed a very small distance between the toxicity level tree and the evolutionary marker-based tree. This indicates that although the toxicity level tree and phylogenetic tree are quite distant when considering topology, when factoring in their branch lengths they are the closest tree combination out of the three possible combinations. However, these two trees have a high degree of dissimilarity based on visual inspection. Previous studies have shown that when considering trees with high dissimilarity, topology-based metrics are preferred as a primary means of comparison (Kuhner & Yamato, 2015). Therefore, although the weighted Robinson-Foulds metric does indicate similarity, the topology-based metric has greater validity. For both the phylogenetic signalling and correlation aspects of the analysis, there were no strong relationships observed. No trait showed an absolute Moran's I coefficient greater than 0.4, and for correlation, the highest Pearson's coefficient was determined to be |-0.32| from the LAAO protein family. Therefore, overall there were no strong relationships found for either phylogenetic signalling of venom composition components and toxicity level, or correlation between venom composition components and toxicity level. In terms of phylogenetic signalling, this could simply be indicative of low conservation of traits between closely related species, however this interpretation should be taken with a grain of salt, as will be explained in a bit. Correlation between the venom proteome and toxicity level was done without considering phylogeny, but since the correlations were quite low regardless, its unlikely that this would have changed the results, especially since there was no strong phylogenetic signalling for the venom proteome traits.

Possibly the most important limitation of this analysis, and the reason as to why there were no strong relationships observed in this data set, was small sample size. The analysis started with the potential of 65 species from the main data set (venom composition), but when joining with the toxicity level data and 12S rRNA marker sequence data, it was reduced to 24 species. This low sample size complicated the downstream analyses quite a bit, as it likely affected all the steps. Phylogenetic reconstruction will be affected by low sample size, especially since we're considering a subfamily which is a relatively high-level taxonomic group for marker-based delineation. For phylogenetic signalling, low sample size likely confounds the results, and combined with the fact that the phylogenetic reconstruction is limited, the results are very difficult to interpret. Due to the small sample size, we simply do not have high confidence in the results. For correlation analysis, the small sample size led to over-fitting when using phylogenetic least squares estimation. Therefore, Pearson's correlation was considered, but it should still be robust as the phylogenetic signalling wasn't strong. The sample size was likely adequate for Pearson's correlation, but the data itself is likely not distributed in a Gaussian fashion, which is an underlying assumption of most correlation techniques (Bewick et al., 2003). This likely confounds the analysis for correlation as well. Overall, in terms of future steps for this type of analysis, sample size must be considered. Perhaps only two of the three data sets can be joined for analysis, which will increase sample size slightly, and provide more validity to the results. The lack of standardized marker data for Crotalinae is also a factor, as many of the sequences for marker regions were of differing length and/or region. If given more time, limiting factors for the data would be analyzed further, such as which species of Crotalinae are missing data from any of the three datasets, and attemps would be made to either retrieve this data from other sources or consider imputation techniques. An important finding was made through continuous trait mapping of toxicity level, where *Crotalus basiliscus* has a very low toxicity level compared to its closely related genus members. Further analysis of this clade could provide some interesting insights into venom toxicity evolution.

------------------------------------------------------------------------

### CITATIONS

Bakken, G. S., Schraft, H. A., Cattell, R. W., Tiu, D. B., & Clark, R. W. (2018). Cooler snakes respond more strongly to infrared stimuli, but we have no idea why. The Journal of Experimental Biology, 221(Pt 17). <https://doi.org/10.1242/jeb.182121>

Bewick, V., Cheek, L., & Ball, J. (2003). Statistics review 7: Correlation and regression. Critical Care, 7(6), 451-459. <https://doi.org/10.1186/cc2401>

Edgar, R.C. (2004). MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Res 32:1792-1797.

Fry, B. G. (2005). From genome to "venome" Molecular origin and evolution of the snake venom proteome inferred from phylogenetic analysis of toxin sequences and related body proteins. Cold Spring Harbor Laboratory Press, 403-420. <https://doi.org/10.1101/gr.3228405>

Kuhner, M. K., & Yamato, J. (2015). Practical performance of tree comparison metrics. Systematic Biology, 64(2), 205-214. <https://doi.org/10.1093/sysbio/syu085>

Mackessy, S. P. (2002). Biochemistry and pharmacology of colubrid snake venoms. Journal of Toxicology - Toxin Reviews, 21(1-2), 43-83. <https://doi.org/10.1081/TXR-120004741>

Ojeda, P. G., Ramírez, D., Alzate-Morales, J., Caballero, J., Kaas, Q., & González, W. (2018). Computational studies of snake venom toxins. Toxins, 10(1), 1-24. <https://doi.org/10.3390/toxins10010008>

Pavoine, S., Ollier, S., Pontier, D., & Chessel, D. (2008). Testing for phylogenetic signal in phenotypic traits: new matrices of phylogenetic proximities. Theoretical Population Biology, 73(1), 79-91. <https://doi.org/10.1016/j.tpb.2007.10.001>

R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL <https://www.R-project.org/>

Revell, L. J. (2012). phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol. 3. 217-223. <doi:10.1111/j.2041-210X.2011.00169.x>

Schliep K.P. (2011). phangorn: phylogenetic analysis in R. Bioinformatics, 27(4) 592-593 Tasoulis, T., & Isbister, G. K. (2017). A review and database of snake venom proteomes. Toxins, 9(9). <https://doi.org/10.3390/toxins9090290>

------------------------------------------------------------------------

### SUPPLEMENTARY MATERIALS

#### Appendix A: Crotalinae Venom Proteome Components

Each component described here is a protein family (Tasoulis & Isbister, 2017).

-   PLA2 - Phospholipase A2
-   SVSP - Snake venom serine protease
-   SVMP - Snake venom metalloprotease
-   LAAO - L-amino acid oxidase
-   CRiSP - Cysteine-rich secretory protein
-   CTL - C-type lectin
-   DIS - Disintegrin
-   NP - Natriuretic peptide
-   DEF - Defensin
-   MPi - Snake venom metalloprotease inhibitor

#### Appendix B: Snake Venom Proteome Data Cleaning

The following steps were done to clean data from the Crotalinae venom proteome data set (Tasoulis & Isbister, 2017):

-   Calculated averages for multiple entries for same species of different geographical location (includes *Bothrops Asper* (2) and *Bothrops Atrox* (4))
-   Calculated averages for multiple subspecies into one entry (includes *Crotalis durissus* (3) and *Lachesis muta* (2))
-   Ranges of composition values were replaced based on average of range
-   Any equalities (e.g. &lt;0.1) were replaced to the number in the equality
