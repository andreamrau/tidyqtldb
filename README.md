# tidyqtldb

This repo contains an R script to format [Animal QTLdb](https://www.animalgenome.org/cgi-bin/QTLdb/index) annotations for use as prior information with the BayesRC+ and BayesRCpi models in the [BayesRCO](https://github.com/fmollandin/BayesRCO) software. A full user's guide for BayesRCO can be found [here](https://github.com/fmollandin/BayesRCO/blob/main/doc/BayesRCO.pdf).

The Animal QTLdb collects publicly available quantitative trait loci mapping and candidate gene and association data for several livestock species. More details about the database can be found at the Animal QTLdb [FAQ](https://www.animalgenome.org/QTLdb/faq) and most recent associated publication:

-   Zhi-Liang Hu, Carissa A. Park, and James M. Reecy (2022). [Bringing the Animal QTLdb and CorrDB into the future: meeting new challenges and providing updated services](https://academic.oup.com/nar/article/50/D1/D956/6437998). Nucleic Acids Research, Volume 50, Issue D1, Pages D956–D961. DOI: [10.1093/nar/gkab1116](https://doi.org/10.1093/nar/gkab1116)

For the BayesRC+ and BayesRCpi models, annotation categories must be constructed for each marker in the genotyping data (e.g., `MAP_PigHeat_LD_GEMMA.txt` in this repo) such that each marker in the user's genotype data is assigned to one or more annotation categories. A final "other" category is used to indicate markers that are unannotated for all other categories. For example,

| seqnames | start | end  | snp   | category_1 | category_2 | other |
|----------|-------|------|-------|------------|------------|-------|
| 1        | 2313  | 2314 | SNP_1 | 0          | 0          | 1     |
| 1        | 7593  | 7594 | SNP_2 | 1          | 0          | 0     |
| 1        | 9872  | 9873 | SNP_3 | 1          | 1          | 0     |
| ...      | ...   | ...  | ...   | ...        | ...        | ...   |

In Animal QTLdb, QTLs are organized according to a trait hierarchy (for an example in pig, see [here](https://www.animalgenome.org/cgi-bin/QTLdb/SS/ontrait?class_ID=1)). In this script, we provide the possibility to use either the highest level trait hierarchy ("Level 1"), or the next sublevel trait hierarchy ("Level 2") for a finer level of detail in creating annotation categories.

## Input files

To use this script, three input files are required:

1.  A `bed` file of providing the base-pair data for QTLs in bed format from the species specific Animal QTLdb page. See `QTLdb_pigSS11.bed` in this repo for an example.
2.  A `gff` file of providing the base-pair data for QTLs in gff format from the species specific Animal QTLdb page. See `QTLdb_pigSS11.gff` in this repo for an example.
3.  A comma-delimited file providing the location of genotype markers in bed format, with column names `snp_name`, `pos`, and `chr`. See `MAP_PigHeat_LD_GEMMA.txt` in this repo for an example.

Note that the bed and gff files described above can be directly downloaded from the species-specific page on Animal QTLdb ("All data by bp") for the desired genome build. For example, for pig, these files can be found [here](https://www.animalgenome.org/cgi-bin/QTLdb/SS/index).

## Script arguments

Beyond providing the names of these input files in the R script, users can modify the following arguments in the header:

-   `chr_include`: Chromosome numbers to be included in constructing annotations. Can be used, for example, to eliminate sex chromosomes if desired.
-   `type`: Specifies the type of annotation window to use for annotating markers. If `"qtldb"`, only markers included in the QTLdb are annotated. If `"qtldb_hard"`, the nearest up- and downstream markers are also included in QTLdb annotations (a so-called "hard window""). If `"qtldb_fuzzy"`, the nearest up- and downstream markers are included in both the QTLdb annotation as well as in the "other" category (a so-called "fuzzy window"").
-   `trait_level`: 1 or 2, according to the desired trait hierarchy levels to be used
-   `min_snps`: Minimum number of annotated QTLs in the Animal QTLdb (used to eliminate traits or subtraits with too few annotated markers)
-   `remove_cM`: If `TRUE`, remove Animal QTLdb associations measured in centimorgans
-   `qtldb_padding`: If non-null, number of basepairs to pad QTLdb positions
-   `species`: Species abbreviation code, one of `c("BT", "GG", "CH", "EC", "SS", "OM", "OA")`, to be used for webscraping of trait subhierarchies.

## Citation

If you make use of this script, please cite this repo.
