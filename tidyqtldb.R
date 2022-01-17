# Input    ---------------------------------------------------------------------

## Header = snp_name, pos, chr // Comma-delimited
marker_data <- "MAP_PigHeat_LD_GEMMA.txt"
## Files downloaded from https://www.animalgenome.org/cgi-bin/QTLdb/index
qtl_bed <- "QTLdb_pigSS11.bed"
qtl_gff <- "QTLdb_pigSS11.gff"
chr_include <- 1:18

type <- "qtldb_fuzzy" ## c("qtldb", "qtldb_hard", "qtldb_fuzzy")
trait_level <- 2 ## c(1, 2)
min_snps <- 400
remove_cM <- TRUE
qtldb_padding <- NULL
species <- "SS" # c("BT", "GG", "CH", "EC", "SS", "OM", "OA")

# Load packages ----------------------------------------------------------------
library(plyranges)
library(rvest)
library(tidyr)
library(dplyr)

# 1. Read marker positions  --------------------------------------------------
snp_pos_nosort <- read.table(marker_data, header=TRUE, sep = ",") %>%
    dplyr::select(seqnames = chr, snp=snp_name) %>%
    filter(seqnames %in% chr_include) %>%
    dplyr::select(snp)
snp_pos <- read.table(marker_data, header=TRUE, sep = ",") %>%
    mutate(end = pos + 1) %>%
    dplyr::select(seqnames = chr, start=pos, end, snp=snp_name) %>%
    mutate(seqnames = factor(seqnames)) %>%
    arrange(seqnames, start, end) %>%
    filter(seqnames %in% chr_include)
snp_pos$seqnames <- droplevels(snp_pos$seqnames)
snp_ranges <- as_granges(snp_pos)

# 2. AnimalQTLdb ---------------------------------------------------------------
## TODO: use html_read to download bed and gff

qtldb <- read.table(qtl_bed, skip=20, sep="\t")
colnames(qtldb)[1:4] <- c("seqnames", "start", "end", "trait_type")
qtldb <- qtldb %>%
    separate(seqnames, into=c("chr", "seqnames")) %>%
    separate(trait_type, into=c("trait_type", "misc"), sep=" QTL") %>%
    select(-chr, -c(V5:V12), -misc) %>%
    filter(seqnames %in% chr_include) %>%
    mutate(seqnames = as.character(seqnames)) %>%
    na.omit() %>%
    unique()

qtldb_traits <- read.table(qtl_gff, skip=21, sep="\t")
colnames(qtldb_traits)[c(3,9)] <- c("trait_cat", "info")
qtldb_traits <- qtldb_traits %>% select(trait_cat, info)
qtldb_traits$trait <- NULL
for(i in 1:nrow(qtldb_traits)) {
    tmp <- unlist(strsplit(qtldb_traits$info[i], split = ";"))
    qtldb_traits$trait[i] <- unlist(strsplit(tmp[grep("trait=", tmp)],
                                             split="="))[2]
}
qtldb_traits <- qtldb_traits %>%
    select(trait_cat, trait_type = trait) %>%
    separate(trait_cat, into = c("trait_cat", "type"), sep = "_QTL",
             extra = "drop", fill = "right") %>%
    separate(trait_cat, into = c("trait_cat", "type2"), sep = "_Association",
             extra = "drop", fill = "right") %>%
    separate(trait_cat, into = c("trait_cat", "type3"), sep = "_eQTL",
             extra = "drop", fill = "right") %>%
    select(-type, -type2, -type3) %>%
    unique()
qtldb_full <- left_join(qtldb, qtldb_traits, by = "trait_type")

## Some end coordinates are before start coordinates
switch_index <- which(qtldb_full$end < qtldb_full$start)
qtldb_full[switch_index, c(2:3)] <- qtldb_full[switch_index, c(3:2)]

## Remove any QTLs measured in cM
if(remove_cM) {
    qtldb_full <- qtldb_full[which(qtldb_full$end - qtldb_full$start < 100),]
}

## Define QTLs using a padding window, if desired
if(!is.null(qtldb_padding)) {
    qtldb_full$start <- qtldb_full$start - qtldb_padding
    qtldb_full$end <- qtldb_full$end + qtldb_padding
    qtldb_full$start <- ifelse(qtldb_full$start < 0, 0, qtldb_full$start)
}

cat_levels <- qtldb_full %>% select(trait_cat) %>% unique() %>%
    unlist() %>% as.character() %>% sort()

# 3. Animal QTLdb subcategories (optional) -------------------------------------
if(trait_level == 2) {
    qtldb_website <-
        read_html(paste0("https://www.animalgenome.org/cgi-bin/QTLdb/",
        species, "/ontrait?class_ID=1"))

    web <- qtldb_website %>%
        html_elements(".collapsibleList") %>%
        html_elements("li") %>%
        html_text2()

    category_indices <- c()
    for(i in cat_levels) {
        category_indices <- c(category_indices,
                              grep(paste0("^",
                              gsub(x=i, pattern="_", replacement=" "),
                              " Traits"), web))
    }
    names(category_indices) <- paste0(gsub(x=cat_levels,
                                    pattern="_", replacement=" "),
                                    " Traits")

    trait_cat_subcat <- data.frame(trait_cat=character(0),
                                   trait_subcat = character(0),
                                   trait_type = character(0))
    for(i in category_indices) {
        tmp <- lapply(list(unlist(strsplit(web[i], "-->"))[-1]),
                      function(x) {strsplit(x, "\n")})
        tmpcat <- lapply(tmp[[1]], function(x) x[-1])
        subcat <-
            lapply(tmpcat, function(x) {
               data.frame(trait_cat = names(category_indices)[which(category_indices==i)],
                          trait_subcat = x[1],
                          trait_type =
                              unlist(lapply(strsplit(x[-1], split=" { ",
                                                     fixed=TRUE),
                                            function(x) x[1])))}) %>%
            bind_rows() %>%
            filter(!grepl("^<!--",trait_type))

        trait_cat_subcat <- bind_rows(trait_cat_subcat, subcat)
    }

    ## => collapsing all QTL/eQTL/association into intermediate trait categories
    ## Two outliers Litter weight, total {birth}, Litter weight, total {weaning}
    qtldb2 <- qtldb %>% separate(trait_type, sep = " {",
                                 into=c("trait_type", "misc"),
                                 extra = "drop", fill = "right") %>%
        select(-misc)
    qtldb_full <- left_join(mutate(qtldb2,
                                   trait_type=tolower(trait_type)),
                            mutate(trait_cat_subcat,
                                   trait_type=tolower(trait_type)),
                            by = "trait_type") %>%
        select(-trait_type, -trait_cat) %>% unique()
    ## Only keep subcats with a minimal number of SNPs
    subcat_choice <- names(table(qtldb_full$trait_subcat))[
        which(table(qtldb_full$trait_subcat) >= min_snps)]
    qtldb_full <- qtldb_full %>%
        filter(trait_subcat %in% subcat_choice)

    ## Some end coordinates are before start coordinates
    switch_index <- which(qtldb_full$end < qtldb_full$start)
    qtldb_full[switch_index, c(2:3)] <- qtldb_full[switch_index, c(3:2)]

    ## Remove any QTLs measured in cM
    qtldb_full <-
        qtldb_full[which(qtldb_full$end - qtldb_full$start < 100),]

    qtldb_full$trait_cat <- qtldb_full$trait_subcat
    cat_levels <- qtldb_full %>% select(trait_cat) %>% unique() %>%
        unlist() %>% as.character() %>% sort()
}

# 4. Match markers to qtldb ----------------------------------------------------

qtldb_level <- olap_level <- vector("list", length(cat_levels))
names(qtldb_level) <- names(olap_level) <- cat_levels
for(i in cat_levels) {
    qtldb_level[[i]] <- qtldb_full %>%
        filter(trait_cat == i) %>%
        as_granges()
    olap_level[[i]] <-
        join_overlap_inner(qtldb_level[[i]], snp_ranges) %>%
        as.data.frame() %>%
        mutate(type = i)
}
olap_all_level <- bind_rows(olap_level[[cat_levels[1]]],
                              olap_level[[cat_levels[2]]])
for(i in 3:length(cat_levels)) {
    olap_all_level <- bind_rows(olap_all_level, olap_level[[i]])
}
olap <- olap_all_level %>%
    select(snp, type) %>%
    mutate(values = 1) %>%
    unique()

annot <- pivot_wider(olap, names_from=type, values_from = values,
                     values_fn = length) %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    left_join(snp_pos, ., by="snp")
tmp_calc <- annot[, cat_levels[1]]
## Create other category and match up with original order of SNPs
annot <- annot %>%
    mutate(other = ifelse(is.na(tmp_calc), 1, 0)) %>%
    mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
    left_join(snp_pos_nosort, ., by = "snp") %>%
    dplyr::select(seqnames, start, end, snp, everything())
colnames(annot) <- gsub(" ", "_", colnames(annot))

if(type != "qtldb") {
    ## Fuzzy windows
    annot_reorder <- annot %>% arrange(seqnames, start, end)
    annot_plus <- annot_reorder[0,] ## Initialize
    ncol1 <- ncol(annot)-1
    for(chr in chr_include) {
        annot_chr <- annot_chr_old <- annot_reorder %>% filter(seqnames == chr)
        for(i in 1:nrow(annot_chr_old)) {
            if(i == 1) {
                annot_chr[i,5:ncol1] <- sign(colSums(annot_chr_old[c(i,i+1),5:ncol1]))
            } else if(i == nrow(annot_chr_old)) {
                annot_chr[i,5:ncol1] <- sign(colSums(annot_chr_old[c(i-1,i),5:ncol1]))
            } else {
                annot_chr[i,5:ncol1] <- sign(colSums(annot_chr_old[c(i-1,i,i+1),5:ncol1]))
            }
        }
        annot_plus <- bind_rows(annot_plus, annot_chr)
    }
    annot_subcat_fuzzywindow <- left_join(snp_pos_nosort, annot_plus, by = "snp") %>%
        dplyr::select(seqnames, start, end, snp, everything())
    colnames(annot_subcat_fuzzywindow) <-
        gsub(" ", "_", colnames(annot_subcat_fuzzywindow))

    ## Hard windows
    nc <- ncol(annot_subcat_fuzzywindow)
    flip_index <- which(annot_subcat_fuzzywindow$other == 1 &
                            rowSums(annot_subcat_fuzzywindow[,-c(1:4,nc)])>0)
    annot_subcat_hardwindow <- annot_subcat_fuzzywindow
    annot_subcat_hardwindow[flip_index,nc] <- 0
    colnames(annot_subcat_hardwindow) <-
        gsub(" ", "_", colnames(annot_subcat_hardwindow))
}

## Write out final results
if(type == "qtldb") annot_final <- annot
if(type == "qtldb_hard") annot_final <- annot_subcat_hardwindow
if(type == "qtldb_fuzzy") annot_final <- annot_subcat_fuzzywindow

colSums(annot_final[-c(1:4)])

write.table(x=annot_final,
            file=paste0("annotation_level-", trait_level, "_",
                        "type-", type, ".txt"),
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")


