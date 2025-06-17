
# Packages ----------------------------------------------------------------

library(tidyverse)
library(phyloseq)
library(janitor)
library(readr)
library(readxl)
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library
library(data.table)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)
install.packages(
  "microViz",
  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
)
library(microViz)
library(RColorBrewer)
library(vegan)
library(ape)
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
library(ggtree)

# Read in Data ------------------------------------------------------------

tax_table(taxonomy1)

table_merged <- read_qza("C:/Users/Robert.Bremer/Documents/seusmbon-merged1thru18-16s-250313/merged_table_16s.qza")

tax_merged <- read_qza("C:/Users/Robert.Bremer/Documents/seusmbon-merged1thru18-16s-250313/merged_taxonomy_16s.qza")
# Read in Metadata --------------------------------------------------------

metadata <- read_excel("C:/Users/Robert.Bremer/Downloads/RB Copy of seusmbon_NOAA_MIMARKS_v1.0.8.survey.water.6.0.xlsx", sheet = "water_sample_data") %>%
  row_to_names(8)

# Combine into Phyloseq Object --------------------------------------------

ps2 <- qza_to_phyloseq(features = "C:/Users/Robert.Bremer/Documents/seusmbon-merged1thru18-16s-250313/merged_table_16s.qza",
                       taxonomy = "C:/Users/Robert.Bremer/Documents/seusmbon-merged1thru18-16s-250313/merged_taxonomy_16s.qza",
                       metadata = "water_sample_metadata.txt",
                       tree = "C:/Users/Robert.Bremer/Documents/phylogeny-align-to-tree-mafft-fasttree/rooted_tree.qza")

ps2

?qza_to_phyloseq
# Reading the tree ---------------------------------------------------------
tree <- read_qza("C:/Users/Robert.Bremer/Documents/phylogeny-align-to-tree-mafft-fasttree/tree.qza")$data
plot_tree(tree)
ggtree::ggtree(tree)
tree

# Preprocessing -----------------------------------------------------------
# What do I do for preprocessing?
sample_sums(ps2)

otu_table(ps2) %>%
  `@`(., ".Data") %>% 
  t() %>% 
  vegan::rarecurve(.,step = 10000)

otu_table(ps2)

mat <- t(otu_table(ps2))
class(mat) <- "matrix"
#> Warning message:
#> In class(mat) <- "matrix" :
#>  Setting class(x) to "matrix" sets attribute to NULL; result will no longer be an S4 object
class(mat)
#> [1] "matrix" "array"

mat <- as(t(otu_table(ps2)), "matrix")
class(mat)
#> [1] "matrix" "array"

raremax <- min(rowSums(mat))
mat

system.time(rarecurve(mat, step = 100, sample = raremax, col = "blue", label = FALSE))


# Visualize Read Depth

sdt = data.table(as(sample_data(ps2), "data.frame"),
                 TotalReads = sample_sums(ps2), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")

seqDepth = ggplot(orderedSamples, aes(y = TotalReads, x = serial_number)) + geom_point() + scale_y_continuous(limits = c(0,30000))
seqDepth

orderedSamples <- sdt %>%
  arrange(TotalReads)

ggplot(orderedSamples, aes(x = as.numeric(row.names(orderedSamples)), y = TotalReads)) +
  geom_point() +
  geom_hline(yintercept = 2006)

ggplot(orderedSamples, aes(x = as.numeric(row.names(orderedSamples)), y = TotalReads)) +
  geom_point() +
  geom_hline(yintercept = 2006)+
  scale_y_continuous(limits = c(0,10000))

# Look at the samples you would be cutting if you rarified here
orderedSamples2 <- orderedSamples %>%
  subset(select = c(SampleID, cruise_id_edna, TotalReads, serial_number)) %>%
  mutate(BlankType = case_when(
    grepl("NTC", SampleID) ~ "NTC",
    grepl("Neg", SampleID) ~ "Negative Control",
    grepl("Blank", SampleID) ~ "Extraction Blank",
    .default = "Not a negative control"))

ggplot(filter(orderedSamples2, TotalReads< 2006), aes(x = cruise_id_edna))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 89 samples are less than 2006 Reads, 20 are NTCs
orderedSamplesBlanks <- orderedSamples2[grep("NTC|Blank|Neg", orderedSamples2$SampleID)] %>%
  mutate(BlankType = case_when(
    grepl("NTC", SampleID) ~ "NTC",
    grepl("Neg", SampleID) ~ "Negative Control",
    grepl("Blank", SampleID) ~ "Extraction Blank",
    .default = "Not a negative control"))

orderedSamples2 <- orderedSamples[grep("NTC|Blank|Neg", orderedSamples$SampleID)] %>%
  mutate(BlankType = case_when(
    grepl("NTC", SampleID) ~ "NTC",
    grepl("Neg", SampleID) ~ "Negative Control",
    grepl("Blank", SampleID) ~ "Extraction Blank",
    .default = "Not a negative control"))

ggplot(orderedSamplesBlanks, aes(x = SampleID, y = TotalReads, color = BlankType))+
  geom_point()+
  geom_hline(yintercept = 2006)

ggplot(orderedSamples2, aes(x = SampleID, y = TotalReads))+
  geom_point(aes(color = factor(BlankType)))+
  geom_hline(yintercept = 2006)+
  scale_color_manual(values =c("red","green","lightgrey","blue"))+
  theme(axis.text.x = element_blank())
  

# Sequencing Depth

pSeqDepth = ggplot(sdt, aes(TotalReads, color = )) + geom_histogram(binwidth = 20) + ggtitle("Sequencing Depth")
pSeqDepth

pSeqDepth + facet_wrap(~cruise_id_edna)

pSeqDepth + facet_wrap(line_name ~ .)

# 
ps2
ntaxa(ps2)
sample_variables(ps2)
rank_names(ps2)
sample_names(ps2)

# Rarefaction
rarecurve(t(otu_table(ps2)), step=5, cex=0.5)
# rarefy without replacement
#ps2.rarefied = rarefy_even_depth(ps2, rngseed=1, sample.size=0.9*min(sample_sums(ps)), replace=F)
ps2.rarefied = rarefy_even_depth(ps2, rngseed=1, sample.size=2006, replace=F)

# Looking at the rarefied object
ps2.rarefied
otu_table(ps2.rarefied)

sdt.rarefied = data.table(as(sample_data(ps2.rarefied), "data.frame"),
                 TotalReads = sample_sums(ps2.rarefied), keep.rownames = TRUE)
pSeqDepth.rarefied = ggplot(sdt.rarefied, aes(TotalReads, color = )) + geom_histogram(binwidth = 20) + ggtitle("Sequencing Depth")
pSeqDepth.rarefied
pSeqDepth


# Visualizing

# Just OTUs
ps2.ord <- ordinate(ps2.rarefied, "NMDS", "bray")
p2_ord = plot_ordination(ps2.rarefied, ps2.ord, type="taxa", color="Phylum", title="taxa")
print(p2_ord)

# By cruise type
p2_ord = plot_ordination(ps2.rarefied, ps2.ord, type="samples", color="cruise_id_edna") 
p2_ord + geom_polygon(aes(fill=cruise_id_edna)) + geom_point(size=5) + ggtitle("samples")

?plot_ordination
??clean_pq

# Alpha diversity  --------
#Alpha diversity (eg, richness, Shannon) by month/year and site/depth https://github.com/aomlomics/seusmbon/issues/2
# Alpha diversity by raw collection date
plot_richness(ps2.rarefied, x="collection_date", measures=c("Observed", "Shannon")) +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_blank())
          
          element_text(angle = -90))axis.text.x = element.blank())

# Alpha diversity by collection month

ps2.1 <- ps2.rarefied %>%
  ps_mutate(collection_date_mod = gsub("T", " ", collection_date)) %>%
  ps_mutate(collection_date_mod = gsub("Z", "", collection_date_mod)) %>%
  ps_mutate(collection_date_mod = as.Date(collection_date_mod, format = ("%Y-%m-%d"))) %>%
  ps_mutate(collection_year = format(as.Date(collection_date_mod, format="%d/%m/%Y"),"%Y")) %>%
  ps_mutate(collection_month = format(as.Date(collection_date_mod, format="%d/%m/%Y"),"%m"))
  
ps1.1@sam_data[["collection_date_mod"]]
ps1.1@sam_data[["collection_year"]]
ps1.1@sam_data[["collection_month"]]

plot_richness(ps2.1, x = "collection_month", measures=c("Observed", "Shannon")) +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

# Alpha diversity of collection year

plot_richness(ps2.1, x = "collection_year", measures=c("Observed", "Shannon")) +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

# Alpha diversity by collection transect

plot_richness(ps2.1, x = "line_name", measures=c("Observed", "Shannon")) +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

# Alpha diversity by collection site (station)

plot_richness(ps2.1, x = "station", measures=c("Observed", "Shannon")) +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

# Alpha diversity by collection depth
plot_richness(ps2.1, x = "ctd_bottle_no", measures=c("Observed", "Shannon")) +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

# Alpha diversity by run
plot_richness(ps2.1, x = "plate_number", measures=c("Observed", "Shannon")) +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

# Beta Diversity ----------------------------------------------------------

#Beta diversity (eg, PCoA) by month/year and site/depth https://github.com/aomlomics/seusmbon/issues/3

# Bray Curtis, PCOA

dist = phyloseq::distance(ps2.1, method="wunifrac")
ordination = ordinate(ps2.1, method="PCoA", distance=dist)

# Beta Diversity by Year
plot_ordination(ps2.1, ordination, color="collection_year") + 
  theme_classic() +
  theme(strip.background = element_blank())

# Beta Diversity by Month
plot_ordination(ps1.1, ordination, color="collection_month") + 
  theme_classic() +
  theme(strip.background = element_blank())

# Beta Diversity by Line
plot_ordination(ps1.1, ordination, color="line_name") + 
  theme_classic() +
  theme(strip.background = element_blank())

# Beta Diversity by Station
plot_ordination(ps1.1, ordination, color="station") + 
  theme_classic() +
  theme(strip.background = element_blank())

# Beta Diversity by Depth
plot_ordination(ps1.1, ordination, color="ctd_bottle_no") + 
  theme_classic() +
  theme(strip.background = element_blank())

# Jaccard NMDS
dist = phyloseq::distance(ps1.1, method="jaccard", binary = TRUE)
ordination = ordinate(ps1.1, method="NMDS", distance=dist)

#Taxa diversity (eg, bubble plots with error bars) by month/year and site/depth https://github.com/aomlomics/seusmbon/issues/4

# Taxa level that will be aggregated and displayed (Person 1 Graph)
tax_aggr <- "Class"
tax_number <- 30 #threshold for how many taxa get displayed; others will be pooled into "other"
tax_col <- "Phylum"

?tax_glom
percentages_glom_phylum <- tax_glom(ps2.rarefied, taxrank = 'Phylum')
View(percentages_glom_phylum@tax_table@.Data)

percentages_glom_genus <- tax_glom(ps2.rarefied, taxrank = 'Genus')
View(percentages_glom_genus@tax_table@.Data)

percentages_glom_order <- tax_glom(ps2.rarefied, taxrank = 'Order')
View(percentages_glom_order@tax_table@.Data)

percentages_df <- psmelt(percentages_glom_phylum)
str(percentages_df)

absolute_glom <- tax_glom(physeq = ps2.rarefied, taxrank = "Phylum")
absolute_df <- psmelt(absolute_glom)
str(absolute_df)

absolute_df$Phylum <- as.factor(absolute_df$Phylum)
phylum_colors_abs<- colorRampPalette(brewer.pal(8,"Dark2")) (length(levels(absolute_df$Phylum)))

absolute_plot <- ggplot(data= absolute_df, aes(x=line_name, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  scale_fill_manual(values = phylum_colors_abs) +
  theme(axis.text.x = element_text(angle = 45, vjust = .5, hjust = .5))

absolute_plot

percentages_df$Phylum <- as.factor(percentages_df$Phylum)
phylum_colors_rel<- colorRampPalette(brewer.pal(8,"Dark2")) (length(levels(percentages_df$Phylum)))
relative_plot <- ggplot(data=percentages_df, aes(x=line_name, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  scale_fill_manual(values = phylum_colors_rel)+
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

# Transforming into relative
ps2.rarefied.relative <- transform_sample_counts(ps2.rarefied, function(x) x /sum(x))
otu_table(ps2.rarefied.relative)

ps2_relative_phylum <- ps2.rarefied.relative %>%
  tax_glom(., taxrank = 'Phylum') 

ps2_relative_genus <- ps2.rarefied.relative %>%
  tax_glom(., taxrank = "Genus")

ps2_relative_species <- ps2.rarefied.relative %>%
  tax_glom(., taxrank = "Species")

# Merge into a different categorical variable like line_id or cruise or season or year
#Year
glom_year_phylum <- merge_samples(ps2_relative_phylum, "collection_date") %>%
  transform_sample_counts(., function(x) x/sum(x)) %>%
  psmelt(.)
glom_year_phylum$Phylum[glom_year_phylum$Abundance < .01] <- "<1% abund"

#Line
glom_line_phylum <- merge_samples(ps2_relative_phylum, "line_id") %>%
  transform_sample_counts(., function(x) x/sum(x)) %>%
  psmelt(.)

glom_line_genus <- merge_samples(ps2_relative_genus, "line_id") %>%
  transform_sample_counts(., function(x) x/sum(x)) %>%
  psmelt(.)
glom_line_genus$Genus[glom_line_genus$Abundance < .01] <- "<1% abund"

glom_line_species <- merge_samples(ps2_relative_species, "line_id") %>%
  transform_sample_counts(., function(x) x/sum(x)) %>%
  psmelt(.)
glom_line_species$Species[glom_line_species$Abundance < .01] <- "<1% abund"

#Cruise
glom_cruise_phylum <- merge_samples(ps2_relative_phylum, "cruise_id") %>%
  transform_sample_counts(., function(x) x/sum(x)) %>%
  psmelt(.)

glom_cruise_genus <- merge_samples(ps2_relative_genus, "cruise_id") %>%
  transform_sample_counts(., function(x) x/sum(x)) %>%
  psmelt(.)
glom_cruise_genus$Genus[glom_cruise_genus$Abundance < .01] <- "<1% abund"
  
glom_cruise_species <- merge_samples(ps2_relative_species, "cruise_id") %>%
  transform_sample_counts(., function(x) x/sum(x)) %>%
  psmelt(.)

#Year
glom_year_phylum <- merge_samples(ps2_relative_phylum, "collection_year") %>%
  transform_sample_counts(., function(x) x/sum(x)) %>%
  psmelt(.)

# Colors
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- append(col_vector, phylum_colors)
col_vector

phylum_colors <- c(
  "grey22", "darkcyan", "orchid1", "green", "orange", "blue", "tomato2", "olivedrab", "grey47",
  "cyan", "coral3", "darkgreen", "magenta", "palegoldenrod", "dodgerblue", "firebrick", "yellow", "purple4",
  "lightblue", "grey77", "mediumpurple1", "tan4", "red", "darkblue", "yellowgreen", "gold")
glom_line_phylum$Phylum <- as.factor(glom_line_phylum$Phylum)
phylum_colors_rel <-colorRampPalette(brewer.pal(8,"Set2")) (length(levels(glom_line_phylum$Phylum)))

#relative_plots
# Year
relative_year_phylum <- ggplot(data=glom_year_phylum, aes(x=Sample, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  #scale_color_brewer(palette = "Set1")+
  scale_fill_manual(values = col_vector)+
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -45, hjust = .25, vjust = .5))

relative_year_phylum

#Line
relative_line_phylum <- ggplot(data=glom_line_phylum, aes(x=Sample, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  #scale_color_brewer(palette = "Set1")+
  scale_fill_manual(values = col_vector)+
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -45, hjust = .25, vjust = .5))

relative_line_phylum

relative_line_genus <- ggplot(data=glom_line_genus, aes(x=Sample, y=Abundance, fill=Genus))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  #scale_color_brewer(palette = "Set1")+
  scale_fill_manual(values = col_vector)+
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -45, hjust = .25, vjust = .5))

relative_line_genus

relative_line_species <- ggplot(data=glom_line_species, aes(x=Sample, y=Abundance, fill=Species))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  #scale_color_brewer(palette = "Set1")+
  scale_fill_manual(values = col_vector)+
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -45, hjust = .25, vjust = .5))

relative_line_species

#Cruise
relative_cruise_phylum <- ggplot(data=glom_cruise_phylum, aes(x=Sample, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  scale_color_brewer(palette = "Set1")+
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

relative_cruise_phylum

relative_cruise_genus <- ggplot(data=glom_cruise_genus, aes(x=Sample, y=Abundance, fill=Genus))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  #scale_color_brewer(palette = "Set1")+
  scale_fill_manual(values = col_vector)+
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

relative_cruise_genus

absolute_plot | relative_plot

plot_bar(ps2, "Family", facet_grid=~year)


#Redundancy analysis with biplots https://github.com/aomlomics/seusmbon/issues/5
#Deicode with biplots https://github.com/aomlomics/seusmbon/issues/6
#General linear models (GLMs) https://github.com/aomlomics/seusmbon/issues/7
#Balance analysis https://github.com/aomlomics/seusmbon/issues/8
#Network analysis https://github.com/aomlomics/seusmbon/issues/9
#Copy number analysis https://github.com/aomlomics/seusmbon/issues/10
#Phylogenetic analysis to detect artifacts or outliers https://github.com/aomlomics/seusmbon/issues/11
#Identify “undefined” sequences (eg, BLAST against GenBank) https://github.com/aomlomics/seusmbon/issues/12
#Volume filtered effects https://github.com/aomlomics/seusmbon/issues/13
#Reference database effects (eg, rCRUX, weighted/regional classifiers) https://github.com/aomlomics/seusmbon/issues/14
#Taxon-specific analyses: Copepods https://github.com/aomlomics/seusmbon/issues/15
#Taxon-specific analyses: Karenia vs. salinity or iron https://github.com/aomlomics/seusmbon/issues/16
#Taxon-specific analyses: Spinning fish https://github.com/aomlomics/seusmbon/issues/17
