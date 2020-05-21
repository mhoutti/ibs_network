setwd("~/Dropbox/ibs_network")

library(dplyr)
library(data.table)

# Load files
# Load files
# Taxa
taxa <- as.data.frame(t(read.table(file="data/taxatable.burst.capitalist.reduced.clr.tsv", sep='\t', header = TRUE, row.names=1, skip=0, check.names=FALSE)))
# DEGs
degs1 <- read.table(file = 'data/NormalizedExpression_DEGs_T1.tsv', sep='\t', header = TRUE, skip=0, quote="", check.names=FALSE)
degs2 <- read.table(file = 'data/NormalizedExpression_DEGs_T2.tsv', sep='\t', header = TRUE, skip=0, quote="", check.names=FALSE)
# Metabolome
metabolome <- read.table(file="data/NMR_metabolite_integrals_IBS_ms.csv", sep=",", header=TRUE, skip=0, quote="")
bile_acids <- read.table(file="data/Bile_acid_TICnorm.tsv", sep="\t", header=TRUE, skip=0, quote="")
# Metadata
meta <- read.csv(file="metadata/Final_Cleaned_Condensed_Metadata.csv", header=TRUE, row.names=1, skip=0, check.names=FALSE)

#simplify taxatable column names (keep only most specific)
colnames(taxa) <- gsub(";NA","", colnames(taxa))
colnames(taxa) <- gsub(".*;","", colnames(taxa))
colnames(taxa) <- gsub("s__", "", colnames(taxa))
colnames(taxa) <- gsub("_", " ", colnames(taxa))

# Collapse taxa to be averages by subject
taxa$Subject <- gsub("\\.T.*", "", row.names(taxa))
taxa_collapsed <- aggregate(taxa, by=list(taxa$Subject), mean)
row.names(taxa_collapsed) <- as.numeric(taxa_collapsed$Group.1)
taxa_collapsed <- taxa_collapsed[,2:(ncol(taxa_collapsed)-1)]
taxa <- taxa_collapsed

# Create data frame for mapping from ID on tube to participant's global study ID
meta <- meta[!is.na(meta$study_id),]
tubeid_to_studyid <- data.frame(tube_id=as.character(meta$ID_on_tube), study_id=meta$study_id)
tubeid_to_studyid <- distinct(tubeid_to_studyid)
row.names(tubeid_to_studyid) <- tubeid_to_studyid$study_id

# Collapse and merge metabolome tables
metabolome <- aggregate(metabolome, by=list(metabolome$ID_on_tube), mean, na.rm=T)
bile_acids <- aggregate(bile_acids, by=list(bile_acids$ID_on_tube), mean, na.rm=T)
row.names(metabolome) <- metabolome$ID_on_tube
row.names(bile_acids) <- bile_acids$ID_on_tube

metabolome <- metabolome[,9:ncol(metabolome)]
bile_acids <- bile_acids[,8:13]
metabolome <- merge(metabolome, bile_acids, by=0, all=T)
row.names(metabolome) <- metabolome$Row.names
metabolome <- metabolome[,2:ncol(metabolome)]

# Clean up and collapse DEGs tables
row.names(degs1) <- degs1$GeneName
degs1 <- degs1[,6:ncol(degs1)]
colnames(degs1) <- substring(colnames(degs1), 2)
degs1 <- as.data.frame(t(degs1))
degs1$ID <- row.names(degs1)

row.names(degs2) <- degs2$GeneName
degs2 <- degs2[,6:ncol(degs2)]
colnames(degs2) <- substring(colnames(degs2), 2)
degs2 <- as.data.frame(t(degs2))
degs2$ID <- row.names(degs2)

degs <- rbindlist(list(degs1, degs2), fill=T)
degs <- aggregate(degs, by=list(degs$ID), mean, na.rm=T)
row.names(degs) <- tubeid_to_studyid[degs$Group.1, "tube_id"]
degs <- degs[,2:ncol(degs)]
degs[is.na(degs)] <- NA
degs <- degs[,c(1:77, 79:ncol(degs))]

create_corr_frame <- function(df1, df2){
  #calculate number of rows we will need in our final data frame
  common_rows <- intersect(row.names(df1), row.names(df2))
  df1 <- df1[common_rows,]
  df2 <- df2[common_rows,]
  rows <- ncol(df1)*ncol(df2)
  #create factor vectors for x and y columns (vertices)
  x <- factor(levels=unique(colnames(df1)))
  y <- factor(levels=unique(colnames(df2)))
  #set levels for x and y ahead of time
  #create numeric vectors for correlation and p-value columns
  correlation <- numeric(rows)
  pval <- numeric(rows)
  index <- 1
  for(i in colnames(df1)){
    for(j in colnames(df2)){
      a <- as.numeric(df1[,i])
      b <- as.numeric(df2[,j])
      cor <- cor.test(a, b, method="spearman")
      x[index] <- i
      y[index] <- j
      correlation[index] <- cor$estimate
      pval[index] <- cor$p.value
      index <- index + 1
    }
  }
  dat <- data.frame(x=x, y=y, correlation=correlation, pval=pval)
  dat
}

met_tax <- create_corr_frame(metabolome, taxa)
met_deg <- create_corr_frame(metabolome, degs)
deg_tax <- create_corr_frame(degs, taxa)

# Add fdr-adjusted p-values column
met_tax$fdr_pval <- p.adjust(met_tax$pval, method="fdr")
met_deg$fdr_pval <- p.adjust(met_deg$pval, method="fdr")
deg_tax$fdr_pval <- p.adjust(deg_tax$pval, method="fdr")

# Keep only significant correlations
# Set desired p-value cut-off here
sig_met_tax <- met_tax[met_tax$fdr_pval < 0.25,]
sig_met_deg <- met_deg[met_deg$fdr_pval < 0.25,]
sig_deg_tax <- deg_tax[deg_tax$fdr_pval < 0.25,]

sig_nodes <- unique(c(as.character(sig_met_tax$x), as.character(sig_met_tax$y), as.character(sig_met_deg$x), as.character(sig_met_deg$y), as.character(sig_deg_tax$x), as.character(sig_deg_tax$y)))

met_nodes <- colnames(metabolome)
deg_nodes <- colnames(degs)
tax_nodes <- colnames(taxa)

# Make sure we know what the category of each node is
met_nodes <- data.frame(node=met_nodes, category=rep("Metabolite",length(met_nodes)))
deg_nodes <- data.frame(node=deg_nodes, category=rep("Gene",length(deg_nodes)))
tax_nodes <- data.frame(node=tax_nodes, category=rep("Microbe",length(tax_nodes)))

nodes <- rbind(met_nodes, deg_nodes)
nodes <- rbind(nodes, tax_nodes)

nodes <- nodes[nodes$node %in% sig_nodes,]

luminalt1_network_table <- rbind(sig_met_tax, rbind(sig_met_deg, sig_deg_tax))
luminalt1_network_table <- luminalt1_network_table[complete.cases(luminalt1_network_table),]
colnames(luminalt1_network_table)[1:3] <- c("from", "to", "correlation")
luminalt1_network_table$direction <- ifelse(luminalt1_network_table$correlation < 0, "neg", "pos")
luminalt1_network_table$weight <- abs(luminalt1_network_table$correlation)

library(igraph)
library(SDMTools)
library(ggraph)

luminalt1_net <- graph_from_data_frame(d=luminalt1_network_table, vertices=nodes, directed=F)

# Set edge and node colors

corr_colors <- c("#CC79A7", "#ffffff", "#0072B2") ## Gradient colors for edges
node_colors <- c("#56B4E9", "#F0E442", "#009E73") ## Colors for nodes

# Set node labels in legend

names(node_colors) <- c("Gene", "Metabolite", "Microbe")

# Set node shapes

node_shapes <- c("circle", "triangle", "square")
names(node_shapes) <- c("Gene", "Metabolite", "Microbe")

luminal <- ggraph(luminalt1_net, layout="fr") + 
  geom_edge_link(aes(color = correlation), edge_width=1) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_color_gradientn(limits = c(-1,1), colors=corr_colors) +
  geom_node_point(aes(color=category, shape=category), size=7) +
  scale_color_manual(values=node_colors) +
  scale_shape_manual(values=node_shapes) +
  geom_node_text(aes(label=name), size=3, repel=TRUE) +
  theme_graph() +
  labs(title = "Luminal (FDR adjusted p < 0.25)", color="Category", shape="Category", edge_color="Spearman Correlation")

ggsave("figures/luminal_aggregated_25.pdf", luminal, device=cairo_pdf, width=20, height=10, dpi=300, scale=1.5)
