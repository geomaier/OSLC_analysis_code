#############################
# Analysis for CGG014_Maier #
#############################

library(arnaq)
library(ggplot2)
source("edgeR_Analysis.R")
library(ape)

do.prelim <- F
do.overview <- F
do.process.gtf <- F
do.main <- F
do.extra.de <- F
do.extra.de.2 <- T
do.gene.list <- F

if(do.prelim){
	# prelim run to build everything
	arnaq(
		treat.groups=c("Group"),
		pca.groups=c("Group", "Position", "Batch", "BRAF", "EscherichiaShigella", "Fusobacterium", "Serratia"),
		svg.export=F
	)
}

# outliers
outlier_mask <- make_outlier_mask("Undetermined") # whatever this is, we do not want it
removed_sample_mask <- sample.metadata$Group!="DELETE"

# Auxiliary groups
sample.metadata$H_ADJ_PooledGroup <- as.character(sample.metadata$Group)
sample.metadata$H_ADJ_PooledGroup[sample.metadata$Group %in% c("H_ADJ_SSL",
   "H_ADJ_HP")] <- "H_ADJ_SSL_HP"
sample.metadata$H_ADJ_PooledGroup <- factor(sample.metadata$H_ADJ_PooledGroup)

sample.metadata$PooledGroup <- as.character(sample.metadata$H_ADJ_PooledGroup)
sample.metadata$PooledGroup[sample.metadata$Group %in% c("SSL", "HP")] <- "SSL_HP"
sample.metadata$PooledGroup <- factor(sample.metadata$PooledGroup)

sample.metadata$SSL_HP_PooledGroup <- as.character(sample.metadata$Group)
sample.metadata$SSL_HP_PooledGroup[sample.metadata$Group %in% c("SSL", "HP")] <- "SSL_HP"
sample.metadata$SSL_HP_PooledGroup <-factor(sample.metadata$SSL_HP_PooledGroup)

sample.metadata$Pairing <- paste(sample.metadata$Position, sample.metadata$Patient, sep="_")
sample.metadata$Pairing[sample.metadata$Group=="DELETE"] <- "UNPAIRED"

if(do.overview){

	# all
	arnaq(model.name="all", sample.mask=outlier_mask & removed_sample_mask,
		treat.groups=c("Group"),
		pca.groups=c("Group", "Position", "Batch", "BRAF", "EscherichiaShigella", "Fusobacterium", "Serratia"),
		svg.export=T)

	# stomach and CA only
	antrum_mask <- sample.metadata$Group=="ANTRUM"
	ca_mask <- sample.metadata$Group=="CA"
	arnaq(model.name="antrum_and_CA", 
		sample.mask=(antrum_mask | ca_mask),
		treat.groups=c("Group"),
		pca.groups=c("Group", "Position", "Batch"),
		svg.export=T
	)

	# Non-stomach
	main_mask <- sample.metadata$Group %in% c("H_COL",
		"H_ADJ_SSL", "H_ADJ_HP", "SSL", "HP", "CA")
	arnaq(model.name="main", sample.mask=main_mask,
		treat.groups=c("Group", "Position"),
		pca.groups=c("Group", "Position", "Batch", "BRAF", "EscherichiaShigella", "Fusobacterium", "Serratia"),
		svg.export=T
	)
}

project.id <- arnaq.run$project.id

if(do.process.gtf){
	# Get gene info from gtf
	cat("Processing gtf file...\n")
	gtf <- read.gff("/well/bsg/resources/RNA_pipeline/genomes/GRCh38.113/gtf/Homo_sapiens.GRCh38.113.gtf", GFF3=T)
	tmp <- gtf[gtf$type=="exon",9]
	tmp2 <- gsub(".*transcript_id \"","",tmp)
	tmp2 <- gsub("\".*","",tmp2)
	tmp3 <- gsub(".*gene_name \"","",tmp)
	tmp3 <- gsub("\".*","",tmp3)
	tmp4 <- gsub(".*gene_id \"","",tmp)
	tmp4 <- gsub("\".*","",tmp4)
	name_map <- data.frame(transcript_id=tmp2, gene_name=tmp3, ensembl_id=tmp4)
	name_map <- unique(name_map)
	name_map$gene_name <- ifelse(name_map$gene_name=="gene_id ",
		name_map$ensembl_id, name_map$gene_name)
	cat("Finished processing gtf file\n")
}

rna.volcano <- function(results.table, p.thresh=0.05, fc.thresh=1, top.n.to.label=10) {
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggrepel))
  tmp <- results.table
	tmp <- tmp[tmp$logCPM>2,]
  tmp$LogSig <- -log10(tmp$FDR)
  tmp$Significance <- ifelse(tmp$FDR<0.05, ifelse(abs(tmp$logFC)>fc.thresh, "HI_FC_SIG",
                                                  "LO_FC_SIG"), "NONSIG")
  tmp$Significance <- factor(tmp$Significance, levels=c("HI_FC_SIG", "LO_FC_SIG",
														 "NONSIG"))
  tmp$Label <- tmp$SYMBOL
  tmp$Label[is.na(tmp$Label)] <- tmp$genes[is.na(tmp$Label)]
  tmp$Label[tmp$Significance=="NONSIG"] <- NA
  tmp$Label[(top.n.to.label+1):nrow(tmp)] <- NA
  out <- ggplot(tmp, aes(x=logFC, y=LogSig, col=Significance, label=Label)) + 
    geom_point(size=0.5, show.legend=T) +
    geom_text_repel(color="black", size=3, na.rm=TRUE) +
    scale_color_manual(values=c("red", "blue", "#606060"), drop=F) +
    ylab("-log10 FDR")
  out
}

if(do.main){

	custom.gene.list <- read.table("custom_gene_list_2.txt", header=F, stringsAsFactor=F)[,1]

	do.group <- function(model.name, sample.mask,
			contrast, normalise="VSN", group="Group"){

		cat(paste("*** Analysis contrast:", model.name, "\n"))

		arnaq(model.name=model.name, sample.mask=sample.mask,
			treat.groups=group,
			pca.groups=c(group, "Position", "Batch", "BRAF", "Gender", "EscherichiaShigella", "Fusobacterium", "Serratia"),
			normalise=normalise,
			#pca.depth=4,
			svg.export=F
		)

		tmp <- analyse(model.name=model.name, data.matrix=count.data,
			col.data=sample.metadata,
			species="Hs",
			sample.mask=sample.mask,
			analysis.groups=group,
			contrast=contrast
		)

		tmp.vol <- rna.volcano(tmp$results[[1]]$table)
		ggsave(paste0("results/volcano/CGG014_", model.name, "_volcano.png"), tmp.vol)

		cat(paste("** GSE analysis\n"))
		gse.name <- paste0("results/gse/CGG014_", gsub(" ", "_", gsub("-", "v", contrast)))
		tmp.gse <- analyse.cP(tmp$results[[1]]$table, gse.name)

		return(tmp)

	}

	do.BRAF <- function(model.name, sample.mask,
			contrast="V600E - WT", normalise="VSN"){

		cat(paste("*** Analysis contrast:", model.name, "\n"))

		arnaq(model.name=model.name, sample.mask=sample.mask,
			treat.groups="BRAF",
			pca.groups=c("Position", "Batch", "BRAF", "Gender", "EscherichiaShigella", "Fusobacterium", "Serratia"),
			normalise=normalise,
			#pca.depth=4,
			svg.export=F
		)

		tmp <- analyse(model.name=model.name, data.matrix=count.data,
			col.data=sample.metadata, species="Hs",
			sample.mask=sample.mask,
			analysis.groups="BRAF",
			contrast=contrast
		)

		tmp.vol <- rna.volcano(tmp$results[[1]]$table)
		ggsave(paste0("results/volcano/CGG014_", model.name, "_volcano.png"), tmp.vol)

		cat(paste("** GSE analysis\n"))
		gse.name <- paste0("results/gse/CGG014_", gsub(" ", "_", gsub("-", "v", contrast)))
		tmp.gse <- analyse.cP(tmp$results[[1]]$table, gse.name)

		return(tmp)

	}

	do.paired <- function(model.name, sample.mask,
			contrast, normalise="VSN", group="Group"){

		cat(paste("*** Paired analysis contrast:", model.name, "\n"))

		arnaq(model.name=model.name, sample.mask=sample.mask,
			treat.groups=group,
			pca.groups=c(group, "Position", "Batch", "BRAF", "Gender", "EscherichiaShigella", "Fusobacterium", "Serratia"),
			normalise=normalise,
			#pca.depth=4,
			svg.export=F
		)

		tmp <- analyse(model.name=model.name, data.matrix=count.data,
			col.data=sample.metadata,
			species="Hs",
			sample.mask=sample.mask,
			analysis.groups=group,
			extra.groups=c("Pairing"),
			contrast=contrast
		)

		tmp.vol <- rna.volcano(tmp$results[[1]]$table)
		ggsave(paste0("results/volcano/CGG014_", model.name, "_volcano.png"), tmp.vol)

		cat(paste("** GSE analysis\n"))
		gse.name <- paste0("results/gse/CGG014_", gsub(" ", "_", gsub("-", "v", contrast)))
		tmp.gse <- analyse.cP(tmp$results[[1]]$table, gse.name)

		return(tmp)

	}

	# Masks
	HCOL.mask <- sample.metadata$Group == "H_COL"
	HADJSSL.mask <- sample.metadata$Group == "H_ADJ_SSL"
	HADJHP.mask <- sample.metadata$Group == "H_ADJ_HP"
	HADJCA.mask <- sample.metadata$Group == "H_ADJ_CA"
	SSL.mask <- sample.metadata$Group == "SSL"
	HP.mask <- sample.metadata$Group == "HP"
	CA.mask <- sample.metadata$Group == "CA"
	ANTRUM.mask <- sample.metadata$Group == "ANTRUM"

	V600E.mask <- sample.metadata$BRAF == "V600E"
	WT.mask <- sample.metadata$BRAF == "WT"

	# Mask for pairing
	tmp <- table(sample.metadata$Pairing)
	tmp <- names(tmp)[tmp==2]
	pairing.mask <- (sample.metadata$Pairing %in% tmp) & 
									 sample.metadata$Pairing != "UNPAIRED"

	# Analyses (annotated with index in list in New_Indications document

	# Unpooled by group
	SSL_HCOL.res <- do.group("SSL_v_H_COL", HCOL.mask | SSL.mask, "SSL - H_COL") ## 1
	SSL_HP.res <- do.group("SSL_v_HP", SSL.mask | HP.mask, "SSL - HP") ## 3
	SSL_HADJSSL.res <- do.group("SSL_v_H_ADJ_SSL", SSL.mask | HADJSSL.mask, "SSL - H_ADJ_SSL") ## 4
	SSL_CA.res <- do.group("SSL_v_CA", SSL.mask | CA.mask, "SSL - CA") ## 5
	SSL_HADJHP.res <- do.group("SSL_v_H_ADJ_HP", SSL.mask | HADJHP.mask, "SSL - H_ADJ_HP") ## 6
	HP_HCOL.res <- do.group("HP_v_H_COL", HCOL.mask | HP.mask, "HP - H_COL") ## 8
	HP_HADJSSL.res <- do.group("HP_v_H_ADJ_SSL", HP.mask | HADJSSL.mask, "HP - H_ADJ_SSL") ## 10
	HP_CA.res <- do.group("HP_v_CA", HP.mask | CA.mask, "HP - CA") ## 11
	HP_HADJHP.res <- do.group("HP_v_H_ADJ_HP", HP.mask | HADJHP.mask, "HP - H_ADJ_HP") ## 12
	CA_HCOL.res <- do.group("CA_v_H_COL", CA.mask | HCOL.mask, "CA - H_COL") ## 20
	CA_HADJSSL.res <- do.group("CA_v_H_ADJ_SSL", HADJSSL.mask | CA.mask, "CA - H_ADJ_SSL") ## 21
	CA_HADJHP.res <- do.group("CA_v_H_ADJ_HP", HADJHP.mask | CA.mask, "CA - H_ADJ_HP") ## 22
	ANTRUM_HCOL.res <- do.group("ANTRUM_v_H_COL", ANTRUM.mask | HCOL.mask, "ANTRUM - H_COL") ## 24
	ANTRUM_HP.res <- do.group("ANTRUM_v_HP", HP.mask | ANTRUM.mask, "ANTRUM - HP") ## 25
	ANTRUM_HADJSSL.res <- do.group("ANTRUM_v_H_ADJ_SSL", HADJSSL.mask | ANTRUM.mask, "ANTRUM - H_ADJ_SSL") ## 26
	ANTRUM_CA.res <- do.group("ANTRUM_v_CA", CA.mask | ANTRUM.mask, "ANTRUM - CA") ## 27
	ANTRUM_HADJHP.res <- do.group("ANTRUM_v_H_ADJ_HP", HADJHP.mask | ANTRUM.mask, "ANTRUM - H_ADJ_HP") ## 28
	ANTRUM_SSL_V600E.res <- do.group("ANTRUM_v_SSL_V600E", ANTRUM.mask | (SSL.mask & V600E.mask), "ANTRUM - SSL") ## 30
	ANTRUM_SSL_WT.res <- do.group("ANTRUM_v_SSL_WT", ANTRUM.mask | (SSL.mask & WT.mask), "ANTRUM - SSL") ## 31
	ANTRUM_HP_V600E.res <- do.group("ANTRUM_v_HP_V600E", ANTRUM.mask | (HP.mask & V600E.mask), "ANTRUM - HP") ## 32
	ANTRUM_HP_WT.res <- do.group("ANTRUM_v_HP_WT", ANTRUM.mask | (HP.mask & WT.mask), "ANTRUM - HP") ## 33
	SSL_V600E_HCOL.res <- do.group("SSL_V600E_v_H_COL",	HCOL.mask | (SSL.mask & V600E.mask), "SSL - H_COL") # 42
	SSL_V600E_HADJSSL.res <- do.group("SSL_V600E_v_H_ADJ_SSL",	HADJSSL.mask | (SSL.mask & V600E.mask), "SSL - H_ADJ_SSL") # 43
	SSL_V600E_HADJHP.res <- do.group("SSL_V600E_v_H_ADJ_HP", HADJHP.mask | (SSL.mask & V600E.mask), "SSL - H_ADJ_HP") # 45
	SSL_V600E_CA.res <- do.group("SSL_V600E_v_CA", CA.mask | (SSL.mask & V600E.mask), "SSL - CA") # 46
	SSL_WT_HCOL.res <- do.group("SSL_WT_v_H_COL",	HCOL.mask | (SSL.mask & WT.mask), "SSL - H_COL") # 47
	SSL_WT_HADJSSL.res <- do.group("SSL_WT_v_H_ADJ_SSL",	HADJSSL.mask | (SSL.mask & WT.mask), "SSL - H_ADJ_SSL") # 48
	SSL_WT_HADJHP.res <- do.group("SSL_WT_v_H_ADJ_HP", HADJHP.mask | (SSL.mask & WT.mask), "SSL - H_ADJ_HP") # 49
	SSL_WT_CA.res <- do.group("SSL_WT_v_CA", CA.mask | (SSL.mask & WT.mask), "SSL - CA") # 51
	HP_V600E_HCOL.res <- do.group("HP_V600E_v_H_COL",	HCOL.mask | (HP.mask & V600E.mask), "HP - H_COL") # 52
	HP_V600E_HADJSSL.res <- do.group("HP_V600E_v_H_ADJ_SSL",	HADJSSL.mask | (HP.mask & V600E.mask), "HP - H_ADJ_SSL") # 53
	HP_V600E_HADJHP.res <- do.group("HP_V600E_v_H_ADJ_HP", HADJHP.mask | (HP.mask & V600E.mask), "HP - H_ADJ_HP") # 54
	HP_V600E_CA.res <- do.group("HP_V600E_v_CA", CA.mask | (HP.mask & V600E.mask), "HP - CA") # 56
	HP_WT_HCOL.res <- do.group("HP_WT_v_H_COL",	HCOL.mask | (HP.mask & WT.mask), "HP - H_COL") # 57
	HP_WT_HADJSSL.res <- do.group("HP_WT_v_H_ADJ_SSL",	HADJSSL.mask | (HP.mask & WT.mask), "HP - H_ADJ_SSL") # 58
	HP_WT_HADJHP.res <- do.group("HP_WT_v_HP", HADJHP.mask | (HP.mask & WT.mask), "HP - H_ADJ_HP") # 59
	HP_WT_CA.res <- do.group("HP_WT_v_CA", CA.mask | (HP.mask & WT.mask), "HP - CA") # 61

	# Pooled as H_ADJ_SSL and H_ADJ_HP
	SSL_HADJSSLHP.res <- do.group("SSL_v_H_ADJ_SSL_HP",
		SSL.mask | HADJSSL.mask | HADJHP.mask, "SSL - H_ADJ_SSL_HP",
		group="H_ADJ_PooledGroup") ## 7
	HP_HADJSSLHP.res <- do.group("HP_v_H_ADJ_SSL_HP",
		HP.mask | HADJSSL.mask | HADJHP.mask, "HP - H_ADJ_SSL_HP",
		group="H_ADJ_PooledGroup") ## 13
	CA_HADJSSLHP.res <- do.group("CA_v_H_ADJ_SSL_HP",
		CA.mask | HADJSSL.mask | HADJHP.mask, "CA - H_ADJ_SSL_HP",
		group="H_ADJ_PooledGroup") ## 23
	ANTRUM_HADJSSLHP.res <- do.group("ANTRUM_v_H_ADJ_SSL_HP",
		ANTRUM.mask | HADJSSL.mask | HADJHP.mask, "ANTRUM - H_ADJ_SSL_HP",
		group="H_ADJ_PooledGroup") ## 29
	HADJSSLHP_HCOL.res <- do.group("H_ADJ_SSL_HP_v_H_COL",
		HADJSSL.mask | HADJHP.mask | HCOL.mask, "H_ADJ_SSL_HP - H_COL",
		group="H_ADJ_PooledGroup") ## 40
	HADJSSLHP_CA.res <- do.group("H_ADJ_SSL_HP_v_CA",
		HADJSSL.mask | HADJHP.mask | CA.mask, "H_ADJ_SSL_HP - CA",
		group="H_ADJ_PooledGroup") ## 41
	SSL_V600E_HADJSSLHP.res <- do.group("SSL_V600E_v_H_ADJ_SSL_HP",
		HADJSSL.mask | HADJHP.mask | (SSL.mask & V600E.mask), "SSL - H_ADJ_SSL_HP", group="H_ADJ_PooledGroup") # 45

	SSL_WT_HADJSSLHP.res <- do.group("SSL_WT_v_H_ADJ_SSL_HP",	HADJSSL.mask | HADJHP.mask | (SSL.mask & WT.mask),
		"SSL - H_ADJ_SSL_HP", group="H_ADJ_PooledGroup") # 50
	HP_V600E_HADJSSLHP.res <- do.group("HP_V600E_v_H_ADJ_SSL_HP",	HADJSSL.mask | HADJHP.mask | (HP.mask & V600E.mask),
		"HP - H_ADJ_SSL_HP", group="H_ADJ_PooledGroup") # 55
	HP_WT_HADJSSLHP.res <- do.group("HP_WT_v_H_ADJ_SSL_HP",	HADJSSL.mask | HADJHP.mask | (HP.mask & WT.mask),
		"HP - H_ADJ_SSL_HP", group="H_ADJ_PooledGroup") # 60
	
	# Pooled as SSL and HP
	SSLHP_HCOL.res <- do.group("SSL_HP_v_H_COL", SSL.mask | HP.mask | HCOL.mask,
		"SSL_HP - H_COL", group="SSL_HP_PooledGroup") ## 14
	SSLHP_HADJSSL.res <- do.group("SSL_HP_v_H_ADJ_SSL", SSL.mask | HP.mask | HADJSSL.mask,
		"SSL_HP - H_ADJ_SSL", group="SSL_HP_PooledGroup") ## 16
	SSLHP_CA.res <- do.group("SSL_HP_v_CA", SSL.mask | HP.mask | CA.mask,
		"SSL_HP - CA", group="SSL_HP_PooledGroup") ## 17
	SSLHP_HADJHP.res <- do.group("SSL_HP_v_H_ADJ_HP", SSL.mask | HP.mask | HADJHP.mask,
		"SSL_HP - H_ADJ_HP", group="SSL_HP_PooledGroup") ## 18
	ANTRUM_SSLHP_WT.res <- do.group("ANTRUM_v_SSL_HP_WT",
		ANTRUM.mask | ((SSL.mask | HP.mask) & WT.mask), "ANTRUM - SSL_HP",
		group="SSL_HP_PooledGroup") ## 34
	ANTRUM_SSLHP_V600E.res <- do.group("ANTRUM_v_SSL_HP_V600E",
		ANTRUM.mask | ((SSL.mask | HP.mask) & V600E.mask), "ANTRUM - SSL_HP",
		group="SSL_HP_PooledGroup") ## 35
	ANTRUM_SSLHP.res <- do.group("ANTRUM_v_SSL_HP",
		ANTRUM.mask | SSL.mask | HP.mask, "ANTRUM - SSL_HP",
		group="SSL_HP_PooledGroup") ## 36
	SSLHP_V600E_HCOL.res <- do.group("SSL_HP_V600E_v_H_COL", ((SSL.mask | HP.mask) & V600E.mask) | HCOL.mask, "SSL_HP - H_COL",
		group="SSL_HP_PooledGroup") ## 63
	SSLHP_WT_HCOL.res <- do.group("SSL_HP_WT_v_H_COL", ((SSL.mask | HP.mask) & WT.mask) | HCOL.mask, "SSL_HP - H_COL",
		group="SSL_HP_PooledGroup") ## 64
	SSLHP_V600E_ANTRUM.res <- do.group("SSL_HP_V600E_v_ANTRUM", ((SSL.mask | HP.mask) & V600E.mask) | ANTRUM.mask, "SSL_HP - ANTRUM",
		group="SSL_HP_PooledGroup") ## 63
	SSLHP_V600E_HADJSSL.res <- do.group("SSL_HP_V600E_v_H_ADJ_SSL", ((SSL.mask | HP.mask) & V600E.mask) | HADJSSL.mask, "SSL_HP - H_ADJ_SSL",
		group="SSL_HP_PooledGroup") ## 68
	SSLHP_WT_HADJSSL.res <- do.group("SSL_HP_WT_v_H_ADJ_SSL", ((SSL.mask | HP.mask) & WT.mask) | HADJSSL.mask, "SSL_HP - H_ADJ_SSL",
		group="SSL_HP_PooledGroup") ## 69
	SSLHP_V600E_CA.res <- do.group("SSL_HP_V600E_v_CA", ((SSL.mask | HP.mask) & V600E.mask) | CA.mask, "SSL_HP - CA",
		group="SSL_HP_PooledGroup") ## 70
	SSLHP_WT_CA.res <- do.group("SSL_HP_WT_v_CA", ((SSL.mask | HP.mask) & WT.mask) | CA.mask, "SSL_HP - CA",
		group="SSL_HP_PooledGroup") ## 71
	SSLHP_V600E_HADJHP.res <- do.group("SSL_HP_V600E_v_H_ADJ_HP", ((SSL.mask | HP.mask) & V600E.mask) | HADJHP.mask, "SSL_HP - H_ADJ_HP",
		group="SSL_HP_PooledGroup") ## 72
	SSLHP_WT_HADJHP.res <- do.group("SSL_HP_WT_v_H_ADJ_HP", ((SSL.mask | HP.mask) & WT.mask) | HADJHP.mask, "SSL_HP - H_ADJ_HP",
		group="SSL_HP_PooledGroup") ## 73


	# Both pooled
	SSLHP_HADJSSLHP.res <- do.group("SSL_HP_v_H_ADJ_SSL_HP", SSL.mask | HP.mask | HADJSSL.mask | HADJHP.mask,
		"SSL_HP - H_ADJ_SSL_HP", group="PooledGroup") ## 19
	SSLHP_V600E_HADJSSLHP.res <- do.group("SSL_HP_V600E_v_H_ADJ_SSL_HP",
		((SSL.mask | HP.mask) & V600E.mask) | HADJSSL.mask | HADJHP.mask,
		"SSL_HP - H_ADJ_SSL_HP", group="PooledGroup") ## 66
	SSLHP_WT_HADJSSLHP.res <- do.group("SSL_HP_WT_v_H_ADJ_SSL_HP",
		((SSL.mask | HP.mask) & WT.mask) | HADJSSL.mask | HADJHP.mask,
		"SSL_HP - H_ADJ_SSL_HP", group="PooledGroup") ## 67

	#	Unpooled by BRAF
	SSL_BRAF.res <- do.BRAF("V600E_v_WT_in_SSL", SSL.mask) ## 2
	HP_BRAF.res <- do.BRAF("V600E_v_WT_in_HP", HP.mask) ## 9
	SSLHP_BRAF.res <- do.BRAF("V600E_v_WT_in_SSL_HP", SSL.mask | HP.mask) ## 15
	HADJSSL_HADJHP.res <- do.group("H_ADJ_SSL_v_H_ADJ_HP",
		HADJHP.mask | HADJSSL.mask, "H_ADJ_SSL - H_ADJ_HP") ## 37

	# Paired

	SSLHP_HADJSSLHP_paired.res <- do.paired("SSL_HP_v_paired_H_ADJ_SSL_HP",
		(SSL.mask | HP.mask | HADJSSL.mask | HADJHP.mask) & pairing.mask,
		"SSL_HP - H_ADJ_SSL_HP", group="PooledGroup") ## 74

	# Done!
	cat("Done!\n\n")
}

if(do.extra.de){
	do.ES <- function(model.name, sample.mask, normalise="VSN", group="EscherichiaShigella"){

		cat(paste("*** Analysis contrast:", model.name, "\n"))

		arnaq(model.name=model.name, sample.mask=sample.mask,
			treat.groups=group,
			pca.groups=group,
			normalise=normalise,
			#pca.depth=4,
			svg.export=F
		)

		tmp <- analyse(model.name=model.name, data.matrix=count.data,
			col.data=sample.metadata,
			species="Hs",
			sample.mask=sample.mask,
			analysis.groups=group,
			contrast="PRESENT - ABSENT"
		)

		tmp.vol <- rna.volcano(tmp$results[[1]]$table)
		ggsave(paste0("results/volcano/CGG014_", model.name, "_volcano.png"), tmp.vol)

		return(tmp)

	}

	ES.mask <- sample.metadata$EscherichiaShigella != "NOT_TESTED"
	Fuso.mask <- sample.metadata$Fusobacterium != "NOT_TESTED"
	Serr.mask <- sample.metadata$Serratia != "NOT_TESTED"

	# EscherichiaShigella	
	ES_SSLHP.res <- do.ES("ES_SSL_HP", ES.mask & (SSL.mask | HP.mask), normalise="VSN")
	ES_SSLHP_WT.res <- do.ES("ES_SSL_HP_WT", WT.mask & ES.mask & (SSL.mask | HP.mask), normalise="VSN")
	ES_SSLHP_V600E.res <- do.ES("ES_SSL_HP_V600E", V600E.mask & ES.mask & (SSL.mask | HP.mask), normalise="VSN")
	ES_H_ADJ_SSLHP.res <- do.ES("ES_H_ADJ_SSL_HP", ES.mask & (HADJSSL.mask | HADJHP.mask), normalise="VSN")
	ES_SSLHP_and_H_ADJ_SSL_HP.res <- do.ES("ES_SSL_HP_and_H_ADJ_SSL_HP", ES.mask & (HADJSSL.mask | HADJHP.mask | SSL.mask | HP.mask), normalise="VSN")
	# Fusobacterium	
	Fuso_SSLHP.res <- do.ES("Fuso_SSL_HP", Fuso.mask & (SSL.mask | HP.mask), normalise="VSN", "Fusobacterium")
	Fuso_H_ADJ_SSLHP.res <- do.ES("Fuso_H_ADJ_SSL_HP", Fuso.mask & (HADJSSL.mask | HADJHP.mask), normalise="VSN", "Fusobacterium")
	Fuso_SSLHP_and_H_ADJ_SSL_HP.res <- do.ES("Fuso_SSL_HP_and_H_ADJ_SSL_HP", Fuso.mask & (HADJSSL.mask | HADJHP.mask | SSL.mask | HP.mask), normalise="VSN", "Fusobacterium")
	# Serratia
	Serr_SSLHP.res <- do.ES("Serr_SSL_HP", Serr.mask & (SSL.mask | HP.mask), normalise="VSN", "Serratia")
	Serr_H_ADJ_SSLHP.res <- do.ES("Serr_H_ADJ_SSL_HP", Serr.mask & (HADJSSL.mask | HADJHP.mask), normalise="VSN", "Serratia")
	Serr_SSLHP_and_H_ADJ_SSL_HP.res <- do.ES("Serr_SSL_HP_and_H_ADJ_SSL_HP", Serr.mask & (HADJSSL.mask | HADJHP.mask | SSL.mask | HP.mask), normalise="VSN", "Serratia")

}

if(do.extra.de.2){

	do.ES.2 <- function(model.name, sample.mask, normalise="VSN",
											groups=c("BRAF", "EscherichiaShigella"),
											contrast){

		cat(paste("*** Analysis contrast:", model.name, "\n"))

		arnaq(model.name=model.name, sample.mask=sample.mask,
			treat.groups=groups,
			pca.groups=groups,
			normalise=normalise,
			#pca.depth=4,
			svg.export=F
		)

		tmp <- analyse(model.name=model.name, data.matrix=count.data,
			col.data=sample.metadata,
			species="Hs",
			sample.mask=sample.mask,
			analysis.groups=groups,
			contrast=contrast
		)

		tmp.vol <- rna.volcano(tmp$results[[1]]$table)
		ggsave(paste0("results/volcano/CGG014_", model.name, "_volcano.png"), tmp.vol)

		return(tmp)

	}

	SSLHP_V600E_ES <- do.ES.2("SSL_HP_V600E_ES_PRES_v_SSL_HP_WT_ES_ABS",
														(SSL.mask | HP.mask) & ES.mask, 
														contrast="V600E_PRESENT - WT_ABSENT")
	SSLHP_V600E_ES_2 <- do.ES.2("SSL_HP_V600E_ES_PRES_v_SSL_HP_WT_ES_PRES",
														(SSL.mask | HP.mask) & ES.mask, 
														contrast="V600E_PRESENT - WT_PRESENT")
	SSLHP_V600E_ES_3 <- do.ES.2("SSL_HP_V600E_ES_ABS_v_SSL_HP_WT_ES_ABS",
														(SSL.mask | HP.mask) & ES.mask, 
														contrast="V600E_ABSENT - WT_ABSENT")
	
}

if(do.gene.list){
	library(arnaq)
	library(dendextend)
	# this mutilates the stored objects, run first step again if you need the full
	# count table after this point

	species.gtf <<- stored.species.gtf
	count.data <<- full.count.data

	custom.gene.list <- read.table("custom_gene_list_2.txt", header=F, stringsAsFactor=F)[,1]
	custom.gene.list.ensembl <- species.gtf$gene_id[
																match(custom.gene.list, species.gtf$gene_name)]

	make_extra_plots <- function(mask, name_frag) {
		cat("Start\n")
		count.data2 <- norm.data
		sample.metadata2 <- sample.metadata[mask,]

		sample.metadata2$Group <- factor(sample.metadata2$Group,
				levels=c("ANTRUM", "SSL", "H_ADJ_SSL_HP", "H_COL", "HP", "CA"))
		count.data2 <- count.data2[,order(sample.metadata2$Group)]
		sample.metadata2 <- sample.metadata2[order(sample.metadata2$Group),]
		cat("Done re-ordering\n")

#		if(!identical(colnames(count.data2),sample.metadata2$Display)){
#			cat(colnames(count.data2))
#			cat("\n")
#			cat(sample.metadata2$Display)
#			cat("\n")
#			stop("WAT")
#		}

		tmp2 <- plot_all_clustering(count.data2)
		tmp2[[1]] <- tmp2[[1]] + theme(axis.text=element_text(size=5,
			colour=sample.metadata2$Group[order(sample.metadata2$Group)]))
		tmp2[[3]] <- tmp2[[3]] + theme(axis.text=element_text(size=5,
			colour=sample.metadata2$Group[order(sample.metadata2$Group)]))

		dend_palette <- c("#000000","#ff0000","#00ff00","#0000ff","#cccc00",
			"#cc00cc","#00cccc","#7f7f7f")
		cat("Done heatmaps\n")

		dend.helper <- function(data, title, ylim=-0.06){
			out <- as.dendrogram(data) %>% set("labels_cex", 0.5)
			out <- out %>% set("labels_colors",
				dend_palette[as.numeric(sample.metadata2$Group[match(labels(out),
					sample.metadata2$Display)])])
			tmp <- get_branches_heights(out)
			cat(paste("Minimum height:",min(tmp),"\n"))
			cat(paste("Maximum height:",max(tmp),"\n"))
			out <- ggplot(as.ggdend(out)) + ylim(ylim, max(get_branches_heights(out))) +
					 labs(title=title)
			return(out)
		}

		tmp2[[2]] <- dend.helper(hclust(dist(t(count.data2))),
			"Distance-based hierarchical clustering")
		tmp2[[4]] <- dend.helper(hclust(as.dist(1-cor(count.data2, 
			method="spearman")), method="ward.D2"), "Spearman's rank-based hierarchical clustering",
			-0.3)
		cat("Done dendrograms\n")

		ggsave(paste0("QC/svg/CGG014_Maier_custom_geneset_",name_frag,"_dist_heatmap.png"),
			tmp2[[1]], width=12, height=12)
		ggsave(paste0("QC/svg/CGG014_Maier_custom_geneset_",name_frag,"_dist_clustering.png"),
			tmp2[[2]], width=12, height=12)
		ggsave(paste0("QC/svg/CGG014_Maier_custom_geneset_",name_frag,"_spearmans_heatmap.png"),
			tmp2[[3]], width=12, height=12)
		ggsave(paste0("QC/svg/CGG014_Maier_custom_geneset_",name_frag,"_spearmans_clustering.png"),
			tmp2[[4]], width=12, height=12)
		ggsave(paste0("QC/svg/CGG014_Maier_custom_geneset_",name_frag,"_PCA.png"),
			all_plots$pca_group_pss[[1]][[1]], width=6, height=6)

	}

	full.count.data <- count.data
	count.data <- count.data[rownames(count.data) %in% custom.gene.list.ensembl,]
	stored.species.gtf <- species.gtf

#	arnaq(resources.file="resources_nobio.yml",
#		model.name="custom_gene_list",
#		treat.groups="Group",
#		sample.mask=outlier_mask,
#		normalise="linear")

#	make_extra_plots(outlier_mask, "all")

	outlier_mask2 <- make_outlier_mask(c("Undetermined","B2_SSL_SIG_9485"))
	sample_mask <- outlier_mask2 & sample.metadata$Group %in% c("Antrum","H_Adj","H_Col","SSL")
	
	arnaq(resources.file="resources_nobio.yml",
		model.name="custom_gene_list_2",
		treat.groups="Group",
		sample.mask=sample_mask,
		normalise="linear")

	make_extra_plots(sample_mask, "prog")


	# Recalc PCA to get loadings
	tmp.pca.loadings <- prcomp(t(norm.data))$rotation
	tmp.pca.loadings <- cbind(ENSEMBL=row.names(tmp.pca.loadings),
			GENE_NAME=stored.species.gtf[
			match(row.names(tmp.pca.loadings), stored.species.gtf$gene_id),"gene_name"],
			tmp.pca.loadings)
	write.table(tmp.pca.loadings, "QC/CGG014_custom_gene_list_2_pca_loadings.txt",
            sep="\t", quote=F, row.names=F)	

	for(PC in c("PC1","PC2","PC3","PC4")){
		
	}
	
	species.gtf <- stored.species.gtf
	count.data <- full.count.data
}
