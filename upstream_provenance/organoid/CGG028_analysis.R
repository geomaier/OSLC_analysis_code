#########################
# Secondary QC and analysis for CGG028
# George Maier's project
#########################

library(arnaq)
library(reshape2)
library(cowplot)
library(ggvenn)

# Requires installation from github
library(ggsankey)

source("edgeR_Analysis.R")

do.prelim <- F
do.infection <- F
do.alt.infection <- F
do.rest.arnaq <- F
do.first.passage.BRAF <- F
do.multipassage.148 <- T
do.passage.plot <- F
do.rest.de <- F
do.rest.venn <- F
do.rest.sankey <- F
do.initial.p.de <- F

if(do.prelim){
	arnaq(normalise="linear", treat.groups=c("BRAF", "Passage"),
		pca.groups=c("Type", "Passage", "Location", "Gender", "BRAF", "Infected", "Patient"))
}	

if(do.infection){
	infected.mask <- !make_outlier_mask(c("GM_GI12181_P4_INF", "GM_GI12181_P4_NONINF",
				"GM_GI12181_P5_1_INF", "GM_GI12181_P5_1_NONINF", "GM_GI12181_P5_2_INF",
				"GM_GI12181_P5_2_NONINF", "GM_GI12181_P6_INF", "GM_GI12181_P6_NONINF"))

	arnaq(treat.groups="Infected", model.name="Infected", sample.mask=infected.mask)

	infected.res <- analyse(model.name="Infected", sample.mask=infected.mask,
													analysis.groups="Infected", contrast="Infected - NonInf")

	infected.cP <- analyse.cP(infected.res$results[[1]]$table, "CGG028_Infected",
														fformat=".png", plot.number=15)
}

if(do.alt.infection){

	save.counts <- function(data, filename) {
		cat(paste("Writing count data to table:", filename, "\n"))
		utils::write.table(data, filename, sep = "\t", quote = FALSE, row.names = TRUE)
	}

	alt_infected.mask <- !make_outlier_mask(c("GM_GI12181_P4_INF", "GM_GI12181_P4_NONINF",
				"GM_GI12181_P6_INF", "GM_GI12181_P6_NONINF"))
	alt_infected.mask <- c(alt_infected.mask, T, T)

	arnaq(treat.groups="Infected", model.name="Alt_Infected", sample.mask=infected.mask)

	alt_count.data <- cbind(count.data,
		GM_GI12181_P5_C_INF=rowSums(count.data[,c("GM_GI12181_P5_1_INF", "GM_GI12181_P5_2_INF")]),
		GM_GI12181_P5_C_NONINF=rowSums(count.data[,c("GM_GI12181_P5_1_NONINF", "GM_GI12181_P5_2_NONINF")]))
	alt_sample.metadata <- rbind(sample.metadata,
													GM_GI12181_P5_C_INF=sample.metadata["GM_GI12181_P5_1_INF",],
													GM_GI12181_P5_C_NONINF=sample.metadata["GM_GI12181_P5_1_NONINF",])
	alt_sample.metadata["GM_GI12181_P5_C_INF","Display"] <- "GM_GI12181_P5_C_INF"
	alt_sample.metadata["GM_GI12181_P5_C_NONINF","Display"] <- "GM_GI12181_P5_C_NONINF"

	tmp.count.data <- alt_count.data[gene.masks[["Genes"]], alt_infected.mask, drop=FALSE]
	alt.cpm <- t(t(tmp.count.data) * 1000000 / colSums(tmp.count.data))
	save.counts(alt.cpm, "QC/CGG028_Alt_Infected_cpm_table.txt")

	alt_infected.res <- analyse(data.matrix=alt_count.data, col.data=alt_sample.metadata,
													model.name="Alt_Infected", sample.mask=alt_infected.mask,
													analysis.groups="Infected", contrast="Infected - NonInf")

	alt_infected.cP <- analyse.cP(alt_infected.res$results[[1]]$table,
														"CGG028_alt_Infected",
														fformat=".png", plot.number=15)

}

if(do.rest.arnaq) {
	arnaq(treat.groups=c("BRAF", "Passage"), model.name="Rest",
					pca.groups=c("BRAF", "Type", "Passage", "Location", "Patient"),
				sample.mask=!infected.mask, svg.export=T)
}

make_arrow_pca <- function(old_plot, sample.mask,
													 new_palette=c("#FF0000", "#000000")){
	tmp.data <- data.frame(PC1=old_plot@data$PC1, PC2=old_plot@data$PC2, BRAF=old_plot@data$BRAF)
  tmp.data$Passage <- sample.metadata$Passage[sample.mask]
  tmp.data$Path <- paste(sample.metadata$Patient[sample.mask],
                         sample.metadata$Type[sample.mask],
                         sample.metadata$Location[sample.mask])
  tmp.data <- cbind(tmp.data[-nrow(tmp.data),], tmp.data[-1, -3])
  colnames(tmp.data)[6:9] <- c("PC1e", "PC2e", "Passage2", "Path2")
	tmp.data <- tmp.data[tmp.data$Path == tmp.data$Path2,]

	newplot <- old_plot + 
								geom_point(size=2.5) +
								scale_colour_discrete(palette=new_palette) +
								ggnewscale::new_scale_colour() +
								geom_segment(data=tmp.data, alpha=0.5,
                             mapping=aes(group=Path, colour=Path, xend=PC1e, yend=PC2e),
                             arrow=grid::arrow(type="closed", length=unit(0.1, "inches"))) +
								scale_colour_discrete(palette=c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78",
																								"#2ca02c", "#98df8a", "#d62728", "#ff9896",
																								"#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
																								"#e377c2", "#f7b6d2", "#bcbd22", "#dbdb8d"))
	# Secondary plot showing average segment directions?

	pca.drift.all <- data.frame(PC1=0, PC2=0)
	pca.drift.V600E <- data.frame(PC1=0, PC2=0)
	pca.drift.WT <-data.frame(PC1=0, PC2=0)
	for(time in 1:8){
		tmp.data.2 <- tmp.data[tmp.data$Passage==paste0("P", time) & 
													 tmp.data$Passage2==paste0("P", time+1),]
		tmp.vec <- data.frame(PC1=tmp.data.2$PC1e-tmp.data.2$PC1,
													PC2=tmp.data.2$PC2e-tmp.data.2$PC2)
		pca.drift.all <- rbind(pca.drift.all, colMeans(tmp.vec))
		pca.drift.V600E <- rbind(pca.drift.V600E, colMeans(tmp.vec[tmp.data.2$BRAF=="V600E",]))
		pca.drift.WT <- rbind(pca.drift.WT, colMeans(tmp.vec[tmp.data.2$BRAF=="WT",]))
	}
	pca.drift.all <- cbind(pca.drift.all[-1,], Passage=paste0("P",2:9), Group="All")
	pca.drift.V600E <- cbind(pca.drift.V600E[-1,], Passage=paste0("P",2:9), Group="V600E")
	pca.drift.WT <- cbind(pca.drift.WT[-1,], Passage=paste0("P",2:9), Group="WT")
	pca.drift <- rbind(pca.drift.all, pca.drift.V600E, pca.drift.WT)
	pca.drift <- cbind(pca.drift, X=0, Y=0)

	pca.drift <- pca.drift[pca.drift$Group !="All",]

	drift_plot <- ggplot(pca.drift, aes(x=PC1, y=PC2, colour=Group)) +
								geom_point(size=2.5) +
								scale_colour_discrete(palette=new_palette) +
								ggnewscale::new_scale_colour() +
								geom_segment(alpha=0.5,
                             mapping=aes(colour=Passage, x=X, y=Y, xend=PC1, yend=PC2),
                             arrow=grid::arrow(type="closed", length=unit(0.1, "inches"))) +
                scale_colour_discrete(palette=c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78",
                                                "#2ca02c", "#98df8a", "#d62728", "#ff9896",
                                                "#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
                                                "#e377c2", "#f7b6d2", "#bcbd22", "#dbdb8d")) +
								ggtitle("Average change in PCA position over time")


	return(list(newplot, drift_plot))
}

if(do.first.passage.BRAF) {
	first_passage.mask <- sample.metadata$Passage=="P1" | 
						(sample.metadata$Passage=="P2" & sample.metadata$Patient=="Pt_12147")
	SSL.mask <- sample.metadata$Type=="SSL"
	HP.mask <- sample.metadata$Type=="HP"

	first_passage_BRAF.res <- analyse(data.matrix=count.data, col.data=sample.metadata,
                                    model.name="First_Passage",
																		sample.mask=first_passage.mask & (SSL.mask | HP.mask),
                                    analysis.groups="BRAF", contrast="V600E - WT")
  first_passage_BRAF.cP <- analyse.cP(first_passage_BRAF.res$results[[1]]$table,
																			"CGG028_first_passage_BRAF",
				                              fformat=".png", plot.number=15)
	
}

if(do.multipassage.148){
	sample.metadata$Type2 <- ifelse(sample.metadata$Type=="H_ADJ", "H_ADJ", "SSL_HP")
	fourth_passage.mask <- sample.metadata$Passage=="P4" & !infected.mask
	eighth_passage.mask <- sample.metadata$Passage=="P8"

	first_passage.res <- analyse(data.matrix=count.data, col.data=sample.metadata,
															 model.name="First_Passage_SSL_HP",
															 sample.mask=first_passage.mask & (SSL_HP.mask | H_ADJ.mask),
															 analysis.groups="Type2", contrast="SSL_HP - H_ADJ",
															 extra.groups="Patient")
	first_passage.cP <- analyse.cP(first_passage.res$results[[1]]$table,
																 "CGG028_first_passage_SSL_HP",
																 fformat=".png", plot.number=15)
	fourth_passage.res <- analyse(data.matrix=count.data, col.data=sample.metadata,
															 model.name="Fourth_Passage_SSL_HP",
															 sample.mask=fourth_passage.mask & (SSL_HP.mask | H_ADJ.mask),
															 analysis.groups="Type2", contrast="SSL_HP - H_ADJ",
															 extra.groups="Patient")
	fourth_passage.cP <- analyse.cP(fourth_passage.res$results[[1]]$table,
																 "CGG028_fourth_passage_SSL_HP",
																 fformat=".png", plot.number=15)
	eighth_passage.res <- analyse(data.matrix=count.data, col.data=sample.metadata,
															 model.name="Eighth_Passage_SSL_HP",
															 sample.mask=eighth_passage.mask & (SSL_HP.mask | H_ADJ.mask),
															 analysis.groups="Type2", contrast="SSL_HP - H_ADJ",
															 extra.groups="Patient")
	eighth_passage.cP <- analyse.cP(eighth_passage.res$results[[1]]$table,
																 "CGG028_eighth_passage_SSL_HP",
																 fformat=".png", plot.number=15)
}

if(do.passage.plot){
	tmp <- make_arrow_pca(all_plots$pca_group_pss$BRAF[[1]], !infected.mask)
	arrow.plot <- tmp[[1]]
	ggsave("CGG028_pca_BRAF_arrows.png", arrow.plot)
	drift.plot <- tmp[[2]]
	ggsave("CGG028_pca_BRAF_drift.png", drift.plot)
}

if(do.rest.de){

	H_ADJ.mask <- sample.metadata$Type=="H_ADJ"
	SSL_HP.mask <- sample.metadata$Type %in% c("SSL", "HP")
	contrast.names <- c("P2-P1", "P3-P1", "P4-P1", "P5-P1", "P6-P1", "P7-P1", "P8-P1",
										 "P9-P1", "P3-P2", "P4-P2", "P5-P2", "P6-P2", "P7-P2", "P8-P2",
										 "P9-P2")

	H_ADJ.res <- analyse(model.name="Passage_H_ADJ", sample.mask=H_ADJ.mask, analysis.groups=c("Passage"),
											 contrast=contrast.names[c(-8,-15)],
											 extra.groups=c("Patient"))
	SSL_HP.res <- analyse(model.name="Passage_SSL_HP", sample.mask=SSL_HP.mask, analysis.groups=c("Passage"),
											 contrast=contrast.names,
											 extra.groups=c("Patient"))
}

if(do.rest.venn){
	# make DE gene table

	H_ADJ.de <- data.frame(do.call(cbind, lapply(H_ADJ.res$results, function(x){x$de})))
	SSL_HP.de <- data.frame(do.call(cbind, lapply(SSL_HP.res$results, function(x){x$de})))
	names(H_ADJ.de) <- contrast.names[c(-8, -15)]
	names(SSL_HP.de) <- contrast.names
	tmp <- melt(H_ADJ.de)
	tmp2 <- melt(SSL_HP.de)
	tmp <- cbind(rownames(H_ADJ.res$results[[1]]$de), tmp)
	tmp2 <- cbind(rownames(SSL_HP.res$results[[1]]$de), tmp2)
	names(tmp) <- names(tmp2) <- c("Gene", "Contrast", "DE")
	tmp <- cbind(tmp, Type="H_ADJ")
	tmp2 <- cbind(tmp2, Type="SSL_HP")
	ALL.de <- rbind(tmp, tmp2)
	
	# Now make venn diagrams

	make.de.venn <- function(melted.dat, contrast){
		tmp <- melted.dat[melted.dat$Contrast == contrast,]
		tmp <- list(
								`H_ADJ Down`=tmp[tmp$DE==-1 & tmp$Type=="H_ADJ", "Gene"],
								`H_ADJ Up`=tmp[tmp$DE==1 & tmp$Type=="H_ADJ", "Gene"],
								`SSL_HP Up`=tmp[tmp$DE==1 & tmp$Type=="SSL_HP", "Gene"],
								`SSL_HP Down`=tmp[tmp$DE==-1 & tmp$Type=="SSL_HP", "Gene"])
		ggvenn(tmp) +
			ggtitle(paste("Significant DE genes in contrast", contrast))
	}

	ALL.de.plot <- lapply(contrast.names, function(x){make.de.venn(ALL.de, x)})
	names(ALL.de.plot) <- contrast.names
	
	save_plot("venn/CGG028_de_venn_1.png", plot_grid(ALL.de.plot[[1]], ALL.de.plot[[2]],
																									 ALL.de.plot[[3]], ALL.de.plot[[4]]),
						base_height=10)
	save_plot("venn/CGG028_de_venn_2.png", plot_grid(ALL.de.plot[[5]], ALL.de.plot[[6]],
																									 ALL.de.plot[[7]], ALL.de.plot[[8]]),
						base_height=10)
	save_plot("venn/CGG028_de_venn_3.png", plot_grid(ALL.de.plot[[9]], ALL.de.plot[[10]],
																									 ALL.de.plot[[11]], ALL.de.plot[[12]]),
						base_height=10)
	save_plot("venn/CGG028_de_venn_4.png", plot_grid(ALL.de.plot[[13]], ALL.de.plot[[14]],
																									 ALL.de.plot[[15]]),
						base_height=10)

}

if(do.rest.sankey) {

	make.de.sankey <- function(dat, type, baseline) {
		ggplot(dat, aes(x=x, next_x=next_x, node = node, next_node=next_node,
										fill=factor(node), label=node)) +
			geom_sankey(flow.alpha = 0.6, node.color ="gray30") +
			geom_sankey_label( color="white", fill="gray70") +
			scale_fill_viridis_d() +
			theme_sankey(base_size=18) +
			labs(x=NULL) +
			theme(legend.position= "none", plot.title=element_text(hjust=0.5)) +
			ggtitle(paste("DE Genes in", type, ", baseline ", baseline))

	}

	tmp <- data.frame(H_ADJ.de[,colnames(H_ADJ.de) %in% contrast.names[1:7]])
	tmp[tmp==1] <- "U"
	tmp[tmp==-1] <- "D"
	tmp[tmp==0] <- "NS"
	names(tmp) <- gsub("\\..*", "", names(tmp))
	tmp <- make_long(tmp, P2, P3, P4, P5, P6, P7, P8)
	H_ADJ.de.sankey1 <- make.de.sankey(tmp, "H_ADJ", "P1")
	ggsave("sankey/CGG028_de_sankey_H_ADJ_1.png", H_ADJ.de.sankey1)

	tmp <- data.frame(SSL_HP.de[,colnames(SSL_HP.de) %in% contrast.names[1:8]])
	tmp[tmp==1] <- "U"
	tmp[tmp==-1] <- "D"
	tmp[tmp==0] <- "NS"
	names(tmp) <- gsub("\\..*", "", names(tmp))
	tmp <- make_long(tmp, P2, P3, P4, P5, P6, P7, P8, P9)
	SSL_HP.de.sankey1 <- make.de.sankey(tmp, "SSL_HP", "P1")
	ggsave("sankey/CGG028_de_sankey_SSL_HP_1.png", SSL_HP.de.sankey1)

	tmp <- data.frame(H_ADJ.de[,colnames(H_ADJ.de) %in% contrast.names[9:14]])
	tmp[tmp==1] <- "U"
	tmp[tmp==-1] <- "D"
	tmp[tmp==0] <- "NS"
	names(tmp) <- gsub("\\..*", "", names(tmp))
	tmp <- make_long(tmp, P3, P4, P5, P6, P7, P8)
	H_ADJ.de.sankey2 <- make.de.sankey(tmp, "H_ADJ", "P2")
	ggsave("sankey/CGG028_de_sankey_H_ADJ_2.png", H_ADJ.de.sankey2)

	tmp <- data.frame(SSL_HP.de[,colnames(SSL_HP.de) %in% contrast.names[9:15]])
	tmp[tmp==1] <- "U"
	tmp[tmp==-1] <- "D"
	tmp[tmp==0] <- "NS"
	names(tmp) <- gsub("\\..*", "", names(tmp))
	tmp <- make_long(tmp, P3, P4, P5, P6, P7, P8, P9)
	SSL_HP.de.sankey2 <- make.de.sankey(tmp, "SSL_HP", "P2")
	ggsave("sankey/CGG028_de_sankey_SSL_HP_2.png", SSL_HP.de.sankey2)

}

if(do.initial.p.de){

	initial.mask <- sample.metadata$Passage=="P1" | (sample.metadata$Passage=="P2" & sample.metadata$Patient=="Pt_12147")

	initial.BRAF.all <- analyse(model.name="Initial_BRAF_all", sample.mask=initial.mask, analysis.groups="BRAF",
															contrast="V600E - WT", extra.groups="Patient")
	initial.BRAF.SSL_HP <- analyse(model.name="Initial_BRAF_SSL_HP", sample.mask=initial.mask & SSL_HP.mask, analysis.groups="BRAF",
															contrast="V600E - WT")
}
