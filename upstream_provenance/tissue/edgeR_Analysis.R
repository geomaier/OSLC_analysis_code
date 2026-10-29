###################
# Uses edgeR for analysis of REX pipeline RNASeq
# Author: bjw78
# Incorporating functions by esh and hel23
###################

# Conda environment: 'arnaq_1.0_analysis.R'
# Based on the current arnaq environment but with extra packages for analysis and
# annotation: edgeR, org.Mm.eg.db and org.Hs.eg.db, stringr
# The latter two are large, so this environment can't sit in user space.
# Available on /well/bsg/conda_envs/arnaq_1.0_analysis.R

# TODO
# replace heatmap with improved heatmap
# See if other columns from the GTF will be useful
# Test!

###################
# Runs the analysis for DE
###################
analyse <- function(proj.name=NULL, model.name=NULL, data.matrix=count.data,
                    col.data=sample.metadata, species="Hs", gene.mask=gene.masks$Genes,
                    sample.mask=NULL, analysis.groups, contrast, extra.groups=c()) {

  # proj.name : Name of projects used for output files. Default of NULL means to look for
  #             an object called project.id, and then for arnaq.run$project.id.
  # model.name : name for this model where multiple models are used for the same project.
  #              Default of NULL means to look for an object called model.name, and then for
  #              arnaq.run$model.name, and then if neither if found use no model.name.
  # data.matrix : matrix of raw counts
  # col.data : Metadata for samples
  # species : Used to annotate the DE results. If "Hs" or "Mm", will use org.XX.eg.db to annotate.
  #           If "GTF", will look for the 'species.gtf' object creatred by ARNAQ and use it to
  #           annotate with gene name only. If a full path to a gtf file, will use that gtf to
  #           annotate with gene name only, and will also store it as 'species.gtf' so future
  #           analyses do not have to load the file again. If "None", perform no annotation.
  # gene.mask : Which genes to use in the analysis. Default value selects
  #             non-zero genes.
  # sample.mask : mask used to select samples for this model. Default is all samples
  # analysis.groups : the name of the column holding the group information.
  #                   If this is a vector, then a new column will be
  #                   appended to the sample.metadata object using concatenation.
  #	                  This allows multi-factor models.
  # contrast : the contrast to test in this group, as a string.
  #            If this is a vector, each contrast will generate its own output.
  #            Multifactor models will need to use the concatenated group
  #            names (using "_").
  # extra.groups : the names of groups that should be added to the model
  #			           but that are considered continuous variables and hence do not
  #				         get included in any concatenation or factor-based analysis.
  #                You can still use indicator variables if desired, they will
  #                just be treated as a simple linear term.

  # Find proj.id
  if (is.null(proj.name)) {
    if ("project.id" %in% ls(1)) {
      proj.name <- project.id
    } else if ("arnaq.run" %in% ls(1)) {
      proj.name <- arnaq.run$project.id
    } else {
      stop("Please give a project id")
    }
  }

  # Find model.name
  if (is.null(model.name)) {
    if ("model.name" %in% ls(1)) {
      model.name <- model.name # Probably redundant but it also makes sure this over-ride arnaq.run
    } else if ("arnaq.run" %in% ls(1)) {
      model.name <- arnaq.run$model.name
    } # Else leave as NULL
  }

  suppressMessages(library(edgeR))

  cat("Starting analysis\n")

  # local copy of count matrix
  d <- data.matrix
  colnames(d) <- col.data$Name
  # Filter to gene list
  d <- d[gene.mask, ]
  # Make sample filter
  if (is.null(sample.mask)) {
    sample.mask <- rep(TRUE, ncol(d))
  }

  cat(paste("Project ID:", proj.name, "\n"))
  cat(paste("Model name:", model.name, "\n"))
  cat(paste("Genes to analyse:", sum(gene.mask), "\n"))
  cat(paste("Samples to analyse:", sum(sample.mask), "\n"))
  cat(paste("Groups with factor terms:", paste(analysis.groups, collapse=", "), "\n"))
  cat(paste("Groups with linear terms:", paste(extra.groups, collapse=", "), "\n"))

  # Create new combined column if necessary
  if (length(analysis.groups)>1) {
    cat("Creating new sample.metadata for compound model\n")
    tmp.col <- col.data[, analysis.groups]
    tmp.col <- apply(tmp.col, 1, function(x) {
                                              paste(x, collapse="_")})
    col.data <- cbind(col.data, tmp.col)
    analysis.groups <- paste(analysis.groups, collapse="_")
    colnames(col.data)[ncol(col.data)] <- analysis.groups
  }

  # Make DGE object
  cat("Making DGE object\n")
  data.y <- prep.DGE(d[, sample.mask], col.data[sample.mask, analysis.groups])

  # Make design object
  cat("Making model design\n")
  model.design <- make.design(data.y$samples, col.data[sample.mask, ], analysis.groups,
                              extra.groups)

  # Generate dispersions
  cat("Calculating dispersions\n")
  data.y <- make.dispersions(data.y, model.design)

  # Create model contrasts
  cat("Making contrasts\n")
  model.contrasts <- makeContrasts(contrasts=contrast, levels=colnames(model.design))

  # Generate results
  cat("Fitting glm model\n")
  model.fit <- glmFit(data.y, model.design)

  # Create filename fragment for output
  base.output.name <- proj.name
  if (!is.null(model.name)) {
    base.output.name <- paste(base.output.name, model.name, sep="_")
  }

  # Handle multiple contrasts
  num.contrasts <- ncol(model.contrasts)
  result.design <- list(y=data.y, sample.mask=sample.mask, gene.mask=gene.mask,
                        design=model.design, fit=model.fit)
  result <- list(design=result.design, results=list())

  # Make sure results directory exists
  out.dir <- "results/"
  if (!dir.exists(out.dir)) {
    dir.create(out.dir, recursive=TRUE)
  }

  # run tests for each contrast
  for (i in 1:num.contrasts){
    cat(paste("Processing contrast:", contrast[i], "\n"))
    model.res <- glmLRT(model.fit, contrast=model.contrasts[, i])
    model.de <- decideTests.DGEExact(model.res, adjust.method="fdr", p.value=0.05)
    model.fdr <- topTags(model.res, sort.by="PValue", n=Inf)

    output.name <- base.output.name
    if (num.contrasts>1) {
      output.name <- paste0(base.output.name, "_", i)
    }

    # Annotate
    model.out <- annotate.table(model.fdr$table, species)

    # Write table
    filename <- paste0(out.dir, output.name, "_DE_results.txt")
    cat(paste("    Writing out table:", filename, "\n"))
    write.table(model.out, filename, sep="\t", quote=FALSE, row.names=FALSE)

    # Return list of objects
    result$results[[i]] <- list(table=model.out, de=model.de, summary=summary(model.de))
  }

  # Plots
  cat("Making summary plots\n")
  pdf(paste0(base.output.name, "_plots.pdf"), height=11, width=8)
  plotMDS(data.y, main=paste(output.name, "MDS"), cex.main=0.7, cex=0.7)
  plotBCV(data.y, main=paste(output.name, "BCV"), cex.main=0.7, cex=0.7)
  dev.off()

  cat("Finished\n")

  # Return list
  return(result)

}

#####################################
# Create violin/boxplot with jitter in GGplot
# Added by Matthieu Miossec 09/03/2022
# Assumes ggplot2 is loaded.
#####################################

make.boxplots <- function(dis) {
  ##takes as its input the object create by analyse function

  library(stringr)
  suppressPackageStartupMessages(library(ggplot2))

  countsPM <- as.matrix(cpm(dis$design$y$counts)) ##do we have CPM already stored somewhere?
  results <- dis$results[[1]]$table
  conds <- str_replace_all(dis$design$fit$samples$group, c("G0"="Control", "G1"="AF"))

  pdf("plots/sigFDR_allplot.pdf", onefile=TRUE)
  # use a loop to iterate over X genes with highest p-value
  options(scipen = 1)
  options(digits = 2)
  i <- 1
  while (results$FDR[i]<0.05) {
    df <- cbind.data.frame(countsPM[match(results$genes[i], row.names(countsPM)), ], conds)
    colnames(df) <- c("counts", "conds")
    p <- ggplot(df, aes(x=conds, y=counts)) +
      geom_violin() +
      geom_boxplot(outlier.shape = NA, width = 0.07) +
      geom_jitter(color="black", size=0.4, alpha=0.9) +
      ggtitle(paste0(results$external_gene_name[i], "  (LogFC: ", formatC(results$logFC[i]),
                     ", P-value: ", formatC(results$PValue[i], digits=3, format="e"), " FDR: ",
                     formatC(results$FDR[i], digits=3, format="e"), ")")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      xlab("Condition") + ylab("Counts per Million")
    print(p) #with boxplot/violinplot/jitter yay!
    i <- i + 1
  }
  dev.off()
}

###################################
# Create heatmap from results set
###################################

rna.heatmap <- function(norm.data, species="human", genes=rownames(norm.data),
                        samples=colnames(norm.data), plot.title="") {
  suppressPackageStartupMessages(library(ggplot2))
  library(reshape2)
  filtered.data <- norm.data[genes, samples]
  tmp <- annotate.table(filtered.data, species)
  rownames(filtered.data) <- tmp$external_gene_name
  filtered.data <- melt(filtered.data)
  colnames(filtered.data) <- c("Gene", "Sample", "Expression")
  out <- ggplot(filtered.data, aes(Gene, Sample, fill=Expression)) + geom_tile()
  out
}

###########################
# Helper functions
###########################

####################
# Prepare DGE object, including group data
prep.DGE <- function(masked.data, group.data) {

  y <- DGEList(counts=masked.data, genes=row.names(masked.data))

  #Add in group information here
  y$samples$group <- group.data

  # Use edgeR's built-in filter, with defaults
  keep <- filterByExpr(y)

  y <- y[keep, ]

  y$samples$lib.size <- colSums(y$counts)

  ## TMM normalisation
  calcNormFactors(y)
}

make.design <- function(samples, masked.col.data, group.name, extra.groups=c()) {
  if (length(extra.groups)>0) {
    contin.dat <- masked.col.data[, extra.groups, drop=FALSE]
  }
  if (length(group.name)>0) {
    levels <- unique(masked.col.data[, group.name])
    group <- factor(masked.col.data[, group.name], levels=levels)
    ff <- "~0+group"
  } else {
    ff <- "~1"
  }
  for (contin in extra.groups){
    ff <- paste0(ff, "+contin.dat$", contin)
    if (is.factor(contin.dat[, contin])) {
      contin.dat[, contin] <- droplevels(contin.dat[, contin])
    }
  }
  design <- model.matrix(formula(ff))
  rownames(design) <- rownames(samples)
  # Fix common names issues to keep them readable
  colnames(design) <- sub("-", "neg", colnames(design))
  colnames(design) <- sub("[+]", "pos", colnames(design))
  colnames(design) <- sub("group", "", colnames(design))
  colnames(design) <- sub("contin.dat\\$", "", colnames(design))

  return(design)
}

make.dispersions <- function(y, design) {
  ## estimate dispersions
  y <- estimateGLMCommonDisp(y, design, verbose=FALSE)
  y <- estimateGLMTrendedDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  y
}

annotate.table <- function(table.fdr, species="Hs") {
  library(AnnotationDbi)

  ids <- row.names(table.fdr)

  if (species=="None" || species=="none") {
    return(table.fdr)
  } else if (species %in% c("Hs", "Mm")) {
    # Use static db
    if (species=="Hs") {
      library(org.Hs.eg.db)
      db <- org.Hs.eg.db
      columns <- c("ENSEMBL", "SYMBOL", "GENENAME", "ENTREZID", "GENETYPE", "MAP", "PATH")
    } else {
      library(org.Mm.eg.db)
      db <- org.Mm.eg.db
      columns <- c("ENSEMBL", "SYMBOL", "GENENAME", "ENTREZID", "GENETYPE", "PATH")
    }
    # dplyr can over-ride select so make the package explicit
    annot <- AnnotationDbi::select(db, keys=ids, keytype="ENSEMBL", columns=columns)

    # This fixes that some Ensembl ids map to multiple
    # other features, including SYMBOL, which is
    # undesirable behaviour from AnnotationDb
    annot <- annot[match(ids, annot$ENSEMBL), ]
  } else {
    # Pull gene names from gtf
    if (species=="GTF" || species=="gtf") {
      if (!("species.gtf" %in% ls())) {
        stop("No 'species.gtf' object to use for annotation")
      }
    } else {
      # Read gtf
      if (!file.exists(species)) {
        stop(paste("No such file:", species, "; 'species' should be Hs, Mm, GTF or the path to a
 .gtf file"))
      }
      species.gtf <- rtracklayer::readGFF(species)
      species.gtf <<- species.gtf[species.gtf$type == "gene", ] # Store in global environment
    }
    annot <- species.gtf[match(ids, species.gtf$ENSEMBL), ]
  }

  if (!identical(ids, annot[, 1])) {
    warning("Error matching gene names in annotate.table(), most likely some identifiers were not 
found")
  }

  return(cbind(table.fdr, annot[, -1]))
}

################################################
# Run standard set of pathway analyses using clusterProfiler
################################################

analyse.cP <- function(results_table, file_prefix, species="human", FDR.thresh=0.05,
                       logCPM.thresh=1, pvalueCutoff=0.1, fformat=".png") {
  # human isn't used currently; at some point adjust the below to handle mouse as well

  library(clusterProfiler)
  library(ReactomePA)
  # Assumes annotated results list

  # Set up entrezID mapping to fix tables that don't convert back to gene symbols
  entrez.mapping <- results_table[, c("SYMBOL", "ENTREZID")]

  # Get significant genes
  sig <- results_table[results_table$logCPM > logCPM.thresh, ]
  sig <- sig[sig$FDR < FDR.thresh, ]

  # Set up output
  out <- list(KEGG=NA, WikiPathways=NA, Reactome=NA, eGO_BP=NA, eGO_CC=NA)

  if (nrow(sig) > 0) {
    # Get Entrez IDs
    entrez_top.id <- as.character(sig$ENTREZID)
    entrez_uni.id <- as.character(results_table$ENTREZID)

    # Run enrichment analyses of various kinds
    out$KEGG <- enrich.helper(enrichKEGG, file_prefix, "clusterProfiler_KEGG", fformat,
                              entrez.mapping=entrez.mapping,
                              gene=entrez_top.id, universe=entrez_uni.id,
                              organism="hsa", pvalueCutoff=pvalueCutoff)
    out$WikiPathways <- enrich.helper(enrichWP, file_prefix, "clusterProfiler_WikiPathways",
                                      fformat, entrez.mapping=entrez.mapping,
                                      entrez_top.id, organism="Homo sapiens",
                                      universe=entrez_uni.id)
    out$Reactome <- enrich.helper(enrichPathway, file_prefix, "clusterProfiler_ReactomePA", fformat,
                                  gene=entrez_top.id, pvalueCutoff=pvalueCutoff,
                                  universe=entrez_uni.id, readable=TRUE)
    out$eGO_BP <- enrich.helper(enrichGO, file_prefix, "clusterProfiler_GO_BP", fformat,
                                gene=entrez_top.id, OrgDb="org.Hs.eg.db", keyType="ENTREZID",
                                ont="BP", pAdjustMethod="BH", pvalueCutoff=pvalueCutoff,
                                universe=entrez_uni.id, readable=TRUE, pool=FALSE)
    out$eGO_CC <- enrich.helper(enrichGO, file_prefix, "clusterProfiler_GO_CC", fformat,
                                gene=entrez_top.id, OrgDb="org.Hs.eg.db", keyType="ENTREZID",
                                ont="CC", pAdjustMethod="BH", pvalueCutoff=pvalueCutoff,
                                universe=entrez_uni.id, readable=TRUE, pool=FALSE)
  }

  return(out)
}

enrich.helper <- function(func, file_prefix, filename_part, fformat, entrez.mapping=data.frame(),
                          ...) {
  out <- func(...)
  if (nrow(entrez.mapping)>0) {
    # Remap entrez ids to gene symbols
    out <- enrich.fix.genes(out, entrez.mapping)
  }
  filename <- paste0(file_prefix, "_", filename_part, ".txt")
  write.table(out, filename, row.names=FALSE, sep="\t", quote=FALSE)
  if (nrow(out)>0) {
    out.plot <- barplot(out, showCategory=25)
    ggsave(sub(".txt", fformat, filename), out.plot)
  }
  out
}

enrich.fix.genes <- function(result.df, entrez.mapping) {
  tmp <- sapply(result.df@result$geneID, enrich.fix.genes.helper,
                entrez.mapping=entrez.mapping, USE.NAMES=FALSE)
  result.df@result$geneID <- tmp
  return(result.df)
}

enrich.fix.genes.helper <- function(geneIDlist, entrez.mapping) {
  entrez.ids <- strsplit(geneIDlist, "/")[[1]]
  gene.symbols <- entrez.mapping$SYMBOL[match(entrez.ids, as.character(entrez.mapping$ENTREZID))]
  out <- paste0(gene.symbols, collapse="/")
  return(out)
}
