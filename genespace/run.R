#setwd("~/bahaha_assembly/synteny/genespace/example/");
.libPaths("/data/projects/rcui/R/x86_64-pc-linux-gnu-library/4.1") 
library(GENESPACE)
runwd <- file.path("/data2/projects/yshao/Goby5Compare/synteny2/rundir2/")
list.files(runwd, recursive = T, full.names = F)

 gpar <- init_genespace(
   genomeIDs = c("TridentigerRadiatus", "TridentigerBarbatus", "TridentigerBifasciatus", "GlossogobiusGiuris", "GobiopsisMacrostoma"),
   speciesIDs = c("Tridentiger_radiatus", "Tridentiger_barbatus", "Tridentiger_bifasciatus", "Glossogobius_giuris", "Gobiopsis_macrostoma"),
   versionIDs = c("1.0", "1.0", "1.0", "1.0", "1.0"),
   ploidy = rep(1,5),
   diamondMode = "default",
   orthofinderMethod = "default",
   wd = runwd,
   nCores = 60,
   minPepLen = 30,
   gffString = "gff",
   pepString = "pep",
   path2orthofinder = "/opt/miniconda3/bin/orthofinder",
   path2diamond = "/opt/miniconda3/bin/diamond",
   path2mcscanx = "/data/software/MCScanX/",
   rawGenomeDir = file.path(runwd, "rawGenomes"))

  
parse_annotations(
  gsParam = gpar,
  gffEntryType = "gene",
  gffIdColumn = "locus",
  gffStripText = "locus=",
  headerEntryIndex = 1,
  headerSep = " ",
  headerStripText = "locus=")

gpar <- run_orthofinder(
  gsParam = gpar)

gpar <- synteny(gsParam = gpar)

arrCol <- rainbow(11);
regs <- data.frame(
  genome = rep('TridentigerRadiatus', 11),
  chr =1:11)

genomes <- c("TridentigerBarbatus", "TridentigerBifasciatus","GlossogobiusGiuris","GobiopsisMacrostoma")
chr_values <- list(c(16, 2, 4, 3, 11, 13, 9, 6, 18), 
                   c(14, 13, 1, 7, 21, 18, 22, 5),
                   c(16, 4, 22, 5,3 ,8 , 11, 6,21, 10,14 ,12 ),
                   c(10, 15, 23, 16, 21,5  ))

genome_col <- rep(genomes, times = sapply(chr_values, length))
chr_col <- unlist(chr_values)

datInvert <- data.frame(genome = genome_col, chr = chr_col)


#datInvert <- data.frame(genome="TridentigerBarbatus", chr=c(16, 2, 4, 3, 11, 13, 9, 6, 18) )

ripdat <- plot_riparianHits(
  gpar,onlyTheseRegions = regs, refChrCols = arrCol,minGenes2plot=50, invertTheseChrs = datInvert,
  blackBg = F, chrFill = "orange",returnSourceData=T,
  chrBorder = "grey", useOrder=F  )#labelTheseGenomes = c('TridentigerRadiatus', 'TridentigerBarbatus')

#dump('ripdat', file="ripdata.R");

pg <- pangenome(gpar)