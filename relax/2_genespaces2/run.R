#setwd("~/bahaha_assembly/synteny/genespace/example/");
.libPaths("/data/projects/rcui/R/x86_64-pc-linux-gnu-library/4.1") 
library(GENESPACE)
runwd <- file.path("/data2/projects/yshao/relax/2_genespace2/rundir/")
#list.files(runwd, recursive = T, full.names = F)

gpar <- init_genespace(
  genomeIDs = c('Boleophthalmuspectinirostris','Glossogobiusgiuris','Gobiopsismacrostoma','Neogobiusmelanostomus','Periophthalmusmodestus','Rhinogobiusformosanus','Rhinogobiussimilis','Siphamiatubifer','Taenioidessp','Tridentigerbarbatus','Tridentigerbifasciatus','Tridentigerradiatus','Oxyeleotrismarmorata','Odontamblyopusrebecca','Periophthalmusmagnuspinnatus'),
  speciesIDs = c('Boleophthalmus_pectinirostris','Glossogobius_giuris','Gobiopsis_macrostoma','Neogobius_melanostomus','Periophthalmus_modestus','Rhinogobius_formosanus','Rhinogobius_similis','Siphamia_tubifer','Taenioides_sp','Tridentiger_barbatus','Tridentiger_bifasciatus','Tridentiger_radiatus','Oxyeleotris_marmorata','Odontamblyopus_rebecca','Periophthalmus_magnuspinnatus'),
  versionIDs = c( "1","1","1","1","1","1","1","1","1","1","1","1","1","1","1"),
  ploidy = rep(1,15),
  diamondMode = "default",
  orthofinderMethod = "default",
  wd = runwd,
  nCores = 80,
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

#arrCol <- rainbow(22);
#regs <- data.frame(
#  genome = rep("Rformfasta", 22),
#  chr =(1:22)

#datInvert <- data.frame(genome="Rformfasta", chr=c(1, 4, 8, 10, 6, 18,19,20, 22) )

#ripdat <- plot_riparianHits(
#  gpar,onlyTheseRegions = regs, refChrCols = arrCol,minGenes2plot=50, invertTheseChrs = datInvert,
#  blackBg = F, chrFill = "orange",returnSourceData=T,
#  chrBorder = "grey", useOrder=F, labelTheseGenomes = c('Rformfasta') )

#dump('ripdat', file="ripdata.R");

pg <- pangenome(gpar)
