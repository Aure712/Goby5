
arrOutgroupSpp='Glossogobius_giuris,Neogobius_melanostomus,Rhinogobius_formosanus,Rhinogobius_similis,Tridentiger_barbatus,Tridentiger_radiatus,Tridentiger_bifasciatus,Taenioides_sp,Boleophthalmus_pectinirostris,Periophthalmus_modestus,Siphamia_tubifer,Odontamblyopus_rebecca,Oxyeleotris_marmorata,Periophthalmus_magnuspinnatus';
arrUnusedClades=''; #c('Aplocheilidae', 'root'); #mark as unused clades "only for hyphy relax"
arrMarkUnusedCladeChildren=''; #c(T , F);

#######################################################################
arrForegroundClades='Gobiopsis_macrostoma'; #mark as foreground
arrMarkForegroundChildren='T';

sTaxonRequirements="min_taxon_requirements.txt";

nMarkStyle='relax'; #either "codeml" or "relax"
sOutDIR="Relax_Gma";

mkdir -p $sOutDIR
Rscript labelbranches.R $sOutDIR  $nMarkStyle  $sTaxonRequirements  $arrOutgroupSpp  $arrForegroundClades  $arrMarkForegroundChildren  $arrUnusedClades  $arrMarkUnusedCladeChildren > $sOutDIR/log.txt &


#sOutDIR="Codeml_all";
#arrMarkForegroundChildren='T';
#nMarkStyle='codeml';
#mkdir -p $sOutDIR
#Rscript labelbranches.R $sOutDIR  $nMarkStyle  $sTaxonRequirements  $arrOutgroupSpp  $arrForegroundClades  $arrMarkForegroundChildren  $arrUnusedClades  $arrMarkUnusedCladeChildren > $sOutDIR/log.txt &

########################################################################

wait
