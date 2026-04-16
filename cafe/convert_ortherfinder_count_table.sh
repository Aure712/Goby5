IN=/data2/projects/yshao/mcmctree/1_genespace/rundir/orthofinder/Results_Jan24/Orthogroups/Orthogroups.GeneCount.tsv
cat $IN | cut -f1-19 | awk '{if ($1=="Orthogroup") {print "Desc\t"$0 } else {print "(null)\t"$0 } }' > orthogroups.genecount.txt
