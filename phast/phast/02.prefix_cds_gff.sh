cat /fast3/group_crf/home/g21shaoy23/goby5/annotations/cds_phase_fix2/Tra.prechange_longest_fixed.gff3 | awk '{if ($3=="CDS") print "Tridentiger_radiatus."$0 }' > cds.gff
