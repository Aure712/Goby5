dustmasker -in scf.fa -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out Tba.asnb

makeblastdb -in scf.fa -input_type fasta -dbtype nucl -parse_seqids -mask_data Tba.asnb -out Tba
