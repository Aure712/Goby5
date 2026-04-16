#!/usr/bin/bash
#SBATCH -p himem
#SBATCH -c 2
#SBATCH --mem=24G

JUICERDIR=/public3/group_crf/software/juicer-1.6/
HIC=/public3/group_crf/home/g21shaoy23/HiC/Juicer/Tradiatus/fastq/

$JUICERDIR/scripts/juicer.sh \
-D /public3/group_crf/software/juicer-1.6/ \
-g ref \
-d `pwd` \
-s DpnII \
-z `pwd`/references/ref.fa \
-p `pwd`/restriction_sites/ref.chrom.sizes \
-y `pwd`/restriction_sites/ref_DpnII.txt \
-q himem \
-l himem \
-t 10 \
> juicer.log 2>&1
