bam=$1
genome=/public3/group_crf/home/g21shaoy23/goby5/vg/vcf/genome/Tra_sort.fa
gatk=/public3/group_crf/home/wanglinlong8405/genome/lty_genome/resequencing/01.SNP_calling/03.gatk/gatk-4.2.1.0/gatk
CPU=10
outstem=$1
SCRATCHDIR=`pwd`/tmp
MAXDEPTH=250

ulimit -n 99999 #change file handler limit!

$gatk CreateSequenceDictionary -R $genome

mkdir -p $SCRATCHDIR

bam=`realpath $bam`
genome=`realpath $genome`

#STEP 1: Sort bam file by coordinate
STEP=01.sortbam

if [ -f $STEP.done ]; then
    echo "$STEP done, skip...";
else
    samtools sort -@ $CPU -o $outstem.sorted.bam $bam >$STEP.log 2>&1 \
    && touch  $STEP.done

    if [  "$?" -ne  "0" ]; then
      exit 1;
    fi
fi


#STEP 2: Add Readgroup tag to bam
STEP=02.addRG

if [ -f $STEP.done ]; then
    echo "$STEP done, skip...";
else
    $gatk --java-options '-Xmx16g'  AddOrReplaceReadGroups I=$outstem.sorted.bam O=$outstem.RG.bam  RGID=$outstem RGLB=$outstem RGPL=mgi RGSM=$outstem RGPU=$outstem CREATE_INDEX=True VALIDATION_STRINGENCY=SILENT TMP_DIR=$SCRATCHDIR  >$STEP.log 2>&1 \
    && touch  $STEP.done

    if [  "$?" -ne  "0" ]; then
      exit 1;
    fi
fi


#STEP 3: Mark PCR Duplicates
STEP=03.dedup

if [ -f $STEP.done ]; then
    echo "$STEP done, skip...";
else
    $gatk --java-options '-Xmx16g'  MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=900  INPUT=$outstem.RG.bam OUTPUT=$outstem.dedup.bam ASSUME_SORTED=true METRICS_FILE=$outstem.dedup.metrics.txt VALIDATION_STRINGENCY=SILENT TMP_DIR=$SCRATCHDIR  >$STEP.log 2>&1 \
    && samtools index $outstem.dedup.bam \
    && touch  $STEP.done

    if [  "$?" -ne  "0" ]; then
      exit 1;
    fi
fi

#GATK 4 no longer needs indel realignment!

#STEP 4: Indel Realignment Create
STEP=04.haplotypecaller

if [ -f $STEP.done ]; then
    echo "$STEP done, skip...";
else
    $gatk --java-options '-Xmx128g' HaplotypeCaller  --native-pair-hmm-threads $CPU  -I $outstem.dedup.bam -O $outstem.g.vcf -R $genome --max-reads-per-alignment-start $MAXDEPTH --minimum-mapping-quality 30 -ERC BP_RESOLUTION --dont-use-soft-clipped-bases true  >$STEP.log 2>&1  \
    && touch  $STEP.done

    if [  "$?" -ne  "0" ]; then
      exit 1;
    fi
fi



#STEP 5: compress
STEP=05.compress

if [ -f $STEP.done ]; then
    echo "$STEP done, skip...";
else
    bgzip -@ $CPU $outstem.g.vcf > $STEP.log 2>&1 && tabix -p vcf $outstem.g.vcf.gz >> $STEP.log 2>&1 \
    && touch  $STEP.done

    if [  "$?" -ne  "0" ]; then
      exit 1;
    fi
fi


#STEP 6: genotype gvcf
STEP=06.genotypegvcf

if [ -f $STEP.done ]; then
    echo "$STEP done, skip...";
else
    $gatk --java-options '-Xmx12g' GenotypeGVCFs -O $outstem.genotyped.g.vcf -R $genome -V $outstem.g.vcf.gz -all-sites true > $STEP.log 2>&1 \
    && touch  $STEP.done

    if [  "$?" -ne  "0" ]; then
      exit 1;
    fi
fi

#STEP 7: compress
STEP=07.compress

if [ -f $STEP.done ]; then
    echo "$STEP done, skip...";
else
    bgzip -@ $CPU $outstem.genotyped.g.vcf > $STEP.log 2>&1 && tabix -p vcf $outstem.genotyped.g.vcf.gz >> $STEP.log 2>&1 \
    && touch  $STEP.done

    if [  "$?" -ne  "0" ]; then
      exit 1;
    fi
fi
