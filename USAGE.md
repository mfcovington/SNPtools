# Prep sample BAM files:

    cd /Users/mfc/git.repos/SNPtools/sample-files/bam/
    samtools view -hb R500.good.bam \
      A01:10001-20000 \
      A02:10001-20000 \
      A03:30001-40000 \
      A04:10001-20000 \
      A05:10001-20000 \
      A06:10001-20000 \
      A07:60001-70000 \
      A08:10001-20000 \
      A09:1-10000 \
      A10:90001-100000 > R500.10kb.bam
    samtools view -hb IMB211.good.bam \
      A01:10001-20000 \
      A02:10001-20000 \
      A03:30001-40000 \
      A04:10001-20000 \
      A05:10001-20000 \
      A06:10001-20000 \
      A07:60001-70000 \
      A08:10001-20000 \
      A09:1-10000 \
      A10:90001-100000 > IMB211.10kb.bam
    samtools index R500.10kb.bam
    samtools index IMB211.10kb.bam

# Identify Polymorphisms

## Find SNPs/indels

    BIN=/Users/mfc/git.repos/SNPtools/bin
    DB_DIR=/Users/mfc/git.repos/SNPtools/sample-files
    OUT_DIR=/Users/mfc/git.repos/SNPtools/sample-files/output

    for ID in R500 IMB211
    do
        $BIN/SNPfinder/snp_finder.pl \
          --id        $ID \
          --bam       $DB_DIR/bam/$ID.10kb.bam \
          --fasta     $DB_DIR/fa/B.rapa_genome_sequence_0830.fa \
          --seq_list  A01,A02,A03,A04,A05,A06,A07,A08,A09,A10 \
          --out_dir   $OUT_DIR \
          --snp_min   0.33 \
          --indel_min 0.33 \
          --threads   3
    done
