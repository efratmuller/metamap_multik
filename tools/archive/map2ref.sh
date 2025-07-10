#!/bin/bash

# Fail on any error
set -euo pipefail

# Commands for testing:

# lab-dir; cd emuller/phages_project/scripts/metamap_multik
# conda activate metamap_test

# (with slurm)
# sbatch -A ALMEIDA-SL2-CPU -J "m2r_v" -o /home/em2035/rds/hpc-work/mapping_tests/ERR6998314/map2ref_v.out -e /home/em2035/rds/hpc-work/mapping_tests/ERR6998314/map2ref_v.err -p cclake --mem=32G --ntasks=1 --cpus-per-task=16 --time=6:00:00 --wrap="./tools/map2ref.sh -t 16 -e 0.5 -b 10 -i /rds/project/rds-aFEMMKDjWlo/rfs_data/metagenomes/fastqs/ERP115476/ERR6998314/ERR6998314_1.fastq.gz -n /rds/project/rds-aFEMMKDjWlo/rfs_data/metagenomes/fastqs/ERP115476/ERR6998314/ERR6998314_2.fastq.gz -r /rds/project/rds-aFEMMKDjWlo/rfs_data/databases/bwa/uhgv_nayfach/filtered_db/votus_mq_plus_no_contam.fna -o /home/em2035/rds/hpc-work/mapping_tests/ERR6998314/ERR6998314_uhgv -c complete"

# sbatch -A ALMEIDA-SL2-CPU -J "m2r_b" -o /home/em2035/rds/hpc-work/mapping_tests/ERR6998314/map2ref_b.out -e /home/em2035/rds/hpc-work/mapping_tests/ERR6998314/map2ref_b.err -p cclake --mem=32G --ntasks=1 --cpus-per-task=16 --time=6:00:00 --wrap="./tools/map2ref.sh -t 16 -e 0.3 -b 5 -i /rds/project/rds-aFEMMKDjWlo/rfs_data/metagenomes/fastqs/ERP115476/ERR6998314/ERR6998314_1.fastq.gz -n /rds/project/rds-aFEMMKDjWlo/rfs_data/metagenomes/fastqs/ERP115476/ERR6998314/ERR6998314_2.fastq.gz -r /rds/project/rds-aFEMMKDjWlo/rfs_data/databases/bwa/uhgg_v1.1/uhgg_v1.1.fa -o /home/em2035/rds/hpc-work/mapping_tests/ERR6998314/ERR6998314_uhgg -c contigs"

# conda deactivate


# SRP128485/SRR6456370
# ERP115476/ERR6998314


usage()
{
cat << EOF
usage: $0 options

Map metagenomic reads against a custom-made BWA reference
and output statistics to infer genome prevalence and abundance.

OPTIONS:
   -t      Number of threads [REQUIRED]
   -i      Input first or only read (.fastq or .fastq.gz format) [REQUIRED]
   -n      Input reverse read (.fastq or fastq.gz [REQUIRED]
   -r      Reference database (.fa or .fasta) [REQUIRED]
   -o      Output files prefix (include path) [REQUIRED]
   -c      Whether entries in ref database are "contigs" or "complete" genomes [REQUIRED]
   -b      Minimal breadth (coverage) for considering a genome present in a sample [Default: 5]
   -e      Minimal ratio of observed/expected breadth for considering a genome present in a sample [Default: 0.5]
EOF
}

# Initialize variables
threads=
ref=
reads=
reads2=
outprefix=
mode=
min_br=5
min_ratio=0.5
clean_intermediate_files=1

# Read command line arguments
while getopts “t:i:n:r:o:c:b:e:” OPTION
do
     case ${OPTION} in
         t)
             threads=${OPTARG}
             ;;
         i)
             reads=${OPTARG}
             ;;
         n)
             reads2=${OPTARG}
             ;;
         r)
             ref=${OPTARG}
             ;;
         o)
             outprefix=${OPTARG}
             ;;
         c)
             mode=${OPTARG}
             ;;
         b)
             min_br=${OPTARG}
             ;;
         e)
             min_ratio=${OPTARG}
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

# Check that all arguments were provided
if [[ -z ${threads} ]] || [[ -z ${reads} ]] || [[ -z ${ref} ]] || [[ -z ${outprefix} ]] || [[ -z ${mode} ]] || [[ -z ${min_br} ]] || [[ -z ${min_ratio} ]]
then
     echo "ERROR : Please supply correct arguments"
     usage
     exit 1
fi

timestamp() {
  date +"%H:%M:%S"
}

# Verify proper "mode" argument
if [[ "${mode}" != "complete" && "${mode}" != "contigs" ]]; then
    echo "Error: The mode variable should be either 'complete' nor 'contigs'." >&2
    exit 1
fi

# Create output directory if it doesnt exist
outdir=$(dirname ${outprefix})
mkdir -p -v ${outdir}

# Set threads for bwa
if [ ${threads} -eq 1 ]
then
    threads_sam=1
else
    threads_sam=$((${threads}-1))
fi

# Start pipeline
echo "$(timestamp) [ map2ref pipeline ] Running map2ref pipeline with arguments:"
echo $@

# Initialize a file to hold statistics about mapped reads at each step
echo -e "category\tkey\tvalue" > ${outprefix}_stats.tab

# Record the total number of qc reads
count=$(zcat ${reads} | wc -l)
echo -e "READS\t1_total_qc_reads\t$((count / 2))" >> ${outprefix}_stats.tab # Dividing by 4 gives the number of read pairs, dividing by 2 gives total number of reads.
if [ ${count} -eq 0 ]
then
    echo "$(timestamp) [ map2ref pipeline ] ERROR: the fastq provided is empty."
    sync # Flush outputs before exiting
    exit 1
fi

# -----------------------------------------------------------------
# (1) Initial mapping agaist the reference database
# -----------------------------------------------------------------

echo "$(timestamp) [ map2ref pipeline ] Running BWA ..."
rm -f ${outprefix}*bam*
bwa mem -M -t ${threads} ${ref} ${reads} ${reads2} \
  | samtools view -@ ${threads_sam} -F 256 -uS - \
  | samtools sort -@ ${threads_sam} - -o ${outprefix}_raw.bam
samtools index -@ ${threads_sam} ${outprefix}_raw.bam

# Record the total number of mapped reads (-F 4 removes unmapped reads)
count=`samtools view -F 4 ${outprefix}_raw.bam | wc -l`
if [ ${count} -eq 0 ]
then
    echo "$(timestamp) [ map2ref pipeline ] ERROR: bwa failed. Verify sufficient resources for bwa mem run, or problems with input fastq files (e.g., issues with paired reads)."
    sync # Flush outputs before exiting
    exit 1
fi
echo -e "READS\t2_mapped_reads\t${count}" >> ${outprefix}_stats.tab

# -----------------------------------------------------------------
# (2) Keep only alignments with >60% alignment fraction & >90% ANI
# -----------------------------------------------------------------

echo "$(timestamp) [ map2ref pipeline ] Filtering by coverage and ANI ..."
tools/bam_ani-filter.py ${outprefix}_raw.bam 90 60 | samtools view - -o ${outprefix}_ani_cov.bam
samtools index -@ ${threads_sam} ${outprefix}_ani_cov.bam

# Record the number of mapped reads following the above filter
count=`samtools view ${outprefix}_ani_cov.bam | wc -l`
echo -e "READS\t3_mapped_reads_90_ANI_60_AF_filters\t${count}" >> ${outprefix}_stats.tab

# -----------------------------------------------------------------
# (3) Extract unique counts (reads mapped to a single reference)
# -----------------------------------------------------------------

echo "$(timestamp) [ map2ref pipeline ] Extracting unique counts ..."
samtools view -@ ${threads_sam} -q 1 ${outprefix}_ani_cov.bam -o ${outprefix}_unique.bam
samtools index -@ ${threads_sam} ${outprefix}_unique.bam

# Record the number of reads with unique mappings
count=`samtools view ${outprefix}_unique.bam | wc -l`
echo -e "READS\t4_mapped_reads_uniques\t${count}" >> ${outprefix}_stats.tab

# -----------------------------------------------------------------
# (4) For the unique counts - also remove paired reads that were 
#     mapped to different genomes
# -----------------------------------------------------------------

samtools view -@ ${threads_sam} -f 2 ${outprefix}_unique.bam -o ${outprefix}_unique_prop_pair.bam
samtools index -@ ${threads_sam} ${outprefix}_unique_prop_pair.bam

# Record the number of uniquely-mapped reads with paired ends properly aligned to the same genome
count=`samtools view ${outprefix}_unique_prop_pair.bam | wc -l`
echo -e "READS\t5_mapped_reads_uniques_proper_pairs\t${count}" >> ${outprefix}_stats.tab

# -----------------------------------------------------------------
# (5) Parse total count output (this time including reads mapped to 
#     multiple references but randomly assigned to one)
# -----------------------------------------------------------------

echo "$(timestamp) [ map2ref pipeline ] Parsing alignment results ..."
samtools idxstats ${outprefix}_ani_cov.bam > ${outprefix}_depth.tab
samtools depth ${outprefix}_ani_cov.bam > ${outprefix}_depth-pos.tab

echo "$(timestamp) [ map2ref pipeline ] Counting total reads per genome ..."
tools/parse_bwa-depth.py ${outprefix}_depth.tab ${outprefix}_depth-pos.tab ${mode} ${outprefix} \
  > ${outprefix}_total.tab

# Record statistics
count=`wc -l < ${outprefix}_total.tab`
echo -e "GENOMES\t1_total_genomes\t$((count - 1))" >> ${outprefix}_stats.tab
count=`awk -F'\t' '$3 > 0' ${outprefix}_total.tab | wc -l`
echo -e "GENOMES\t2_genomes_with_mapped_reads\t${count}" >> ${outprefix}_stats.tab

# -----------------------------------------------------------------
# (6) Parse unique count output
# -----------------------------------------------------------------

echo "$(timestamp) [ map2ref pipeline ] Counting unique reads per genome ..."
samtools idxstats ${outprefix}_unique_prop_pair.bam > ${outprefix}_unique_depth.tab
samtools depth ${outprefix}_unique_prop_pair.bam > ${outprefix}_unique_depth-pos.tab
tools/parse_bwa-depth.py ${outprefix}_unique_depth.tab ${outprefix}_unique_depth-pos.tab ${mode} ${outprefix} \
  > ${outprefix}_unique.tab
count=`awk -F'\t' '$3 > 0' ${outprefix}_unique.tab | wc -l`
echo -e "GENOMES\t3_genomes_with_uniquely_mapped_reads\t${count}" >> ${outprefix}_stats.tab

# -----------------------------------------------------------------
# (7) Infer which genomes are present based on the observed breadth
#     and observed/expected breadth ratio.
# -----------------------------------------------------------------

echo "$(timestamp) [ map2ref pipeline ] Inferring presence/absence of species ..."

# Note: The filtering of genomes is based on coverage and expected coverage as calculated based 
#  on the *total* number of reads assigned to each genome, and not only uniquely-assigned reads.
awk -F'\t' -v min_br=${min_br} -v min_ratio=${min_ratio} \
  '$3 > 0 && $5 > min_br && $8 > min_ratio { print $1 }' \
  ${outprefix}_total.tab > ${outprefix}_present_genomes.txt
count=`wc -l < ${outprefix}_present_genomes.txt`
echo -e "GENOMES\t4_genomes_present_after_breadth_filters\t${count}" >> ${outprefix}_stats.tab

# Check if no genomes were inferred as present at all
if [ ${count} -eq 0 ]
then
    echo "$(timestamp) [ map2ref pipeline ] No genomes inferred as present in this sample when using the provided coverage thresholds."
    echo -e "READS\t6_mapped_reads_uniques_proper_pairs_genomes_present\t0" >> ${outprefix}_stats.tab
    echo -e "GENOMES\t5_genomes_present_with_uniquely_mapped_reads\t0" >> ${outprefix}_stats.tab
    
    # Create empty files to support downstream commands / rules
    echo "" > ${outprefix}_aligned_reads_fwd.txt
    gzip -f ${outprefix}_aligned_reads_fwd.txt
    echo "" > ${outprefix}_aligned_reads_rev.txt
    gzip -f ${outprefix}_aligned_reads_rev.txt
    head -n 1 ${outprefix}_unique.tab > ${outprefix}_unique_filtered.tab
else

    # -------------------------------------------------------------
    # (8) Save a list of aligned reads (uniquely mapped, after 
    #     filtering genomes inferred as present in the sample)
    # -------------------------------------------------------------

    # This list will be used to calculate the number of reads that were mapped to multiple reference catalogs.
    echo "$(timestamp) [ map2ref pipeline ] Filtering read counts by inferred species presence ..."

    # Create a temporary file with patterns for fast filtering with grep 
    if [[ ${mode} == "complete" ]]; then
        sed ':a;N;$!ba;s/\n/\t|/g' ${outprefix}_present_genomes.txt \
          | head -c -1 > ${outprefix}_present_genomes_temp.txt
        echo -e "\t" >> ${outprefix}_present_genomes_temp.txt
    elif [[ ${mode} == "contigs" ]]; then
        sed ':a;N;$!ba;s/\n/_|/g' ${outprefix}_present_genomes.txt | \
          head -c -1 > ${outprefix}_present_genomes_temp.txt
        echo "_" >> ${outprefix}_present_genomes_temp.txt
    fi

    # Filter read alignments using the above patterns
    # (We record forward and reverse reads seperatly, as even though reads were properly paired, 
    #  in some cases one of the strands did not pass the ANI/coverage thresholds)
    samtools view -F 16 ${outprefix}_unique_prop_pair.bam \
      | grep -E -f ${outprefix}_present_genomes_temp.txt \
      | cut -f 1 > ${outprefix}_aligned_reads_fwd.txt

    samtools view -f 16 ${outprefix}_unique_prop_pair.bam \
      | grep -E -f ${outprefix}_present_genomes_temp.txt \
      | cut -f 1 > ${outprefix}_aligned_reads_rev.txt

    # Record statistics (note that each read pair is listed once here, so the line count is multiplied by 2)
    count=$(cat ${outprefix}_aligned_reads_fwd.txt ${outprefix}_aligned_reads_rev.txt | wc -l)
    echo -e "READS\t6_mapped_reads_uniques_proper_pairs_genomes_present\t${count}" >> ${outprefix}_stats.tab

    # Compress
    gzip -f ${outprefix}_aligned_reads_fwd.txt
    gzip -f ${outprefix}_aligned_reads_rev.txt

    # -------------------------------------------------------------
    # (9) Filter the unique/total read counts to only include 
    #     genomes inferred as present
    # -------------------------------------------------------------

    # (Again, for fast filtering, create a grep pattern and then use grep)
    sed ':a;N;$!ba;s/\n/\t|/g' ${outprefix}_present_genomes.txt \
      | head -c -1 > ${outprefix}_present_genomes_temp.txt
    echo -e "\t" >> ${outprefix}_present_genomes_temp.txt
    head -n 1 ${outprefix}_unique.tab > ${outprefix}_unique_filtered.tab
    grep -E -f ${outprefix}_present_genomes_temp.txt ${outprefix}_unique.tab \
      >> ${outprefix}_unique_filtered.tab

    # Also remove species with no uniquely-mapped reads mapped to them
    awk -F'\t' 'NR==1 || $3 != 0' ${outprefix}_unique_filtered.tab \
      > ${outprefix}_unique_filtered_tmp.tab \
      && mv ${outprefix}_unique_filtered_tmp.tab ${outprefix}_unique_filtered.tab

    # How many species are we left with after excluding species with no uniquely mapped reads?
    count=`wc -l < ${outprefix}_unique_filtered.tab`
    echo -e "GENOMES\t5_genomes_present_with_uniquely_mapped_reads\t$((count-1))" >> ${outprefix}_stats.tab
fi

# -----------------------------------------------------------------
# (10) Clean intrmediate files
# -----------------------------------------------------------------

if [ ${clean_intermediate_files} -eq 1 ]
then
    echo "$(timestamp) [ map2ref pipeline ] Cleaning tmp files ..."
    rm -rf ${outprefix}_unique.ba*
    rm -rf ${outprefix}_unique_depth*
    rm -rf ${outprefix}_unique_prop_pair.ba*
    rm -rf ${outprefix}_raw.ba*
    rm -rf ${outprefix}_ani_cov.ba*
    rm -rf ${outprefix}_depth*
    rm -rf ${outprefix}_present_genomes*
fi

echo "$(timestamp) [ map2ref pipeline ] Done."
