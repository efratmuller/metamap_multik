#!/bin/bash

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
EOF
}

# variables
threads=
ref=
reads=
reads2=
outprefix=
mode=

while getopts “t:m:i:n:r:o:c:” OPTION
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
         ?)
             usage
             exit
             ;;
     esac
done

# check arguments
if [[ -z ${threads} ]] || [[ -z ${reads} ]] || [[ -z ${ref} ]] || [[ -z ${outprefix} ]] || [[ -z ${mode} ]]
then
     echo "ERROR : Please supply correct arguments"
     usage
     exit 1
fi


timestamp() {
  date +"%H:%M:%S"
}

readname=$(basename ${outprefix})
refix=$(basename ${ref%%.fa*})
if [ ${threads} -eq 1 ]
then
    threads_sam=1
else
    threads_sam=$((${threads}-1))
fi

# initial mapping and sorting
echo "$(timestamp) [ map2ref pipeline ] Running BWA ..."
rm -f ${outprefix}_${refix}*bam*
bwa mem -M -t ${threads} ${ref} ${reads} ${reads2} | samtools view -@ ${threads_sam} -F 256 -uS - | samtools sort -@ ${threads_sam} - -o ${outprefix}_${refix}_raw.bam
samtools index -@ ${threads_sam} ${outprefix}_${refix}_raw.bam

# keep only alignments where >60% of the read had >90% identity to the reference
tools/bam_ani-filter.py ${outprefix}_${refix}_raw.bam 90 60 | samtools view -u - -o ${outprefix}_${refix}_sorted.bam

# extract unique counts (reads mapped to a single reference)
echo "$(timestamp) [ map2ref pipeline ] Extracting unique counts ..."
samtools view -@ ${threads_sam} -q 1 -f 2 -u ${outprefix}_${refix}_sorted.bam -o ${outprefix}_${refix}_unique_sorted.bam
samtools index -@ ${threads_sam} ${outprefix}_${refix}_unique_sorted.bam

# parse unique count output
echo "$(timestamp) [ map2ref pipeline ] Parsing results ..."
samtools idxstats ${outprefix}_${refix}_unique_sorted.bam > ${outprefix}_${refix}_unique_depth.tab
samtools depth ${outprefix}_${refix}_unique_sorted.bam > ${outprefix}_${refix}_unique_depth-pos.tab
tools/parse_bwa-depth.py ${outprefix}_${refix}_unique_depth.tab ${outprefix}_${refix}_unique_depth-pos.tab ${mode} _${refix} > ${outprefix}_${refix}_unique.tab

# parse total count output (this time including reads mapped to multiple references but randomly assigned to one)
samtools index -@ ${threads_sam} ${outprefix}_${refix}_sorted.bam
samtools idxstats ${outprefix}_${refix}_sorted.bam > ${outprefix}_${refix}_depth.tab
samtools depth ${outprefix}_${refix}_sorted.bam > ${outprefix}_${refix}_depth-pos.tab
tools/parse_bwa-depth.py ${outprefix}_${refix}_depth.tab ${outprefix}_${refix}_depth-pos.tab ${mode} _${refix} > ${outprefix}_${refix}_total.tab

# save a list of aligned reads
# (to reduce space, each read-pair is noted once and unique mappings are saved seperatly from non-unique mappings)
# samtools view ${outprefix}_${refix}_sorted.bam | cut -f 1,2,5 | gzip -c > ${outprefix}_${refix}_aligned_reads.tab.gz
samtools view ${outprefix}_${refix}_unique_sorted.bam | cut -f 1 | sort | uniq > ${outprefix}_${refix}_aligned_reads_uniq.txt
samtools view ${outprefix}_${refix}_sorted.bam | cut -f 1 | sort | uniq | comm -13 ${outprefix}_${refix}_aligned_reads_uniq.txt - > ${outprefix}_${refix}_aligned_reads_non_uniq.txt
gzip ${outprefix}_${refix}_aligned_reads_uniq.txt
gzip ${outprefix}_${refix}_aligned_reads_non_uniq.txt

# clean tmp files
echo "$(timestamp) [ map2ref pipeline ] Cleaning tmp files ..."
rm -rf ${outprefix}_${refix}_unique_sorted.ba* ${outprefix}_${refix}_unique_depth*
rm -rf ${outprefix}_${refix}_raw.ba* ${outprefix}_${refix}_unsorted.bam ${outprefix}_${refix}_sorted.ba* ${outprefix}_${refix}_depth*
