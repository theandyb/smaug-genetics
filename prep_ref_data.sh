#!/bin/bash

#############################################################################
# Script does computationally heavier formatting of the reference data
#############################################################################
source ./parse_yaml.sh

# read yaml file
eval $(parse_yaml _config.yaml "config_")
refdir="$config_parentdir/reference_data"
curdir=${PWD}
cd $refdir

#############################################################################
# hg19 chromosome lengths
#############################################################################

# Make fixed-width windows
bedtools makewindows -g "$refdir/hg19.genome" -w 1000000 | grep -Ev "_|X|Y|M" | sort -k 1,1 -k2,2n > "$refdir/genome.1000kb.sorted.bed"

bedtools makewindows -g "$refdir/hg19.genome" -w 5000000 | grep -Ev "_|X|Y|M" | sort -k 1,1 -k2,2n > "$refdir/genome.5000kb.sorted.bed"

bedtools makewindows -g "$refdir/hg19.genome" -w 100000 | grep -Ev "_|X|Y|M" | sort -k 1,1 -k2,2n > "$refdir/genome.100kb.sorted.bed"

bedtools makewindows -g "$refdir/hg19.genome" -w 10000 | grep -Ev "_|X|Y|M" | sort -k 1,1 -k2,2n > "$refdir/genome.10kb.sorted.bed"

bedtools makewindows -g "$refdir/hg19.genome" -w 3000000000 | grep -Ev "_|X|Y|M" | sort -k 1,1 -k2,2n > "$refdir/genome.full.sorted.bed"

#############################################################################
# Compress and index reference genomes
#############################################################################
# v 37
for i in `seq 1 22`; do
	samtools faidx "$refdir/human_g1k_v37/human_g1k_v37.fasta" $i | bgzip -c > "$refdir/human_g1k_v37/chr$i.fasta.gz"
	samtools faidx "$refdir/human_g1k_v37/chr$i.fasta.gz"
done

# mask v 37
for i in `seq 1 22`; do
	samtools faidx "$refdir/human_g1k_v37_mask/human_g1k_v37.mask.fasta" $i | bgzip -c > "$refdir/human_g1k_v37_mask/chr$i.fasta.gz"
	samtools faidx "$refdir/human_g1k_v37_mask/chr$i.fasta.gz"
done

# ancestral genome
for i in `seq 1 22`; do
	echo "$refdir/human_ancestor_GRCh37_e59/human_ancestor_$i.fa"
	cat "$refdir/human_ancestor_GRCh37_e59/human_ancestor_$i.fa" | sed "s,^>.*,>$i," | \
		bgzip -c > "$refdir/human_ancestor_GRCh37_e59/human_ancestor_$i.fa.gz"
	samtools faidx "$refdir/human_ancestor_GRCh37_e59/human_ancestor_$i.fa.gz"
done
