#!/bin/bash

# Usage: ./prepare_treesort_dataset.sh [--segments "..." --fast] fasta_path reference_segment outdir
# Using --fast will make all trees to be inferred with FastTree.
# By default (without --fast) the reference tree is inferred with IQ-Tree, which is recommended for better accuracy.
# Example usage: ./prepare_treesort_dataset.sh --segments "HA,NA" segments.fasta HA myoutdir
# Example with default segments:  ./prepare_treesort_dataset.sh segments.fasta HA myoutdir

# These are the default segment names
declare -a segments=("PB2" "PB1" "PA" "HA" "NP" "NA" "MP" "NS")
FAST=0

POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
	case $1 in
		--segments)
			SEGMENTS_STR="$2"
			segments=(${SEGMENTS_STR//,/ })
			shift  # past argument
			shift  # past value
			;;
		--fast)
			FAST=1
			shift
			;;
		-*|--*)
			echo "Unrecognized option $1"
			exit 1
			;;
		*)
			POSITIONAL_ARGS+=("$1")  # save positional arg
			shift  # past argument
			;;
	esac
done

set -- "${POSITIONAL_ARGS[@]}"

# Required arguments:
main_fasta="$1"  # Provide a path to a fasta file with all segments
ref_seg="$2"  # Name of the segment to use as the reference (typically - HA)
outdir="$3"  # Path to the directory to store the results

# Make sure the FASTA input file exists and is not empty.
if [ ! -e "$main_fasta" ] || [ ! -s "$main_fasta" ]; then
   echo -e "The FASTA input file does not exist or is empty.\n"
   exit 1
fi

rm -r $outdir  # Clear out the directory
mkdir $outdir  # Re-create the directory

name=${main_fasta##*/}

# Maintain an array of segments found in the FASTA file. 
declare -a found_segments=()

# Split out the segments and align them
for seg in "${segments[@]}"
do
   # Copy all this segment's sequences into a segment-specific FASTA file.
   cat $main_fasta | smof grep "|${seg}|" > "${outdir}/${seg}-${name}"

   # Was the FASTA file created and is it non-empty?
   if [[ -s "${outdir}/${seg}-${name}" ]]; then
      
      found_segments+=("$seg")

      echo -e "Aligning ${seg}...\n"
      mafft --thread 6 "${outdir}/${seg}-${name}" | sed "s/|${seg}|/|/g"> "${outdir}/${seg}-${name}.aln"

      if [ $? -ne 0 ]; then
         echo "MAFFT alignment failed"
         exit 1
      fi
   fi 

   rm "${outdir}/${seg}-${name}"
done

if [ ${#found_segments[@]} -eq 0 ]; then
   echo -e "No segments found in the input file.\n"
   exit 1
fi

# Calculate the total number of non-empty alignment files.
aln_count=$(find ${outdir} -maxdepth 1 -type f -path "*.aln" | wc -l)
empty_aln_count=$(find ${outdir} -maxdepth 1 -type f -path "*.aln" -empty | wc -l)

# Were alignment files generated?
if [ $((aln_count - empty_aln_count)) -lt 1 ]; then
   echo -e "No alignment files were generated.\n"
   exit 1
fi

if [ $FAST -eq 0 ]; then
	# Build fasttree trees in parallel for non-reference segments
	echo -e "Building non-reference trees in parallel with FastTree...\n"
	for seg in "${found_segments[@]}"
	do
		if [ $seg != $ref_seg ]; then
			fasttree -nt -gtr -gamma ${outdir}/${seg}-${name}.aln > ${outdir}/${seg}-${name}.tre &
		fi
	done
	wait  # Wait to finish.

	# Build an IQ-Tree tree for the reference segment. We use the GTR+F+R5 model by default which can be changed
	echo -e "Building the reference tree with IQ-Tree...\n"
	iqtree2 -s ${outdir}/${ref_seg}-${name}.aln -T 6 --prefix "${outdir}/${ref_seg}-${name}" -m GTR+G+R5
	mv ${outdir}/${ref_seg}-${name}.treefile ${outdir}/${ref_seg}-${name}.tre

else
	# Build all trees with FastTree in parallel.
	echo -e "Building trees in parallel with FastTree...\n"
	for seg in "${found_segments[@]}"
	do
		fasttree -nt -gtr -gamma "${outdir}/${seg}-${name}.aln" > "${outdir}/${seg}-${name}.tre" &
	done
	wait  # Wait to finish.
fi

# Calculate the total number of non-empty tree files.
tre_count=$(find ${outdir} -maxdepth 1 -type f -path "*.tre" | wc -l)
empty_tre_count=$(find ${outdir} -maxdepth 1 -type f -path "*.tre" -empty | wc -l)

# Were tree files generated?
if [ $((tre_count - empty_tre_count)) -lt 1 ]; then
   echo -e "No tree files were generated.\n"
   exit 1
fi

# Root the trees with a custom rooting script (in parallel)
echo -e "Rooting trees with TreeTime...\n"
for seg in "${found_segments[@]}"
do
	treetime-root "${outdir}/${seg}-${name}.tre" "${outdir}/${seg}-${name}.aln" &
done
wait

# Calculate the total number of non-empty rooted tree files.
rooted_count=$(find ${outdir} -maxdepth 1 -type f -path "*.aln.rooted.tre" | wc -l)
empty_rooted_count=$(find ${outdir} -maxdepth 1 -type f -path "*.aln.rooted.tre" -empty | wc -l)

# Were tree files generated?
if [ $((rooted_count - empty_rooted_count)) -lt 1 ]; then
   echo -e "No rooted tree files were generated.\n"
   exit 1
fi

# Create a descriptor file
descriptor=${outdir}/descriptor.csv
for seg in "${found_segments[@]}"
do
	if [ $seg == $ref_seg ]; then
		echo -n "*" >> $descriptor
	fi
	echo "${seg},${outdir}/${seg}-${name}.aln,${outdir}/${seg}-${name}.aln.rooted.tre" >> $descriptor
done
echo -e "The descriptor file was written to ${descriptor}\n"

# Send found_segments to stdout so it can be retrieved by run_treesort.py.
echo -e "FOUND_SEGMENTS: ${found_segments[*]}"

exit 0