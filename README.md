# pt_genome_remapping

Here contains all the files and scripts used to perform remapping of Phatr3 annoation data from the GCA_000150955.2 (Diatom Consortium) genome assembly to the GCA_914521175.1 (Telomere-to-telomere) genome assembly of Phaeodactylum tricornutum.

## Step 1: Clean up data files

### 1. Clean up metadata

- Inputs
  - Uncleaned up metadata in `.gff` format
    - `genome_data/old/uncleaned_metadata_file/phatr3_gene_models_with_href.gff3`
- Script
  - used `sed`

```bash
cat genome_data/old/uncleaned_metadata_file/phatr3_gene_models_with_href.gff3 | sed 's=\\==g; s=\/==g; s=%==g; s=(==g; s=)==g; s=|==g; s=&==g; s={==g; s=}==g; s=<==g; s=>==g; s=>==g; s=\*==g; s=?==g; s=\$==g; s=!==g; s=@==g; s="==g; s= ==g; s=:=_=g; s/=//g; s=,==g; s=;=_=g' > tmp.gff
cat tmp.gff | sed "s='==g" > genome_data/old/cleaned_up_metadata_file/old_annotation_data_fully_cleaned_and_w_unique_lines.gff
```

- Outputs
  - cleaned up metadata in `.gff` format
    - `genome_data/old/cleaned_up_metadata_file/old_annotation_data_fully_cleaned_and_w_unique_lines.gff`

## Step 2: Extract Sequence Data

- Inputs
  - old genome assembly in `.fasta` format
    - `assemblies/old_nuclear_assembly.fasta`
  - cleaned up metadata in `.gff` format
    - `genome_data/old/cleaned_up_metadata_file/old_annotation_data_fully_cleaned_and_w_unique_lines.gff`
- Script
  - use `bedtools getfasta`
  - Need the `–name` and `–s` flags to make sure everything downstream works properly

```bash
bedtools getfasta -name -s -fi assemblies/old_nuclear_assembly.fasta -bed genome_data/old/cleaned_up_metadata_file/old_annotation_data_fully_cleaned_and_w_unique_lines.gff > gff_sequences.fa
```

- Outputs
  - File containing all sequence  in `fasta` format
    - `gff_sequences.fa`

## Step 3: Filter out based off Sequence Data

### 1. Filter out sequences if N's are present

- Inputs
  - Names of `.fasta` files in a list (strip the `.fa` or `.fasta` for this script to work)
    - `gff_sequences`
- Script
  - `./utils/n_check/countNdetailPerFeature.sh`

```bash
# cd into the folder with gff_sequences
source utils/n_check/countNdetailPerFeature.sh gff_sequences
```

- Outputs
  - List of feature names without N's in the sequence (`.txt` format)
    - `gff_sequences_features_with_no_Ns.txt`
  - List of feature names with N's in the sequence (`.txt` format)
    - `gff_sequences_features_with_Ns.txt`
  - List of feature names with N's in the sequence and how many N's they posess (`.txt` format)
    - `gff_sequences_features_withNs.txt`

### 2. Filter `.fasta` file to exclude features with N's

- Inputs
  - file containing all the sequences of genomic features from `.gff` file in `.fasta` format
    - `gff_sequences.fa`
  - list of feature names without N's in the sequence (`.txt` format)
    - `gff_sequences_features_with_no_Ns.txt`
- Script
  - used `seqtk subseq`

```bash
seqtk subseq gff_sequences.fa gff_sequences_features_with_no_Ns.txt > gff_sequences_no_ns.fa
```

- Outputs
  - File containing all the sequences of genomic features from `.gff` file excluding those with N's in `.fasta` format
    - `gff_sequences_no_ns.fa`

### 3. Filter out sequences smaller than 20 bp

- Inputs
  - Names of `.fasta` files in a list (strip the `.fa` or `.fasta` for this script to work)
    - `gff_sequences_no_ns`
- Script
  - `./utils/sequence_length_check/check_seq_length.sh`
    - This script requires `bioawk` to work. Make sure you have `bioawk` installed in your conda environment and modify the script to source that conda environment from within the script.

```bash
# cd into the folder with gff_sequences
source utils/sequence_length_check/check_seq_length.sh gff_sequences_no_ns
```

- Outputs
  - List of feature names without N's in the sequence and that are bigger or equal to 20bp (`.txt` format)
    - `gff_sequences_no_ns_bigger_than_or_eq_to_20_bp.txt`
  - List of feature names without N's in the sequence and that are smaller than 20bp (`.txt` format)
    - `gff_sequences_no_ns_than_20_bp.txt`
  - List of feature names without N's in the sequence and that are smaller than 20bp with their sizes(`.txt` format)
    - `gff_sequences_no_ns_than_20_bp_with_their_sizes.txt`

### 4. Filter `.fasta` file to exclude features that are smaller than 20bp 

- Inputs
  - File containing all the sequences of genomic features from `.gff` file in without N's `.fasta` format
    - `gff_sequences_no_ns.fa`
  - List of feature names without N's in the sequence and that are bigger or equal to 20bp (`.txt` format)
    - `gff_sequences_no_ns_bigger_than_or_eq_to_20_bp.txt`
- Script
  - used `seqtk subseq`

```bash
seqtk subseq gff_sequences_no_ns.fa gff_sequences_no_ns_bigger_than_or_eq_to_20_bp.txt > gff_sequences_no_ns_bigger_than_or_eq_to_20_bp.fa
```

- Outputs
  - File containing all the sequences of genomic features from `.gff` file excluding those with N's and those that are bigger or equal to 20bp in `.fasta` format
    - `gff_sequences_no_ns_bigger_than_or_eq_to_20_bp.fa`

## Step 4: Perform Local Alignments (BLAST)

### 1. Create BLAST database of new genome assembly

- Inputs
  - New genome assembly in `.fasta` format
    - `assemblies/new_nuclear_assembly.fasta`
- Script
  - use `makeblastdb`

```bash
makeblastdb -in assemblies/new_nuclear_assembly.fasta -dbtype nucl
```

- Outputs
  - BLAST database index files

### 2. BLAST filtered features to new assembly

- Inputs
  - New genome assembly in `.fasta` format
    - `assemblies/new_nuclear_assembly.fasta`
  - File containing all the sequences of genomic features from `.gff` file excluding those with N's and those that are bigger or equal to 20bp in `.fasta` format
    - `gff_sequences_no_ns_bigger_than_or_eq_to_20_bp.fa`
- Script
  - used `blastn`
    - Make sure to have `sstrand` included as a column because that takes in stand info

```bash
blastn -db assemblies/new_nuclear_assembly.fasta -query gff_sequences_no_ns_bigger_than_or_eq_to_20_bp.fa -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand” -out blast_outputs/features_to_new_filtered_blast.tsv
```

- Outputs
  - Table with BLAST outputs in `.tsv` format
    - `blast_outputs/features_to_new_filtered_blast.tsv`

### 3. Append Row Numbers to BLAST output

- Inputs
  - Table with BLAST outputs in `.tsv` format
    - `blast_outputs/features_to_new_filtered_blast.tsv`
- Script
  - used `awk`

```bash
awk 'BEGIN{OFS="\t"} {print $0, NR-1}' blast_outputs/features_to_new_filtered_blast.tsv > blast_outputs/features_to_new_filtered_blast_rows_appended.tsv
```

- Outputs
  - Table with BLAST outputs and row numbers appended in `.tsv` format
    - `blast_outputs/features_to_new_filtered_blast_rows_appended.tsv`
    - Note: Rows are 0-indexed

## Step 5: Perform Global Alignment (Needle)

### 1. Create a `.bed` file from the BLAST results to generate a `.fasta` file for the global alignment

- Inputs
  - Table with BLAST outputs and row numbers appended in `.tsv` format
    - `blast_outputs/features_to_new_filtered_blast_rows_appended.tsv`
- Script
  - `utils/format/blastToBedToolsInput.py`

```bash
python3 blastToBedToolsInput.py -bo blast_outputs/features_to_new_filtered_blast_rows_appended.tsv -of oldGenesBedToolsInput.txt
```

- Outputs
  - BLAST outputs in `.bed` format
    - `features_to_new_filtered_blast.bed`

### 2. Create a `.fasta` file from the BLAST results for the global alignment

- Inputs
  - New genome assembly in `.fasta` format
    - `assemblies/new_nuclear_assembly.fasta`
  - BLAST outputs in `.bed` format
    - `features_to_new_filtered_blast.bed`
- Script
  - used `bedtools getfasta`

```bash
bedtools getfasta -fi assemblies/new_nuclear_assembly.fasta -bed features_to_new_filtered_blast.bed -name -s -fo features_to_new_filtered_blast.fasta
```

- Outputs
  - BLAST outputs in `.fasta` format
    - `features_to_new_filtered_blast.fasta`

### 3. Create an intermediate feature file (needed for global alignment script to work)

- Inputs
  - BLAST outputs in `.fasta` format
    - `features_to_new_filtered_blast.fasta`
- Script
  - used `awk`

```bash
awk -F "\n" 'BEGIN{RS=">"} {print $1}' features_to_new_filtered_blast.fasta > features_to_new_filtered_blast_list.txt
```

- Outputs
  - List of header names from `.fasta` file (`.txt` file)
    - `features_to_new_filtered_blast_list.txt`
  
### 4. Perform Global Alignment

- Inputs
  - File with BLAST queries sequences in `.fasta` format
    - `gff_sequences_no_ns_bigger_than_or_eq_to_20_bp.fa`
  - File with BLAST subject sequences in `.fasta` format
    - `features_to_new_filtered_blast.fasta`
  - List of header names from `.fasta` file (`.txt` file)
    - `features_to_new_filtered_blast_list.txt`
  - Table with BLAST outputs and row numbers appended in `.tsv` format
    - `blast_outputs/features_to_new_filtered_blast_rows_appended.tsv`
- Script
  - `utils/global_alignment/needleAlignBlastUpdaterV2.py`

```bash
python3 utils/global_alignment/needleAlignBlastUpdaterV2.py -nf features_to_new_filtered_blast_list.txt -gf gff_sequences_no_ns_bigger_than_or_eq_to_20_bp.fa -sf features_to_new_filtered_blast.fasta -bf blast_outputs/features_to_new_filtered_blast_rows_appended.tsv -obf blast_outputs/all_blasted_to_new_assembly_POST_GLOBAL_ALIGNMENT.tsv
```

- Outputs
  - a `.tsv` file with all the blast info post global alignment
    - `blast_outputs/all_blasted_to_new_assembly_POST_GLOBAL_ALIGNMENT.tsv`

## Step 6: Perform Assembly to Assembly Alignment

### 1. Alignment

- Inputs
  - new assembly without chloroplast and mitochondria in `.fasta` format
    - `assemblies/new_nuclear_assembly.fasta`
  - old assembly without chloroplast and mitochondria in `.fasta` format
    - `assemblies/old_nuclear_assembly.fasta`
- Script
  - `minimap2` alignment of old to new alignment
  - `-x asm5` flag specified to optimize aligner for full genome/assembly alignment

```bash
minimap2 -x asm5 assemblies/new_nuclear_assembly.fasta assemblies/old_nuclear_assembly.fasta > minimap_outputs/oldToNewMap.paf
```

- Outputs
  - old assembly to new assembly alignment in `.paf` format
    - `minimap_outputs/oldToNewMap.paf`

### 2. Pre-Filter the Old Assembly to New Assembly Alignment

- Remove all the small and low scoring alignments.
- A cut off of **450bp** was chosen as the minimum alignment length. This is because the smallest scaffold in the old assembly with annotated data is 450bp.
- A cutoff score of **20** was used for my alignments. `.paf` scores use Phred-scaled quality scores. Generally speaking, a Phred score of 20 or above is acceptable, because this means that whatever it qualifies is >99% accurate, with a <1% chance of error.

- Inputs
  - old assembly to new assembly alignment in in `.paf` format
    - `minimap_outputs/oldToNewMap.paf`
- Script
  - `utils/pre-filter/remove_by_size_and_score.sh`
    - Take in 2 positional arguments
      - 1. The minimum size of mapping you want to allow
      - 2. The minimum quality score you want to allow

```bash
./utils/pre-filter/remove_by_size_and_score.sh  minimap_outputs/oldToNewMap.paf 450 20 > filtered_alignments/pre-filtered/oldToNewMap-450size-20score.paf
```

- Outputs
  - old assembly to new assembly alignment mappings >450bp and >20 Q-score in `.paf` format
    - `filtered_alignments/pre-filtered/oldToNewMap-450size-20score.paf`

### 3. Get chromosome sizes of new assembly nuclear genome

Needed for the following filtering script

- Inputs
  - new assembly without chloroplast and mitochondria in `.fasta` format
    - `assemblies/new_nuclear_assembly.fasta`
- Script
  - `utils/chr_sizes/get_sizes.sh`
    - Make sure you have `bioawk` in your environment before running this script

```bash
source ./utils/chr_sizes/get_sizes.sh assemblies/new_nuclear_assembly.fasta > ./chr_sizes/new_nuclear_assembly_chr_sizes.txt
```

- Outputs
  - Sequences sizes of all the chromosomes in the new assembly (and `chromosome_` stripped from chromosome name); in a key, value format.
    - Col 1 = chromosome number
    - Col 2 = chromosome size
    - `chr_sizes/new_nuclear_assembly_chr_sizes.txt`

### 4. Remove Gaps from alignment

- Inputs
  - old assembly to new assembly alignment mappings >450bp and >20 Q-score in `.paf` format
    - `filtered_alignments/pre-filtered/oldToNewMap-450size-20score.paf`
  - Sequences sizes of all the chromosomes in the new assembly (and `chromosome_` stripped from chromosome name); in a key, value format.
    - Col 1 = chromosome number
    - Col 2 = chromosome size
    - `chr_sizes/new_nuclear_assembly_chr_sizes.txt`
  - Your output file name in `.paf` format
  - The percentage of the size of the mapping that you're willing to allow overlap with another plasmid `a number ranging from 0-1`
  - `--noSave` if you don't want to save intermediate tables
  - `--noLoad` if you don't want to load intermediate tables
- Script
  - `utils/greedy_lenient_nogaps/filterPaf_lenient_nogaps.py`

```bash
python3 utils/greedy_lenient_nogaps/filterPaf_lenient_nogaps.py -f1 filtered_alignments/pre-filtered/oldToNewMap-450size-20score.paf -gs chr_sizes/new_nuclear_assembly_chr_sizes.txt -of filtered_alignments/greedy-lenient-nogaps/oldToNewMap-450size-20score_20pct-overlap_nogaps.paf --sensitivity 0.20 --noLoad
```

- Outputs
  - old assembly to new assembly alignment mappings >450bp and >20 Q-score and updated mapping coordiantes to ensure there are no gaps in `.paf` format
    - `filtered_alignments/greedy-lenient-nogaps/oldToNewMap-450size-20score_20pct-overlap_nogaps.paf`

## Step 7: Figure out how many remaps you'd expect from the old assembly to new assembly remapping

### 1. Update the old GFF file with the information from the assembly to assembly mapping

For each entry in the old GFF, check if its coordinates fall within a region that got remapped to the filtered old to new assembly mapping. Append which mapping it got mapped to and where it should fall on the new assembly

- Inputs
  - old assembly to new assembly alignment mappings >450bp and >20 Q-score and updated mapping coordiantes to ensure there are no gaps in `.paf` format
    - `filtered_alignments/greedy-lenient-nogaps/oldToNewMap-450size-20score_20pct-overlap_nogaps.paf`
  - genome data file in `.gff` format with the length of the feature calculated and appended to `column 13`
    - `genome_data/old/cleaned_up_metadata_file/old_annotation_data_fully_cleaned_and_w_unique_lines_and_feature_sizes_in_last_col.gff`
- Script
  - `utils/extract_gffs_from_mappings/extract_gffs_from_mappings.py`

```bash
python3 utils/extract_gffs_from_mappings/extract_gffs_from_mappings.py --pafFile filtered_alignments/greedy-lenient-nogaps/oldToNewMap-450size-20score_20pct-overlap_nogaps.paf --gffFile genome_data/old/cleaned_up_metadata_file/old_annotation_data_fully_cleaned_and_w_unique_lines_and_feature_sizes_in_last_col.gff --outputFile genome_data/filtered_by_mapping/old_annotation_data_fully_cleaned_and_w_unique_lines_and_mapping_info_added-450size-20score_20pct-overlap_nogaps.gff
```

- Outputs
  - genome data file in `.gff` format with the length of the feature calculated and appended to `column 13`
  - columns list is `["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes", "featureLength", "line_id_on_mapping_file", "target_chrom", "mapping_start_on_old", "mapping_end_on_old"]`
    - `genome_data/filtered_by_mapping/old_annotation_data_fully_cleaned_and_w_unique_lines_and_mapping_info_added-450size-20score_20pct-overlap_nogaps.gff`

### 2. Calculate the number of expected remaps expected in each chromosome

- Inputs

  - genome data file in `.gff` format with the length of the feature calculated and appended to `column 13`
  - columns list is `["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes", "featureLength", "line_id_on_mapping_file", "target_chrom", "mapping_start_on_old", "mapping_end_on_old"]`
    - `genome_data/filtered_by_mapping/old_annotation_data_fully_cleaned_and_w_unique_lines_and_mapping_info_added-450size-20score_20pct-overlap_nogaps.gff`

- Script

```bash
awk '{print $12}' genome_data/filtered_by_mapping/old_annotation_data_fully_cleaned_and_w_unique_lines_and_mapping_info_added-450size-20score_20pct-overlap_nogaps.gff | sort | uniq -c | awk 'BEGIN {OFS="\t"} {if ($2 ~ /^chromosome_*/) print $2, $1}' | sort -n -t "_" -k 2 | awk 'BEGIN {OFS="\t"} {total += $2; print $0} END {print "total", total}' > genome_data/new/expected_nuclear_remaps_summary.txt
```

- Outputs
  - Expected number of features aligned to each new chromosome simply based off where that feature falls within the old assembly to new assembly alignment mappings; in a key, value format.
    - Col 1 = New chromosome
    - Col 2 = Expected number of features that remapped
    - `genome_data/new/expected_nuclear_remaps_summary.txt`

## Step 8: Update the BLAST file with the all the necessary info for filtering

Want to create a file with all the BLAST info plus the global alignment info plus the original info from the gff plus the info from the assembly to assembly mapping in order to begin performing reciprical best hit and figuring out which BLAST hits are accurate remaps of old features

- Inputs
  - a `.tsv` file with all the blast info post global alignment
    - `blast_outputs/all_blasted_to_new_assembly_POST_GLOBAL_ALIGNMENT.tsv`

```bash
# columns list is ["qseqid", #1, "sseqid", #2, "%_identity", #3, "alignment_length", #4, "mismatch", #5, "gapopen", #6, "qstart", #7, "qend", #8, "sstart", #9, "send", #10, "evalue", #11, "bitscore", #12, "subject_strand", #13, "line_in_og_BLAST", #14, "Needle_score", #15, "First_10_bp_matching", #16, "Last_10_bp_matching", #17, "First_9_bp_matching", #18, "Last_9_bp_matching", #19, "First_8_bp_matching", #20, "Last_8_bp_matching", #21, "First_7_bp_matching", #22, "Last_7_bp_matching", #23, "First_6_bp_matching", #24, "Last_6_bp_matching", #25, "First_5_bp_matching", #26, "Last_5_bp_matching", #27, "First_4_bp_matching", #28, "Last_4_bp_matching", #29, "First_3_bp_matching", #30, "Last_3_bp_matching", #31, "First_2_bp_matching", #32, "Last_2_bp_matching", #33, "First_1_bp_matching", #34, "Last_1_bp_matching", #35]
```

- .
  - genome data file in `.gff` format with the length of the feature calculated and appended to `column 13`
    - `genome_data/filtered_by_mapping/old_annotation_data_fully_cleaned_and_w_unique_lines_and_mapping_info_added-450size-20score_20pct-overlap_nogaps.gff`

```bash
# columns list is ["seqid" #1, "source" #2, "type" #3, "start" #4, "end" #5, "score" #6, "strand" #7, "phase" #8, "attributes" #9, "featureLength" #10, "line_id_on_mapping_file" #11, "target_chrom" #12, "mapping_start_on_old" #13, "mapping_end_on_old" #14, "if_query_and_target_on_same_strand" #15; ‘+’ if query/target on the same strand; ‘-’ if opposite, "mapping_start_on_new" #16, "mapping_end_on_new" #17, "mapping_q_score" #18]
```

- Script
  - `utils/add_gff_info_to_blast/pull_info_from_gff_into_blast.py`
  - `--noLoad` if you don't want to load intermediate tables

```bash
nohup python3 -u pull_info_from_gff_into_blast.py --blastFile blast_outputs/all_blasted_to_new_assembly_POST_GLOBAL_ALIGNMENT.tsv --gffFile old_annotation_data_fully_cleaned_and_w_unique_lines_and_mapping_info_added-450size-20score_20pct-overlap_nogaps.gff --outBlastFile blast_outputs/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL.tsv --noLoad > all.log 2> all.errlog &
```

- Outputs
  - A `.tsv` file with all the BLAST info plus the global alignment info plus the original info from the gff plus the info from the assembly to assembly mapping
    - `blast_outputs/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL.tsv`

```bash
# columns list is ["qseqid", #1, "sseqid", #2, "%_identity", #3, "alignment_length", #4, "mismatch", #5, "gapopen", #6, "qstart", #7, "qend", #8, "sstart", #9, "send", #10, "evalue", #11, "bitscore", #12, "subject_strand", #13, "line_in_og_BLAST", #14, "Needle_score", #15, "First_10_bp_matching", #16, "Last_10_bp_matching", #17, "First_9_bp_matching", #18, "Last_9_bp_matching", #19, "First_8_bp_matching", #20, "Last_8_bp_matching", #21, "First_7_bp_matching", #22, "Last_7_bp_matching", #23, "First_6_bp_matching", #24, "Last_6_bp_matching", #25, "First_5_bp_matching", #26, "Last_5_bp_matching", #27, "First_4_bp_matching", #28, "Last_4_bp_matching", #29, "First_3_bp_matching", #30, "Last_3_bp_matching", #31, "First_2_bp_matching", #32, "Last_2_bp_matching", #33, "First_1_bp_matching", #34, "Last_1_bp_matching", #35, "query_feature_length", #36, "query_feature_chromosome", #37, "query_feature_start", #38, "query_feature_end", #39, "query_feature_strand", #40, "query_feature_source", #41, "query_feature_type", #42, "query_feature_score", #43, "query_feature_phase", #44, "query_feature_attributes_cleaned", #45, "mapping_row_in_non-sorted_paf_file_0_based", #46, "expected_chromosome_to_map_to", #47, "expected_start_position_to_map_to", #48, "expected_end_position_to_map_to", #49, "subject_length" #50]
```

## Step 9: Filter the BLAST file to perform reciprical best hit and figuring out which BLAST hits are accurate remaps of old features

### 1. Prefilter the BLAST file

- Inputs
  - A `.tsv` file with all the BLAST info plus the global alignment info plus the original info from the gff plus the info from the assembly to assembly mapping
    - `blast_outputs/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL.tsv`

```bash
# columns list is ["qseqid", #1, "sseqid", #2, "%_identity", #3, "alignment_length", #4, "mismatch", #5, "gapopen", #6, "qstart", #7, "qend", #8, "sstart", #9, "send", #10, "evalue", #11, "bitscore", #12, "subject_strand", #13, "line_in_og_BLAST", #14, "Needle_score", #15, "First_10_bp_matching", #16, "Last_10_bp_matching", #17, "First_9_bp_matching", #18, "Last_9_bp_matching", #19, "First_8_bp_matching", #20, "Last_8_bp_matching", #21, "First_7_bp_matching", #22, "Last_7_bp_matching", #23, "First_6_bp_matching", #24, "Last_6_bp_matching", #25, "First_5_bp_matching", #26, "Last_5_bp_matching", #27, "First_4_bp_matching", #28, "Last_4_bp_matching", #29, "First_3_bp_matching", #30, "Last_3_bp_matching", #31, "First_2_bp_matching", #32, "Last_2_bp_matching", #33, "First_1_bp_matching", #34, "Last_1_bp_matching", #35, "query_feature_length", #36, "query_feature_chromosome", #37, "query_feature_start", #38, "query_feature_end", #39, "query_feature_strand", #40, "query_feature_source", #41, "query_feature_type", #42, "query_feature_score", #43, "query_feature_phase", #44, "query_feature_attributes_cleaned", #45, "mapping_row_in_sorted_paf_file_0_based", #46, "expected_chromosome_to_map_to", #47, "expected_start_position_to_map_to", #48, "expected_end_position_to_map_to", #49, "subject_length" #50]
```

- Script
  - `utils/filter_blast_remaps/pre_filter_blast_remaps.sh`

```bash
./utils/filter_blast_remaps/pre_filter_blast_remaps.sh blast_outputs/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL.tsv
```

- Outputs
  - A `.tsv` file with all the BLAST info plus the global alignment info plus the original info from the gff plus the info from the assembly to assembly mapping, filtered according to the following pre-filters
    - Pre-Filter 1: Remove all entries where the first base pair of the subject sequence is not the same as the last base pair of the query sequence
    - Pre-Filter 2: Remove all entries where the last base pair of the subject sequence is not the same as the last base pair of the query sequence
    - Pre-Filter 3: Remove all entries where the first 3 base pairs of the subject sequence is not the same as the first base pair of the query sequence if the sequence is a CDS or exon
    - Pre-Filter 4: Remove all entries where the last 3 base pairs of the subject sequence is not the same as the last base pair of the query sequence if the sequence is a CDS or exon
    - Pre-Filter 5: Remove all entries where the first 8 base pairs of the subject sequence are not the same as the first 8 base pairs of the query sequence
    - Pre-Filter 6: Remove all entries where the last 8/10 base pairs of the subject sequence are not the same as the last 8 base pairs of the query sequence
    - Pre-Filter 7: Remove all entries where the length of the subject sequence is 50% smaller than the length of the query sequence
    - Pre-Filter 8: Remove all entries where the length of the subject sequence is 50% larger than the length of the query sequence
    - Pre-Filter 9: Remove all entries where the length of the subject sequence is 10% smaller than the length of the query sequence if the query is a CDS or exon
    - Pre-Filter 10: Remove all entries where the length of the subject sequence is 10% larger than the length of the query sequence if the query is a CDS or exon
    - `blast_outputs/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_and_PRE_FILTERED.tsv`

### 2. Filter the BLAST file to make sure the BLAST mappings fall within a region deemed acceptable by the old assembly to new assembly alignment

- Inputs

  - A prefiltered`.tsv` file with all the BLAST info plus the global alignment info plus the original info from the gff plus the info from the assembly to assembly mapping
    - `blast_outputs/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_and_PRE_FILTERED.tsv`
  - old assembly to new assembly alignment mappings >450bp and >20 Q-score and updated mapping coordiantes to ensure there are no gaps in `.paf` format
    - `filtered_alignments/greedy-lenient-nogaps/oldToNewMap-450size-20score_20pct-overlap_nogaps.paf`
  - `--noLoad` if you don't want to load intermediate tables

- Script
  - `utils/filter_blast_remaps/filter_blast_remaps.py`

```bash
python3 utils/filter_blast_remaps/filter_blast_remaps.py --blastFile blast_outputs/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_and_PRE_FILTERED.tsv --pafFile filtered_alignments/greedy-lenient-nogaps/oldToNewMap-450size-20score_20pct-overlap_nogaps.paf --outFile ./blast_outputs/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_and_FILTERED_POST_SYNTENY_CHECK.tsv --noLoad
```

- Outputs
  - A filtered`.tsv` file with all the BLAST info plus the global alignment info plus the original info from the gff plus the info from the assembly to assembly mapping, filtered to remove blast hits that don't align with old assembly to new assembly alignment synteny
    - `./blast_outputs/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_and_FILTERED_POST_SYNTENY_CHECK.tsv`

### 5. Filter the BLAST file to only hits with the highest Global Alignment Score (Reciprical best hit)

<!-- TODO make this script have proper inputs and outputs. At the moment it is just run by itself -->

- Inputs
  - A filtered`.tsv` file with all the BLAST info plus the global alignment info plus the original info from the gff plus the info from the assembly to assembly mapping, filtered to remove blast hits that don't align with old assembly to new assembly alignment synteny
    - `./blast_outputs/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_and_FILTERED_POST_SYNTENY_CHECK.tsv`
  - .
    - ``
- Script

  - `utils/filter_blast_remaps/filter_blast_remaps.Rmd`
    - Run in R Studio

- Outputs
  - A `.gff` file of all the features that remapped (without their coordinates updated)
    - `genome_data/new/remapped_features_COORDS_NOT_UPDATED.gff`
  - A `.gff` file of all the features that were found in the blast outputs but that did not successfully remap (without their coordinates updated)
    - `genome_data/new/genome_data/new/non_remapped_features_FROM_BLAST_OUTPUT_ONLY.gff`
  - A `.gff` file of all the features that remapped, merged with all the data at from the blast table (merged by the id column)
    - `genome_data/new/merged_remapped_table_with_BLAST_info.gff`

```bash
# columns list is ["seqid", #1, "source", #2, "type", #3, "start", #4, "end", #5, "score", #6, "strand", #7, "phase", #8, "q_metadata", #9, "qseqid", #10, "sseqid", #11, "pident", #12, "length", #13, "mismatch", #14, "gapopen", #15, "qstart", #16, "qend", #17, "sstart", #18, "send", #19, "evalue", #20, "bitscore", #21, "sstrand", #22, "order", #23, "needle_score", #24, "first_10", #25, "last_10", #26, "first_9", #27, "last_9", #28, "first_8", #29, "last_8", #30, "first_7", #31, "last_7", #32, "first_6", #33, "last_6", #34, "first_5", #35, "last_5", #36, "first_4", #37, "last_4", #38, "first_3", #39, "last_3", #40, "first_2", #41, "last_2", #42, "first_1", #43, "last_1", #44, "q_length", #45, "q_chr", #46, "q_old_start", #47, "q_old_end", #48, "q_old_strand", #49, "q_old_source", #50, "q_old_type", #51, "q_old_score", #52, "q_old_phase", #53, "mapping_row_in_paf_0i", #54, "expected_chromosome", #55, "mapping_start_on_old_chromosome", #56, "mapping_end_on_old_chromosome", #57, "if_query_and_target_on_same_strand", #58, "mapping_start_on_new_chromosome", #59, "mapping_end_on_new_chromosome", #60, "mapping_q_score", #61, "s_length", #62, "in_acceptable_region", #63, "orientation_check" #64]
```

## Step 10: Update and clean up the remapped GFF file

Update and clean up the remapped GFF file so that each remapped feature has its new chromosome and new start and end coordinates corresponding with the new assembly. Also order the file first by the chromosome number, and then by the start coordinates of the features (ascending order). Then append all the chloroplast and mitochondrial features as well and then finally replace the cleaned up metadata with the original ones.

### 1. Update the chromosome name and the coordinates of the features in the new gff file

- Inputs
  - A `.gff` file of all the features that remapped, merged with all the data at from the blast table (merged by the id column)
    - `genome_data/new/merged_remapped_table_with_BLAST_info.gff`
- Script
  - `utils/update_remapping_coordinates/update_remapping_deails.sh`

```bash
./utils/update_remapping_coordinates/update_remapping_deails.sh genome_data/new/merged_remapped_table_with_BLAST_info.gff
```

- Outputs
  - A `.gff` file with the successfully remapped features updated with their new chromosomes and new start and end coordinates corresponding with the new assembly.
    - `genome_data/new/final_remapped_data_updated_coordinates_minus_cp_and_mt_minus_original_metadata.gff`
  - Also outputs an intermediate `.gff` file with all the blast info still appended
    - `genome_data/new/final_remapped_data_updated_coordinates_plus_all_blast_data_minus_cp_and_mt_minus_original_metadata.gff`

### 2. Sort the file first by the chromosome number, and then by the start coordinates of the features (ascending order)

- Inputs
  - A `.gff` file with the successfully remapped features updated with their new chromosomes and new start and end coordinates corresponding with the new assembly.
    - `genome_data/new/final_remapped_data_updated_coordinates_minus_cp_and_mt_minus_original_metadata.gff`
- Script
  - `utils/update_remapping_coordinates/sort_gff.R`

```bash
./utils/update_remapping_coordinates/sort_gff.R genome_data/new/final_remapped_data_updated_coordinates_minus_cp_and_mt_minus_original_metadata.gff
```

- Outputs
  - A `.gff` file with the successfully remapped features updated with their new chromosomes and new start and end coordinates corresponding with the new assembly, sorted by the chromosome number, and then by the start coordinates of the features in ascending order
    - `genome_data/new/final_remapped_data_updated_coordinates_sorted_minus_cp_and_mt_minus_original_metadata.gff`

### 3. Add Mitochondria and Chloroplast data

- Inputs

  - A `.gff` file with the successfully remapped features updated with their new chromosomes and new start and end coordinates corresponding with the new assembly, sorted by the chromosome number, and then by the start coordinates of the features in ascending order
    - `genome_data/new/final_remapped_data_updated_coordinates_sorted_minus_cp_and_mt_minus_original_metadata.gff`
  - A `.gff` file containing all the features from the mitochondrial genome (needs to have cleaned up metadata)
    - `genome_data/old/cleaned_up_metadata_file/mt_annotation_data_fully_cleaned_and_w_unique_lines.gff`
  - A `.gff` file containing all the features from the chloroplast genome (needs to have cleaned up metadata)
    - `genome_data/old/cleaned_up_metadata_file/cp_annotation_data_fully_cleaned_and_w_unique_lines.gff`

- Script
  - `utils/update_remapping_coordinates/add_mt_and_cp_details.sh`
    - Takes in 3 positional arguments
      - 1. The `.gff` file to be updated
      - 2. The `.gff` file with the mitochondria annotation data
      - 3. The `.gff` file with the chloroplast annotation data

```bash
./utils/update_remapping_coordinates/add_mt_and_cp_details.sh genome_data/new/final_remapped_data_updated_coordinates_sorted_minus_cp_and_mt_minus_original_metadata.gff genome_data/old/cleaned_up_metadata_file/mt_annotation_data_fully_cleaned_and_w_unique_lines.gff genome_data/old/cleaned_up_metadata_file/cp_annotation_data_fully_cleaned_and_w_unique_lines.gff
```

- Outputs
  - A `.gff` file with the successfully remapped features updated with their new chromosomes and new start and end coordinates corresponding with the new assembly, sorted by the chromosome number, and then by the start coordinates of the features in ascending order, updated with the mitochondria and chloroplast annotation data
    - `genome_data/new/final_remapped_data_updated_coordinates_sorted_plus_cp_and_mt_minus_original_metadata.gff`

### 4. Update the Metadata

#### 4.a. Create a key value pair file for old metadata (unique IDs to what they should be replaced with)

```bash
# Remove headers from the old features file
awk 'BEGIN {FS = "\t"} {if (NF == 9) print $0}' genome_data/old/uncleaned_metadata_file/phatr3_gene_models_with_href.gff3 > genome_data/old/uncleaned_metadata_file/phatr3_gene_models_with_href_HEADERS_REMOVED.gff3

# Create key value pairs of cleaned up IDS to what they should be replaced with at the end
paste cleaned_up_metadata_file/old_annotation_data_fully_cleaned_and_w_unique_lines.gff genome_data/old/uncleaned_metadata_file/phatr3_gene_models_with_href_HEADERS_REMOVED.gff3 | awk 'BEGIN {FS = "\t"; OFS = "\t"} {print $9, $19}' > txts/new_id_to_old_id_kv_pairs.txt
```

#### 4.b. Replace the cleaned up metadata with the original ones

- Inputs
  - A `.gff` file with the successfully remapped features updated with their new chromosomes and new start and end coordinates corresponding with the new assembly, sorted by the chromosome number, and then by the start coordinates of the features in ascending order, updated with the mitochondria and chloroplast annotation data
    - `genome_data/new/final_remapped_data_updated_coordinates_sorted_plus_cp_and_mt_minus_original_metadata.gff`
  - A key value pair file for old metadata
    - Keys: unique IDs
    - Values: original metadata that the unique ids should be replaced with
    - `txts/new_id_to_old_id_kv_pairs.txt`
- Script
  - `utils/update_remapping_coordinates/revert_to_original_metadata.py`

```bash
python3 utils/update_remapping_coordinates/revert_to_original_metadata.py --gffFile genome_data/new/final_remapped_data_updated_coordinates_sorted_plus_cp_and_mt_minus_original_metadata.gff --metaTable txts/new_id_to_old_id_kv_pairs.txt --outputGffFile genome_data/new/final_remapped_data_updated_coordinates_sorted_plus_cp_and_mt_plus_original_metadata.gff
```

- Outputs
  - A `.gff` file with the successfully remapped features updated with their new chromosomes and new start and end coordinates corresponding with the new assembly, sorted by the chromosome number, and then by the start coordinates of the features in ascending order, updated with the mitochondria and chloroplast annotation data, and with original metadata
    - `genome_data/new/final_remapped_data_updated_coordinates_sorted_plus_cp_and_mt_plus_original_metadata.gff`

### 5. Add whole chromosome data

Manually add whole chromosome data to the correct lines of the remapped `.gff` (before the features start for the given chromsome). The start coordinate is 1 and the end coordinate is the length of the chromosome, and the strand is "."

- Inputs

  - A `.gff` file with the successfully remapped features updated with their new chromosomes and new start and end coordinates corresponding with the new assembly, sorted by the chromosome number, and then by the start coordinates of the features in ascending order, updated with the mitochondria and chloroplast annotation data, and with original metadata
    - `genome_data/new/final_remapped_data_updated_coordinates_sorted_plus_cp_and_mt_plus_original_metadata.gff`

- Script
  - `utils/update_remapping_coordinates/revert_to_original_metadata.py`

```bash
# For gff with original metadata
python3 utils/update_remapping_coordinates/add_whole_chromosome_data.py --gffFile genome_data/new/final_remapped_data_updated_coordinates_sorted_plus_cp_and_mt_plus_original_metadata.gff --outputFile genome_data/new/final_remapped_data.gff
# For gff with cleaned up metadata
python3 utils/update_remapping_coordinates/add_whole_chromosome_data.py --gffFile genome_data/new/final_remapped_data_updated_coordinates_sorted_plus_cp_and_mt_minus_original_metadata.gff --outputFile genome_data/new/final_remapped_data_minus_original_metadata.gff
```

- Outputs
  - A `.gff` file with the successfully remapped features updated with their new chromosomes and new start and end coordinates corresponding with the new assembly, sorted by the chromosome number, and then by the start coordinates of the features in ascending order, updated with the mitochondria and chloroplast annotation data, with original metadata, and with chromosome information
    - `genome_data/new/final_remapped_data.gff`
  - If you choose to keep the cleaned up metadata, can be found here
    - `genome_data/new/final_remapped_data_minus_original_metadata.gff`
