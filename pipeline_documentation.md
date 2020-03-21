# Rhea workflow documentation
This document describes the individual rules of the pipeline for information purposes. For instructions on installation and usage please refer to the [README](README.md).    

## Overview
### General
* read samples table
* create log directories
* **create_index_star**
* **extract_transcriptome**
* **create_index_salmon**
* **create_index_kallisto**
* **extract_transcripts_as_bed12**
* **index_genomic_alignment_samtools**
* **star_rpm**
* **rename_star_rpm_for_alfa**
* **calculate_TIN_scores**
* **merge_TIN_scores**
* **plot_TIN_scores**
* **salmon_quantmerge_genes**
* **salmon_quantmerge_transcripts**
* **generate_alfa_index**
* **alfa_qc**
* **alfa_qc_all_samples**
* **alfa_concat_results**
* **prepare_files_for_report**
* **prepare_MultiQC_config**
* **MULTIQC_report**

### Sequencing mode specific
* **(pe_)fastqc**
* **(pe_)remove_adapters_cutadapt**
* **(pe_)remove_polya_cutadapt**
* **(pe_)map_genome_star**
* **(pe_)quantification_salmon**
* **(pe_)genome_quantification_kallisto**




## Detailed description of steps
The pipeline consists of three snakefiles: A main Snakefile and an individual Snakefile for each sequencing mode (single-end and paired-end), as parameters to individual tools differ between the sequencing modes. The main Snakefile contains some general rules for the creation of indices, rules that are applicable to both sequencing modes, and rules that deal with summary steps and combining results across samples of the run.     
Individual rules of the pipeline are described briefly, and links to the respective software manuals are given. If parameters can be customised by the user (via the samples table) they are also described.
Description of paired and single-end rules are combined, only differences are highlighted.


### General
#### read samples table
 
Requirements: 
* tab separated file
* first row has to contain parameter names as in  [samples.tsv](tests/input_files/samples.tsv)
* First column will be used as indices (sample identifiers)

Parameter name | Description
--- | ---
sample | descriptive sample name (type=STRING)
seqmode | "paired_end" or "single_end" (type=STRING)
fq1 | PATH/TO/INPUT_FILE.mate_1.fastq.gz (type=STRING)
fq2 | PATH/TO/INPUT_FILE.mate_2.fastq.gz (type=STRING)
fq1_3p | 3' adapter of mate 1 (use "XXXXXXXXXXXXXXX" if none); for cutadapt (type=STRING)
fq1_5p | 5' adapter of mate 1 (use "XXXXXXXXXXXXXXX" if none); for cutadapt (type=STRING)
fq2_3p | 3' adapter of mate 2 (use "XXXXXXXXXXXXXXX" if none); for cutadapt (type=STRING)
fq2_5p | 5' adapter of mate 2 (use "XXXXXXXXXXXXXXX" if none); for cutadapt (type=STRING)
fq1_polya3p | stretch of As or Ts, depending on read orientation (use "XXXXXXXXXXXXXXX" if none), trimmed from the 3' end of the read; for cutadapt (type=STRING)
fq1_polya5p | stretch of As or Ts, depending on read orientation (use "XXXXXXXXXXXXXXX" if none), trimmed from the 5' end of the read; for cutadapt (type=STRING)
fq2_polya3p| stretch of As or Ts, depending on read orientation (use "XXXXXXXXXXXXXXX" if none), trimmed from the 3' end of the read; for cutadapt (type=STRING)
fq2_polya5p| stretch of As or Ts, depending on read orientation (use "XXXXXXXXXXXXXXX" if none), trimmed from the 5' end of the read; for cutadapt (type=STRING)
index_size | ideally `max(ReadLength)-1`. Required for STAR mapping and index creation
kmer | 31; for salmon index creation. A default value of 31 seems to work for reads of 75 bp or longer. If you get poor mappings, consider smaller values for kmer. (type=STRING or type=INT)
organism | format has to correspond to the naming of provided genome files and directories, like "ORGANISM" in the path below. Use e.g. "homo_sapiens" (type=STRING)
gtf | PATH/TO/ORGANISM/annotation.gtf; for star index (type=STRING)
gtf_filtered | PATH/TO/ORGANISM/annotation.gtf; for salmon quantification (type=STRING)
genome | PATH/TO/ORGANISM/genome.fa; for star index (type=STRING)
tr_fasta_filtered | PATH/TO/ORGANISM/transcriptome.fa, for salmon and kallisto index creation (type=STRING)
sd | Estimated standard deviation of fragment length; for single-end kallisto quantification (type=STRING or type=INT)
mean | Estimated average fragment length; for single-end kallisto quantification (type=STRING or type=INT)
multimappers | max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped; for star mapping (type=STRING or type=INT)
soft_clip | "Local": standard local alignment with soft-clipping allowed. "EndToEnd": force end-to-end read alignment, do not soft-clip; for star mapping (type=STRING)
pass_mode | "None": 1-pass mapping; "Basic": basic 2-pass mapping, with all 1st pass junctions inserted into the genome indices on the fly; for star mapping (type=STRING)
libtype | "A": automatically infer. For more info see [salmon manual](https://salmon.readthedocs.io/en/latest/salmon.html) (type=STRING)
kallisto_directionality | "--fr-stranded":Strand specific reads, first read forward. "--rf-stranded": Strand specific reads, first read reverse; for kallisto (type=STRING)

#### create log directories
Currently not implemented as Snakemake rule, but general statement.


#### create_index_star
Create index for STAR alignments. Supply the reference genome sequences (FASTA files) and annotations (GTF file), from which STAR generates genome indexes that are utilized in the 2nd (mapping) step. The genome indexes are saved to disk and are only be generated once for each genome/annotation/index size combination. [STAR manual](http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf#section.2)    

**Input:** genome fasta file, gtf file    
**Parameters:** sjdbOverhang (This is the `index_size` specified in the samples table).    
**Output:** chrNameLength.txt will be used for STAR mapping; chrName.txt

#### extract_transcriptome
Create transcriptome from genome and gene annotations using [gffread](https://github.com/gpertea/gffread).


**Input:** `genome` and `gtf` of the input samples table    
**Output:** transcriptome fasta file.    

 
#### create_index_salmon
Create index for [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) quantification. Salmon index of transcriptome, required for mapping-based mode of Salmon. The index is created via an auxiliary k-mer hash over k-mers of length 31. While mapping algorithms will make use of arbitrarily long matches between the query and reference, the k-mer size selected here will act as the minimum acceptable length for a valid match.  A k-mer size of 31 seems to work well for reads of 75bp or longer, although smaller size might improve sensitivity. A smaller k-mer size is suggested when working with shorter reads.

**Input:** transcriptome fasta file for transcripts to be quantified    
**Parameters:** kmer length (`kmer` in the input samples table).    
**Output:** salmon index, used for quantification.


#### create_index_kallisto
Create index for [Kallisto](https://pachterlab.github.io/kallisto/manual) quantification. Similar to salmon index described above. A default kmer-size of 31 is used in this workflow.

**Input:** transcriptome fasta file for transcripts to be quantified    
**Output:** kallisto index, used for kallisto quantification.    


#### extract_transcripts_as_bed12
Convert transcripts from gtf to bed12 format. Required for the TIN score calculation. No user customised parameters. [GitLab repository](https://git.scicore.unibas.ch/zavolan_group/tools/gtf_transcript_type_to_bed12/)    

**Input:** gtf file    
**Output:** "full_transcripts_protein_coding.bed"    


#### index_genomic_alignment_samtools
Index the genomic alignment with [samtools index](http://quinlanlab.org/tutorials/samtools/samtools.html#samtools-index). Indexing a genome sorted BAM file enables quick extraction of alignments overlapping particular genomic regions. It is also required by genome viewers such as IGV allowing for quick display of read coverages in specific genomic regions chosen by the user.    
Required for TIN score calculation and bedgraph coverage calculation.    

**Input:** bam file    
**Output:** bam.bai index file    


#### star_rpm
Create stranded bedgraph coverage based on RPM normalisation provided by STAR. There are two ways of counting the coverage: using *Unique* reads alone, or using *UniqueMultiple* (unique and multi-mapping) reads.
Description [here](https://ycl6.gitbooks.io/rna-seq-data-analysis/visualization.html)
STAR RPM uses SAM flags to correctly infer where each read and its mate are mapped. The assignment to either plus or minus strand is based on mate1. If mate1 is mapped to the plus strand, then mate2 is mapped to the minus strand but both reads will be considered for the plus strand by STAR.

In `bedtools genomecov -bg -split`, each reads is assigned to its respective strand irrespective of its mate.

**Input:** .bam, .bam.bai index
**Output:** coverage bedGraphs
**Non-customisable arguments:**
--outWigStrans “Stranded”
--outWigNorm “RPM”


#### rename_star_rpm_for_alfa
Local rule to rename and copy the stranded bedgraph coverage tracks to comply with [ALFA](https://github.com/biocompibens/ALFA).
The renaming to `plus.bg` and `minus.bg` depends on the library orientation, as specified in `kallisto_directionality`. 

**Input:** .bg coverage tracks    
**Output:** renamed and copied bedgraph files    


#### calculate_TIN_scores
Calculation of Transcript Integrity Number (TIN) for each transcript [GitLab repository](https://git.scicore.unibas.ch/zavolan_group/tools/tin_score_calculation). Requires a set of BAM files and a BED file containing the gene annotation. TIN is conceptually similar to RIN (RNA integrity number) but provides transcript level measurement of RNA quality and is more sensitive in measuring low quality RNA samples:

* TIN score of a transcript measures the RNA integrity of the transcript.
* Median TIN score across all transcripts measures RNA integrity of that "RNA sample".
* TIN ranges from 0 (the worst) to 100 (the best). TIN = 60 means: 60% of the transcript would be covered if the read coverage was uniform.
* TIN is 0 if the transcript has no coverage or it is lower than a default cutoff.

**Input:** aligned reads.bam.bai, "full_transcripts_protein_coding.bed"         
**Output:** TIN score tsv file



#### merge_TIN_scores
Concatenate the tsv files of all samples into one table.

**Input:** TIN score tsv files per sample    
**Output:** TIN score tsv file for all samples



#### plot_TIN_scores
Generate sample-wise [boxplots](https://en.wikipedia.org/wiki/Box_plot) of TIN scores.

**Input:** TIN score tsv file for all samples  
**Output:** .pdf and .png files with boxplots



#### salmon_quantmerge_genes
Merge the salmon quantification *gene* results for all samples into a single file. Gene expression metrics are TPM (Transcripts Per Million) and number of reads. For both of these metrics a separate table is created.

**Input:** All salmon quant genes files of same seqmode   
**Output:** Two tsv files for gene quantifications, one for tpm and one for number of reads.   

#### salmon_quantmerge_transcripts
Merge the salmon quantification *transcript* results for all samples into a single file. Gene expression metrics are TPM (Transcripts Per Million) and number of reads. For both of these metrics a separate table is created.

**Input:** All salmon quant transcript files of same seqmode     
**Output:** Two tsv files for transcript quantifications, one for tpm and one for number of reads.   


#### generate_alfa_index
Create ALFA index required by [ALFA](https://github.com/biocompibens/ALFA), for a given organism.    

**Input:** .gtf genome annotation, chrNameLength.txt file containing chromosome names and lengths    
**Output:** two ALFA index files, one stranded and one unstranded    


#### alfa_qc
Run [ALFA](https://github.com/biocompibens/ALFA) from stranded bedgraph tracks.
The library orientation is required as *fr-firststrand* and *fr-secondstrand*. Currently, the values from `kallisto_directionality` are re-used.
ALFA counts features in the bedgraph coverage tracks, which were generated in `rename_star_rpm_for_alfa`. It uses the library orientation and the ALFA index files to count features. The resulting counts are stored in `ALFA_feature_counts.tsv`.
The main output of ALFA are two plots, `ALFA_Biotypes.pdf` and `ALFA_Categories.pdf`. They display the nucleotide distributions among the different features and their enrichment. For details see [ALFA documentation](https://github.com/biocompibens/ALFA).


**Input:** the renamed .bg files (suffixed with `out.plus.bg` and `out.minus.bg`), library orientation, the stranded ALFA index file
**Output:** ALFA_Biotypes.pdf and ALFA_Categories.pdf; ALFA_feature_counts.tsv containing table for the plots



#### alfa_qc_all_samples
Combine the output of all samples into one plot generated by [ALFA](https://github.com/biocompibens/ALFA).
The rule uses the output of `alfa_qc`, for both bedgraph tracks *Unique* and *UniqueMultiple*. These tracks count unique and unique and multi-mapping reads respectively. See `star_rpm` for more information.

**Input:** ALFA_feature_counts.tsv from each sample in `samples.tsv`
**Output:** Unique/ALFA_Biotypes.pdf, Unique/ALFA_Categories.pdf, UniqueMultiple/ALFA_Biotypes.pdf and UniqueMultiple/ALFA_Categories.pdf for all samples together


#### alfa_concat_results
Concatenate and convert ALFA output plots into one `.png`.

**Input:** Unique/ALFA_Biotypes.pdf, Unique/ALFA_Categories.pdf, UniqueMultiple/ALFA_Biotypes.pdf and UniqueMultiple/ALFA_Categories.pdf
**Output:** ALFA_plots.concat.png
**Parameters** density; ImageMagick’s density parameter for image output


#### prepare_files_for_report
This is an internal rule with `run` directive. It gathers all the output files, restructures the `log` and `results` directories and modifies some `stdout` and `stderr` streams of previous rules for proper parsing of sample names in the final report.


#### prepare_MultiQC_config
Prepares a dedicated config file for [MultiQC](https://multiqc.info/).

**Input:** Currently directories created during `prepare_files_for_report` serve as input.
**Output:** Config file in .yaml format


#### MULTIQC_report
Interactive report of the various workflow steps. [MultiQC](https://multiqc.info/) gathers results and logs after each distinct workflow step, parses them and presents the output graphically in an HTML file.

**Input:** Config file fort MultiQC in .yaml format
**Output:** Directory with automatically generated HTML report



### Sequencing mode specific rules
#### (pe_)fastqc
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) enables quality control checks on raw sequence data coming from high throughput sequencing workflows. It produces a modular set of analyses which provide insights on data quality and other issues affecting further analysis steps.   

**Input:** raw fastq file(s)    
**Output:** fastqc report (.txt) and several figures (.png)



#### (pe_)remove_adapters_cutadapt
[Cutadapt](https://cutadapt.readthedocs.io/en/stable/) detects and removes adapter sequences, primers, and other types of unwanted sequence contaminants from high-throughput sequencing reads.   

**Input:** fastq reads    
**Parameters:** Adapters to be removed, specified by user in the columns `fq1_3p`, `fq1_5p`, `fq2_3p`, `fq2_5p` respectively.    
**Output:** fastq files with adapters removed  


**Non-customisable arguments:**        
-e 0.1  maximum error-rate of 10%    
-j 8    use 8 threads    
-m 10   Discard processed reads that are shorter than 10    
-n 2    search for all the given adapter sequences repeatedly, either until no adapter match was found or until 2 rounds have been performed.    

*paired end:*    
--pair-filter=any      filtering criteria must apply to any of the two reads in order for a read pair to be discarded

#### (pe_)remove_polya_cutadapt
Here, [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) is used to remove poly(A) tails. 

**Input:** fastq reads    
**Parameters:** Adapters to be removed, specified by user in the columns 'fq1_polya', 'fq2_polya', respectively.    
**Output:** fastq files with poly(A) tails removed, reads shorter than 10nt will be discarded. 

**Arguments similar to remove_adapters_cutadapt and additionally:**    
-n 1    search for all the given adapter sequences repeatedly, either until no adapter match was found or until 1 round has been performed.    
*paired end:*
--pair-filter=any      filtering criteria must apply to both reads in order for a read pair to be discarded

*single end:*    
-O 1    minimal overlap of 1

#### (pe_)map_genome_star
Spliced Transcripts Alignment to a Reference; Read the [Publication](https://www.ncbi.nlm.nih.gov/pubmed/23104886) or check out the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf).    

**Input:** STAR_index, reads as .fastq.gz    
**Parameters:**    
* index size, specified by user in column `index_size`
* multimappers, specified by user in column `multimappers`   
* soft_clip, specified by user in column `soft_clip`
* pass_mode, specified by user in column `pass_mode`     

**Output:** aligned reads as .bam, STAR logfile

**Non-customisable arguments:**    
--outSAMunmapped None: do not output unmapped reads in SAM file    
--outFilterMultimapScoreRange 0: the score range below the maximum score for multimapping alignments    
--outSAMattributes All: NH HI AS nM NM MD jM jI MC ch    
--outStd BAM_SortedByCoordinate:  which output will be directed to stdout (standard out)    
--outSAMtype BAM SortedByCoordinate:   type of SAM/BAM output    
--outFilterMismatchNoverLmax 0.04: alignment will be output only if its ratio of mismatches to *mapped* length is less than or equal to this value    
--outFilterScoreMinOverLread 0.3: same as outFilterScoreMin, but normalized to read length (sum of mates’ lengths for paired-end reads)    
--outFilterMatchNminOverLread 0.3:  Minimal fraction of aligned bases    
--outFilterType BySJout: reduces the number of ”spurious” junctions     
--outReadsUnmapped None: do not output unmapped reads.

*Same for single- and paired-end.*


#### (pe_)quantification_salmon
[Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) is a tool for wicked-fast transcript quantification from RNA-seq data.

**Input:** 
* .fastq.gz reads, adapters and poly(A)tails removed. 
* filtered annotation .gtf
* salmon index, from **create_index_salmon**    

**Parameters:**  libType, specified by user as `libtype`    
**Output:** gene and transcript quantifications

**Non-customisable arguments:**    
--seqBias: [Correct for sequence specific biases](https://salmon.readthedocs.io/en/latest/salmon.html#seqbias)    
--validateMappings: Enables selective alignment of the sequencing reads when mapping them to the transcriptome. This can improve both the sensitivity and specificity of mapping and, as a result, can [improve quantification accuracy](https://salmon.readthedocs.io/en/latest/salmon.html#validatemappings).    
--writeUnmappedNames: write out the names of reads (or mates in paired-end reads) that do not map to the transcriptome.

*additionally for single end:*    
**Parameters:** 
* --fldMean: fragment length, user specified as `mean`
* --fldSD: fragment length SD, user specified as `sd`   



#### (pe_)genome_quantification_kallisto
[kallisto](http://pachterlab.github.io/kallisto/manual.html) is a program for quantifying abundances of transcripts from RNA-Seq data, or more generally of target sequences using high-throughput sequencing reads. It is based on the novel idea of pseudoalignment for rapidly determining the compatibility of reads with targets, without the need for alignment.    

**Input:** 
* .fastq.gz reads, adapters and poly(A)tails removed.
* kallisto index, from **create_index_kallisto** 

**Parameters:** directionality, which is `kallisto_directionality` from samples table    

**Output:** Pseudoalignment .sam file

*additionally for single end:*    
* -l: fragment length, user specified as `mean`
* -s: fragment length SD, user specified as `sd` 

