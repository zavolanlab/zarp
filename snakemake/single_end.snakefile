
rule fastqc:
	''' A quality control tool for high throughput sequence data. '''
	input:
		reads = os.path.join(get_fq_1("{sample}"))
	output:
		outdir = directory(os.path.join(config["output_dir"], "single_end", "{sample}", "fastqc"))
	singularity:
		"docker://zavolab/fastqc:0.11.8"
	log:
		os.path.join(config["local_log"], "single_end", "{sample}", "fastqc.log")
	shell:
		"(mkdir -p {output.outdir}; \
		fastqc \
		--outdir {output.outdir} \
		{input.reads}) &> {log}"


rule htseq_qa:
	''' Assess the technical quality of a run. '''
	input:
		reads = os.path.join(get_fq_1("{sample}"))
	output:
		qual_pdf = os.path.join(config["output_dir"], "single_end", "{sample}", "htseq_quality.pdf")
	singularity:
		"docker://zavolab/python_htseq:3.6.5_0.10.0"
	log:
		os.path.join(config["local_log"], "single_end", "{sample}", "htseq_qa.log")
	shell:
		"(htseq-qa \
		-t fastq \
		-o {output.qual_pdf} \
		{input.reads} ) &> {log}"


rule remove_adapters_cutadapt:
	''' Remove adapters '''
	input:
		reads = os.path.join(get_fq_1("{sample}"))
	output:
		reads = os.path.join(config["output_dir"], "single_end", "{sample}", "{sample}.remove_adapters.fastq.gz")
	params:
		adapters_3 = lambda wildcards: 
			sample_table.loc[wildcards.sample, 'fq1_3p']
		adapters_5 = lambda wildcards: 
			sample_table.loc[wildcards.sample, 'fq1_5p']

	singularity:
		"docker://zavolab/cutadapt:1.16"
	threads: 8
	log:
		os.path.join(config["local_log"], "single_end", "{sample}", "remove_adapters_cutadapt.log")
	shell:
		"cutadapt \
		-e 0.1 \
		-O 1 \
		-j {threads} \
		-m 10 \
		-n 3 \
		-a {params.adapters_3} \
		-g {params.adapters_5} \
		-o {output.reads} \
		{input.reads}) &> {log}"


rule remove_polya_cutadapt:
	''' Remove ployA  tails'''
	input:
		reads = os.path.join(get_fq_1("{sample}"))
	output:
		reads = os.path.join(config["output_dir"], "single_end", "{sample}", "{sample}.remove_polya.fastq.gz")
	params:
		polya_3 = 
		adapters_3 = lambda wildcards: 
			sample_table.loc[wildcards.sample, "fq1_polya"]
	singularity:
		"docker://zavolab/cutadapt:1.16"
	threads: 8
	log:
		os.path.join(config["local_log"], "single_end", "{sample}", "remove_polya_cutadapt.log")
	shell:
		"(cutadapt \
		--match-read-wildcards \
		-j {threads} \
		-n 2 \
		-e 0.1 \
		-O 1 \
		-q 6 \
		-m 10  \
		-a {params.polya_3} \
		-o {output.reads} \
		{input.reads}) &> {log}"


rule create_index_star:
	''' Create index using STAR. '''
	input:
		genome = sample_table.loc["{sample}", "genome"],
		gtf = sample_table.loc["{sample}", "gtf"]
	output:
		chromosome_info = os.path.join(
					config["star_indexes"],
					sample_table.loc["{sample}", "organism"],
					sample_table.loc["{sample}", "index_size"], 
					"STAR_index","chrNameLength.txt"),
		chromosome_names = os.path.join(
					config["star_indexes"],
					sample_table.loc["{sample}", "organism"], 
					sample_table.loc["{sample}", "index_size"], 
					"STAR_index", "chrName.txt")
	params:
		output_dir = lambda wildcards:
				os.path.join(
					config["star_indexes"],
					sample_table.loc["{sample}", "organism"],
					sample_table.loc[wildcards.sample, "index_size"], 
					"STAR_index"),
		outFileNamePrefix = lambda wildcards:
				os.path.join(
					config["star_indexes"],
					sample_table.loc["{sample}", "organism"],
					sample_table.loc[wildcards.sample, "index_size"], 
					"STAR_index/STAR_"),
		sjdbOverhang = lambda wildcards:
			sample_table.loc[wildcards.sample, "index_size"]
	singularity:
		"docker://zavolab/star:2.6.0a"
	threads: 12
	log:
		os.path.join(config["local_log"], "single_end", "{sample}", "create_index_star.log")
	shell:
		"(mkdir -p {params.output_dir}; \
		chmod -R 777 {params.output_dir}; \
		STAR \
		--runMode genomeGenerate \
		--sjdbOverhang {params.sjdbOverhang} \
		--genomeDir {params.output_dir} \
		--genomeFastaFiles {input.genome} \
		--runThreadN {threads} \
		--outFileNamePrefix {params.outFileNamePrefix} \
		--sjdbGTFfile {input.gtf}) &> {log}"


rule map_genome_star:
	''' Map to genome using STAR. '''
	input:
		index = os.path.join(
					config["star_indexes"],
					sample_table.loc["{sample}", "organism"],
					sample_table.loc["{sample}", "index_size"], 
					"STAR_index","chrNameLength.txt"),
		reads = os.path.join(config["output_dir"], "single_end", "{sample}", "{sample}.remove_polya.fastq.gz")
	output:
		bam = os.path.join(config["output_dir"], "single_end", 
			"{sample}", 
			"map_genome", 
			"{sample}_Aligned.sortedByCoord.out.bam"),
		logfile = os.path.join(config["output_dir"], "single_end", 
			"{sample}", 
			"map_genome", 
			"{sample}_Log.final.out")
	params:
		sample_id = "{sample}"
		index = lambda wildcards:
				os.path.join(
					config["star_indexes"],
					sample_table.loc["{sample}", "organism"],
					sample_table.loc[wildcards.sample, "index_size"], 
					"STAR_index"),
		outFileNamePrefix = lambda wildcards:
				os.path.join(
					config["output_dir"], 
					"single_end", 
					"{sample}", "map_genome", "{sample}_"),
		multimappers = lambda wildcards:
				sample_table.loc[wildcards.sample, "multimappers"],
		soft_clip = lambda wildcards:
				sample_table.loc[wildcards.sample, "soft_clip"],
		pass_mode = lambda wildcards:
				sample_table.loc[wildcards.sample, "pass_mode"],		
	singularity:
		"docker://zavolab/star:2.6.0a"
	threads: 12
	log:
		os.path.join(config["local_log"], "single_end", "{sample}", "map_genome_star.log")
	shell:
		"(STAR \
		--runMode alignReads \
		-- twopassMode {params.pass_mode} \
		--runThreadN {threads} \
		--genomeDir {params.index} \
		--readFilesIn {input.reads} \
		--readFilesCommand zcat \
		--outSAMunmapped None  \
		--outFilterMultimapNmax {params.multimappers} \
		--outFilterMultimapScoreRange 1 \
		--outFileNamePrefix {params.outFileNamePrefix} \
		--outSAMattributes All \
		--outStd BAM_SortedByCoordinate \
		--outSAMtype BAM SortedByCoordinate \
		--outFilterMismatchNoverLmax 0.04 \
		--outFilterScoreMinOverLread 0.3 \
		--outFilterMatchNminOverLread 0.3 \
		--outFilterType BySJout \
		--outReadsUnmapped None \
		--outSAMattrRGline ID:rcrunch SM:{params.sample_id} \
		--alignEndsType {params.soft_clip}} > {output.bam};) &> {log}"


rule index_genomic_alignment_samtools:
    '''Index genome bamfile using samtools.'''
    input:
        bam = os.path.join(config["output_dir"], 
        	"single_end", 
        	"{sample}", 
			"map_genome", 
			"{sample}_Aligned.sortedByCoord.out.bam")
    output:
        bai = os.path.join(config["output_dir"], 
        	"single_end", 
			"{sample}", 
			"map_genome", 
			"{sample}_Aligned.sortedByCoord.out.bam.bai")
    singularity:
        "docker://zavolab/samtools:1.8"
    threads: 1
    log:
		os.path.join(config["local_log"], "single_end", "{sample}", "index_genomic_alignment_samtools.log")
    shell:
        "(samtools index {input.bam} {output.bai};) &> {log}"


rule create_index_salmon:
	''' Create index for Salmon quantification.  '''
	input:
		transcriptome = sample_table.loc[wildcards.sample, 'tr_fasta_filtered']
	output:
		index = os.path.join(
			config["salmon_indexes"],
			sample_table["{sample}", 'organism'],
			"salmon.idx")
	params:
		kmerLen = lambda wildcards:
			sample_table.loc[wildcards.sample, 'kmer']
	singularity:
		"docker://zavolab/salmon:0.11.0"
	log:
		os.path.join(config["local_log"], "single_end", "{sample}", "create_index_salmon.log")
	threads:	8
	shell:
		"(salmon index \
		--t {input.transcriptome} \
		--i {output.index} \
		--k {params.kmerLen} \
		--threads {threads}) &> {log}"


rule quantification_salmon:
	''' Quantification at transcript and gene level using Salmon. '''
    input:
    	reads = os.path.join(
			config["output_dir"], 
			"single_end", 
			"{sample}", 
			"{sample}.remove_polya.fastq.gz"),
        index = os.path.join(
			config["salmon_indexes"],
			sample_table["{sample}", 'organism'],
			"salmon.idx"),
       	gtf = sample_table.loc["{sample}", "gtf_filtered"]
    output:
    	gn_estimates = os.path.join(
			config["output_dir"], 
			"single_end", 
			"{sample}", 
			"salmon_quant",
			"quant.genes.sf"),
    	tr_estimates = os.path.join(
			config["output_dir"], 
			"single_end", 
			"{sample}", 
			"salmon_quant",
			"quant.sf")
    params:
        output_dir = os.path.join(
			config["output_dir"], 
			"single_end", 
			"{sample}", 
			"salmon_quant"),
        libType = lambda wildcards:
     			sample_table.loc[wildcards.sample, "libtype"]
    log:
		os.path.join(config["local_log"], "single_end", "{sample}", "quantification_salmon.log")
    threads:    12
    conda:
        "envs/salmon.yaml"
    shell:
        "(salmon quant \
        --libType {params.libType} \
        --seqBias \
        --validateMappings \
        --threads {threads} \
        --writeUnmappedNames \
        --index {input.index} \
        --geneMap {input.gtf} \
        --unmatedReads {input.reads} \
        -o {params.output_dir}) &> {log}"



rule create_index_kallisto:
	''' Create index for running Kallisto. '''
	input:
		transcriptome = sample_table.loc["{sample}", 'tr_fasta_filtered']
	output:
		index = os.path.join(
			config["kallisto_indexes"],
			sample_table.loc["{sample}", "organism"],
		 	"kallisto.idx")
	params:
		output_dir = lambda wildcards:
			os.path.join(
				config["kallisto_indexes"], 
				sample_table.loc["{sample}", "organism"]),
	log:
		os.path.join(config["local_log"], "single_end", "{sample}", "create_index_kallisto.log")
	singularity:
		"docker://zavolab/kallisto:0.9"
	shell:
		"(mkdir -p {params.output_dir}; \
		chmod -R 777 {params.output_dir}; \
		kallisto index -i {output.index} {input.transcriptome}) &> {log}"

rule genome_quantification_kallisto:
	''' Quantification at transcript and gene level using Kallisto. '''
	input:
		reads = os.path.join(
			config["output_dir"], 
			"single_end", 
			"{sample}", 
			"{sample}.remove_polya.fastq.gz"),
		index = os.path.join(
			config["kallisto_indexes"],
			sample_table.loc["{sample}", "organism"],
		 	"kallisto.idx")
	output:
		pseudoalignment = os.path.join(
			config["output_dir"], 
			"single_end", 
			"{sample}", 
			"{sample}.kallisto.pseudo.sam")
	params:
		output_dir = lambda wildcards:
			os.path.join(
				config["output_dir"], 
				"single_end", 
				"{sample}", 
				"quant_kallisto"),
		fraglen = lambda wildcards: sample_table.loc[wildcards.sample, 'mean'],
		fragsd = lambda wildcards: sample_table.loc[wildcards.sample, 'sd'],
		directionality = lambda wildcards: sample_table.loc[wildcards.sample, 'kallisto_directionality']
	threads:	    8
	log:
		os.path.join(config["local_log"],"kallisto_align_{sample}.log")
	singularity:
		"docker://zavolab/kallisto:0.9"
	shell:
		"(kallisto quant \
		-i {input.index} \
		-o {params.output_dir} \
		--single \
		-l {params.fraglen} \
		-s {params.fragsd} \
		--pseudobam \
		--{params.directionality}-stranded \
		{input.reads} > {output.pseudoalignment}) &> {log}"

		