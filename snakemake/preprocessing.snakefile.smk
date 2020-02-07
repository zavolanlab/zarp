

rule index_genome_STAR:
    '''
    Create Star index
    '''
	input:
		genome = os.path.join(config["output_dir"], "genome.fa"),
		annotation = os.path.join(config["output_dir"], "annotation.gtf")
	output:
		output = os.path.join(config["database_path"], config['organism'], config['STAR_idx_folder], "STAR_index" + {sjdb})
	params:
		outputdir = os.path.join(config["output_dir"],"STAR_index"),
		sjdb = lambda wildcards: samples.loc['sjdb']
	threads:	8
	singularity:
		"docker://zavolab/star:2.6.0a"
	log:
		os.path.join(config["local_log"],"index_genome_STAR.log")
	shell:
		"mkdir -p {output.output}; \
		chmod -R 777 {output.output}; \
		(STAR --runMode genomeGenerate \
		--sjdbOverhang {params.sjdbOverhang} \
		--genomeDir {params.outputdir} \
		--genomeFastaFiles {input.genome} \
		--runThreadN {threads} \
		--sjdbGTFfile {input.annotation}) &> {log}"