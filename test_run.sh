snakemake --snakefile workflow/rules/sra_download.smk \
	--config samples=tests/input_files/sra_samples.tsv \
		outdir=sra_downloads \
		samples_out=sra_downloads/sra_samples.out.tsv \
		cluster_log_dir=logs/cluster_log \
	--cores 8 --verbose \
	--profile profiles/local-conda
