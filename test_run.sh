snakemake --snakefile workflow/rules/sra_download.smk \
	--config samples=tests/input_files/sra_samples.tsv \
		outdir=sra_downloads \
		samples_out=sra_downloads/sra_samples.out.tsv \
	--cores 8 --use-conda --verbose 