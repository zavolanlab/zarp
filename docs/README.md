<figure markdown>
  ![ZARP logo](./images/zarp_logo.384px.png){ width="384" }
</figure>

# ZARP

**Welcome to the _ZARP_ documentation pages!**

_ZARP_ is a generic RNA-Seq analysis workflow that allows users to process and analyze Illumina short-read sequencing libraries with minimum effort. Better yet: With our companion [**ZARP-cli**](https://github.com/zavolanlab/zarp-cli) command line interface, you can start ZARP runs with the simplest and most intuitive commands.

_RNA-seq analysis doesn't get simpler than that!_

The workflow is developed in [Snakemake][snakemake], a widely used workflow management system in the bioinformatics community. ZARP will pre-process, align and quantify your single- or paired-end stranded bulk RNA-seq sequencing libraries with publicly available state-of-the-art bioinformatics tools. ZARP's browser-based rich reports and visualizations will give you meaningful initial insights in the quality and composition of your sequencing experiments - fast and simple. Whether you are an experimentalist struggling with large scale data analysis or an experienced bioinformatician, when there's RNA-seq data to analyze, just _zarp 'em_!

## How does it work?

ZARP requires conda or mamba to install the basic dependencies. Each individual step of the workflow run either in its own Apptainer (Singularity) container or in its own Conda virtual environment.

Once the installation is complete, you fill in a [config.yaml](https://github.com/zavolanlab/zarp/blob/dev/tests/input_files/config.yaml) file with parameters and a [samples.tsv](https://github.com/zavolanlab/zarp/blob/dev/tests/input_files/samples.tsv) file with sample specific information. You can easily trigger ZARP by making a call to snakemake with the appropriate parameters.

The pipeline can be executed in different systems or HPC clusters. ZARP generates multiple output files that help you Quality Control (QC) your data and proceed with downstream analyses. Apart from running the main ZARP workflow, you can also run a second pipeline that pulls sequencing sample data from the Sequence Read Archive (SRA), and a third pipeline that populates a file with the samples and infers missing metadata.

## How to cite

If you use _ZARP_ in your work, please kindly cite the following article:

**ZARP: A user-friendly and versatile RNA-seq analysis workflow**  
_Maria Katsantoni, Foivos Gypas, Christina J. Herrmann, Dominik Burri, Maciej
Bak, Paula Iborra, Krish Agarwal, Meric Ataman, Máté Balajti, Noè Pozzan, Niels
Schlusser, Youngbin Moon, Aleksei Mironov, Anastasiya Börsch, Mihaela Zavolan,
Alexander Kanitz_  
F1000Research 2024, 13:533  
<https://doi.org/10.12688/f1000research.149237.1>

## Info materials

### Poster

<p float="left">
  <a href="https://f1000research.com/posters/13-968"><img alt="ZARP poster" src="./images/poster_ZARP_latest.jpg" width="100" /></a>
</p>

## Reach out

There are several ways to get in touch with us:

- For ZARP usage questions, please use the [_ZARP_ Q&A forum][zarp-qa] (requires
  [GitHub registration][github-signup])
- For feature suggestions and bug reports, please use either the
  [ZARP][zarp-issue-tracker] or [ZARP-cli][zarp-cli-issue-tracker] issue
  tracker (requires [GitHub registration][github-signup])
- For any other requests, please reach out to us via [email][contact]

!!! info "Contributors welcome!"

    Open source contributors are always welcome, for [_ZARP_][zarp],
    [_ZARP-cli_][zarp-cli] or any other of the [Zavolab
    projects][zavolab-gh]. Simply reach out by [email][contact] to schedule
    an onboarding call.

## Acknowledgements

[![Zavolab](images/zavolab_logo.200px.png)](https://www.biozentrum.unibas.ch/research/research-groups/research-groups-a-z/overview/unit/research-group-mihaela-zavolan)
[![Biozentrum, University of Basel](images/biozentrum_logo.200px.png)](https://www.biozentrum.unibas.ch/)
[![Swiss Institute of Bioinformatics](images/sib_logo.200px.png)](https://www.sib.swiss/)
