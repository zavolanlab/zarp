FROM continuumio/miniconda3:24.7.1-0

COPY install/environment.yml /environment.yml
COPY workflow /workflow
COPY resources /resources
COPY tests/input_files/config.yaml /config.yaml
COPY tests/input_files/samples.tsv /samples.tsv
COPY tests/input_files/rule_config.yaml /rule_config.yaml
COPY tests/input_files/project1/synthetic.mate_1.fastq.gz /project1/synthetic.mate_1.fastq.gz
COPY tests/input_files/project1/synthetic.mate_2.fastq.gz /project1/synthetic.mate_2.fastq.gz
COPY tests/input_files/project2/synthetic.mate_1.fastq.gz /project2/synthetic.mate_1.fastq.gz
COPY tests/input_files/homo_sapiens/annotation.gtf /annotation.gtf
COPY tests/input_files/homo_sapiens/genome.fa /genome.fa

RUN sed -i 's#  - conda-forge##' workflow/envs/STAR.yaml && \
  sed -i 's#2.7.11#2.7.10#' workflow/envs/STAR.yaml && \
  sed -i 's#../input_files/project1/#/project1/#g' /samples.tsv && \
  sed -i 's#../input_files/project2/#/project2/#g' /samples.tsv && \
  sed -i 's#../input_files/homo_sapiens/##g' /samples.tsv && \
  sed -i 's#../input_files/##' /config.yaml

RUN conda install -c conda-forge mamba --yes && \
  mamba env create -f /environment.yml && \
  conda clean --all --yes

RUN echo "source activate zarp" > ~/.bashrc

ENV SNAKEMAKE_CONDA_PREFIX="/conda_envs"
ENV PATH=/opt/conda/envs/zarp/bin:$PATH

RUN snakemake -p --snakefile /workflow/Snakefile --configfile /config.yaml --cores 4 --use-conda --conda-create-envs-only --verbose && \
  conda clean --all --yes

RUN rm /config.yaml /samples.tsv /rule_config.yaml /project1/synthetic.mate_1.fastq.gz  /project1/synthetic.mate_2.fastq.gz /project2/synthetic.mate_1.fastq.gz

RUN mkdir -p /data
