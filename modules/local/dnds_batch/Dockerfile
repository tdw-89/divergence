FROM mambaorg/micromamba:1.5.8

RUN micromamba install -y -n base -c conda-forge -c bioconda \
    biopython=1.84 \
    paml=4.10.7 \
    procps-ng \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}"
