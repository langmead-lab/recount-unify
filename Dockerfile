FROM continuumio/miniconda3:4.5.12

USER root

RUN apt-get update -y && apt-get install -y gcc g++ libz-dev make unzip zip

# install conda environment
COPY env.yml /env-base.yml
RUN chmod a+r /env-base.yml

RUN conda env create -q -f /env-base.yml && conda clean -y -p -t

# install Snakefile
COPY Snakefile /Snakefile
RUN chmod a+r /Snakefile

# install generic workflow wrapper script
COPY workflow.bash /workflow.bash
RUN chmod a+rx /workflow.bash

USER recount

WORKDIR /data

CMD ["bash", "-c", "source activate recount-unify && /workflow.bash"]
