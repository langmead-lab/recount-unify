FROM continuumio/miniconda3:4.5.12

RUN apt-get update -y && apt-get install -y gcc g++ make git python2.7 python-pip python-dev ant default-jdk sqlite3 tabix

RUN /usr/bin/pip install numpy

RUN mkdir -p /recount-unify
COPY rejoin/ /recount-unify/rejoin/
COPY merge/ /recount-unify/merge/
COPY scripts/ /recount-unify/scripts/
COPY log_qc/ /recount-unify/log_qc/
COPY annotate/ /recount-unify/annotate/
COPY sample_ids/ /recount-unify/sample_ids/

WORKDIR /recount-unify/rejoin/
RUN make all

WORKDIR /recount-unify/merge/
RUN make all

WORKDIR /recount-unify/scripts/
RUN make all

WORKDIR /

# install conda environment
COPY env.yml /recount-unify/env-base.yml
RUN chmod a+r /recount-unify/env-base.yml

RUN conda env create -q -f /recount-unify/env-base.yml && conda clean -y -p -t

RUN apt-get remove --purge -y gcc g++ make

# install Snakefile
COPY Snakefile /recount-unify/Snakefile
COPY Snakefile.study_jxs /recount-unify/Snakefile.study_jxs
RUN chmod a+r /recount-unify/Snakefile
RUN chmod a+r /recount-unify/Snakefile.study_jxs

# install generic workflow wrapper script
COPY workflow.bash /recount-unify/workflow.bash
RUN chmod a+rx /recount-unify/workflow.bash

COPY list_of_zeros.gz /recount-unify/list_of_zeros.gz

CMD ["bash", "-c", "source activate recount-unify && /recount-unify/workflow.bash"]
