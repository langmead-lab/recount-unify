FROM quay.io/broadsword/recount-unify-base:1.1.0
 

WORKDIR /

COPY rejoin/ /recount-unify/rejoin/
COPY merge/ /recount-unify/merge/
COPY scripts/ /recount-unify/scripts/
COPY log_qc/ /recount-unify/log_qc/
COPY annotate/ /recount-unify/annotate/
COPY sample_ids/ /recount-unify/sample_ids/
COPY metadata/ /recount-unify/metadata/
COPY recount-pump/metadata/scripts /recount-unify/recount-pump/metadata/scripts/

# install Snakefile
COPY Snakefile /recount-unify/Snakefile
COPY Snakefile.study_jxs /recount-unify/Snakefile.study_jxs
RUN chmod a+r /recount-unify/Snakefile
RUN chmod a+r /recount-unify/Snakefile.study_jxs

# install generic workflow wrapper script
COPY workflow.bash /recount-unify/workflow.bash
RUN chmod a+rx /recount-unify/workflow.bash

RUN echo 'pigz --stdout -p2 -d $1' > /usr/bin/pcat
RUN chmod a+rx /usr/bin/pcat

WORKDIR /recount-unify/rejoin/
RUN make clean
RUN make all

#need latest megadepth to do fast summing across samples BigWigs for AUC check against what was rejoined
RUN wget https://github.com/ChristopherWilks/megadepth/releases/download/1.2.0/megadepth -O /recount-unify/rejoin/megadepth120
RUN chmod a+rx /recount-unify/rejoin/megadepth120

WORKDIR /recount-unify/merge/
RUN make clean
RUN make all

WORKDIR /recount-unify/scripts/
RUN make clean
RUN make all

CMD ["bash", "-c", "export PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin && source activate recount-unify && /recount-unify/workflow.bash"]
