FROM continuumio/miniconda3:4.5.12
 
RUN apt-get update -y && apt-get install -y gcc g++ make git python2.7 python-pip python-dev ant default-jdk sqlite3 wget curl parallel libz-dev

## set Snaptron related up here
RUN mkdir -p /recount-unify
# clone master branch of Snaptron
RUN git clone https://github.com/ChristopherWilks/snaptron.git /recount-unify/snaptron
# cribbed from https://bitbucket.org/coady/docker/src/tip/pylucene/Dockerfile
#RUN conda deactivate
#RUN export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
RUN mkdir -p /usr/src/pylucene
WORKDIR /usr/src/pylucene
RUN which python
RUN curl http://snaptron.cs.jhu.edu/data/pylucene-6.5.0-src.tar.gz | tar -xz --strip-components=1
RUN cd jcc && JCC_JDK=/usr/lib/jvm/default-java /usr/bin/python setup.py install
RUN mkdir -p /root/.ant/lib/
RUN wget https://repo1.maven.org/maven2/org/apache/ivy/ivy/2.3.0/ivy-2.3.0.jar -O /root/.ant/lib/ivy-2.3.0.jar
RUN make all install JCC='/usr/bin/python -m jcc' ANT=ant PYTHON=/usr/bin/python NUM_FILES=8

WORKDIR ..
RUN rm -rf pylucene
WORKDIR /

#get no compression version of bgzip/tabix
#doesn't matter it's old, it's still comptabile with later versions
#and we don't need multi-thread support in bgzip because we aren't compressing
RUN wget http://snaptron.cs.jhu.edu/data/htslib-1.2.1_nocomp.tar.gz
RUN tar -zxvf htslib-1.2.1_nocomp.tar.gz
WORKDIR /htslib-1.2.1_nocomp
RUN make
RUN mv tabix /bin/tabix_nocomp
RUN mv bgzip /bin/bgzip_nocomp
#now get zstd versions of bgzip/tabix 1.2.1
WORKDIR /bin/
RUN wget http://snaptron.cs.jhu.edu/data/bgzip_zstd
RUN wget http://snaptron.cs.jhu.edu/data/tabix_zstd
RUN chmod a+x bgzip_zstd tabix_zstd tabix_nocomp bgzip_nocomp
WORKDIR /

#RUN /usr/bin/pip install numpy

RUN /usr/bin/pip install -r /recount-unify/snaptron/requirements.txt

#now back to recount-unify
#RUN export PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
#RUN conda activate
WORKDIR /
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

#RUN apt-get remove --purge -y gcc g++ make

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

CMD ["bash", "-c", "export PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin && source activate recount-unify && /recount-unify/workflow.bash"]
