FROM continuumio/anaconda3:5.0.1

USER root

WORKDIR /home/lael/software

#RUN conda install python=3.6.9

# Update all conda to avoid inconsistencies.
RUN conda update --all

RUN conda install python=3.6.9

# This fails:
#RUN conda install -c anaconda gcc

# Install a C compiler.
#RUN apt install build-essential

# Install biopython for parsing various bioinformatics software output files,
# etc.
RUN conda install -c anaconda biopython 

# Install NCBI BLAST for similarity searching.
#RUN conda install -c bioconda blast
RUN echo '***********************************' && \
    echo '*******Installing ncbi-blast+ *****' && \
    echo '***********************************' && \
    wget "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.10.0+-x64-linux.tar.gz" &&\
    tar zxvpf ncbi-blast-2.10.0+-x64-linux.tar.gz && \
    rm ncbi-blast-2.10.0+-x64-linux.tar.gz

ENV PATH "$PATH:/home/lael/software/ncbi-blast-2.10.0+/bin"


# Install HMMer3 for similarity searching.
#RUN conda install -c bioconda hmmer
#RUN echo '***********************************' && \
#    echo '*******Installing IQ-TREE     *****' && \
#    echo '***********************************' && \
#    wget "http://eddylab.org/software/hmmer/hmmer-3.3.tar.gz" &&\
#    tar zxvpf hmmer-3.3.tar.gz &&\
#    rm hmmer-3.3.tar.gz &&\
#    cd hmmer-3.3 &&\
#    ./configure &&\
#    make &&\
#    make check &&\
#    make install &&\
#    cd easel; make install
#
#WORKDIR /home/lael/software
#
#ENV PATH "$PATH:/home/lael/software/hmmer-3.3/binaries"

RUN conda install -c bioconda hmmer

# Install iqtree for phylogenetic analysis.
RUN echo '***********************************' && \
    echo '*******  Installing IQ-TREE   *****' && \
    echo '***********************************' && \
    wget "https://github.com/Cibiv/IQ-TREE/releases/download/v1.6.12/iqtree-1.6.12-Linux.tar.gz" && \
    tar zxvpf iqtree-1.6.12-Linux.tar.gz &&\
    rm iqtree-1.6.12-Linux.tar.gz

ENV PATH "$PATH:/home/lael/software/iqtree-1.6.12-Linux/bin"
#RUN conda install -c bioconda iqtree

# Install GCC for compiling mrbayes.
#RUN conda install -c anaconda gcc # Didn't work because of missing dependancy.
# Install clang instead for compiling mrbayes.
#RUN conda install -c anaconda clang # Also doesn't work!

# Install MrBayes for phylogenetic analysis.
#RUN conda install -c bioconda mrbayes
#RUN echo '***********************************' && \
#    echo '*******        MrBayes        *****' && \
#    echo '***********************************' && \
#    git clone --depth=1 https://github.com/NBISweden/MrBayes.git  --branch v3.2.7a &&\
#    cd MrBayes &&\
#    ./configure &&\
#    make && sudo make install 
#
#WORKDIR /home/lael/software

# Install MUSCLE for multiple sequence alignment.
#RUN conda install -c bioconda muscle
RUN echo '***********************************' && \
    echo '*******  Installing MUSCLE    *****' && \
    echo '***********************************' && \
    wget "https://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz" && \
    tar zxvpf muscle3.8.31_i86linux64.tar.gz && \
    rm muscle3.8.31_i86linux64.tar.gz && \
    mv muscle3.8.31_i86linux64 muscle

ENV PATH "$PATH:/home/lael/software"

# Install additional python packages.
RUN pip install --upgrade pip
RUN pip install gffutils
RUN pip install PyPDF2
RUN pip install reportlab

# Install ete3 for working with phylogenetic trees.
#RUN conda install -c etetoolkit ete3
RUN pip install ete3 

# Install exonerate for exon prediction (at least better than TBLASTN).
#RUN conda install -c bioconda exonerate
RUN echo '***********************************' && \
    echo '******* Installing exonerate  *****' && \
    echo '***********************************' && \
    wget "http://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/exonerate-2.2.0-x86_64.tar.gz" &&\
    tar zxvpf exonerate-2.2.0-x86_64.tar.gz && \
    rm exonerate-2.2.0-x86_64.tar.gz

ENV PATH "$PATH:/home/lael/software/exonerate-2.2.0-x86_64/bin"

# Install latex environments for jupyter markdown cells (all these don't work).
#RUN conda install -c conda-forge jupyter_latex_envs
#RUN conda install -c conda-forge texlive-core
#RUN conda install -c anaconda nbconvert
#RUN apt-get install texlive-xetex texlive-fonts-recommended texlive-generic-recommended


# Add amoebae repository directory to path so that amoebae scripts can be run
# from the docker container.
ENV PATH "$PATH:/opt/notebooks"






