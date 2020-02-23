FROM continuumio/anaconda3

# miniconda may be a more efficient way to set up amoebae in the future.
#FROM continuumio/miniconda3

USER root

WORKDIR /home/amoebae_user/software

# Update all conda to avoid inconsistencies.
#RUN conda update --all

# These packages are necessary for installing ete3:
# ete3 3.1.1 does not work with Python 3.7 or above.
RUN conda install python=3.6.9
RUN conda install lxml
# libgl1-mesa-dev is a dependency of PyQt, which is apparently not provided via
# conda.
RUN apt-get update
RUN apt-get install -y libgl1-mesa-dev
# PyQt5 is necessary for ete3 to generate graphics.
RUN conda install pyqt
# Check that PyQt5 can actually be imported via python3.
RUN python3 -c "from PyQt5 import QtGui"

# Install ete3 (the ete3 external tools are not available for install via conda
# on linux).
RUN conda install -c etetoolkit ete3 ete_toolchain
# Check that ete3 was installed correctly.
RUN ete3 build check



# This fails:
#RUN conda install -c anaconda gcc

# Install a C compiler.
#RUN apt install build-essential

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
#WORKDIR /home/amoebae_user/software




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

ENV PATH "$PATH:/home/amoebae_user/software/ncbi-blast-2.10.0+/bin"


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
#WORKDIR /home/amoebae_user/software
#
#ENV PATH "$PATH:/home/amoebae_user/software/hmmer-3.3/binaries"
# Install via conda.
RUN conda install -c bioconda hmmer

# Install iqtree for phylogenetic analysis.
RUN echo '***********************************' && \
    echo '*******  Installing IQ-TREE   *****' && \
    echo '***********************************' && \
    wget "https://github.com/Cibiv/IQ-TREE/releases/download/v1.6.12/iqtree-1.6.12-Linux.tar.gz" && \
    tar zxvpf iqtree-1.6.12-Linux.tar.gz &&\
    rm iqtree-1.6.12-Linux.tar.gz

ENV PATH "$PATH:/home/amoebae_user/software/iqtree-1.6.12-Linux/bin"
#RUN conda install -c bioconda iqtree


# Install MUSCLE for multiple sequence alignment.
#RUN conda install -c bioconda muscle
RUN echo '***********************************' && \
    echo '*******  Installing MUSCLE    *****' && \
    echo '***********************************' && \
    wget "https://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz" && \
    tar zxvpf muscle3.8.31_i86linux64.tar.gz && \
    rm muscle3.8.31_i86linux64.tar.gz && \
    mv muscle3.8.31_i86linux64 muscle

ENV PATH "$PATH:/home/amoebae_user/software"

# Install additional python packages.
RUN pip install --upgrade pip
RUN pip install gffutils
RUN pip install PyPDF2
RUN pip install reportlab


# Install exonerate for exon prediction (at least better than TBLASTN).
#RUN echo '***********************************' && \
#    echo '******* Installing exonerate  *****' && \
#    echo '***********************************' && \
#    wget "http://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/exonerate-2.2.0-x86_64.tar.gz" &&\
#    tar zxvpf exonerate-2.2.0-x86_64.tar.gz && \
#    rm exonerate-2.2.0-x86_64.tar.gz
#
#ENV PATH "$PATH:/home/amoebae_user/software/exonerate-2.2.0-x86_64/bin"
# Install via conda.
RUN conda install -c bioconda exonerate

# Install latex environments for jupyter markdown cells (all these don't work).
#RUN conda install -c conda-forge jupyter_latex_envs
#RUN conda install -c conda-forge texlive-core
#RUN conda install -c anaconda nbconvert
#RUN apt-get install texlive-xetex texlive-fonts-recommended texlive-generic-recommended

# Install runipy for running Jupyter notebooks as scripts from command line. 
RUN apt-get install python3-pip
RUN pip3 install runipy

# Install jupyter notebook extensions.
# This includes the TableOfContents extension for navigating large jupyter
# notebooks: https://jupyter-contrib-nbextensions.readthedocs.io/en/latest/nbextensions/toc2/README.html
# documentation: https://github.com/ipython-contrib/jupyter_contrib_nbextensions
RUN conda install -c conda-forge jupyter_contrib_nbextensions

# Install program for displaying directory contents in a tree structure.
RUN apt-get install tree


# Add amoebae repository directory to path so that amoebae scripts can be run
# from the docker container.
ENV PATH "$PATH:/opt/notebooks"






