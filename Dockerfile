FROM ubuntu:18.04

WORKDIR /opt

RUN apt-get update && apt-get install --yes \
    build-essential \
    gcc-multilib \
    gfortran \
    apt-utils \
    zlib1g-dev \
    liblzma-dev \
    xorg-dev \
    libreadline-dev \
    libpcre++-dev \
    libcurl4 \
    libcurl4-openssl-dev \
    libpango1.0-dev \
    openjdk-8-jdk \
    vim-common \
    git \
    g++ \
    python3-pip \
    libbz2-dev \
    wget && \
    export JAVA_HOME="/opt/java/jre/" && export PATH=$PATH:$JAVA_HOME/bin && \
    pip3 install snakemake && \
    wget https://cran.r-project.org/bin/linux/ubuntu/bionic-cran35/r-base_3.6.3.orig.tar.gz && \
    tar -xvzf /opt/r-base_3.6.3.orig.tar.gz && cd R-3.6.3/ &&  ./configure && make 

RUN /opt/R-3.6.3/bin/R -e "install.packages(c('BiocManager', 'stringr', 'optparse', 'foreach', 'doParallel', 'jsonlite'), dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
    /opt/R-3.6.3/bin/R -e "BiocManager::install(c('DNAcopy','GenomicRanges'))"

RUN mkdir -p /mnt/DATA6/mouseData/ && mkdir -p /usr/bin/ && \
    pip3 install pandas && \
    cp /opt/R-3.6.3/bin/R /usr/bin/ && \
    cp /opt/R-3.6.3/bin/Rscript /usr/bin/ 
    

COPY 20200908mouseCnaScript.R 20210517segmentationScript.R 20210521snakeFileInputs.R cnSnakefileV2 /mnt/DATA6/mouseData/

WORKDIR /mnt/DATA6/mouseData/

CMD /usr/bin/Rscript --vanilla 20210521snakeFileInputs.R && \
    snakemake --snakefile cnSnakefileV2 --jobs 3 
