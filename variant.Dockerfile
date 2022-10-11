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
    pip3 install snakemake && pip3 install pandas && \
    wget https://cran.r-project.org/bin/linux/ubuntu/bionic-cran35/r-base_3.6.3.orig.tar.gz && \
    tar -xvzf /opt/r-base_3.6.3.orig.tar.gz && cd R-3.6.3/ &&  ./configure && make && cd /opt/ && \
    mkdir -p /mnt/DATA6/mouseData/ && mkdir -p /usr/bin/ && mkdir -p /mnt/tmp/ && \
    cp /opt/R-3.6.3/bin/R /usr/bin/ && cp /opt/R-3.6.3/bin/Rscript /usr/bin/ && \
    wget https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2 -O bcftools.tar.bz2 && \
    tar -xjvf bcftools.tar.bz2 && cd bcftools-1.3.1/ && make && make prefix=/usr/local/bin install && \
    ln -s /usr/local/bin/bin/bcftools /usr/bin/bcftools && \
    cd /opt && wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz && \
    tar -xvzf annovar.latest.tar.gz

# to create mousedbs for annovar & required R libraries
RUN perl /opt/annovar/annotate_variation.pl -downdb -buildver mm10 -webfrom annovar refGene mousedb/ && \
    perl /opt/annovar/annotate_variation.pl --buildver mm10 --downdb seq mousedb/mm10_seq && \
    perl /opt/annovar/retrieve_seq_from_fasta.pl mousedb/mm10_refGene.txt -seqdir mousedb/mm10_seq -format refGene -outfile mousedb/mm10_refGeneMrna.fa && \
    /usr/bin/R -e "install.packages(c('stringr', 'optparse'), dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
    cd /mnt/DATA6/mouseData/ && apt-get update && apt-get install --yes tabix

COPY 20210520vcfInputs.R annotationSnakefileV2 20210524processingVarAnno.R combinedAnno /mnt/DATA6/mouseData/
COPY mm10_mgp.v6.combinedMouseFilt.txt /opt/mousedb/

WORKDIR /mnt/DATA6/mouseData/

# when I add genotype steps make sure I have snakemake.ignore so combined vcf isn't processed
# loop may ensue i.e thinks there is new vcf from combined -> new combined

CMD /usr/bin/Rscript --vanilla 20210520vcfInputs.R && \
    snakemake --snakefile annotationSnakefileV2 -k --jobs 10 || true && \
    snakemake --snakefile combinedAnno -k --jobs 5 
