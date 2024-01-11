FROM base_pipeline

WORKDIR /opt

RUN /usr/bin/R -e "install.packages(c('BiocManager', 'stringr', 'optparse', 'foreach', 'doParallel', 'jsonlite'), dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
    /usr/bin/R -e "BiocManager::install(c('DNAcopy','GenomicRanges'))"    

COPY 20200908mouseCnaScript.R 20210517segmentationScript.R 20210521snakeFileInputs.R cnSnakefileV2 /mnt/DATA6/mouseData/

WORKDIR /mnt/DATA6/mouseData/

#CMD /usr/bin/Rscript --vanilla 20210521snakeFileInputs.R && \
#    snakemake --snakefile cnSnakefileV2 -k --jobs 3 
