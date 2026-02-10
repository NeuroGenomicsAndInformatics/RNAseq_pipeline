################################################################################
# Dockerfile for DCC
################################################################################
FROM python:3.8

ENV PYTHONPATH "${PYTHONPATH}:/usr/local/bin"
#ENV PATH /usr/local/bin:$PATH
# Get DCC pacakge dependencies
RUN pip install numpy==1.21.4 pandas==1.3.4 pysam==0.18.0 HTSeq==1.99.2 pybedtools


# Get DCC 
RUN cd opt && wget --no-check-certificate https://github.com/dieterich-lab/DCC/archive/refs/tags/v0.5.0.tar.gz && \
    tar -xvf v0.5.0.tar.gz && rm v0.5.0.tar.gz 

# Copy in the changed script for memory problem
COPY circFilter.py /opt/DCC-0.5.0/DCC
COPY genecount.py /opt/DCC-0.5.0/DCC

# then run the setup script 
RUN cd opt && cd DCC-0.5.0 && python setup.py install 

# chmod the script to get the seq of the circ coords 
RUN chmod +rx /opt/DCC-0.5.0/scripts/*

# get bedtools 
RUN cd /opt && \
    wget --no-check-certificate https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz && \
    tar -xf bedtools-2.30.0.tar.gz && rm bedtools-2.30.0.tar.gz && \
    cd bedtools2 && make && make install && make clean


ENV PATH /usr/local/bin:$PATH
