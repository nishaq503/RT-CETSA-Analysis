FROM r-base:4.4.0

	
ARG EXEC_DIR="/opt/executables"
ARG DATA_DIR="/data"

RUN apt-get update && apt-get install -y python3

#Create folders
RUN mkdir -p /{EXEC_DIR} \
    && mkdir -p /{DATA_DIR}/csv \
    && mkdir -p /{DATA_DIR}/outputs

# Install R packages external dependencies
RUN apt-get install -y libssl-dev \
    && apt-get install -y libcurl4-openssl-dev libfontconfig1-dev \
    && apt-get install -y libharfbuzz-dev libfribidi-dev \
    && apt-get install -y libxml2-dev \
    && apt-get install -y libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    && apt-get install -y cmake 

# Install required R packages
ADD requirements/Requirements_1.R /{EXEC_DIR}/Requirements_1.R
RUN Rscript /{EXEC_DIR}/Requirements_1.R

ADD requirements/Requirements_2.R /{EXEC_DIR}/Requirements_2.R
RUN Rscript /{EXEC_DIR}/Requirements_2.R

ADD requirements/Requirements_3.R /{EXEC_DIR}/Requirements_3.R
RUN Rscript /{EXEC_DIR}/Requirements_3.R

#Copy executable
#TODO FIX extra / in front of EXCEC_DIR (will need to rebuild the whole container!)
COPY src/ /{EXEC_DIR}/

WORKDIR /{EXEC_DIR}

CMD ["Rscript", "main.R"]