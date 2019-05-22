#!/usr/bin/env bash

R -e "library(devtools); devtools::install_github('hms-dbmi/UpSetR',   dependencies=FALSE, args=c('--library=${PREFIX}/lib/R/library'))"
R -e "library(devtools); devtools::install_github('kassambara/ggpubr', dependencies=FALSE, args=c('--library=${PREFIX}/lib/R/library'))"

mkdir -p ${PREFIX}/bin
chmod +x ${SRC_DIR}/*.py
mv ${SRC_DIR}/*.py ${PREFIX}/bin/
