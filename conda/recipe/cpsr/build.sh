#!/bin/bash

export DISABLE_AUTOBREW=1
${R} -e "install.packages('devtools', repos = 'https://cloud.r-project.org/', lib = '${PREFIX}/lib/R/library')"
${R} -e "devtools::install_github(repo = 'JanMarvin/openxlsx2', ref = 'v1.6', lib = '${PREFIX}/lib/R/library')"
${R} CMD INSTALL --build . ${R_ARGS}
