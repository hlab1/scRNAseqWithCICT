#!/bin/bash

BASEDIR=${PWD}
SIF_DIR=${BASEDIR}/images_singularity
EXT3_DIR=${BASEDIR}/images_singularity_overlay
SIF_FILE=centos-8.2.2004.sif
EXT3_FILE=overlay-5GB-200K.ext3.gz
ALG=CICT

singularity exec --overlay ${EXT3_DIR}/${ALG}.ext3 ${SIF_DIR}/${ALG}.sif /bin/bash
