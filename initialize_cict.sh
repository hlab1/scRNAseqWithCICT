#!/bin/bash

# Building docker for the different algorithms 
echo "Initialize bare singularity image and overlay for CICT"

BASEDIR=$(pwd)
SIF_DIR_ARCHIVE=/archive/ch153/envs/singularity
MC_DIR_ARCHIVE=/archive/ch153/envs/miniconda
EXT3_DIR_ARCHIVE=/archive/ch153/envs/overlay-fs-ext3
CONDA_DIR=${BASEDIR}/conda_greene
SIF_DIR=${BASEDIR}/images_singularity
EXT3_DIR=${BASEDIR}/images_singularity_overlay
SIF_FILE=centos-8.2.2004.sif
EXT3_FILE=overlay-5GB-200K.ext3.gz
MC3_INSTALL_SCRIPT=Miniconda3-py37_4.10.3-Linux-x86_64.sh
MC3_OVERLAY_SCRIPT=overlay_ext3_mc3.sh


mkdir -p ${SIF_DIR} ${EXT3_DIR} ${CONDA_DIR}
cp -u ${MC_DIR_ARCHIVE}/${MC3_INSTALL_SCRIPT} ${CONDA_DIR}
cp -u ${EXT3_DIR_ARCHIVE}/${MC3_OVERLAY_SCRIPT} ${CONDA_DIR}
cp -u ${SIF_DIR_ARCHIVE}/${SIF_FILE} ${CONDA_DIR}
cp -u ${EXT3_DIR_ARCHIVE}/${EXT3_FILE} ${CONDA_DIR}

ALG=CICT
cp -u ${CONDA_DIR}/${SIF_FILE} ${SIF_DIR}/${ALG}.sif
cp -rp ${CONDA_DIR}/${EXT3_FILE} ${EXT3_DIR}/${ALG}.ext3.gz
gunzip ${EXT3_DIR}/${ALG}.ext3.gz
cd ${BASEDIR}; singularity exec --overlay ${EXT3_DIR}/${ALG}.ext3 ${SIF_DIR}/${ALG}.sif \
			   /bin/sh -c "
sh ${CONDA_DIR}/Miniconda3-py37_4.10.3-Linux-x86_64.sh -b -p /ext3/miniconda3
cp ${CONDA_DIR}/overlay_ext3_mc3.sh /ext3/env.sh
source /ext3/env.sh
conda update -n base conda -y
conda clean --all --yes
conda install -n base --override-channels -c conda-forge mamba=0.27.0 'python_abi=*=*cp*' #conda install -c conda-forge mamba 
mamba create -y --name ${ALG} -c conda-forge r-base=4.1.2
conda activate ${ALG}
mamba install -y -c conda-forge libgit2 gmp time 
mamba install -c conda-forge r-ragg
R -e \"install.packages('remotes',repos='https://cloud.r-project.org')\"
"
echo "Singularity files for ${ALG}: image is ${SIF_DIR}/${ALG}.sif, overlay is ${EXT3_DIR}/${ALG}.ext3"
# CICT_RENV; setting up conda environment with renv
cd ${BASEDIR}; singularity exec --overlay ${EXT3_DIR}/${ALG}.ext3 ${SIF_DIR}/${ALG}.sif \
			   /bin/sh -c "
mamba create -y --name CICT_RENV -c conda-forge r-base=4.1.2
conda activate CICT_RENV
mamba install -y -c conda-forge libgit2 gmp time
R -e \"install.packages('renv',repos='https://cloud.r-project.org')\"
"
###To initialize lock file
#PKG_CONFIG_PATH=/ext3/miniconda3/envs/CICT_RENV/lib/pkgconfig R --no-restore
#In R: renv::init('.',bioconductor=TRUE); renv::install(c("Bioc::RCy3","Bioc::WGCNA","Bioc::minet")); renv::install("kmlShape@0.9.5"); renv::snapshot()
###To restore from lock file, simplly run
# PKG_CONFIG_PATH=/ext3/miniconda3/envs/CICT_RENV/lib/pkgconfig R --no-restore
# In R: renv::restore()