#!/bin/bash

LIGAND_NAME="lig1"
OUTPUT_DIR="output/${LIGAND_NAME}"
mkdir -p "${OUTPUT_DIR}"

# Предполагается, что lig1_bcc_manual.itp уже создан с GAFF2
LIGAND_ITP_SOURCE="${OUTPUT_DIR}/${LIGAND_NAME}_bcc_manual.itp"
LIGAND_GRO="${OUTPUT_DIR}/${LIGAND_NAME}_GMX.gro"
LIGAND_POSRES="${OUTPUT_DIR}/posre_${LIGAND_NAME}.itp"

# Очистка топологии лиганда
awk 'BEGIN {skip=0} /^\[ atomtypes \]/ {skip=1} skip && /^\[/ {skip=0} !skip' "${LIGAND_ITP_SOURCE}" | sed '/#include.*atomtypes/d' > "${OUTPUT_DIR}/clean_${LIGAND_NAME}_bcc_manual.itp"

echo "Топология лиганда очищена в ${OUTPUT_DIR}"
