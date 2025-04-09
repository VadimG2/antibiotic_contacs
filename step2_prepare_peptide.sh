#!/bin/bash

PEPTIDE_NAME="pept1"
OUTPUT_DIR="output/${PEPTIDE_NAME}"
mkdir -p "${OUTPUT_DIR}"

gmx_mpi pdb2gmx -f "${PEPTIDE_NAME}.pdb" -o "${OUTPUT_DIR}/${PEPTIDE_NAME}_optimized.gro" -p "${OUTPUT_DIR}/${PEPTIDE_NAME}.top" -ignh -ff amber14sb -water tip3p

# Извлечение чистой топологии молекулы
awk 'BEGIN {printit=0} /^\[ moleculetype \]/ {printit=1} printit && /^\[ system \]|^#include.*tip3p.itp|^#include.*ions.itp/ {printit=0} printit' "${OUTPUT_DIR}/${PEPTIDE_NAME}.top" > "${OUTPUT_DIR}/${PEPTIDE_NAME}_mol.itp"

echo "Топология пептида подготовлена в ${OUTPUT_DIR}"
