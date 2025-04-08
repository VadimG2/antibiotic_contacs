#!/bin/bash

# --- Директории ---
OUTPUT_BASE_DIR="output"
SIMULATIONS_DIR="simulations"

# --- Параметры сборки ---
BOX_DISTANCE=2.0
SALT_CONCENTRATION=0.15
SOLVENT_MODEL="spc216.gro"
FORCEFIELD_DIR="amber14sb.ff"

echo "--- НАЧАЛО: Сборка всех систем Пептид-Лиганд ---"
mkdir -p "${SIMULATIONS_DIR}" || { echo "ОШИБКА: Не удалось создать ${SIMULATIONS_DIR}"; exit 1; }

# --- Поиск директорий пептидов и лигандов ---
peptide_dirs=()
ligand_dirs=()
for d in "${OUTPUT_BASE_DIR}"/*/; do
    if [ -n "$(find "${d}" -maxdepth 1 -name '*_optimized.gro' -print -quit)" ]; then
        peptide_dirs+=("${d%/}")
    fi
done
for d in "${OUTPUT_BASE_DIR}"/*/; do
    if [ -n "$(find "${d}" -maxdepth 1 -name '*_GMX.gro' -print -quit)" ]; then
        ITP_FILE=$(find "${d}" -maxdepth 1 -name '*_bcc_manual.itp' -print -quit)
        POSRES_FILE=$(find "${d}" -maxdepth 1 -name 'posre_*.itp' -print -quit)
        if [ -n "${ITP_FILE}" ] && [ -n "${POSRES_FILE}" ]; then
            ligand_dirs+=("${d%/}")
        else
            echo "ВНИМАНИЕ: В ${d} отсутствуют .itp или posre файлы. Пропускаем."
        fi
    fi
done
if [ ${#peptide_dirs[@]} -eq 0 ]; then echo "ОШИБКА: Не найдено пептидов"; exit 1; fi
if [ ${#ligand_dirs[@]} -eq 0 ]; then echo "ОШИБКА: Не найдено лигандов"; exit 1; fi
echo "Найдено пептидов: ${#peptide_dirs[@]}, Лигандов: ${#ligand_dirs[@]}"

# --- Основной цикл ---
for PEPTIDE_DIR in "${peptide_dirs[@]}"; do
    PEPTIDE_BASENAME=$(basename "${PEPTIDE_DIR}")
    PEPTIDE_OPT_GRO="${PEPTIDE_DIR}/${PEPTIDE_BASENAME}_optimized.gro"
    PEPTIDE_TOP_SOURCE="${PEPTIDE_DIR}/${PEPTIDE_BASENAME}.top"
    PEPTIDE_POSRES_SOURCE=$(find "${PEPTIDE_DIR}" -maxdepth 1 -name 'posre_*.itp' -print -quit)
    PEPTIDE_ITP_SOURCE="${PEPTIDE_DIR}/${PEPTIDE_BASENAME}.itp"

    if [ ! -f "${PEPTIDE_ITP_SOURCE}" ]; then
        echo "--> Создание ${PEPTIDE_ITP_SOURCE} из ${PEPTIDE_TOP_SOURCE}"
        awk '/^\[ moleculetype \]/,/^\[ system \]/' "${PEPTIDE_TOP_SOURCE}" | sed '$d' > "${PEPTIDE_ITP_SOURCE}"
        if [ $? -ne 0 ] || [ ! -s "${PEPTIDE_ITP_SOURCE}" ]; then
            echo "ОШИБКА извлечения moleculetype пептида из ${PEPTIDE_TOP_SOURCE}"
            continue
        fi
    fi
    PEPTIDE_MOL_NAME=$(awk '/^\[ moleculetype \]/{getline; print $1; exit}' "${PEPTIDE_ITP_SOURCE}" || echo "Protein_${PEPTIDE_BASENAME}")

    if [ ! -f "${PEPTIDE_OPT_GRO}" ] || [ ! -f "${PEPTIDE_ITP_SOURCE}" ] || [ ! -f "${PEPTIDE_POSRES_SOURCE}" ]; then
        echo "ВНИМАНИЕ: Отсутствуют .gro, .itp или posre файлы для пептида ${PEPTIDE_BASENAME}. Пропускаем."
        continue
    fi

    for LIGAND_DIR in "${ligand_dirs[@]}"; do
        LIGAND_BASENAME=$(basename "${LIGAND_DIR}")
        LIGAND_GRO_SOURCE=$(find "${LIGAND_DIR}" -maxdepth 1 -name '*_GMX.gro' -print -quit)
        LIGAND_ITP_SOURCE=$(find "${LIGAND_DIR}" -maxdepth 1 -name '*_bcc_manual.itp' -print -quit)
        LIGAND_POSRES_SOURCE=$(find "${LIGAND_DIR}" -maxdepth 1 -name 'posre_*.itp' -print -quit)
        LIGAND_MOL_NAME=$(awk '/^\[ moleculetype \]/{getline; print $1; exit}' "${LIGAND_ITP_SOURCE}" || echo "${LIGAND_BASENAME}")

        if [ ! -f "${LIGAND_GRO_SOURCE}" ] || [ ! -f "${LIGAND_ITP_SOURCE}" ] || [ ! -f "${LIGAND_POSRES_SOURCE}" ]; then
            echo "ВНИМАНИЕ: Отсутствуют файлы для лиганда ${LIGAND_BASENAME}. Пропускаем комбинацию."
            continue
        fi

        SYSTEM_NAME="${PEPTIDE_BASENAME}_${LIGAND_BASENAME}"
        SIMULATION_OUTPUT_DIR="${SIMULATIONS_DIR}/${SYSTEM_NAME}"

        echo ""
        echo "--- Сборка системы: ${SYSTEM_NAME} в ${SIMULATION_OUTPUT_DIR} ---"
        mkdir -p "${SIMULATION_OUTPUT_DIR}" || { echo "ОШИБКА: Не удалось создать ${SIMULATION_OUTPUT_DIR}"; continue; }

        FINAL_TOPOLOGY="${SIMULATION_OUTPUT_DIR}/${SYSTEM_NAME}.top"
        FINAL_COORDS="${SIMULATION_OUTPUT_DIR}/${SYSTEM_NAME}_solv_ions.gro"
        FINAL_INDEX="${SIMULATION_OUTPUT_DIR}/index.ndx"
        PEPTIDE_ITP_LOCAL="${SIMULATION_OUTPUT_DIR}/$(basename ${PEPTIDE_ITP_SOURCE})"
        PEPTIDE_POSRES_LOCAL="${SIMULATION_OUTPUT_DIR}/$(basename ${PEPTIDE_POSRES_SOURCE})"
        LIGAND_ITP_LOCAL="${SIMULATION_OUTPUT_DIR}/$(basename ${LIGAND_ITP_SOURCE})"
        LIGAND_POSRES_LOCAL="${SIMULATION_OUTPUT_DIR}/$(basename ${LIGAND_POSRES_SOURCE})"
        LIGAND_GRO_LOCAL="${SIMULATION_OUTPUT_DIR}/$(basename ${LIGAND_GRO_SOURCE})"
        PEPTIDE_CENTER_GRO="${SIMULATION_OUTPUT_DIR}/${PEPTIDE_BASENAME}_center.gro"
        COMBINED_GRO="${SIMULATION_OUTPUT_DIR}/${SYSTEM_NAME}_combined.gro"
        BOX_GRO="${SIMULATION_OUTPUT_DIR}/${SYSTEM_NAME}_box.gro"
        SOLV_GRO="${SIMULATION_OUTPUT_DIR}/${SYSTEM_NAME}_solv.gro"
        IONS_TPR="${SIMULATION_OUTPUT_DIR}/${SYSTEM_NAME}_ions.tpr"
        IONS_MDP="${SIMULATION_OUTPUT_DIR}/ions.mdp"

        # --- 5.A: Копирование файлов ---
        echo "--> 5.A: Копирование файлов в ${SIMULATION_OUTPUT_DIR}"
        cp "${PEPTIDE_ITP_SOURCE}" "${PEPTIDE_ITP_LOCAL}" || { echo "ОШИБКА копирования ${PEPTIDE_ITP_SOURCE}"; continue; }
        cp "${PEPTIDE_POSRES_SOURCE}" "${PEPTIDE_POSRES_LOCAL}" || { echo "ОШИБКА копирования ${PEPTIDE_POSRES_SOURCE}"; continue; }
        cp "${LIGAND_ITP_SOURCE}" "${LIGAND_ITP_LOCAL}" || { echo "ОШИБКА копирования ${LIGAND_ITP_SOURCE}"; continue; }
        cp "${LIGAND_POSRES_SOURCE}" "${LIGAND_POSRES_LOCAL}" || { echo "ОШИБКА копирования ${LIGAND_POSRES_SOURCE}"; continue; }
        cp "${LIGAND_GRO_SOURCE}" "${LIGAND_GRO_LOCAL}" || { echo "ОШИБКА копирования ${LIGAND_GRO_SOURCE}"; continue; }

        # --- 5.B: Создание топологии ---
        echo "--> 5.B: Создание ${FINAL_TOPOLOGY}"
        cat << EOF > "${FINAL_TOPOLOGY}"
; Основной файл топологии для системы ${SYSTEM_NAME}

; Включаем основное силовое поле (AMBER14SB)
#include "${FORCEFIELD_DIR}/forcefield.itp"

; Включаем топологию лиганда
#include "./$(basename ${LIGAND_ITP_LOCAL})"

; Включаем топологию пептида
#include "./$(basename ${PEPTIDE_ITP_LOCAL})"

; Включаем модель воды TIP3P
#include "${FORCEFIELD_DIR}/tip3p.itp"

[ system ]
${SYSTEM_NAME} in water

[ molecules ]
; Compound        nmols
${PEPTIDE_MOL_NAME}     1
${LIGAND_MOL_NAME}      1
; Вода и ионы будут добавлены автоматически
EOF
        if [ $? -ne 0 ] || [ ! -s "${FINAL_TOPOLOGY}" ]; then echo "ОШИБКА создания ${FINAL_TOPOLOGY}"; continue; fi

        # --- 5.C: Сборка координат ---
        echo "--> 5.C: Сборка координат системы..."
        gmx_mpi editconf -f "${PEPTIDE_OPT_GRO}" -o "${PEPTIDE_CENTER_GRO}" -c -box 6.0 6.0 6.0 || { echo "ОШИБКА editconf"; continue; }
        gmx_mpi insert-molecules -f "${PEPTIDE_CENTER_GRO}" -ci "${LIGAND_GRO_LOCAL}" -nmol 1 -o "${COMBINED_GRO}" -try 50 -radius 0.2 -dr 1.0 1.0 1.0 || { echo "ОШИБКА insert-molecules"; continue; }
        gmx_mpi editconf -f "${COMBINED_GRO}" -o "${BOX_GRO}" -c -d ${BOX_DISTANCE} -bt cubic || { echo "ОШИБКА editconf"; continue; }

        # --- 5.D: Сольватация и ионизация ---
        echo "--> 5.D: Сольватация и ионизация..."
        cd "${SIMULATION_OUTPUT_DIR}" || { echo "ОШИБКА: Не удалось перейти в ${SIMULATION_OUTPUT_DIR}"; continue; }

        gmx_mpi solvate -cp "$(basename ${BOX_GRO})" -cs ${SOLVENT_MODEL} -o "$(basename ${SOLV_GRO})" -p "$(basename ${FINAL_TOPOLOGY})" || { echo "ОШИБКА solvate"; cd -; continue; }
        NUM_SOL=$(grep "^SOL" "$(basename ${FINAL_TOPOLOGY})" | awk '{print $2}')
        echo "Добавлено ${NUM_SOL} молекул воды."

        printf "integrator = md\nnsteps = 0\n" > ions.mdp
        gmx_mpi grompp -f ions.mdp -c "$(basename ${SOLV_GRO})" -p "$(basename ${FINAL_TOPOLOGY})" -o "$(basename ${IONS_TPR})" -maxwarn 1 || { echo "ОШИБКА grompp"; cd -; continue; }
        echo "SOL" | gmx_mpi genion -s "$(basename ${IONS_TPR})" -o "$(basename ${FINAL_COORDS})" -p "$(basename ${FINAL_TOPOLOGY})" -pname NA -nname CL -neutral -conc ${SALT_CONCENTRATION} || { echo "ОШИБКА genion"; cd -; continue; }
        echo "--> Ионы добавлены."

        # --- 5.E: Создание индексного файла ---
        echo "--> 5.E: Создание ${FINAL_INDEX}"
        NDX_LIST_OUTPUT=$(echo "q" | gmx_mpi make_ndx -f "$(basename ${FINAL_COORDS})" -o junk_temp.ndx 2>&1)
        rm -f junk_temp.ndx
        PROTEIN_GROUP_NUM=$(echo "${NDX_LIST_OUTPUT}" | grep -E '^[[:space:]]*[0-9]+[[:space:]]+Protein[[:space:]]+\(' | awk '{print $1}')
        LIGAND_GROUP_NUM=$(echo "${NDX_LIST_OUTPUT}" | grep -E '^[[:space:]]*[0-9]+[[:space:]]+'"${LIGAND_MOL_NAME}"'[[:space:]]+\(' | awk '{print $1}')
        SOL_GROUP_NUM=$(echo "${NDX_LIST_OUTPUT}" | grep -E '^[[:space:]]*[0-9]+[[:space:]]+SOL[[:space:]]+\(' | awk '{print $1}')
        NA_GROUP_NUM=$(echo "${NDX_LIST_OUTPUT}" | grep -E '^[[:space:]]*[0-9]+[[:space:]]+NA[[:space:]]+\(' | awk '{print $1}')
        CL_GROUP_NUM=$(echo "${NDX_LIST_OUTPUT}" | grep -E '^[[:space:]]*[0-9]+[[:space:]]+CL[[:space:]]+\(' | awk '{print $1}')
        [[ -z "${PROTEIN_GROUP_NUM}" ]] && { echo "ОШИБКА: Не найдена группа Protein"; cd -; continue; }
        [[ -z "${LIGAND_GROUP_NUM}" ]] && { echo "ОШИБКА: Не найдена группа ${LIGAND_MOL_NAME}"; cd -; continue; }
        MAKE_NDX_COMMANDS="${PROTEIN_GROUP_NUM} | ${LIGAND_GROUP_NUM}\nname Protein_Ligand\n${SOL_GROUP_NUM} | ${NA_GROUP_NUM} | ${CL_GROUP_NUM}\nname Water_and_ions\nq\n"
        echo -e "${MAKE_NDX_COMMANDS}" | gmx_mpi make_ndx -f "$(basename ${FINAL_COORDS})" -o "$(basename ${FINAL_INDEX})" || { echo "ОШИБКА make_ndx"; cd -; continue; }

        # --- 5.F: Очистка ---
        echo "--> 5.F: Очистка временных файлов..."
        rm -f "$(basename ${PEPTIDE_CENTER_GRO})" "$(basename ${COMBINED_GRO})" "$(basename ${BOX_GRO})" "$(basename ${SOLV_GRO})" "$(basename ${IONS_TPR})" ions.mdp mdout*.mdp \#*#
        cd - || exit 1

        echo "--- Сборка системы ${SYSTEM_NAME} УСПЕШНО ЗАВЕРШЕНА ---"
    done
done

echo "--- СБОРКА ВСЕХ СИСТЕМ завершена ---"
echo "Результаты в ${SIMULATIONS_DIR}"
