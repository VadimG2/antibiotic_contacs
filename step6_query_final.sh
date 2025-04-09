#!/bin/bash


mv ./output/*+* ./simulations
# --- Базовые директории ---
BASE_DIR=$(pwd)
OUTPUT_DIR="${BASE_DIR}/output" # Папка с ИСХОДНЫМИ itp/gro пептидов/лигандов
SIMULATIONS_DIR="${BASE_DIR}/simulations" # Папка с папками для ЗАПУСКА симуляций

# --- Имена MDP файлов ---
MDP_MIN="${BASE_DIR}/step6.1_minimization.mdp"
MDP_NVT="${BASE_DIR}/step6.2_equilibration_nvt.mdp"
MDP_NPT="${BASE_DIR}/step6.3_equilibration_npt.mdp"
MDP_PROD="${BASE_DIR}/step6.4_production.mdp"

echo "============================================="
echo "    ЗАПУСК КОНВЕЙЕРА МД СИМУЛЯЦИЙ (Шаг 6) v18 - POSRES в MDP"
echo "============================================="

# --- Проверка наличия MDP файлов ---
for mdp_file in "${MDP_MIN}" "${MDP_NVT}" "${MDP_NPT}" "${MDP_PROD}"; do
    if [ ! -f "${mdp_file}" ]; then echo "ОШИБКА: MDP файл ${mdp_file} не найден!"; exit 1; fi
     # Проверяем наличие define = -DPOSRES в NVT/NPT
     if [[ "${mdp_file}" == *nvt.mdp || "${mdp_file}" == *npt.mdp ]]; then
          if ! grep -q "^define *= *-DPOSRES" "${mdp_file}"; then
               echo "ОШИБКА: В файле ${mdp_file} отсутствует или закомментирована строка 'define = -DPOSRES'. Пожалуйста, добавьте/раскомментируйте ее."
               exit 1
          fi
     fi
     # Проверяем отсутствие define = -DPOSRES в Prod
     if [[ "${mdp_file}" == *production.mdp ]]; then
          if grep -q "^define *= *-DPOSRES" "${mdp_file}"; then
               echo "ОШИБКА: В файле ${mdp_file} присутствует активная строка 'define = -DPOSRES'. Пожалуйста, закомментируйте/удалите ее."
               exit 1
          fi
     fi
done

if [ ! -d "${SIMULATIONS_DIR}" ]; then echo "ОШИБКА: Директория ${SIMULATIONS_DIR} не найдена!"; exit 1; fi

# --- Итерация по папкам с системами в SIMULATIONS_DIR ---
for RUN_DIR in ${SIMULATIONS_DIR}/*; do # Ищем папки внутри simulations/
    if [ ! -d "${RUN_DIR}" ]; then continue; fi

    # --- ИСПРАВЛЕНО: Извлечение имен пептида и лиганда ---
    SYS_NAME_FULL=$(basename "${RUN_DIR}") # Получаем "pept1+lig1_Sim" или "pept1+lig1"
    SYS_NAME=${SYS_NAME_FULL%_Sim} # Удаляем "_Sim" в конце, если есть -> "pept1+lig1"
    PROTEIN_NAME=$(echo "${SYS_NAME}" | awk -F'+' '{print $1}') # -> "pept1"
    LIGAND_NAME=$(echo "${SYS_NAME}" | awk -F'+' '{print $2}') # -> "lig1"

    if [ -z "${PROTEIN_NAME}" ] || [ -z "${LIGAND_NAME}" ]; then
        echo "ПРЕДУПРЕЖДЕНИЕ: Не удалось извлечь имена пептида/лиганда из ${SYS_NAME_FULL}. Пропускаем."
        continue
    fi

    PROTEIN_PREP_DIR="${OUTPUT_DIR}/${PROTEIN_NAME}"
    LIGAND_PREP_DIR="${OUTPUT_DIR}/${LIGAND_NAME}"

    echo ""
    echo "--- Обработка системы: ${SYS_NAME} ---"
    cd "${RUN_DIR}" || { echo "ОШИБКА: Не удалось перейти в ${RUN_DIR}"; cd ${BASE_DIR}; continue; }
    echo "Рабочая директория: $(pwd)"

    # --- Определяем имена файлов ---
    GRO_START_FILE="system_solv_ions.gro"
    TOP_FILE="system.top"
    LIGAND_ITP_TARGET="ligand.itp"
    ATOMTYPES_ITP_TARGET="gaff2_atomtypes.itp"
    POSRE_LIGAND_TARGET="posre_${LIGAND_NAME}.itp"
    POSRE_PROTEIN_TARGET="posre_protein.itp"
    LIGAND_ITP_SOURCE="${LIGAND_PREP_DIR}/${LIGAND_NAME}.itp"
    LIGAND_ATOMTYPES_SOURCE="${LIGAND_PREP_DIR}/atomtypes_${LIGAND_NAME}.itp"
    POSRE_LIGAND_SOURCE="${LIGAND_PREP_DIR}/posre_${LIGAND_NAME}.itp"
    POSRE_PROTEIN_SOURCE="${PROTEIN_PREP_DIR}/posre_protein.itp"

    # --- Проверка и копирование недостающих ITP файлов ---
    # ... (Блок копирования ITP как в v16) ...
    echo "--> Проверка и копирование ITP файлов в $(pwd)..."
    FILES_OK=1
    declare -A FILES_TO_COPY=(
        ["${LIGAND_ITP_SOURCE}"]="${LIGAND_ITP_TARGET}"
        ["${LIGAND_ATOMTYPES_SOURCE}"]="${ATOMTYPES_ITP_TARGET}"
        ["${POSRE_LIGAND_SOURCE}"]="${POSRE_LIGAND_TARGET}"
        ["${POSRE_PROTEIN_SOURCE}"]="${POSRE_PROTEIN_TARGET}"
    )
    for source_file in "${!FILES_TO_COPY[@]}"; do
        target_file="${FILES_TO_COPY[$source_file]}"
        if [ ! -f "${target_file}" ]; then
             if [ -f "${source_file}" ]; then
                 echo "Копирование ${source_file} -> ${target_file}"
                 cp "${source_file}" "${target_file}"
                 if [ $? -ne 0 ]; then echo "ОШИБКА: Не удалось скопировать ${source_file}"; FILES_OK=0; break; fi
             else
                 echo "ОШИБКА: Не найден исходный файл ${source_file}"
                 FILES_OK=0; break
             fi
        fi
    done
    if [ ${FILES_OK} -eq 0 ]; then echo "ОШИБКА: Отсутствуют или не удалось скопировать ITP файлы. Пропускаем."; cd ${BASE_DIR}; continue; fi
    echo "Все необходимые ITP файлы присутствуют."

    # --- Проверка наличия стартовых файлов симуляции ---
    if [ ! -f "${GRO_START_FILE}" ] || [ ! -f "${TOP_FILE}" ]; then
        echo "ОШИБКА: В директории ${RUN_DIR} отсутствуют файлы ${GRO_START_FILE} или ${TOP_FILE}. Пропускаем."
        cd ${BASE_DIR}; continue;
    fi

    # --- Создание индексного файла ---
    # ... (Блок создания index.ndx как в v16) ...
    echo "--> Создание индексного файла index.ndx..."
    if [ ! -f "index.ndx" ]; then
        MAKE_NDX_CMDS=$(cat <<EOF
5 | 2
name 21 Protein_Ligand
17 | 3 | 4
name 22 SOL_Ion
q
EOF
)
        echo "Попытка создать группы 'Protein_Ligand' (5|2) и 'SOL_Ion' (17|3|4)"
        gmx_mpi make_ndx -f "${GRO_START_FILE}" -o index.ndx <<< "${MAKE_NDX_CMDS}" > make_ndx.log 2>&1
        NDX_EXIT_CODE=$?
        if [ ${NDX_EXIT_CODE} -ne 0 ] || ! grep -q "\[ Protein_Ligand \]" index.ndx || ! grep -q "\[ SOL_Ion \]" index.ndx; then
             echo "ОШИБКА: Не удалось автоматически создать группы в index.ndx. Пожалуйста, создайте их вручную в $(pwd)."
             cat make_ndx.log
             cd ${BASE_DIR}; continue
        fi
        echo "Индексный файл index.ndx создан."
        rm make_ndx.log
    else
        echo "Индексный файл index.ndx уже существует."
    fi
    # Определяем имена групп для использования в MDP файлах (если понадобится адаптация)
    # SOLUTE_GROUP="Protein_Ligand"
    # SOLVENT_GROUP="SOL_Ion"

    # --- БЛОК КОПИРОВАНИЯ/АДАПТАЦИИ MDP УДАЛЕН ---

    # --- Шаг 6.1: Минимизация ---
    echo "--> Шаг 6.1: Минимизация энергии для ${SYS_NAME}..."
    MIN_TPR="${SYS_NAME}_min.tpr"
    MIN_GRO="${SYS_NAME}_min.gro"
    if [ ! -f "${MIN_GRO}" ]; then
        # Используем ПРЯМОЙ путь к MDP
        gmx_mpi grompp -f "${MDP_MIN}" -c "${GRO_START_FILE}" -p "${TOP_FILE}" -o "${MIN_TPR}" -maxwarn 5
        if [ $? -ne 0 ]; then echo "ОШИБКА: grompp минимизации"; cd ${BASE_DIR}; continue; fi
        gmx_mpi mdrun -deffnm "${SYS_NAME}_min" -v -nb gpu
        if [ $? -ne 0 ]; then echo "ОШИБКА: mdrun минимизации"; cd ${BASE_DIR}; continue; fi
        echo "Минимизация завершена."
    else
        echo "Минимизация уже выполнена, пропускаем."
        MIN_TPR="${SYS_NAME}_min.tpr"
    fi

    # --- Шаг 6.2: NVT Уравновешивание ---
    echo "--> Шаг 6.2: NVT Уравновешивание для ${SYS_NAME}..."
    NVT_TPR="${SYS_NAME}_nvt.tpr"
    NVT_GRO="${SYS_NAME}_nvt.gro"
    NVT_CPT="${SYS_NAME}_nvt.cpt"
    if [ ! -f "${NVT_GRO}" ]; then
        # Используем ПРЯМОЙ путь к MDP. Grompp сам прочитает 'define = -DPOSRES' из файла MDP_NVT.
        gmx_mpi grompp -f "${MDP_NVT}" -c "${MIN_GRO}" -r "${MIN_GRO}" -p "${TOP_FILE}" -n index.ndx -o "${NVT_TPR}" -maxwarn 5
        if [ $? -ne 0 ]; then echo "ОШИБКА: grompp NVT"; cd ${BASE_DIR}; continue; fi
        gmx_mpi mdrun -deffnm "${SYS_NAME}_nvt" -v -nb gpu -pme gpu -update gpu -bonded gpu
        if [ $? -ne 0 ]; then echo "ОШИБКА: mdrun NVT"; cd ${BASE_DIR}; continue; fi
        echo "NVT уравновешивание завершено."
    else
        echo "NVT уравновешивание уже выполнено, пропускаем."
        NVT_TPR="${SYS_NAME}_nvt.tpr"
        NVT_CPT="${SYS_NAME}_nvt.cpt"
    fi

    # --- Шаг 6.3: NPT Уравновешивание ---
    echo "--> Шаг 6.3: NPT Уравновешивание для ${SYS_NAME}..."
    NPT_TPR="${SYS_NAME}_npt.tpr"
    NPT_GRO="${SYS_NAME}_npt.gro"
    NPT_CPT="${SYS_NAME}_npt.cpt"
     if [ ! -f "${NPT_GRO}" ]; then
        if [ ! -f "${NVT_CPT}" ]; then echo "ОШИБКА: Отсутствует ${NVT_CPT}"; cd ${BASE_DIR}; continue; fi
        # Используем ПРЯМОЙ путь к MDP. Grompp сам прочитает 'define = -DPOSRES' из файла MDP_NPT.
        gmx_mpi grompp -f "${MDP_NPT}" -c "${NVT_GRO}" -r "${NVT_GRO}" -t "${NVT_CPT}" -p "${TOP_FILE}" -n index.ndx -o "${NPT_TPR}" -maxwarn 5
        if [ $? -ne 0 ]; then echo "ОШИБКА: grompp NPT"; cd ${BASE_DIR}; continue; fi
        gmx_mpi mdrun -deffnm "${SYS_NAME}_npt" -v -nb gpu -pme gpu -update gpu -bonded gpu
        if [ $? -ne 0 ]; then echo "ОШИБКА: mdrun NPT"; cd ${BASE_DIR}; continue; fi
        echo "NPT уравновешивание завершено."
    else
        echo "NPT уравновешивание уже выполнено, пропускаем."
        NPT_TPR="${SYS_NAME}_npt.tpr"
        NPT_CPT="${SYS_NAME}_npt.cpt"
    fi

    # --- Шаг 6.4: Продуктивная МД ---
    echo "--> Шаг 6.4: Продуктивная МД для ${SYS_NAME}..."
    PROD_TPR="${SYS_NAME}_prod.tpr"
    PROD_XTC="${SYS_NAME}_prod.xtc"
    if [ ! -f "${PROD_XTC}" ]; then
        if [ ! -f "${NPT_CPT}" ]; then echo "ОШИБКА: Отсутствует ${NPT_CPT}"; cd ${BASE_DIR}; continue; fi
        # Используем ПРЯМОЙ путь к MDP. Grompp прочитает MDP_PROD (где define закомментирован).
        gmx_mpi grompp -f "${MDP_PROD}" -c "${NPT_GRO}" -t "${NPT_CPT}" -p "${TOP_FILE}" -n index.ndx -o "${PROD_TPR}" -maxwarn 5
        if [ $? -ne 0 ]; then echo "ОШИБКА: grompp Production"; cd ${BASE_DIR}; continue; fi
        echo "Запуск продуктивной МД... Это займет много времени!"
        gmx_mpi mdrun -deffnm "${SYS_NAME}_prod" -v -nb gpu -pme gpu -update gpu -bonded gpu
        if [ $? -ne 0 ]; then echo "ОШИБКА: mdrun Production"; cd ${BASE_DIR}; continue; fi
        echo "Продуктивная МД завершена."
    else
        echo "Продуктивная МД уже выполнена, пропускаем."
    fi

    echo "--- Обработка системы ${SYS_NAME} завершена ---"
    cd ${BASE_DIR} || exit 1

done # Конец итерации по папкам симуляций

echo ""
echo "============================================="
echo "    КОНВЕЙЕР СИМУЛЯЦИЙ ЗАВЕРШЕН"
echo "============================================="
