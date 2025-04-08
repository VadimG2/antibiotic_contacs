#!/bin/bash

# --- Директории ---
OUTPUT_BASE_DIR="output" # Директория, где лежат папки с результатами Шага 2

# --- Параметры отжига ---
ANNEAL_TIME_PS=500
ANNEAL_TEMP_LOW=5
ANNEAL_TEMP_HIGH=400
ANNEAL_TIME_RAMP_PS=250

echo "--- НАЧАЛО: Отжиг для всех подготовленных пептидов ---"

# Проверяем наличие базовой выходной директории
if [ ! -d "${OUTPUT_BASE_DIR}" ]; then
    echo "ОШИБКА: Базовая директория ${OUTPUT_BASE_DIR} не найдена! Запустите сначала prepare_all_peptides.sh."
    exit 1
fi

# --- Итерация по всем поддиректориям в output ---
for PEPTIDE_DIR in "${OUTPUT_BASE_DIR}"/*/; do
    # Убираем слэш в конце имени директории
    PEPTIDE_DIR=${PEPTIDE_DIR%/}
    # Получаем базовое имя пептида из имени директории
    PEPTIDE_BASENAME=$(basename "${PEPTIDE_DIR}")

    echo ""
    echo "--- Обработка пептида: ${PEPTIDE_BASENAME} (в директории ${PEPTIDE_DIR}) ---"

    # Определяем пути к входным файлам
    PEPTIDE_GRO_INPUT="${PEPTIDE_DIR}/${PEPTIDE_BASENAME}.gro"
    PEPTIDE_TOP_INPUT="${PEPTIDE_DIR}/${PEPTIDE_BASENAME}.top"

    # Проверяем наличие входных файлов
    if [ ! -f "${PEPTIDE_GRO_INPUT}" ] || [ ! -f "${PEPTIDE_TOP_INPUT}" ]; then
        echo "ВНИМАНИЕ: Отсутствуют файлы ${PEPTIDE_BASENAME}.gro или ${PEPTIDE_BASENAME}.top в ${PEPTIDE_DIR}. Пропускаем."
        continue
    fi

    # Определяем имена промежуточных и выходных файлов внутри директории пептида
    PEPTIDE_BOX_GRO="${PEPTIDE_DIR}/${PEPTIDE_BASENAME}_box.gro"
    ANNEAL_MDP="${PEPTIDE_DIR}/anneal_${PEPTIDE_BASENAME}.mdp" # Уникальное имя MDP
    ANNEAL_TPR="${PEPTIDE_DIR}/${PEPTIDE_BASENAME}_anneal.tpr"
    ANNEAL_DEFFNM="${PEPTIDE_DIR}/${PEPTIDE_BASENAME}_anneal" # Базовое имя для mdrun
    ANNEAL_XTC="${ANNEAL_DEFFNM}.xtc" # Имя файла траектории
    OPTIMIZED_GRO="${PEPTIDE_DIR}/${PEPTIDE_BASENAME}_optimized.gro" # Финальный выход

    # --- Шаг 3.1: Создание бокса ---
    echo "--> 3.1: Создание бокса ${PEPTIDE_BOX_GRO}"
    gmx_mpi editconf -f "${PEPTIDE_GRO_INPUT}" -o "${PEPTIDE_BOX_GRO}" -c -d 1.0 -bt cubic
    if [ $? -ne 0 ] || [ ! -f "${PEPTIDE_BOX_GRO}" ]; then
        echo "ОШИБКА на этапе editconf для ${PEPTIDE_BASENAME}."
        continue
    fi

    # --- Шаг 3.2: Создание MDP файла ---
    echo "--> 3.2: Создание ${ANNEAL_MDP}"
    cat << EOF > "${ANNEAL_MDP}"
; --- Параметры отжига для ${PEPTIDE_BASENAME} (NPT) ---
title           = ${PEPTIDE_BASENAME} Annealing
integrator      = md
nsteps          = $(($ANNEAL_TIME_PS * 1000))
dt              = 0.001
nstxout-compressed = 5000
nstlog          = 5000
nstenergy       = 5000
constraints     = h-bonds
constraint_algorithm = LINCS
cutoff-scheme   = Verlet
nstlist         = 10
rcoulomb        = 1.0
rvdw            = 1.0
DispCorr        = EnerPres
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.16
tcoupl          = V-rescale
tc-grps         = Protein       ; Ожидаем группу Protein
tau_t           = 0.1
ref_t           = $ANNEAL_TEMP_HIGH
pcoupl          = Berendsen
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = 1.0
compressibility = 4.5e-5
annealing       = single
annealing-npoints = 3
annealing-time  = 0 $ANNEAL_TIME_RAMP_PS $ANNEAL_TIME_PS
annealing-temp  = $ANNEAL_TEMP_LOW $ANNEAL_TEMP_HIGH $ANNEAL_TEMP_HIGH
gen_vel         = yes
gen_temp        = $ANNEAL_TEMP_LOW
gen_seed        = -1
pbc             = xyz
EOF

    # --- Шаг 3.3: Запуск grompp ---
    echo "--> 3.3: Запуск grompp для ${PEPTIDE_BASENAME}"
    gmx_mpi grompp -f "${ANNEAL_MDP}" -c "${PEPTIDE_BOX_GRO}" -p "${PEPTIDE_TOP_INPUT}" -o "${ANNEAL_TPR}" -maxwarn 2
    if [ $? -ne 0 ] || [ ! -f "${ANNEAL_TPR}" ]; then
        echo "ОШИБКА на этапе grompp для ${PEPTIDE_BASENAME}."
        continue
    fi

    # --- Шаг 3.4: Запуск mdrun ---
    echo "--> 3.4: Запуск mdrun для отжига ${PEPTIDE_BASENAME} (${ANNEAL_TIME_PS} ps)..."
    # Запускаем mdrun, используя -deffnm с полным путем
    gmx_mpi mdrun -deffnm "${ANNEAL_DEFFNM}" -v
    if [ $? -ne 0 ] || [ ! -f "${ANNEAL_XTC}" ]; then
        echo "ОШИБКА на этапе mdrun для ${PEPTIDE_BASENAME}."
        continue
    fi

    # --- Шаг 3.5: Извлечение кадра ---
    echo "--> 3.5: Извлечение финального кадра для ${PEPTIDE_BASENAME}"
    # Получаем номер группы Protein (или запасной вариант)
    NDX_OUTPUT=$(echo "q" | gmx_mpi make_ndx -f "${PEPTIDE_BOX_GRO}" -o junk_index.ndx 2>&1)
    rm -f junk_index.ndx
    PROTEIN_GROUP_NUM=$(echo "${NDX_OUTPUT}" | grep "Protein" | head -n 1 | awk '{print $1}')
    if [[ -z "${PROTEIN_GROUP_NUM}" || ! "${PROTEIN_GROUP_NUM}" =~ ^[0-9]+$ ]]; then
        echo "ВНИМАНИЕ: Не удалось найти группу 'Protein' для ${PEPTIDE_BASENAME}. Используем группу 1."
        PROTEIN_GROUP_NUM=1
    fi
    echo "Для trjconv ${PEPTIDE_BASENAME} будет выбрана группа: ${PROTEIN_GROUP_NUM}"

    # Запускаем trjconv
    echo "${PROTEIN_GROUP_NUM}" | gmx_mpi trjconv -s "${ANNEAL_TPR}" -f "${ANNEAL_XTC}" -o "${OPTIMIZED_GRO}" -dump $ANNEAL_TIME_PS
    if [ $? -ne 0 ] || [ ! -f "${OPTIMIZED_GRO}" ]; then
        echo "ОШИБКА на этапе trjconv для ${PEPTIDE_BASENAME} или файл ${OPTIMIZED_GRO} не создан."
        continue
    fi
    echo "--> Оптимизированные координаты сохранены в: ${OPTIMIZED_GRO}"

    # --- Шаг 3.6: Очистка ---
    echo "--> 3.6: Очистка промежуточных файлов для ${PEPTIDE_BASENAME}..."
    rm -f "${PEPTIDE_BOX_GRO}" "${ANNEAL_MDP}" "${ANNEAL_DEFFNM}.log" "${ANNEAL_DEFFNM}.xtc" "${ANNEAL_DEFFNM}.edr" "${ANNEAL_DEFFNM}.tpr" "${ANNEAL_DEFFNM}.gro" "${ANNEAL_DEFFNM}.cpt" "${PEPTIDE_DIR}/mdout.mdp" "${PEPTIDE_DIR}/\#*#"
    echo "Очистка для ${PEPTIDE_BASENAME} завершена."

    echo "--- Обработка пептида ${PEPTIDE_BASENAME} УСПЕШНО ЗАВЕРШЕНА ---"

done

echo ""
echo "--- ОТЖИГ ВСЕХ ПЕПТИДОВ завершен ---"
