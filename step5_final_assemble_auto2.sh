#!/bin/bash

# --- Базовые директории ---
BASE_DIR=$(pwd) # Директория, где лежит скрипт
INPUT_DIR_UNUSED="${BASE_DIR}/input"
OUTPUT_DIR="${BASE_DIR}/output"
PROTEIN_PREP_DIR_BASE="${OUTPUT_DIR}"
LIGAND_PREP_DIR_BASE="${OUTPUT_DIR}"

# --- Параметры ---
FORCEFIELD="amber14sb"
WATERMODEL="tip3p"
BOX_DISTANCE=2.0
ION_CONC=0.15
LIGAND_CHARGE=0
GAFF_VERSION="gaff2"

mkdir -p "${OUTPUT_DIR}"

echo "============================================="
echo "    ЗАПУСК КОНВЕЙЕРА МД СИМУЛЯЦИЙ (v13 - Рабочая логика)"
echo "============================================="

# --- Итерация по папкам ГОТОВЫХ пептидов ---
for PROTEIN_PREP_DIR in ${PROTEIN_PREP_DIR_BASE}/pept*; do
    if [ ! -d "${PROTEIN_PREP_DIR}" ]; then continue; fi

    PROTEIN_NAME=$(basename "${PROTEIN_PREP_DIR}")
    PROTEIN_OPT_GRO="${PROTEIN_PREP_DIR}/${PROTEIN_NAME}_optimized.gro"
    PROTEIN_TOP_FILE="${PROTEIN_PREP_DIR}/${PROTEIN_NAME}.top"
    POSRE_PROTEIN_FILE="${PROTEIN_PREP_DIR}/posre_protein.itp" # Ожидаемое имя PosRe пептида

    # Проверка файлов пептида
    if [ ! -f "${PROTEIN_OPT_GRO}" ] || [ ! -f "${PROTEIN_TOP_FILE}" ] || [ ! -f "${POSRE_PROTEIN_FILE}" ]; then
        echo "ПРЕДУПРЕЖДЕНИЕ: Пропускаем пептид ${PROTEIN_NAME}, отсутствуют подготовленные файлы."
        continue
    fi
    # Определяем имя молекулы пептида
    PROTEIN_MOL_NAME_FROM_TOP=$(awk '/^\[ *moleculetype *\]/{f=1; next} f && !/^;/ && NF>=2{print $1; exit}' "${PROTEIN_TOP_FILE}")
    if [ -z "${PROTEIN_MOL_NAME_FROM_TOP}" ]; then echo "ОШИБКА: Не удалось определить имя пептида из ${PROTEIN_TOP_FILE}"; continue; fi
    echo ""
    echo "--- Найден пептид: ${PROTEIN_NAME} (Имя молекулы: ${PROTEIN_MOL_NAME_FROM_TOP}) ---"

    # --- Итерация по папкам ГОТОВЫХ лигандов ---
    for LIGAND_PREP_DIR in ${LIGAND_PREP_DIR_BASE}/lig*; do
        if [ ! -d "${LIGAND_PREP_DIR}" ]; then continue; fi

        LIGAND_NAME=$(basename "${LIGAND_PREP_DIR}")
        LIGAND_ITP_SOURCE="${LIGAND_PREP_DIR}/${LIGAND_NAME}.itp"
        LIGAND_GRO_SOURCE="${LIGAND_PREP_DIR}/${LIGAND_NAME}_GMX.gro"
        LIGAND_ATOMTYPES_SOURCE="${LIGAND_PREP_DIR}/atomtypes_${LIGAND_NAME}.itp"
        POSRE_LIGAND_SOURCE="${LIGAND_PREP_DIR}/posre_${LIGAND_NAME}.itp"

        echo "  --- Поиск лиганда: ${LIGAND_NAME} ---"
        # Проверка файлов лиганда
        if [ ! -f "${LIGAND_ITP_SOURCE}" ] || [ ! -f "${LIGAND_GRO_SOURCE}" ] || [ ! -f "${LIGAND_ATOMTYPES_SOURCE}" ] || [ ! -f "${POSRE_LIGAND_SOURCE}" ]; then
             echo "  ПРЕДУПРЕЖДЕНИЕ: Пропускаем лиганд ${LIGAND_NAME}, отсутствуют все необходимые файлы в ${LIGAND_PREP_DIR}"
             continue
        fi
        # --- Шаг 4Б: Коррекция имени молекулы и PosRe в ITP лиганда ---
        LIGAND_MOL_NAME_TARGET="${LIGAND_NAME}"
        LIGAND_MOL_NAME_ORIG=$(awk '/^\[ *moleculetype *\]/{f=1; next} f && !/^;/ && NF>=2{print $1; exit}' "${LIGAND_ITP_SOURCE}")
        POSRE_LIGAND_BASENAME=$(basename "${POSRE_LIGAND_SOURCE}")
        # Исправляем только если имя не совпадает или #include для PosRe отсутствует/неправильный
        if [ "${LIGAND_MOL_NAME_ORIG}" != "${LIGAND_MOL_NAME_TARGET}" ] || ! grep -q "#include \+\"${POSRE_LIGAND_BASENAME}\"" "${LIGAND_ITP_SOURCE}"; then
             echo "  Коррекция файла ${LIGAND_ITP_SOURCE}..."
             cp "${LIGAND_ITP_SOURCE}" "${LIGAND_ITP_SOURCE}.bak"
             # Заменяем старое имя молекулы на новое (более точная замена)
             sed -i "/^\[ *moleculetype *\]/{n; s/^ *${LIGAND_MOL_NAME_ORIG}[[:space:]]\+/ ${LIGAND_MOL_NAME_TARGET} /}" "${LIGAND_ITP_SOURCE}"
             # Удаляем старый блок PosRes, если он есть
             sed -i '/; Include Position restraint file/,/#endif/d' "${LIGAND_ITP_SOURCE}"
             # Добавляем новый правильный блок PosRes
             echo -e "\n; Include Position restraint file\n#ifdef POSRES\n#include \"${POSRE_LIGAND_BASENAME}\"\n#endif" >> "${LIGAND_ITP_SOURCE}"
        fi
        LIGAND_MOL_NAME_FROM_ITP="${LIGAND_MOL_NAME_TARGET}" # Используем исправленное имя
        echo "  --- Найден лиганд: ${LIGAND_NAME} (Имя молекулы: ${LIGAND_MOL_NAME_FROM_ITP}) ---"

        # --- ШАГ 5: Сборка системы ---
        SIMULATION_DIR="${OUTPUT_DIR}/${PROTEIN_NAME}+${LIGAND_NAME}_Sim"
        SIMULATION_TOP_FILE="${SIMULATION_DIR}/system.top"
        SIMULATION_GRO_FILE="${SIMULATION_DIR}/system_solv_ions.gro"

        echo "  --> Шаг 5: Сборка системы ${PROTEIN_NAME}+${LIGAND_NAME} в ${SIMULATION_DIR}"
        mkdir -p "${SIMULATION_DIR}"
        cd "${SIMULATION_DIR}" || { echo "ОШИБКА: Не могу перейти в ${SIMULATION_DIR}"; cd ${BASE_DIR}; continue; }

        if [ ! -f "${SIMULATION_GRO_FILE}" ]; then
            echo "  Выполняем сборку системы..."
            # Копируем файлы
            cp "${PROTEIN_OPT_GRO}" peptide.gro
            cp "${PROTEIN_TOP_FILE}" system.top # Копируем исходный top пептида
            cp "${LIGAND_GRO_SOURCE}" ligand.gro
            cp "${LIGAND_ITP_SOURCE}" ligand.itp # Копируем ИСПРАВЛЕННЫЙ itp лиганда
            cp "${LIGAND_ATOMTYPES_SOURCE}" gaff2_atomtypes.itp
            cp "${POSRE_PROTEIN_FILE}" posre_protein.itp
            cp "${POSRE_LIGAND_SOURCE}" posre_ligand.itp
            if [ $? -ne 0 ]; then echo "ОШИБКА: Копирование файлов"; cd ${BASE_DIR}; continue; fi

            # Объединяем и создаем бокс
            gmx_mpi editconf -f peptide.gro -o peptide_center.gro -c -box 6.0 6.0 6.0
            gmx_mpi insert-molecules -f peptide_center.gro -ci ligand.gro -nmol 1 -o complex_combined.gro -try 50 -radius 0.2 -dr 1.0 1.0 1.0
            if [ $? -ne 0 ]; then echo "ОШИБКА: insert-molecules"; cd ${BASE_DIR}; continue; fi
            gmx_mpi editconf -f complex_combined.gro -o complex_box.gro -c -d ${BOX_DISTANCE} -bt cubic
            if [ $? -ne 0 ]; then echo "ОШИБКА: editconf (финальный бокс)"; cd ${BASE_DIR}; continue; fi

            # --- ПОДГОТОВКА system.top (v13 - Модификация копии pept1.top) ---
            echo "  Подготовка system.top..."

            # 1. Вставляем #include для GAFF2 типов и лиганда ПОСЛЕ forcefield.itp
            # Используем временный файл для sed
            sed '/#include.*forcefield\.itp/a \
; Include GAFF2 atom types\
#include "gaff2_atomtypes.itp"\
; Include ligand topology\
#include "ligand.itp"
' system.top > system.top.temp && mv system.top.temp system.top

            # 2. Удаляем встроенное определение SOL (если оно есть)
            # Эта команда удалит блок, начинающийся с '[ moleculetype ]' и содержащий 'SOL' в третьем поле следующей строки
            awk 'BEGIN{p=1} /^\[ *moleculetype *\]/{moltype_line=$0; getline; if ($1 == "SOL" || $3 == "SOL") p=0; else {print moltype_line; print $0; p=1}} p' system.top > system.top.temp && mv system.top.temp system.top

            # 3. ДОБАВЛЯЕМ лиганд в секцию [ molecules ] (перед SOL, который добавит solvate)
            # Используем awk для вставки перед строкой, начинающейся с SOL (или просто в конец, если SOL нет)
            awk -v lig_name="${LIGAND_MOL_NAME_FROM_ITP}" '
            /^\[ *molecules *\]/ { print; print " " lig_name "              1"; next } # Вставляем сразу после [ molecules ]
            { print }
            ' system.top > system.top.temp && mv system.top.temp system.top
            # --- ИСПРАВЛЕНИЕ: Проверяем и исправляем #include PosRe в system.top (копии pept1.top) ---
            echo "  Проверка/исправление PosRe include в system.top..."
            # Заменяем ЛЮБОЕ имя posre файла на правильное posre_protein.itp
            sed -i 's/#include "posre.*\.itp"/#include "posre_protein.itp"/' system.top
            # Убеждаемся, что блок #ifdef POSRES существует вокруг include
            if grep -q '#include "posre_protein.itp"' system.top && ! grep -B 1 '#include "posre_protein.itp"' system.top | grep -q '#ifdef POSRES'; then
                echo "  Добавляем #ifdef POSRES вокруг include posre_protein.itp..."
                sed -i '/#include "posre_protein.itp"/i ; Include Position restraint file\n#ifdef POSRES' system.top
                sed -i '/#include "posre_protein.itp"/a #endif' system.top
            fi
            # --- Конец подготовки system.top ---

            # Сольватация
            gmx_mpi solvate -cp complex_box.gro -cs spc216.gro -o complex_solv.gro -p system.top
            if [ $? -ne 0 ]; then echo "ОШИБКА: solvate"; cd ${BASE_DIR}; continue; fi

            # Ионы
            printf "integrator = md\nnsteps = 0\n" > ions.mdp
            echo "  Запуск grompp для ионов..."
            gmx_mpi grompp -f ions.mdp -c complex_solv.gro -p system.top -o ions.tpr -maxwarn 5
            if [ $? -ne 0 ]; then echo "ОШИБКА: grompp ионов"; cd ${BASE_DIR}; continue; fi
            echo "  Запуск genion (выбираем SOL)..."
            echo "SOL" | gmx_mpi genion -s ions.tpr -o system_solv_ions.gro -p system.top -pname NA -nname CL -neutral -conc ${ION_CONC}
            if [ $? -ne 0 ]; then echo "ОШИБКА: genion"; cd ${BASE_DIR}; continue; fi

            # Очистка
            rm -f peptide.gro peptide_center.gro ligand.gro complex_combined.gro complex_box.gro complex_solv.gro ions.mdp ions.tpr \#*# *~ mdout.mdp gaff2_atomtypes.itp ligand.itp posre_*.itp

            echo "  Сборка системы ${PROTEIN_NAME}+${LIGAND_NAME} завершена."
        else
            echo "  Собранная система ${PROTEIN_NAME}+${LIGAND_NAME} уже существует, пропускаем сборку."
        fi
        cd ${BASE_DIR} || exit 1

        # --- Следующий шаг: Запуск симуляций ---
        # ... (остается как было) ...

    done # Конец цикла по лигандам

done # Конец цикла по пептидам

echo ""
echo "============================================="
echo "    КОНВЕЙЕР СБОРКИ ЗАВЕРШЕН"
echo "============================================="
