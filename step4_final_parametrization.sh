#!/bin/bash

# --- Директории ---
INPUT_LIGAND_DIR="input/ligand"
OUTPUT_BASE_DIR="output"

# --- Параметры ---
LIGAND_CHARGE=0      # Чистый заряд лиганда (можно сделать параметром скрипта!)
ATOM_TYPE="gaff2"
AMBER_COMPATIBILITY="-k amber" # Флаг для acpype

# --- Рабочая директория для временных файлов AmberTools ---
# Создаем временную папку для файлов каждого лиганда, чтобы избежать конфликтов имен
WORK_DIR_BASE="ligand_work"

echo "--- НАЧАЛО: Параметризация всех лигандов (Обходной путь) ---"

# Проверяем наличие входной директории
if [ ! -d "${INPUT_LIGAND_DIR}" ]; then
    echo "ОШИБКА: Входная директория ${INPUT_LIGAND_DIR} не найдена!"
    exit 1
fi

# Создаем базовую выходную директорию, если ее нет
mkdir -p "${OUTPUT_BASE_DIR}"
if [ $? -ne 0 ]; then echo "ОШИБКА: Не удалось создать директорию ${OUTPUT_BASE_DIR}"; exit 1; fi

# --- Итерация по всем PDB файлам в input/ligand ---
shopt -s extglob nocaseglob
for LIGAND_PDB_PATH in "${INPUT_LIGAND_DIR}"/*.@(pdb); do
    shopt -u nocaseglob

    if [ -f "${LIGAND_PDB_PATH}" ]; then
        LIGAND_BASENAME=$(basename "${LIGAND_PDB_PATH}" .pdb)
        echo ""
        echo "--- Обработка лиганда: ${LIGAND_BASENAME} ---"

        # Создаем уникальную выходную директорию для лиганда
        LIGAND_OUTPUT_DIR="${OUTPUT_BASE_DIR}/${LIGAND_BASENAME}"
        mkdir -p "${LIGAND_OUTPUT_DIR}"
        if [ $? -ne 0 ]; then echo "ОШИБКА: Не удалось создать директорию ${LIGAND_OUTPUT_DIR}"; continue; fi

        # Создаем и переходим во временную рабочую директорию
        CURRENT_WORK_DIR="${WORK_DIR_BASE}_${LIGAND_BASENAME}"
        rm -rf "${CURRENT_WORK_DIR}" # Удаляем старую, если есть
        mkdir -p "${CURRENT_WORK_DIR}"
        if [ $? -ne 0 ]; then echo "ОШИБКА: Не удалось создать рабочую директорию ${CURRENT_WORK_DIR}"; continue; fi
        cd "${CURRENT_WORK_DIR}"
        echo "--> Рабочая директория: $(pwd)"

        # Копируем исходный PDB во временную директорию
        cp "../${LIGAND_PDB_PATH}" ./ligand.pdb # Используем стандартное имя внутри папки
        LIGAND_PDB="ligand.pdb" # Переопределяем для удобства

        # --- Определяем имена файлов внутри рабочей директории ---
        LIGAND_MOL2="ligand.mol2"
        LIGAND_GAFF2_MOL2="ligand_gaff2.mol2"
        LIGAND_FRCMOD="ligand.frcmod"
        LIGAND_PRMTOP="ligand.prmtop"
        LIGAND_INPCRD="ligand.inpcrd"
        TLEAP_SCRIPT="tleap.in"
        ACPYPE_BASENAME="${LIGAND_BASENAME}_bcc_manual" # Уникальное базовое имя
        ACPYPE_OUTDIR="${ACPYPE_BASENAME}.amb2gmx"

        # --- Шаги параметризации (как в предыдущем скрипте, но внутри CURRENT_WORK_DIR) ---

        # 4.1: Конвертация PDB в MOL2
        echo "--> 4.1: Конвертация PDB -> MOL2..."
        obabel ${LIGAND_PDB} -O ${LIGAND_MOL2}
        if [ $? -ne 0 ] || [ ! -f "${LIGAND_MOL2}" ]; then echo "ОШИБКА: Open Babel."; cd ..; rm -rf "${CURRENT_WORK_DIR}"; continue; fi

        # 4.2: Запуск antechamber
        echo "--> 4.2: Запуск antechamber (типы ${ATOM_TYPE}, заряды bcc)..."
        antechamber -i ${LIGAND_MOL2} -fi mol2 -o ${LIGAND_GAFF2_MOL2} -fo mol2 -c bcc -nc ${LIGAND_CHARGE} -at ${ATOM_TYPE} -s 2
        if [ ! -f "${LIGAND_GAFF2_MOL2}" ]; then echo "ОШИБКА: antechamber."; cd ..; rm -rf "${CURRENT_WORK_DIR}"; continue; fi

        # 4.3: Запуск parmchk2
        echo "--> 4.3: Запуск parmchk2..."
        parmchk2 -i ${LIGAND_GAFF2_MOL2} -f mol2 -o ${LIGAND_FRCMOD}
        if [ $? -ne 0 ] || [ ! -f "${LIGAND_FRCMOD}" ]; then echo "ОШИБКА: parmchk2."; cd ..; rm -rf "${CURRENT_WORK_DIR}"; continue; fi

        # 4.4: Создание скрипта tleap.in
        echo "--> 4.4: Создание ${TLEAP_SCRIPT}..."
        cat << EOF > ${TLEAP_SCRIPT}
source leaprc.${ATOM_TYPE}
loadamberparams ${LIGAND_FRCMOD}
lig = loadmol2 ${LIGAND_GAFF2_MOL2}
saveamberparm lig ${LIGAND_PRMTOP} ${LIGAND_INPCRD}
quit
EOF

        # 4.5: Запуск tleap
        echo "--> 4.5: Запуск tleap..."
        tleap -f ${TLEAP_SCRIPT} > leap.log 2>&1
        if [ $? -ne 0 ] || [ ! -f "${LIGAND_PRMTOP}" ] || [ ! -f "${LIGAND_INPCRD}" ]; then echo "ОШИБКА: tleap. См. leap.log."; cd ..; rm -rf "${CURRENT_WORK_DIR}"; continue; fi
        echo "--> Файлы AMBER созданы."

        # 4.6: Конвертация AMBER -> GROMACS (через acpype)
        echo "--> 4.6: Запуск acpype для конвертации..."
        acpype -p ${LIGAND_PRMTOP} -x ${LIGAND_INPCRD} -b ${ACPYPE_BASENAME} ${AMBER_COMPATIBILITY}
        if [ $? -ne 0 ] || [ ! -d "${ACPYPE_OUTDIR}" ]; then echo "ОШИБКА: acpype."; cd ..; rm -rf "${CURRENT_WORK_DIR}"; continue; fi
        echo "--> Директория ${ACPYPE_OUTDIR} создана."

        # 4.7: Извлечение и подготовка финальных файлов GROMACS в выходную папку
        echo "--> 4.7: Извлечение и подготовка финальных файлов GROMACS в ${LIGAND_OUTPUT_DIR}"

        # Определяем пути к файлам внутри папки acpype
        ACPYPE_GRO_SOURCE="${ACPYPE_OUTDIR}/${ACPYPE_BASENAME}_GMX.gro"
        ACPYPE_TOP_SOURCE="${ACPYPE_OUTDIR}/${ACPYPE_BASENAME}_GMX.top"
        ACPYPE_POSRES_SOURCE="${ACPYPE_OUTDIR}/posre_${ACPYPE_BASENAME}.itp"

        # Определяем финальные пути в папке output/ligand_name
        FINAL_LIGAND_ITP_DEST="../${LIGAND_OUTPUT_DIR}/${ACPYPE_BASENAME}.itp"
        FINAL_LIGAND_GRO_DEST="../${LIGAND_OUTPUT_DIR}/${ACPYPE_BASENAME}_GMX.gro"
        FINAL_LIGAND_POSRES_DEST="../${LIGAND_OUTPUT_DIR}/posre_${ACPYPE_BASENAME}.itp"
        FINAL_ATOMTYPES_DEST="../${LIGAND_OUTPUT_DIR}/atomtypes_${ACPYPE_BASENAME}.itp" # Файл для atomtypes

        if [ ! -f "${ACPYPE_GRO_SOURCE}" ] || [ ! -f "${ACPYPE_TOP_SOURCE}" ] || [ ! -f "${ACPYPE_POSRES_SOURCE}" ]; then
             echo "ОШИБКА: Не найдены .gro, .top или posre*.itp в ${ACPYPE_OUTDIR}."
             cd ..; rm -rf "${CURRENT_WORK_DIR}"; continue
        fi

        # Копируем .gro
        cp "${ACPYPE_GRO_SOURCE}" "${FINAL_LIGAND_GRO_DEST}" || { echo "ОШИБКА копирования gro"; cd ..; rm -rf "${CURRENT_WORK_DIR}"; continue; }
        # Копируем PosRes
        cp "${ACPYPE_POSRES_SOURCE}" "${FINAL_LIGAND_POSRES_DEST}" || { echo "ОШИБКА копирования PosRes"; cd ..; rm -rf "${CURRENT_WORK_DIR}"; continue; }

        # Извлекаем [ moleculetype ]
        awk '/^\[ moleculetype \]/,/^\[ system \]/' "${ACPYPE_TOP_SOURCE}" | sed '$d' > "${FINAL_LIGAND_ITP_DEST}"
        if [ $? -ne 0 ] || [ ! -s "${FINAL_LIGAND_ITP_DEST}" ]; then echo "ОШИБКА извлечения moleculetype"; cd ..; rm -rf "${CURRENT_WORK_DIR}"; continue; fi

        # Добавляем #ifdef POSRES в .itp лиганда
        printf "\n; Include Position restraint file for ligand\n#ifdef POSRES\n#include \"%s\"\n#endif\n" "$(basename ${FINAL_LIGAND_POSRES_DEST})" >> "${FINAL_LIGAND_ITP_DEST}"
        if [ $? -ne 0 ]; then echo "ОШИБКА добавления #ifdef POSRES"; cd ..; rm -rf "${CURRENT_WORK_DIR}"; continue; fi

        # Извлекаем [ atomtypes ]
        awk '/^\[ atomtypes \]/ { printing=1; next } /^\[.*\]/ && printing { printing=0 } printing' "${ACPYPE_TOP_SOURCE}" > "${FINAL_ATOMTYPES_DEST}"
        if [ -s "${FINAL_ATOMTYPES_DEST}" ]; then
            sed -i '1i\[ atomtypes ]' "${FINAL_ATOMTYPES_DEST}"
            echo "--> Типы атомов сохранены в ${FINAL_ATOMTYPES_DEST}"
        else
            echo "ВНИМАНИЕ: Секция [ atomtypes ] не найдена или пуста для ${LIGAND_BASENAME}."
            rm -f "${FINAL_ATOMTYPES_DEST}"
        fi

        # Возвращаемся из временной директории и удаляем ее
        cd ..
        echo "--> Очистка временной директории ${CURRENT_WORK_DIR}..."
        rm -rf "${CURRENT_WORK_DIR}"

        echo "--- Обработка лиганда ${LIGAND_BASENAME} УСПЕШНО ЗАВЕРШЕНА ---"
        echo "Результаты в директории: ${LIGAND_OUTPUT_DIR}"

    fi
    shopt -s nocaseglob
done

shopt -u extglob nocaseglob

echo ""
echo "--- ПАРАМЕТРИЗАЦИЯ ВСЕХ ЛИГАНДОВ завершена ---"
