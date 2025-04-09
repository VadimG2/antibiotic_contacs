#!/bin/bash

# --- Пути к директориям ---
INPUT_LIGAND_DIR="input/ligand"
OUTPUT_DIR="output"
WORK_DIR_BASE="ligand_work" # Используем временную папку в BASE_DIR

# --- Параметры ---
# FORCEFIELD="amber14sb" # Не нужен здесь
LIGAND_CHARGE=0
ATOM_TYPE="gaff2"

echo "=== НАЧАЛО ПАРАМЕТРИЗАЦИИ ЛИГАНДОВ (v2) ==="

# Проверка инструментов
for tool in antechamber parmchk2 tleap acpype obabel awk sed; do
    if ! command -v "$tool" &> /dev/null; then echo "ОШИБКА: $tool не найден."; exit 1; fi
done

mkdir -p "${OUTPUT_DIR}"
BASE_DIR=$(pwd) # Запоминаем базовую директорию

# Итерация по всем файлам .pdb в input/ligand
for ligand_pdb in "${INPUT_LIGAND_DIR}"/*.pdb; do
    if [ ! -f "${ligand_pdb}" ]; then continue; fi # Проверяем, что это файл

    ligand_base=$(basename "${ligand_pdb}" .pdb)
    ligand_output_dir="${OUTPUT_DIR}/${ligand_base}"
    mkdir -p "${ligand_output_dir}"

    # Создаем временную рабочую директорию в BASE_DIR
    work_dir="${BASE_DIR}/${WORK_DIR_BASE}_${ligand_base}" # Полный путь к временной папке
    rm -rf "${work_dir}"
    mkdir -p "${work_dir}"
    cd "${work_dir}" || { echo "ОШИБКА: Не могу перейти в ${work_dir}"; continue; }

    echo ""
    echo "--- ШАГ 4: Параметризация лиганда ${ligand_base} в ${work_dir} ---"

    # Копируем PDB во временную директорию
    cp "${BASE_DIR}/${ligand_pdb}" "ligand.pdb" # Используем полный путь к исходному PDB

    # Шаг 4.1: PDB -> MOL2
    echo "--> 4.1: Конвертация PDB -> MOL2"
    obabel "ligand.pdb" -O "ligand.mol2" -h
    if [ $? -ne 0 ] || [ ! -s "ligand.mol2" ]; then echo "ОШИБКА: obabel"; cd "${BASE_DIR}"; rm -rf "${work_dir}"; continue; fi

    # Шаг 4.2: antechamber
    echo "--> 4.2: Запуск antechamber"
    antechamber -i "ligand.mol2" -fi mol2 -o "ligand_gaff.mol2" -fo mol2 -c bcc -nc "${LIGAND_CHARGE}" -at "${ATOM_TYPE}" -s 2
    if [ $? -ne 0 ] || [ ! -s "ligand_gaff.mol2" ]; then echo "ОШИБКА: antechamber"; cd "${BASE_DIR}"; rm -rf "${work_dir}"; continue; fi

    # Шаг 4.3: parmchk2
    echo "--> 4.3: Запуск parmchk2"
    parmchk2 -i "ligand_gaff.mol2" -f mol2 -o "ligand.frcmod" -a "${ATOM_TYPE}"
    if [ $? -ne 0 ] || [ ! -s "ligand.frcmod" ]; then echo "ОШИБКА: parmchk2"; cd "${BASE_DIR}"; rm -rf "${work_dir}"; continue; fi

    # Шаг 4.4: tleap
    echo "--> 4.4: Создание и запуск tleap"
    cat << EOF > "tleap.in"
source leaprc.${ATOM_TYPE}
loadamberparams ligand.frcmod
lig = loadmol2 ligand_gaff.mol2
saveamberparm lig ligand.prmtop ligand.inpcrd
quit
EOF
    tleap -f "tleap.in" > "tleap.log" 2>&1
    if [ $? -ne 0 ] || [ ! -s "ligand.prmtop" ] || [ ! -s "ligand.inpcrd" ]; then echo "ОШИБКА: tleap. См. tleap.log"; cd "${BASE_DIR}"; rm -rf "${work_dir}"; continue; fi

    # Шаг 4.5: acpype AMBER -> GROMACS
    echo "--> 4.5: Конвертация AMBER -> GROMACS с acpype"
    acpype -p "ligand.prmtop" -x "ligand.inpcrd" -b "${ligand_base}" -k amber
    if [ $? -ne 0 ] || [ ! -d "${ligand_base}.amb2gmx" ]; then echo "ОШИБКА: acpype"; cd "${BASE_DIR}"; rm -rf "${work_dir}"; continue; fi

    # --- ШАГ 4.6: Извлечение и сохранение нужных файлов в OUTPUT директорию ---
    echo "--> 4.6: Извлечение и сохранение файлов в ${ligand_output_dir}"
    ACPYPEDIR="${work_dir}/${ligand_base}.amb2gmx" # Путь к папке acpype
    TARGETDIR="${BASE_DIR}/${ligand_output_dir}" # Путь к финальной папке output/ligX

    # Извлекаем [ moleculetype ] в ligX.itp
    awk '/^\[ *moleculetype *\]/{f=1} f{print} /^\[ *system *\]/{f=0}' "${ACPYPEDIR}/${ligand_base}_GMX.top" > "${TARGETDIR}/${ligand_base}.itp"
    if [ ! -s "${TARGETDIR}/${ligand_base}.itp" ]; then echo "ОШИБКА: Не удалось извлечь moleculetype"; cd "${BASE_DIR}"; rm -rf "${work_dir}"; continue; fi

    # Извлекаем [ atomtypes ] в atomtypes_ligX.itp
    awk '/^\[ *atomtypes *\]/{f=1} f{print} /^\[ *moleculetype *\]/{f=0}' "${ACPYPEDIR}/${ligand_base}_GMX.top" > "${TARGETDIR}/atomtypes_${ligand_base}.itp"
     if [ ! -s "${TARGETDIR}/atomtypes_${ligand_base}.itp" ]; then echo "ОШИБКА: Не удалось извлечь atomtypes"; cd "${BASE_DIR}"; rm -rf "${work_dir}"; continue; fi

    # Копируем GRO
    cp "${ACPYPEDIR}/${ligand_base}_GMX.gro" "${TARGETDIR}/${ligand_base}_GMX.gro"
    if [ ! -f "${TARGETDIR}/${ligand_base}_GMX.gro" ]; then echo "ОШИБКА: Не удалось скопировать _GMX.gro"; cd "${BASE_DIR}"; rm -rf "${work_dir}"; continue; fi

    # Копируем PosRe
    cp "${ACPYPEDIR}/posre_${ligand_base}.itp" "${TARGETDIR}/posre_${ligand_base}.itp"
    if [ ! -f "${TARGETDIR}/posre_${ligand_base}.itp" ]; then echo "ОШИБКА: Не удалось скопировать posre_*.itp"; cd "${BASE_DIR}"; rm -rf "${work_dir}"; continue; fi

    # Шаг 4.7: Коррекция ligX.itp (имя молекулы и PosRe #include)
    echo "--> 4.7: Коррекция ${ligand_base}.itp"
    LIGAND_ITP_FINAL="${TARGETDIR}/${ligand_base}.itp"
    POSRE_LIGAND_FINAL_BASENAME="posre_${ligand_base}.itp" # Имя PosRe файла в финальной папке

    LIGAND_MOL_NAME_ORIG=$(awk '/^\[ *moleculetype *\]/{f=1; next} f && !/^;/ && NF>=2{print $1; exit}' "${LIGAND_ITP_FINAL}")
    LIGAND_MOL_NAME_TARGET="${ligand_base}" # Используем имя файла как целевое имя молекулы

    # Исправляем, если нужно
    if [ "${LIGAND_MOL_NAME_ORIG}" != "${LIGAND_MOL_NAME_TARGET}" ] || ! grep -q "#include \+\"${POSRE_LIGAND_FINAL_BASENAME}\"" "${LIGAND_ITP_FINAL}"; then
         cp "${LIGAND_ITP_FINAL}" "${LIGAND_ITP_FINAL}.bak"
         # Заменяем имя молекулы
         sed -i "/^\[ *moleculetype *\]/{n; s/^ *${LIGAND_MOL_NAME_ORIG}[[:space:]]\+/ ${LIGAND_MOL_NAME_TARGET} /}" "${LIGAND_ITP_FINAL}"
         # Удаляем старый блок PosRes
         sed -i '/; Include Position restraint file/,/#endif/d' "${LIGAND_ITP_FINAL}"
         # Добавляем новый правильный блок PosRes
         echo -e "\n; Include Position restraint file\n#ifdef POSRES\n#include \"${POSRE_LIGAND_FINAL_BASENAME}\"\n#endif" >> "${LIGAND_ITP_FINAL}"
         echo "Файл ${LIGAND_ITP_FINAL} скорректирован."
    fi

    # Шаг 4.8: Очистка временной директории
    echo "--> 4.8: Очистка временной директории ${work_dir}"
    cd "${BASE_DIR}"; rm -rf "${work_dir}"

    echo "--- ШАГ 4 ЗАВЕРШЕН: Параметризация лиганда ${ligand_base} завершена. Файлы в ${ligand_output_dir} ---"
done

echo ""
echo "=== ПАРАМЕТРИЗАЦИЯ ВСЕХ ЛИГАНДОВ ЗАВЕРШЕНА ==="
