#!/bin/bash

# --- Входные/Выходные данные ---
LIGAND_PDB="lig1.pdb"
LIGAND_MOL2="lig1.mol2"               # Промежуточный MOL2
LIGAND_GAFF2_MOL2="lig1_gaff2.mol2"     # Выход antechamber
LIGAND_FRCMOD="lig1.frcmod"           # Выход parmchk2
LIGAND_PRMTOP="lig1.prmtop"           # Выход tleap
LIGAND_INPCRD="lig1.inpcrd"           # Выход tleap
TLEAP_SCRIPT="tleap.in"
ACPYPE_BASENAME="lig1_bcc_manual"     # Базовое имя для acpype
ACPYPE_OUTDIR="${ACPYPE_BASENAME}.amb2gmx" # Директория вывода acpype

# Имена финальных файлов GROMACS в основной директории
FINAL_LIGAND_ITP="${ACPYPE_BASENAME}.itp"
FINAL_LIGAND_GRO="${ACPYPE_BASENAME}_GMX.gro"
FINAL_LIGAND_POSRES="${ACPYPE_BASENAME}.itp" # Имя файла PosRes, которое создаст acpype
FINAL_LIGAND_POSRES_TARGET="posre_${ACPYPE_BASENAME}.itp" # Новое имя для PosRes

# --- Параметры ---
LIGAND_CHARGE=0      # Чистый заряд лиганда
ATOM_TYPE="gaff2"

echo "--- ШАГ 4 (Обходной): Параметризация лиганда ${LIGAND_PDB} ---"

# 4.1: Конвертация PDB в MOL2
echo "--> 4.1: Конвертация PDB -> MOL2..."
obabel ${LIGAND_PDB} -O ${LIGAND_MOL2}
if [ $? -ne 0 ] || [ ! -f "${LIGAND_MOL2}" ]; then
    echo "ОШИБКА: Open Babel не смог конвертировать PDB в MOL2."
    exit 1
fi

# 4.2: Запуск antechamber
echo "--> 4.2: Запуск antechamber (типы ${ATOM_TYPE}, заряды bcc)..."
antechamber -i ${LIGAND_MOL2} -fi mol2 -o ${LIGAND_GAFF2_MOL2} -fo mol2 -c bcc -nc ${LIGAND_CHARGE} -at ${ATOM_TYPE} -s 2
# Проверяем, создался ли выходной mol2. Ручной запуск не всегда выдает код ошибки при проблемах с sqm.
if [ ! -f "${LIGAND_GAFF2_MOL2}" ]; then
    echo "ОШИБКА: antechamber не создал файл ${LIGAND_GAFF2_MOL2}. Проверьте вывод antechamber и файлы sqm."
    exit 1
fi

# 4.3: Запуск parmchk2
echo "--> 4.3: Запуск parmchk2..."
parmchk2 -i ${LIGAND_GAFF2_MOL2} -f mol2 -o ${LIGAND_FRCMOD}
if [ $? -ne 0 ] || [ ! -f "${LIGAND_FRCMOD}" ]; then
    echo "ОШИБКА: parmchk2 не завершился успешно или не создал ${LIGAND_FRCMOD}."
    exit 1
fi

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
tleap -f ${TLEAP_SCRIPT} > leap.log 2>&1 # Сохраняем вывод в лог
if [ $? -ne 0 ] || [ ! -f "${LIGAND_PRMTOP}" ] || [ ! -f "${LIGAND_INPCRD}" ]; then
    echo "ОШИБКА: tleap не завершился успешно или не создал ${LIGAND_PRMTOP}/${LIGAND_INPCRD}. См. leap.log."
    exit 1
fi
echo "--> Файлы AMBER (${LIGAND_PRMTOP}, ${LIGAND_INPCRD}) созданы."

# 4.6: Конвертация AMBER -> GROMACS (через acpype)
echo "--> 4.6: Запуск acpype для конвертации AMBER -> GROMACS..."
acpype -p ${LIGAND_PRMTOP} -x ${LIGAND_INPCRD} -b ${ACPYPE_BASENAME} -k amber
if [ $? -ne 0 ] || [ ! -d "${ACPYPE_OUTDIR}" ]; then
    echo "ОШИБКА: acpype не завершился успешно или не создал директорию ${ACPYPE_OUTDIR}."
    exit 1
fi
echo "--> Директория ${ACPYPE_OUTDIR} с файлами GROMACS создана."

# 4.7: Извлечение и подготовка финальных файлов GROMACS
echo "--> 4.7: Извлечение и подготовка финальных файлов GROMACS..."

# Проверяем наличие нужных файлов в директории acpype
ACPYPE_GRO="${ACPYPE_OUTDIR}/${ACPYPE_BASENAME}_GMX.gro"
ACPYPE_TOP="${ACPYPE_OUTDIR}/${ACPYPE_BASENAME}_GMX.top"
ACPYPE_POSRES="${ACPYPE_OUTDIR}/posre_${ACPYPE_BASENAME}.itp" # Правильное имя PosRes из acpype

if [ ! -f "${ACPYPE_GRO}" ] || [ ! -f "${ACPYPE_TOP}" ] || [ ! -f "${ACPYPE_POSRES}" ]; then
     echo "ОШИБКА: Не найдены необходимые файлы .gro, .top или posre*.itp в директории ${ACPYPE_OUTDIR}."
     exit 1
fi

# Копируем .gro файл лиганда
echo "--> Копирование ${ACPYPE_GRO} -> ${FINAL_LIGAND_GRO}"
cp "${ACPYPE_GRO}" "${FINAL_LIGAND_GRO}"
if [ $? -ne 0 ]; then echo "ОШИБКА копирования .gro файла"; exit 1; fi

# Копируем и переименовываем .itp файл с позиционными ограничениями
echo "--> Копирование ${ACPYPE_POSRES} -> ${FINAL_LIGAND_POSRES_TARGET}"
cp "${ACPYPE_POSRES}" "${FINAL_LIGAND_POSRES_TARGET}"
if [ $? -ne 0 ]; then echo "ОШИБКА копирования файла PosRes"; exit 1; fi

# Извлекаем секцию [ moleculetype ] из ACPYPE_TOP в FINAL_LIGAND_ITP
echo "--> Извлечение [ moleculetype ] из ${ACPYPE_TOP} в ${FINAL_LIGAND_ITP}"
# Используем awk: печатаем строки начиная с '[ moleculetype ]' и до строки перед '[ system ]'
awk '/^\[ moleculetype \]/,/^\[ system \]/' "${ACPYPE_TOP}" | sed '$d' > "${FINAL_LIGAND_ITP}"
# sed '$d' удаляет последнюю строку (которая будет '[ system ]')
if [ $? -ne 0 ] || [ ! -s "${FINAL_LIGAND_ITP}" ]; then # -s проверяет, что файл не пустой
    echo "ОШИБКА: Не удалось извлечь [ moleculetype ] или результат пуст."
    exit 1
fi

# Добавляем блок #ifdef POSRES в конец извлеченного FINAL_LIGAND_ITP
echo "--> Добавление #ifdef POSRES в ${FINAL_LIGAND_ITP}"
# \n используется для добавления новой строки перед блоком
printf "\n; Include Position restraint file for ligand\n#ifdef POSRES\n#include \"%s\"\n#endif\n" "${FINAL_LIGAND_POSRES_TARGET}" >> "${FINAL_LIGAND_ITP}"
if [ $? -ne 0 ]; then echo "ОШИБКА добавления блока #ifdef POSRES"; exit 1; fi

echo "--> 4.8: Извлечение секции [ atomtypes ] из ${ACPYPE_TOP}"
ATOMTYPES_FILE="atomtypes_${ACPYPE_BASENAME}.itp" # Имя временного файла
# Используем awk: печатаем строки начиная с '[ atomtypes ]' и до следующей секции '['
awk '/^\[ atomtypes \]/ { printing=1; next } /^\[.*\]/ && printing { printing=0 } printing' "${ACPYPE_TOP}" > "${ATOMTYPES_FILE}"

# Проверяем, что что-то извлеклось
if [ -s "${ATOMTYPES_FILE}" ]; then # -s проверяет, что файл не пустой
    # Добавляем заголовок [ atomtypes ] обратно, так как awk его пропускает
    sed -i '1i\[ atomtypes ]' "${ATOMTYPES_FILE}"
    echo "--> Секция [ atomtypes ] сохранена в ${ATOMTYPES_FILE}"
else
    echo "ВНИМАНИЕ: Секция [ atomtypes ] не найдена в ${ACPYPE_TOP} или пуста. Файл ${ATOMTYPES_FILE} не создан или пуст."
    rm -f "${ATOMTYPES_FILE}" # Удаляем пустой файл
fi

echo "--- ШАГ 4 (Обходной) ЗАВЕРШЕН: Топология лиганда подготовлена ---"
echo "Финальные файлы:"
echo "  Топология лиганда: ${FINAL_LIGAND_ITP}"
echo "  Координаты лиганда: ${FINAL_LIGAND_GRO}"
echo "  Ограничения лиганда: ${FINAL_LIGAND_POSRES_TARGET}"
if [ -f "${ATOMTYPES_FILE}" ]; then
    echo "  Типы атомов GAFF2:  ${ATOMTYPES_FILE}"
fi

# --- Очистка промежуточных файлов (можно закомментировать для отладки) ---
echo "--> Очистка промежуточных файлов..."
rm -f ${LIGAND_MOL2} ${LIGAND_GAFF2_MOL2} ${LIGAND_FRCMOD} ${LIGAND_PRMTOP} ${LIGAND_INPCRD} ${TLEAP_SCRIPT} leap.log
rm -rf ${ACPYPE_OUTDIR} # Удаляем директорию acpype
rm -f ANTECHAMBER_* ATOMTYPE.INF sqm.in sqm.out sqm.pdb # Файлы antechamber/sqm
# НЕ УДАЛЯЕМ ${ATOMTYPES_FILE}, он нам нужен на следующем шаге
echo "Очистка завершена (кроме файла atomtypes)."
