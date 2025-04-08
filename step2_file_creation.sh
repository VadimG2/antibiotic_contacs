#!/bin/bash

# --- Входные данные ---
PEPTIDE_PDB="pept1.pdb"
PEPTIDE_GRO="pept1.gro"
PEPTIDE_TOP="pept1.top"
PEPTIDE_POSRES_OLD="posre.itp"          # Имя файла, создаваемого pdb2gmx
PEPTIDE_POSRES_NEW="posre_protein.itp" # Новое имя файла ограничений

# --- Параметры pdb2gmx ---
FORCEFIELD="amber14sb"
WATERMODEL="tip3p"

echo "--- ШАГ 2: Подготовка топологии пептида ---"

# 2.1: Генерация топологии GROMACS
echo "--> 2.1: Запуск pdb2gmx для ${PEPTIDE_PDB}..."
gmx_mpi pdb2gmx -f ${PEPTIDE_PDB} -o ${PEPTIDE_GRO} -p ${PEPTIDE_TOP} -ignh -ff ${FORCEFIELD} -water ${WATERMODEL}
# Проверка успешности выполнения pdb2gmx
if [ $? -ne 0 ] || [ ! -f "${PEPTIDE_TOP}" ] || [ ! -f "${PEPTIDE_GRO}" ] || [ ! -f "${PEPTIDE_POSRES_OLD}" ]; then
    echo "ОШИБКА: gmx pdb2gmx не завершился успешно или не создал необходимые файлы."
    exit 1
fi
echo "--> Топология (${PEPTIDE_TOP}), координаты (${PEPTIDE_GRO}) и ограничения (${PEPTIDE_POSRES_OLD}) успешно сгенерированы."

# 2.2: Переименование файла ограничений
echo "--> 2.2: Переименование ${PEPTIDE_POSRES_OLD} -> ${PEPTIDE_POSRES_NEW}"
if [ -f "${PEPTIDE_POSRES_OLD}" ]; then
    mv ${PEPTIDE_POSRES_OLD} ${PEPTIDE_POSRES_NEW}
    if [ $? -ne 0 ]; then
        echo "ОШИБКА: Не удалось переименовать файл ограничений."
        exit 1
    fi
else
    echo "ВНИМАНИЕ: Файл ${PEPTIDE_POSRES_OLD} не найден. Возможно, pdb2gmx не создал его."
    # Можно добавить exit 1, если этот файл обязателен
fi

# 2.3: Замена include в файле топологии
echo "--> 2.3: Обновление #include в файле ${PEPTIDE_TOP}..."
# Используем sed для поиска и замены строки.
# s/шаблон/замена/ - команда замены
# \#include \"posre\.itp\" - искомый шаблон. Кавычки и точка экранированы '\'.
# \#include \"posre_protein\.itp\" - строка для замены.
# g - флаг глобальной замены (хотя здесь ожидается одно вхождение)
sed -i.bak 's/#include "posre\.itp"/#include "posre_protein.itp"/g' ${PEPTIDE_TOP}

# Проверка, произошла ли замена (сравнением с бэкапом)
if cmp -s "${PEPTIDE_TOP}" "${PEPTIDE_TOP}.bak"; then
    echo "ВНИМАНИЕ: Строка '#include \"posre.itp\"' не найдена или не заменена в ${PEPTIDE_TOP}."
    # Здесь можно решить, критична ли ошибка. Если PosRes нужен, то exit 1.
    # rm "${PEPTIDE_TOP}.bak" # Удалить бэкап, если замена не нужна/не прошла
else
    echo "--> #include успешно обновлен на '${PEPTIDE_POSRES_NEW}'."
    rm "${PEPTIDE_TOP}.bak" # Удаляем бэкап, так как замена прошла успешно
fi

echo "--- ШАГ 2 ЗАВЕРШЕН: Топология пептида подготовлена ---"
