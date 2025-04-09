#!/bin/bash

# --- Пути к директориям ---
INPUT_PROTEIN_DIR="input/protein"
OUTPUT_DIR="output"

# --- Параметры pdb2gmx ---
FORCEFIELD="amber14sb"
WATERMODEL="tip3p"

# Создаем основную выходную директорию, если ее нет
mkdir -p "${OUTPUT_DIR}"

echo "=== НАЧАЛО ОБРАБОТКИ ВСЕХ КОМБИНАЦИЙ ПЕПТИДОВ И ЛИГАНДОВ ==="

# Итерация по всем файлам .pdb в input/protein
for protein_pdb in "${INPUT_PROTEIN_DIR}"/*.pdb; do
    # Проверяем, есть ли файлы в директории
    if [ ! -e "${protein_pdb}" ]; then
        echo "ОШИБКА: Не найдено .pdb файлов в ${INPUT_PROTEIN_DIR}"
        exit 1
    fi

    # Извлекаем имя файла без пути и расширения
    protein_base=$(basename "${protein_pdb}" .pdb)
    protein_output_dir="${OUTPUT_DIR}/${protein_base}"
    mkdir -p "${protein_output_dir}"

    # --- Входные и выходные файлы для пептида ---
    PEPTIDE_PDB="${protein_pdb}"
    PEPTIDE_GRO="${protein_output_dir}/${protein_base}.gro"
    PEPTIDE_TOP="${protein_output_dir}/${protein_base}.top"
    PEPTIDE_POSRES_OLD="${protein_output_dir}/posre.itp"          # Куда мы переместим posre.itp
    PEPTIDE_POSRES_NEW="${protein_output_dir}/posre_protein.itp"  # Новое имя файла ограничений

    echo "--- ШАГ 2: Подготовка топологии пептида ${protein_base} ---"

    # 2.1: Генерация топологии GROMACS
    echo "--> 2.1: Запуск pdb2gmx для ${protein_base}..."
    gmx_mpi pdb2gmx -f "${PEPTIDE_PDB}" -o "${PEPTIDE_GRO}" -p "${PEPTIDE_TOP}" -ignh -ff "${FORCEFIELD}" -water "${WATERMODEL}"
    
    # Перемещаем posre.itp из текущей директории, если он там создался
    if [ -f "posre.itp" ]; then
        mv "posre.itp" "${PEPTIDE_POSRES_OLD}"
        if [ $? -ne 0 ]; then
            echo "ОШИБКА: Не удалось переместить posre.itp в ${PEPTIDE_POSRES_OLD}"
            continue
        fi
    fi

    # Проверка успешности выполнения pdb2gmx
    if [ $? -ne 0 ] || [ ! -f "${PEPTIDE_TOP}" ] || [ ! -f "${PEPTIDE_GRO}" ] || [ ! -f "${PEPTIDE_POSRES_OLD}" ]; then
        echo "ОШИБКА: gmx pdb2gmx не завершился успешно или не создал необходимые файлы для ${protein_base}."
        continue
    fi
    echo "--> Топология (${PEPTIDE_TOP}), координаты (${PEPTIDE_GRO}) и ограничения (${PEPTIDE_POSRES_OLD}) успешно сгенерированы."

    # 2.2: Переименование файла ограничений
    echo "--> 2.2: Переименование ${PEPTIDE_POSRES_OLD} -> ${PEPTIDE_POSRES_NEW}"
    if [ -f "${PEPTIDE_POSRES_OLD}" ]; then
        mv "${PEPTIDE_POSRES_OLD}" "${PEPTIDE_POSRES_NEW}"
        if [ $? -ne 0 ]; then
            echo "ОШИБКА: Не удалось переименовать файл ограничений для ${protein_base}."
            continue
        fi
    else
        echo "ВНИМАНИЕ: Файл ${PEPTIDE_POSRES_OLD} не найден для ${protein_base}. Возможно, pdb2gmx не создал его."
    fi

    # 2.3: Замена include в файле топологии
    echo "--> 2.3: Обновление #include в файле ${PEPTIDE_TOP}..."
    sed -i.bak 's/#include "posre\.itp"/#include "posre_protein.itp"/g' "${PEPTIDE_TOP}"
    if cmp -s "${PEPTIDE_TOP}" "${PEPTIDE_TOP}.bak"; then
        echo "ВНИМАНИЕ: Строка '#include \"posre.itp\"' не найдена или не заменена в ${PEPTIDE_TOP}."
    else
        echo "--> #include успешно обновлен на 'posre_protein.itp' для ${protein_base}."
        rm "${PEPTIDE_TOP}.bak"
    fi

    echo "--- ШАГ 2 ЗАВЕРШЕН: Топология пептида ${protein_base} подготовлена ---"
done

echo "=== ОБРАБОТКА ШАГА 2 ДЛЯ ВСЕХ ПЕПТИДОВ ЗАВЕРШЕНА ==="
