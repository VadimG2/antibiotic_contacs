#!/bin/bash

# --- Директории ---
INPUT_PROTEIN_DIR="input/protein"
OUTPUT_BASE_DIR="output"

# --- Параметры pdb2gmx ---
FORCEFIELD="amber14sb"
WATERMODEL="tip3p"

echo "--- НАЧАЛО: Подготовка топологий для всех пептидов ---"

# Проверяем наличие входной директории
if [ ! -d "${INPUT_PROTEIN_DIR}" ]; then
    echo "ОШИБКА: Входная директория ${INPUT_PROTEIN_DIR} не найдена!"
    exit 1
fi

# Создаем базовую выходную директорию, если ее нет
mkdir -p "${OUTPUT_BASE_DIR}"
if [ $? -ne 0 ]; then echo "ОШИБКА: Не удалось создать директорию ${OUTPUT_BASE_DIR}"; exit 1; fi

# --- Итерация по всем PDB файлам в input/protein ---
# Игнорируем регистр расширения (.pdb, .PDB)
shopt -s extglob nocaseglob
for PEPTIDE_PDB_PATH in "${INPUT_PROTEIN_DIR}"/*.@(pdb); do
    # Восстанавливаем чувствительность к регистру для внутренних операций
    shopt -u nocaseglob

    # Проверяем, является ли найденный элемент файлом
    if [ -f "${PEPTIDE_PDB_PATH}" ]; then
        # Получаем базовое имя файла без расширения
        PEPTIDE_BASENAME=$(basename "${PEPTIDE_PDB_PATH}" .pdb)
        echo "" # Пустая строка для разделения логов
        echo "--- Обработка пептида: ${PEPTIDE_BASENAME} ---"

        # Создаем директорию для вывода этого пептида
        PEPTIDE_OUTPUT_DIR="${OUTPUT_BASE_DIR}/${PEPTIDE_BASENAME}"
        mkdir -p "${PEPTIDE_OUTPUT_DIR}"
        if [ $? -ne 0 ]; then
            echo "ОШИБКА: Не удалось создать директорию ${PEPTIDE_OUTPUT_DIR}"
            continue # Переходим к следующему файлу
        fi

        # Определяем имена выходных файлов внутри директории пептида
        PEPTIDE_GRO="${PEPTIDE_OUTPUT_DIR}/${PEPTIDE_BASENAME}.gro"
        PEPTIDE_TOP="${PEPTIDE_OUTPUT_DIR}/${PEPTIDE_BASENAME}.top"
        PEPTIDE_POSRES_OLD="${PEPTIDE_OUTPUT_DIR}/posre.itp" # pdb2gmx создаст его здесь
        PEPTIDE_POSRES_NEW="${PEPTIDE_OUTPUT_DIR}/posre_protein.itp"

        # --- Шаг 2.1: Генерация топологии GROMACS ---
        echo "--> 2.1: Запуск pdb2gmx для ${PEPTIDE_PDB_PATH}..."
        # Запускаем pdb2gmx, указывая полные пути для вывода
        gmx_mpi pdb2gmx -f "${PEPTIDE_PDB_PATH}" -o "${PEPTIDE_GRO}" -p "${PEPTIDE_TOP}" -ignh -ff ${FORCEFIELD} -water ${WATERMODEL}
        # pdb2gmx создаст posre.itp в ТЕКУЩЕЙ директории, если не указать -po
        # Перемещаем его сразу, если он создался в текущей директории
        if [ -f "posre.itp" ]; then
            mv "posre.itp" "${PEPTIDE_POSRES_OLD}"
        fi

        # Проверка успешности и наличия файлов
        if [ $? -ne 0 ] || [ ! -f "${PEPTIDE_TOP}" ] || [ ! -f "${PEPTIDE_GRO}" ] || [ ! -f "${PEPTIDE_POSRES_OLD}" ]; then
            echo "ОШИБКА: gmx pdb2gmx не завершился успешно для ${PEPTIDE_BASENAME} или не создал необходимые файлы в ${PEPTIDE_OUTPUT_DIR}."
            continue # Переходим к следующему файлу
        fi
        echo "--> Топология, координаты и ограничения для ${PEPTIDE_BASENAME} сгенерированы."

        # --- Шаг 2.2: Переименование файла ограничений ---
        echo "--> 2.2: Переименование posre.itp -> ${PEPTIDE_POSRES_NEW}"
        if [ -f "${PEPTIDE_POSRES_OLD}" ]; then
            mv "${PEPTIDE_POSRES_OLD}" "${PEPTIDE_POSRES_NEW}"
            if [ $? -ne 0 ]; then
                echo "ОШИБКА: Не удалось переименовать файл ограничений для ${PEPTIDE_BASENAME}."
                continue
            fi
        else
            echo "ВНИМАНИЕ: Файл ${PEPTIDE_POSRES_OLD} не найден для ${PEPTIDE_BASENAME}. Возможно, pdb2gmx не создал его."
        fi

        # --- Шаг 2.3: Замена include в файле топологии ---
        echo "--> 2.3: Обновление #include в файле ${PEPTIDE_TOP}..."
        # Заменяем старое имя файла PosRes на новое, с привязкой к имени пептида
        sed -i.bak "s/#include \"posre\.itp\"/#include \"posre_${PEPTIDE_BASENAME}\.itp\"/g" "${PEPTIDE_TOP}"

        # Проверка замены
        if cmp -s "${PEPTIDE_TOP}" "${PEPTIDE_TOP}.bak"; then
            # Проверяем, может быть там уже было правильное имя (если pdb2gmx сам его подставил?)
             if grep -q "#include \"posre_${PEPTIDE_BASENAME}\.itp\"" "${PEPTIDE_TOP}"; then
                 echo "--> #include уже содержит правильное имя '${PEPTIDE_POSRES_NEW}'."
                 rm "${PEPTIDE_TOP}.bak"
             else
                 echo "ВНИМАНИЕ: Строка '#include \"posre.itp\"' не найдена или не заменена в ${PEPTIDE_TOP}."
                 # rm "${PEPTIDE_TOP}.bak"
             fi
        else
            echo "--> #include успешно обновлен на '${PEPTIDE_POSRES_NEW}'."
            rm "${PEPTIDE_TOP}.bak" # Удаляем бэкап
        fi

        echo "--- Обработка пептида ${PEPTIDE_BASENAME} ЗАВЕРШЕНА ---"
        echo "Результаты в директории: ${PEPTIDE_OUTPUT_DIR}"
        echo "ВАЖНО: Проверьте ${PEPTIDE_TOP} на предмет корректности циклизации!"

    fi
    # Включаем обратно чувствительность к регистру для следующей итерации
    shopt -s nocaseglob
done

# Выключаем extglob в конце
shopt -u extglob nocaseglob

echo ""
echo "--- ВСЕ ПЕПТИДЫ обработаны ---"
