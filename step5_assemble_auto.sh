#!/bin/bash

# --- Имена файлов пептида ---
PEPTIDE_OPT_GRO="pept1_optimized.gro" # Из Шага 3 (отжиг)
PEPTIDE_TOP_BASE="pept1.top"         # Исходный .top от Шага 2
PEPTIDE_NAME="pept1"
PEPTIDE_MOL_NAME="Protein_chain_A"    # Имя из [ moleculetype ] в PEPTIDE_TOP_BASE

# --- Имена файлов лиганда ---
LIGAND_GRO="lig1_bcc_manual_GMX.gro" # Из Шага 4
LIGAND_ITP="lig1_bcc_manual.itp"    # Из Шага 4
# Временный файл с [ atomtypes ] лиганда (ожидается от Шага 4)
LIGAND_ATOMTYPES_ITP="atomtypes_lig1_bcc_manual.itp"
# Имя молекулы лиганда из LIGAND_ITP
LIGAND_MOL_NAME="lig1_bcc_manual"

# --- Имена файлов для вывода ---
SYSTEM_NAME="${PEPTIDE_NAME}_${LIGAND_MOL_NAME}"
FINAL_TOPOLOGY="${SYSTEM_NAME}.top" # Итоговый файл топологии системы
FINAL_COORDS="${SYSTEM_NAME}_solv_ions.gro"
FINAL_INDEX="index.ndx"

echo "--- ШАГ 5 (Автоматизированный): Сборка системы и индексирование ---"

# --- 5.A: Подготовка финального файла топологии ---
echo "--> 5.A: Создание и модификация ${FINAL_TOPOLOGY} из ${PEPTIDE_TOP_BASE}"

# Проверяем наличие входных файлов
if [ ! -f "${PEPTIDE_TOP_BASE}" ] || [ ! -f "${LIGAND_ITP}" ]; then
    echo "ОШИБКА: Не найдены входные файлы топологии ${PEPTIDE_TOP_BASE} или ${LIGAND_ITP}"
    exit 1
fi
# Проверяем наличие файла с типами атомов, если он ожидается
ATOMTYPES_EXPECTED=true # Поставьте false, если типы атомов уже в основном силовом поле
if ${ATOMTYPES_EXPECTED} && [ ! -f "${LIGAND_ATOMTYPES_ITP}" ]; then
     echo "ОШИБКА: Ожидался файл с типами атомов ${LIGAND_ATOMTYPES_ITP}, но он не найден."
     exit 1
fi


# Копируем основу - топологию пептида
cp "${PEPTIDE_TOP_BASE}" "${FINAL_TOPOLOGY}"
if [ $? -ne 0 ]; then echo "ОШИБКА копирования ${PEPTIDE_TOP_BASE}"; exit 1; fi

# --- Вставка #include для лиганда и его типов атомов ---
FF_INCLUDE_LINE=$(grep '#include.*forcefield\.itp"' "${FINAL_TOPOLOGY}")
if [ -z "${FF_INCLUDE_LINE}" ]; then
    echo "ВНИМАНИЕ: Не найдена строка #include для forcefield.itp в ${FINAL_TOPOLOGY}. Попытка вставки в начало."
    # Создаем текст для вставки
    INCLUDE_TEXT=""
    if ${ATOMTYPES_EXPECTED} && [ -f "${LIGAND_ATOMTYPES_ITP}" ]; then
        INCLUDE_TEXT=$(cat "${LIGAND_ATOMTYPES_ITP}")
        INCLUDE_TEXT="${INCLUDE_TEXT}\n" # Добавляем перенос строки
    fi
    INCLUDE_TEXT="${INCLUDE_TEXT}; Include ligand topology\n#include \"${LIGAND_ITP}\""
    # Вставляем в начало (менее надежно)
    echo -e "$INCLUDE_TEXT\n$(cat ${FINAL_TOPOLOGY})" > tmp_topology.txt && mv tmp_topology.txt "${FINAL_TOPOLOGY}"
else
    echo "--> Вставка #include и atomtypes после ${FF_INCLUDE_LINE}"
    # Создаем текст для вставки
    INCLUDE_TEXT=""
    if ${ATOMTYPES_EXPECTED} && [ -f "${LIGAND_ATOMTYPES_ITP}" ]; then
        INCLUDE_TEXT=$(cat "${LIGAND_ATOMTYPES_ITP}")
        INCLUDE_TEXT="${INCLUDE_TEXT}\n"
    fi
    INCLUDE_TEXT="${INCLUDE_TEXT}; Include ligand topology\n#include \"${LIGAND_ITP}\""

    # Используем awk для вставки ПОСЛЕ найденной строки
    awk -v include_block="$INCLUDE_TEXT" -v ff_line="$FF_INCLUDE_LINE" '
    $0 == ff_line { print; print include_block; next }
    { print }
    ' "${FINAL_TOPOLOGY}" > tmp_topology.txt && mv tmp_topology.txt "${FINAL_TOPOLOGY}"
fi

# Проверяем, что include лиганда добавлен
if ! grep -q "#include \"${LIGAND_ITP}\"" "${FINAL_TOPOLOGY}"; then
    echo "ОШИБКА: Не удалось добавить #include для ${LIGAND_ITP} в ${FINAL_TOPOLOGY}"
    exit 1
fi
# Проверяем, что atomtypes добавлены (если ожидались)
if ${ATOMTYPES_EXPECTED} && ! grep -q '\[ atomtypes \]' "${FINAL_TOPOLOGY}" && ! grep -q "#include \".*${LIGAND_ATOMTYPES_ITP}\"" "${FINAL_TOPOLOGY}"; then
     echo "ОШИБКА: Не удалось добавить определения [ atomtypes ] в ${FINAL_TOPOLOGY}"
     # Можно добавить exit 1, если типы критичны
fi

# --- Удаление и добавление секции [ molecules ] ---
echo "--> Подготовка секции [ molecules ] в ${FINAL_TOPOLOGY}"
# Удаляем все со строки [ molecules ] до конца, если она есть
if grep -q '^\[ molecules \]' "${FINAL_TOPOLOGY}"; then
    sed -i '/^\[ molecules \]/,$d' "${FINAL_TOPOLOGY}"
fi
# Удаляем секцию [ system ] тоже, чтобы добавить ее перед molecules
if grep -q '^\[ system \]' "${FINAL_TOPOLOGY}"; then
    sed -i '/^\[ system \]/,$d' "${FINAL_TOPOLOGY}"
fi

# Добавляем правильные секции (без воды и ионов пока)
cat << EOF >> "${FINAL_TOPOLOGY}"

[ system ]
; Name - будет обновлено позже
${SYSTEM_NAME} Initial Structure

[ molecules ]
; Compound        nmols
${PEPTIDE_MOL_NAME}     1
${LIGAND_MOL_NAME}      1
EOF
echo "--> Финальная топология ${FINAL_TOPOLOGY} подготовлена."

# --- 5.B: Сборка системы в боксе с водой и ионами ---
echo "--> 5.B: Сборка системы в боксе..."

# 5.B.1: Комбинирование и создание бокса
echo "--> Позиционирование пептида и лиганда..."
gmx_mpi editconf -f ${PEPTIDE_OPT_GRO} -o ${PEPTIDE_NAME}_center.gro -c -box 6.0 6.0 6.0
gmx_mpi insert-molecules -f ${PEPTIDE_NAME}_center.gro -ci ${LIGAND_GRO} -nmol 1 -o ${SYSTEM_NAME}_combined.gro -try 50 -radius 0.2 -dr 1.0 1.0 1.0
if [ ! -f "${SYSTEM_NAME}_combined.gro" ]; then echo "Ошибка: Не удалось вставить молекулу лиганда."; exit 1; fi
echo "--> Создание финального бокса..."
gmx_mpi editconf -f ${SYSTEM_NAME}_combined.gro -o ${SYSTEM_NAME}_box.gro -c -d 2.0 -bt cubic
if [ $? -ne 0 ]; then echo "Ошибка на этапе editconf (финальный бокс)"; exit 1; fi

# 5.B.2: Сольватация
echo "--> Сольватация системы..."
gmx_mpi solvate -cp ${SYSTEM_NAME}_box.gro -cs spc216.gro -o ${SYSTEM_NAME}_solv.gro -p "${FINAL_TOPOLOGY}"
if [ $? -ne 0 ]; then echo "Ошибка на этапе solvate"; exit 1; fi
NUM_SOL=$(grep "^SOL" "${FINAL_TOPOLOGY}" | awk '{print $2}')
echo "Добавлено ${NUM_SOL} молекул воды."

# 5.B.3: Добавление ионов
echo "--> Добавление ионов..."
printf "integrator = md\nnsteps = 0\n" > ions.mdp
echo "--> Запуск grompp для ионов..."
gmx_mpi grompp -f ions.mdp -c ${SYSTEM_NAME}_solv.gro -p "${FINAL_TOPOLOGY}" -o ${SYSTEM_NAME}_ions.tpr -maxwarn 5
if [ $? -ne 0 ]; then echo "Ошибка на этапе grompp для ионов"; exit 1; fi
echo "--> Запуск genion..."
echo "SOL" | gmx_mpi genion -s ${SYSTEM_NAME}_ions.tpr -o "${FINAL_COORDS}" -p "${FINAL_TOPOLOGY}" -pname NA -nname CL -neutral -conc 0.15
if [ $? -ne 0 ]; then echo "Ошибка на этапе genion"; exit 1; fi
echo "--> Ионы добавлены. Финальные координаты: ${FINAL_COORDS}"

# --- 5.C: Автоматическое создание индексного файла ---
echo "--> 5.C: Создание индексного файла ${FINAL_INDEX}"

# Получаем список групп, запустив make_ndx и передав q
NDX_LIST_OUTPUT=$(echo "q" | gmx_mpi make_ndx -f "${FINAL_COORDS}" -o junk_temp.ndx 2>&1)
rm -f junk_temp.ndx

# --- Извлекаем номера групп из вывода ---
# Используем grep -m 1 для надежности (первое совпадение) и awk
PROTEIN_GROUP_NUM=$(echo "${NDX_LIST_OUTPUT}" | grep -E '^[[:space:]]*[0-9]+[[:space:]]+Protein[[:space:]]+\(' | grep -m 1 Protein | awk '{print $1}')
LIGAND_GROUP_NUM=$(echo "${NDX_LIST_OUTPUT}" | grep -E '^[[:space:]]*[0-9]+[[:space:]]+AAY[[:space:]]+\(' | grep -m 1 AAY | awk '{print $1}') # Используем AAY как имя лиганда
SOL_GROUP_NUM=$(echo "${NDX_LIST_OUTPUT}" | grep -E '^[[:space:]]*[0-9]+[[:space:]]+SOL[[:space:]]+\(' | grep -m 1 SOL | awk '{print $1}')
NA_GROUP_NUM=$(echo "${NDX_LIST_OUTPUT}" | grep -E '^[[:space:]]*[0-9]+[[:space:]]+NA[[:space:]]+\(' | grep -m 1 NA | awk '{print $1}')
CL_GROUP_NUM=$(echo "${NDX_LIST_OUTPUT}" | grep -E '^[[:space:]]*[0-9]+[[:space:]]+CL[[:space:]]+\(' | grep -m 1 CL | awk '{print $1}')
WATER_IONS_GROUP_NUM=$(echo "${NDX_LIST_OUTPUT}" | grep -E '^[[:space:]]*[0-9]+[[:space:]]+Water_and_ions[[:space:]]+\(' | grep -m 1 Water_and_ions | awk '{print $1}')

# Проверка и запасные варианты (можно улучшить, но для начала так)
[[ -z "${PROTEIN_GROUP_NUM}" ]] && PROTEIN_GROUP_NUM=1 # Предполагаем 1, если не найдено
[[ -z "${LIGAND_GROUP_NUM}" ]] && LIGAND_GROUP_NUM=13 # Предполагаем 13 (AAY), если не найдено
[[ -z "${SOL_GROUP_NUM}" ]] && SOL_GROUP_NUM=17     # Предполагаем 17
[[ -z "${NA_GROUP_NUM}" ]] && NA_GROUP_NUM=14       # Предполагаем 14
[[ -z "${CL_GROUP_NUM}" ]] && CL_GROUP_NUM=15       # Предполагаем 15
[[ -z "${WATER_IONS_GROUP_NUM}" ]] && WATER_IONS_GROUP_NUM="" # Оставляем пустым, если готовой нет

echo "Обнаружены/Предположены группы: Protein=${PROTEIN_GROUP_NUM}, Ligand=${LIGAND_GROUP_NUM}, SOL=${SOL_GROUP_NUM}, NA=${NA_GROUP_NUM}, CL=${CL_GROUP_NUM}, Water_and_ions=${WATER_IONS_GROUP_NUM}"

# --- Создаем команды для make_ndx ---
MAKE_NDX_COMMANDS="${PROTEIN_GROUP_NUM} | ${LIGAND_GROUP_NUM}\nname Protein_Ligand\n" # Объединяем белок и лиганд

# Проверяем, есть ли готовая группа Water_and_ions
if [ -z "${WATER_IONS_GROUP_NUM}" ]; then
    echo "Группа Water_and_ions не найдена, создаем вручную..."
    MAKE_NDX_COMMANDS+="${SOL_GROUP_NUM} | ${NA_GROUP_NUM} | ${CL_GROUP_NUM}\nname Water_and_ions\n" # Объединяем воду и ионы
else
    echo "Используем существующую группу Water_and_ions (номер ${WATER_IONS_GROUP_NUM})."
    # Ничего не добавляем, она уже есть
fi
MAKE_NDX_COMMANDS+="q\n" # Команда выхода

# --- Запускаем make_ndx с командами ---
echo "--> Запуск make_ndx для создания финального ${FINAL_INDEX}..."
echo -e "${MAKE_NDX_COMMANDS}" | gmx_mpi make_ndx -f "${FINAL_COORDS}" -o "${FINAL_INDEX}"

# Финальная проверка
if [ $? -ne 0 ] || [ ! -f "${FINAL_INDEX}" ]; then
    echo "ОШИБКА: Не удалось создать индексный файл ${FINAL_INDEX} автоматически."
    echo "Попробуйте создать его вручную: gmx_mpi make_ndx -f ${FINAL_COORDS} -o ${FINAL_INDEX}"
    exit 1
fi
# Проверяем наличие созданных групп (хотя бы одной)
if ! grep -q '\[ Protein_Ligand \]' "${FINAL_INDEX}" || ! grep -q '\[ Water_and_ions \]' "${FINAL_INDEX}"; then
    echo "ВНИМАНИЕ: Ожидаемые группы Protein_Ligand или Water_and_ions могут отсутствовать в ${FINAL_INDEX}. Проверьте файл вручную!"
    # Не выходим с ошибкой, но предупреждаем
fi

echo "--> Индексный файл ${FINAL_INDEX} создан."

# --- 5.D: Очистка временных файлов ---
echo "--> 5.D: Очистка временных файлов..."
rm -f ${PEPTIDE_NAME}_center.gro ${SYSTEM_NAME}_combined.gro ${SYSTEM_NAME}_box.gro ${SYSTEM_NAME}_solv.gro ${SYSTEM_NAME}_ions.tpr ions.mdp \#*# *~ mdout.mdp
# Удаляем временный файл atomtypes, если он был создан и ожидался
if ${ATOMTYPES_EXPECTED} && [ -f "${LIGAND_ATOMTYPES_ITP}" ]; then
   echo "--> Удаление ${LIGAND_ATOMTYPES_ITP}"
   rm -f "${LIGAND_ATOMTYPES_ITP}"
fi

echo "=========================================================="
echo " ШАГ 5 (АВТОМАТИЧЕСКИЙ) ЗАВЕРШЕН."
echo "=========================================================="
echo "Готовые файлы для Шага 6 (Минимизация):"
echo "  Координаты: ${FINAL_COORDS}"
echo "  Топология:  ${FINAL_TOPOLOGY}"

