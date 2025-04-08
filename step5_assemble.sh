#!/bin/bash

# --- Имена файлов пептида (результат Шага 3) ---
PEPTIDE_OPT_GRO="pept1_optimized.gro" # <--- Оптимизированные координаты
PEPTIDE_TOP="pept1.top"               # <--- Топология пептида
PEPTIDE_NAME="pept1"
# -------------------------------------------------

# --- Имена файлов лиганда (антибиотика) (результат Шага 4 + конвертации) ---
ANTIBIOTIC_GRO="lig1_bcc_manual_gmx.gro" # Координаты лиганда
ANTIBIOTIC_ITP="lig1_bcc_manual.itp"    # Топология лиганда (ITP)
# Убедитесь, что это имя совпадает с именем в [ moleculetype ] внутри ANTIBIOTIC_ITP
ANTIBIOTIC_MOL_NAME="lig1_bcc_manual" # <--- ПРОВЕРЬТЕ ИМЯ ВНУТРИ ITP!
# ---------------------------------------------

SYSTEM_NAME="${PEPTIDE_NAME}_${ANTIBIOTIC_MOL_NAME}" # Имя системы

# --- Шаг 5: Сборка системы ---

echo "=========================================================="
echo " НАЧИНАЕМ ШАГ 5: СБОРКА СИСТЕМЫ ${SYSTEM_NAME}"
echo "=========================================================="

# 5.1: Комбинирование и создание бокса
echo "--> Шаг 5.1: Позиционирование пептида и лиганда..."
# Создаем временный бокс для пептида
gmx_mpi editconf -f ${PEPTIDE_OPT_GRO} -o ${PEPTIDE_NAME}_center.gro -c -box 6.0 6.0 6.0 #-bt cubic
# Вставляем лиганд со случайным смещением
gmx_mpi insert-molecules -f ${PEPTIDE_NAME}_center.gro -ci ${ANTIBIOTIC_GRO} -nmol 1 -o ${SYSTEM_NAME}_combined.gro -try 50 -radius 0.2 -dr 1.0 1.0 1.0
if [ ! -f "${SYSTEM_NAME}_combined.gro" ]; then echo "Ошибка: Не удалось вставить молекулу лиганда."; exit 1; fi

# Создаем финальный бокс симуляции (2.0 нм от края)
echo "--> Шаг 5.2: Создание финального бокса симуляции..."
gmx_mpi editconf -f ${SYSTEM_NAME}_combined.gro -o ${SYSTEM_NAME}_box.gro -c -d 2.0 -bt cubic
if [ $? -ne 0 ]; then echo "Ошибка на этапе editconf (финальный бокс)"; exit 1; fi

# 5.3: Редактирование топологии пептида (${PEPTIDE_TOP})
echo "--------------------------------------------------------------------"
echo "ДЕЙСТВИЕ: Отредактируйте файл ${PEPTIDE_TOP}"
echo "1. Добавьте строку ПЕРЕД секцией [ system ]:"
echo "   ;; Include ligand topology"
echo "   #include \"./${ANTIBIOTIC_ITP}\""
echo "2. Обновите секцию [ molecules ] в КОНЦЕ файла (ЗАМЕНИТЕ существующую секцию):"
echo "   [ molecules ]"
echo "   ; Compound        nmols"
echo "   Protein_chain_A   1      ; (или как называется ваш пептид в .top)"
echo "   ${ANTIBIOTIC_MOL_NAME}         1      ; (Используйте имя из ANTIBIOTIC_MOL_NAME)"
echo "--------------------------------------------------------------------"
read -p "Нажмите [Enter] после редактирования файла ${PEPTIDE_TOP}..."

# 5.4: Сольватация (добавление воды)
echo "--> Шаг 5.4: Сольватация системы..."
# Сначала проверим, существует ли уже запись SOL в molecules. Если да - удалим перед сольватацией.
if grep -q "^SOL" ${PEPTIDE_TOP}; then
    echo "Обнаружена строка SOL в топологии. Удаляем перед повторной сольватацией..."
    sed -i '/^SOL/d' ${PEPTIDE_TOP}
fi
gmx_mpi solvate -cp ${SYSTEM_NAME}_box.gro -cs spc216.gro -o ${SYSTEM_NAME}_solv.gro -p ${PEPTIDE_TOP}
if [ $? -ne 0 ]; then echo "Ошибка на этапе solvate"; exit 1; fi
# Узнаем, сколько молекул воды было добавлено
NUM_SOL=$(grep "^SOL" ${PEPTIDE_TOP} | awk '{print $2}')
echo "Добавлено ${NUM_SOL} молекул воды."

# 5.5: Добавление ионов
echo "--> Шаг 5.5: Добавление ионов..."
# Создаем корректный ions.mdp
echo "--> Запуск grompp для ионов..."
gmx_mpi grompp -f ions.mdp -c ${SYSTEM_NAME}_solv.gro -p ${PEPTIDE_TOP} -o ${SYSTEM_NAME}_ions.tpr -maxwarn 5
if [ $? -ne 0 ]; then echo "Ошибка на этапе grompp для ионов"; exit 1; fi

echo "--> Запуск genion... Выберите группу 'SOL' для замены на ионы."
# Используем 'SOL' для замены молекул воды ионами
echo "SOL" | gmx_mpi genion -s ${SYSTEM_NAME}_ions.tpr -o ${SYSTEM_NAME}_solv_ions.gro -p ${PEPTIDE_TOP} -pname NA -nname CL -neutral -conc 0.15
if [ $? -ne 0 ]; then echo "Ошибка на этапе genion"; exit 1; fi

echo "--> Создание индексных файлов"
echo "--> Сейчас вам нужно ввести цифры, которым соответствует группа Protein и группа AYY в формате 1 | 2 из списка, идущего далее, в консоль"
echo "--> Введите цифры групп, соответствующих пептиду и лиганду, через '|' и нажмите q для ВЫХОДА"
gmx_mpi make_ndx -f pept1_lig1_bcc_manual_solv_ions.gro -o index_updated.ndx
mv index_updated.ndx index.ndx

# 5.6: Очистка временных файлов
echo "--> Шаг 5.6: Очистка временных файлов..."
rm -f ${PEPTIDE_NAME}_center.gro ${SYSTEM_NAME}_combined.gro ${SYSTEM_NAME}_ions.tpr ions.mdp \#*# *~ mdout.mdp

echo "=========================================================="
echo " ШАГ 5 (СБОРКА СИСТЕМЫ) ЗАВЕРШЕН."
echo "=========================================================="
echo "Готовые файлы для следующего шага (Шаг 6 - Минимизация):"
echo "  Координаты: ${SYSTEM_NAME}_solv_ions.gro"
echo "  Топология:  ${PEPTIDE_TOP} (обновленная)"
echo "----------------------------------------------------------"
