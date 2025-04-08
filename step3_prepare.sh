#!/bin/bash

# --- Имена файлов ---
PEPTIDE_GRO="pept1.gro"         # Входные координаты от pdb2gmx_mpi
PEPTIDE_TOP="pept1.top"         # Входная топология от pdb2gmx_mpi (с ручными правками для циклизации, если нужно!)
PEPTIDE_NAME="pept1"            # Базовое имя
OPTIMIZED_GRO="${PEPTIDE_NAME}_optimized.gro" # Финальный выход этого скрипта

# --- Параметры ---
ANNEAL_TIME_PS=500   # Общее время отжига в пикосекундах (как в PDF)
ANNEAL_TEMP_LOW=5    # Начальная температура (К)
ANNEAL_TEMP_HIGH=400 # Максимальная температура (К)
ANNEAL_TIME_RAMP_PS=250 # Время нагрева (половина общего времени, как в PDF)

echo "--- ШАГ 3: Оптимизация геометрии пептида (Отжиг) ---"

# Шаг 3.1: Создание бокса для симуляции
echo "--> Шаг 3.1: Создание бокса для ${PEPTIDE_NAME}"
gmx_mpi editconf -f ${PEPTIDE_GRO} -o ${PEPTIDE_NAME}_box.gro -c -d 1.0 -bt cubic
if [ $? -ne 0 ]; then echo "ОШИБКА на этапе editconf"; exit 1; fi

# Шаг 3.2: Создание файла параметров для отжига (anneal.mdp)
echo "--> Шаг 3.2: Создание anneal.mdp"
cat << EOF > anneal.mdp
; --- Параметры для отжига пептида в "вакууме" (NPT) ---
title           = Peptide Annealing
integrator      = md            ; Molecular dynamics
nsteps          = $(($ANNEAL_TIME_PS * 1000)) ; Время = $ANNEAL_TIME_PS ps, шаг = 1 fs (1000 steps/ps)
dt              = 0.001         ; 1 fs time step
nstxout-compressed  = 5000
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


# Шаг 3.3: Запуск grompp (подготовка бинарного файла .tpr)
echo "--> Шаг 3.3: Запуск grompp для отжига"
# Проверяем, существует ли группа Protein в топологии
# Эта проверка теперь менее важна, так как tc-grps ниже все равно будет проверена grompp
if ! grep -q 'Protein' ${PEPTIDE_TOP}; then
  echo "ВНИМАНИЕ: Группа 'Protein' не найдена в ${PEPTIDE_TOP}. Убедитесь, что она используется в anneal.mdp (tc-grps)!"
fi

gmx_mpi grompp -f anneal.mdp -c ${PEPTIDE_NAME}_box.gro -p ${PEPTIDE_TOP} -o ${PEPTIDE_NAME}_anneal.tpr -maxwarn 2
if [ $? -ne 0 ]; then echo "ОШИБКА на этапе grompp"; exit 1; fi

# Шаг 3.4: Запуск mdrun (выполнение симуляции отжига)
echo "--> Шаг 3.4: Запуск mdrun для отжига (${ANNEAL_TIME_PS} ps)..."
# Запускаем без GPU, так как это быстрая симуляция пептида
gmx_mpi mdrun -deffnm ${PEPTIDE_NAME}_anneal -v
if [ $? -ne 0 ]; then echo "ОШИБКА на этапе mdrun"; exit 1; fi

#Шаг 3.5: Извлечение последнего кадра как оптимизированной структуры
echo "--> Шаг 3.5: Извлечение финального кадра отжига"

# --- АВТОМАТИЧЕСКИЙ ВЫБОР ГРУППЫ (Исправленный) ---
# Запускаем make_ndx, передаем 'q' через stdin, ловим stdout в переменную
NDX_OUTPUT=$(echo "q" | gmx_mpi make_ndx -f ${PEPTIDE_NAME}_anneal.gro -o junk_index.ndx 2>&1) # Перенаправляем stderr в stdout
rm -f junk_index.ndx # Удаляем ненужный файл индекса

# Ищем номер группы 'Protein' в пойманном выводе
PROTEIN_GROUP_NUM=$(echo "${NDX_OUTPUT}" | grep "Protein" | head -n 1 | awk '{print $1}')

# Проверяем, нашли ли номер. Если нет, используем 1 (Protein) или 0 (System) как запасной вариант.
if [[ -z "${PROTEIN_GROUP_NUM}" || ! "${PROTEIN_GROUP_NUM}" =~ ^[0-9]+$ ]]; then
    echo "ВНИМАНИЕ: Не удалось автоматически определить номер группы 'Protein'. Используем группу 1."
    PROTEIN_GROUP_NUM=1
fi
echo "Для trjconv будет выбрана группа с номером: ${PROTEIN_GROUP_NUM}"

# Передаем номер группы через echo в trjconv
echo "${PROTEIN_GROUP_NUM}" | gmx_mpi trjconv -s ${PEPTIDE_NAME}_anneal.tpr -f ${PEPTIDE_NAME}_anneal.xtc -o ${OPTIMIZED_GRO} -dump $ANNEAL_TIME_PS
# ----------------------------------

if [ $? -ne 0 ] || [ ! -f "${OPTIMIZED_GRO}" ]; then
    echo "ОШИБКА на этапе trjconv или файл ${OPTIMIZED_GRO} не создан."
    exit 1
fi

# Шаг 3.6: Очистка (опционально)
echo "--> Шаг 3.6: Очистка промежуточных файлов отжига..."
rm -f ${PEPTIDE_NAME}_box.gro anneal.mdp ${PEPTIDE_NAME}_anneal.log ${PEPTIDE_NAME}_anneal.xtc ${PEPTIDE_NAME}_anneal.edr ${PEPTIDE_NAME}_anneal.tpr ${PEPTIDE_NAME}_anneal.gro ${PEPTIDE_NAME}_anneal.cpt \#*# *~ mdout.mdp
echo "Очистка завершена."

echo "--------------------------------------------------------------------"
echo "Шаг 3 (Отжиг пептида) ЗАВЕРШЕН."
echo "Оптимизированные координаты сохранены в: ${OPTIMIZED_GRO}"
echo "Топология (проверьте правки для циклизации!): ${PEPTIDE_TOP}"
echo "--------------------------------------------------------------------"
