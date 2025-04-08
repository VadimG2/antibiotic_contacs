#!/bin/bash

# --- Имена ---
SYS_NAME="pept1_lig1_bcc_manual" # Базовое имя системы
TOP_FILE="pept1.top"

# --- Шаг 6.1: Минимизация ---
echo "--> Шаг 6.1: Минимизация энергии..."
gmx_mpi grompp -f step6.1_minimization.mdp -c ${SYS_NAME}_solv_ions.gro -p ${TOP_FILE} -o ${SYS_NAME}_min.tpr -maxwarn 5
if [ $? -ne 0 ]; then echo "Ошибка grompp на этапе минимизации"; exit 1; fi
gmx_mpi mdrun -deffnm ${SYS_NAME}_min -v -nb gpu
if [ $? -ne 0 ]; then echo "Ошибка mdrun на этапе минимизации"; exit 1; fi
echo "Минимизация завершена."

# --- Шаг 6.2: NVT Уравновешивание (Требует index.ndx и posre*.itp!) ---
echo "--> Шаг 6.2: NVT Уравновешивание (100 ps)..."
# Убедитесь, что в nvt.mdp указаны ПРАВИЛЬНЫЕ tc-grps из вашего index.ndx!
gmx_mpi grompp -f step6.2_equilibration_nvt.mdp -c ${SYS_NAME}_min.gro -r ${SYS_NAME}_min.gro -p ${TOP_FILE} -n index.ndx -o ${SYS_NAME}_nvt.tpr -maxwarn 5
# Флаг -r нужен для позиционных ограничений (reference structure)
# Флаг -n указывает на индексный файл с группами температуры
if [ $? -ne 0 ]; then echo "Ошибка grompp на этапе NVT"; exit 1; fi
gmx_mpi mdrun -deffnm ${SYS_NAME}_nvt -v -nb gpu -pme gpu -update gpu -bonded gpu
if [ $? -ne 0 ]; then echo "Ошибка mdrun на этапе NVT"; exit 1; fi
echo "NVT уравновешивание завершено."

# --- Шаг 6.3: NPT Уравновешивание (Требует index.ndx и posre*.itp!) ---
echo "--> Шаг 6.3: NPT Уравновешивание (500 ps)..."
# Убедитесь, что в npt.mdp указаны ПРАВИЛЬНЫЕ tc-grps и pcoupl-grps (если нужно) из index.ndx!
gmx_mpi grompp -f step6.3_equilibration_npt.mdp -c ${SYS_NAME}_nvt.gro -r ${SYS_NAME}_nvt.gro -t ${SYS_NAME}_nvt.cpt -p ${TOP_FILE} -n index.ndx -o ${SYS_NAME}_npt.tpr -maxwarn 5
# Флаг -t передает checkpoint файл с информацией о скоростях и состоянии с NVT этапа
if [ $? -ne 0 ]; then echo "Ошибка grompp на этапе NPT"; exit 1; fi
gmx_mpi mdrun -deffnm ${SYS_NAME}_npt -v -nb gpu -pme gpu -update gpu -bonded gpu
if [ $? -ne 0 ]; then echo "Ошибка mdrun на этапе NPT"; exit 1; fi
echo "NPT уравновешивание завершено."

# --- Шаг 6.4: Продуктивная МД (200 ns) ---
echo "--> Шаг 6.4: Продуктивная МД (200 ns)..."
# Убедитесь, что в prod.mdp ОТКЛЮЧЕНЫ posres и указаны ПРАВИЛЬНЫЕ группы!
gmx_mpi grompp -f step6.4_production.mdp -c ${SYS_NAME}_npt.gro -t ${SYS_NAME}_npt.cpt -p ${TOP_FILE} -n index.ndx -o ${SYS_NAME}_prod.tpr -maxwarn 5
if [ $? -ne 0 ]; then echo "Ошибка grompp на этапе Production"; exit 1; fi
echo "Запуск продуктивной МД... Это займет много времени!"
gmx_mpi mdrun -deffnm ${SYS_NAME}_prod -v -nb gpu -pme gpu -update gpu -bonded gpu
# Здесь может потребоваться запуск на кластере или в фоне
lig1_bcc_manual.amb2gmxif [ $? -ne 0 ]; then echo "Ошибка mdrun на этапе Production"; exit 1; fi
echo "Продуктивная МД завершена."

echo "--> Удаление временных файлов..."

rm *#*

echo "=========================================================="
echo " ВСЕ ЭТАПЫ СИМУЛЯЦИИ ЗАВЕРШЕНЫ "
echo "=========================================================="
echo "Основной результат - траектория: ${SYS_NAME}_prod.xtc"
echo "Используйте ее для Шага 7 (Анализ)."
