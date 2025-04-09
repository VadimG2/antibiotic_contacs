#!/bin/bash

# --- Пути к директориям ---
INPUT_DIR="output"  # Папка, где лежат папки с подготовленными пептидами (pept1, pept2...)
OUTPUT_DIR="output" # Сюда же будем класть optimized.gro

echo "=== НАЧАЛО ОПТИМИЗАЦИИ ГЕОМЕТРИИ ПЕПТИДОВ (ОТЖИГ - ТЕСТ 5ps) ==="

# --- Параметры ---
# --- ДЛЯ ТЕСТА ---
TEST_STEPS=5000           # 5000 шагов
TEST_TIME_PS=5             # Соответствует 5 ps при dt=0.001
TEST_DUMP_TIME=5           # Время для извлечения кадра
# --- Финальные параметры (закомментированы) ---
# ANNEAL_STEPS=500000        # 500 ps
# ANNEAL_TIME_PS=500
# ANNEAL_DUMP_TIME=500
# --- Общие параметры отжига ---
ANNEAL_TEMP_LOW=5
ANNEAL_TEMP_HIGH=400
ANNEAL_TIME_RAMP_PS=2.5 # Для 5ps теста делаем нагрев за половину времени

# --- Определяем параметры для текущего запуска (тест или полный) ---
NSTEPS=${TEST_STEPS}
DUMP_TIME=${TEST_DUMP_TIME}
ANNEAL_TOTAL_TIME=${TEST_TIME_PS} # Общее время для mdp (в ps)
# ANNEAL_TOTAL_TIME=${ANNEAL_TIME_PS} # Для полного запуска
ANNEAL_RAMP_TIME=$(echo "$ANNEAL_TOTAL_TIME / 2" | bc -l) # Время нагрева - половина

# Итерация по всем папкам пептидов в output/
for peptide_dir in "${INPUT_DIR}"/*; do
    if [ -d "${peptide_dir}" ]; then
        peptide_base=$(basename "${peptide_dir}")
        PEPTIDE_GRO="${peptide_dir}/${peptide_base}.gro"  # Входной файл от шага 2
        PEPTIDE_TOP="${peptide_dir}/${peptide_base}.top"  # Входной файл от шага 2
        OPTIMIZED_GRO="${peptide_dir}/${peptide_base}_optimized.gro"  # Выходной файл

        # Проверка наличия входных файлов
        if [ ! -f "${PEPTIDE_GRO}" ] || [ ! -f "${PEPTIDE_TOP}" ]; then
            echo "Пропуск: ${peptide_dir} не содержит файлов пептида (${peptide_base}.gro или ${peptide_base}.top)"
            continue
        fi

        echo ""
        echo "--- ШАГ 3: Оптимизация геометрии пептида ${peptide_base} (Отжиг - Тест ${ANNEAL_TOTAL_TIME} ps) ---"

        # Удаляем старые файлы отжига перед запуском
        rm -f "${peptide_dir}/${peptide_base}_box.gro" "${peptide_dir}/anneal.mdp" "${peptide_dir}/${peptide_base}_anneal."* "${peptide_dir}/mdout.mdp" "${OPTIMIZED_GRO}"

        # Шаг 3.1: Создание бокса
        echo "--> Шаг 3.1: Создание бокса для ${peptide_base}"
        gmx_mpi editconf -f "${PEPTIDE_GRO}" -o "${peptide_dir}/${peptide_base}_box.gro" -c -d 1.0 -bt cubic
        if [ $? -ne 0 ]; then echo "ОШИБКА на этапе editconf для ${peptide_base}"; continue; fi

        # Шаг 3.2: Создание anneal.mdp
        echo "--> Шаг 3.2: Создание anneal.mdp для ${peptide_base}"
        cat << EOF > "${peptide_dir}/anneal.mdp"
; --- Параметры для отжига пептида (Тест ${ANNEAL_TOTAL_TIME} ps) ---
integrator      = md
nsteps          = ${NSTEPS} ; <--- Используем переменную для числа шагов
dt              = 0.001
nstxout-compressed  = 1000 ; Сохраняем чаще для короткого теста (каждые 1 ps)
nstlog          = 1000
nstenergy       = 1000
xtc-precision   = 1000  ; <-- Добавлено для явного создания XTC
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
tc-grps         = System ; Используем System, т.к. только пептид
tau_t           = 0.1
ref_t           = $ANNEAL_TEMP_HIGH
pcoupl          = Berendsen ; Используем Berendsen, как было в anneal.mdp ранее
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = 1.0
compressibility = 4.5e-5
annealing       = single
annealing-npoints = 3
annealing-time  = 0 ${ANNEAL_RAMP_TIME} ${ANNEAL_TOTAL_TIME} ; Время отжига
annealing-temp  = $ANNEAL_TEMP_LOW $ANNEAL_TEMP_HIGH $ANNEAL_TEMP_HIGH ; Температуры отжига
gen_vel         = yes
gen_temp        = $ANNEAL_TEMP_LOW
gen_seed        = -1
pbc             = xyz
EOF

        # Шаг 3.3: Запуск grompp
        echo "--> Шаг 3.3: Запуск grompp для отжига ${peptide_base}"
        gmx_mpi grompp -f "${peptide_dir}/anneal.mdp" -c "${peptide_dir}/${peptide_base}_box.gro" -p "${PEPTIDE_TOP}" -o "${peptide_dir}/${peptide_base}_anneal.tpr" -maxwarn 2
        if [ $? -ne 0 ]; then echo "ОШИБКА на этапе grompp для ${peptide_base}"; continue; fi

        # Шаг 3.4: Запуск mdrun (без явного указания GPU/CPU, как в рабочем варианте)
        echo "--> Шаг 3.4: Запуск mdrun для отжига (${ANNEAL_TOTAL_TIME} ps) ${peptide_base}"
        gmx_mpi mdrun -deffnm "${peptide_dir}/${peptide_base}_anneal" -v
        if [ $? -ne 0 ]; then echo "ОШИБКА на этапе mdrun для ${peptide_base}"; continue; fi

        # Проверка, дошла ли симуляция до конца
        LAST_STEP=$(grep "Statistics over " "${peptide_dir}/${peptide_base}_anneal.log" | tail -n 1 | awk '{print $4}')
        if [ "${LAST_STEP}" != "${NSTEPS}" ]; then
             echo "ОШИБКА: Отжиг завершился на шаге ${LAST_STEP}, а не ${NSTEPS}."
             continue
        fi

        # Шаг 3.5: Извлечение последнего кадра
        echo "--> Шаг 3.5: Извлечение финального кадра (${DUMP_TIME} ps) отжига для ${peptide_base}"
        # Пытаемся извлечь из .xtc
        echo "0" | gmx_mpi trjconv -s "${peptide_dir}/${peptide_base}_anneal.tpr" -f "${peptide_dir}/${peptide_base}_anneal.xtc" -o "${OPTIMIZED_GRO}" -dump ${DUMP_TIME} &> /dev/null # Скрываем вывод trjconv
        # Проверяем, создался ли файл и не пустой ли он
        if [ ! -s "${OPTIMIZED_GRO}" ]; then
            echo "Предупреждение: Не удалось извлечь из .xtc или файл пуст. Попытка извлечь из .trr..."
            if [ -f "${peptide_dir}/${peptide_base}_anneal.trr" ]; then
                echo "0" | gmx_mpi trjconv -s "${peptide_dir}/${peptide_base}_anneal.tpr" -f "${peptide_dir}/${peptide_base}_anneal.trr" -o "${OPTIMIZED_GRO}" -dump ${DUMP_TIME} &> /dev/null
                if [ ! -s "${OPTIMIZED_GRO}" ]; then
                    echo "ОШИБКА: Не удалось извлечь кадр и из .trr."; continue;
                fi
                 echo "Успешно извлечено из .trr."
            else
                echo "ОШИБКА: Файл .trr не найден. Не удалось извлечь финальный кадр."; continue;
            fi
        else
             echo "Успешно извлечено из .xtc."
        fi


        # Шаг 3.6: Очистка (опционально, можно закомментировать)
        # echo "--> Шаг 3.6: Очистка промежуточных файлов для ${peptide_base}"
        # rm -f "${peptide_dir}/${peptide_base}_box.gro" "${peptide_dir}/anneal.mdp" "${peptide_dir}/${peptide_base}_anneal.log" "${peptide_dir}/${peptide_base}_anneal.xtc" "${peptide_dir}/${peptide_base}_anneal.trr" "${peptide_dir}/${peptide_base}_anneal.edr" "${peptide_dir}/${peptide_base}_anneal.tpr" "${peptide_dir}/${peptide_base}_anneal.gro" "${peptide_dir}/${peptide_base}_anneal.cpt" "${peptide_dir}/#*" "${peptide_dir}/*~"

        echo "--- ШАГ 3 ЗАВЕРШЕН: Оптимизация пептида ${peptide_base} завершена (Тест ${ANNEAL_TOTAL_TIME} ps) ---"
    fi
done

echo ""
echo "=== ОПТИМИЗАЦИЯ ВСЕХ ПЕПТИДОВ ЗАВЕРШЕНА (Тест ${ANNEAL_TOTAL_TIME} ps) ==="
