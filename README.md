# Моделирование взаимодействия пептидов с антибиотиками (GROMACS)

Этот репозиторий содержит набор скриптов для автоматизированного моделирования молекулярной динамики (МД) комплексов пептид-антибиотик с использованием GROMACS и анализа результатов.

## Обзор пайплайна

Пайплайн выполняет следующие шаги для каждой пары пептид-лиганд, найденной во входных директориях:

1.  **Подготовка пептида:**
    *   Генерация GROMACS топологии (`.top`), координат (`.gro`) и файла позиционных ограничений (`posre_*.itp`) из PDB файла пептида с помощью `gmx pdb2gmx`.
    *   *Требуется ручная проверка/правка `.top` файла для обеспечения цикличности пептида (если применимо).*
    *   Добавление директивы `#ifdef POSRES` в `.top` файл для использования файла позиционных ограничений.
2.  **Оптимизация пептида (Отжиг):**
    *   Проведение короткой МД симуляции (отжига) в вакууме для релаксации начальной структуры пептида.
    *   Извлечение последней структуры как оптимизированной (`_optimized.gro`).
3.  **Параметризация лиганда:**
    *   Конвертация PDB файла лиганда в MOL2 с помощью `Open Babel`.
    *   Назначение типов атомов GAFF2 и расчет BCC зарядов с помощью `antechamber` (из AmberTools).
    *   Поиск недостающих параметров с помощью `parmchk2`.
    *   Сборка AMBER топологии (`.prmtop`) и координат (`.inpcrd`) с помощью `tleap`.
    *   Конвертация AMBER файлов в формат GROMACS (`.itp`, `_GMX.gro`) с помощью `ACPYPE`.
    *   Извлечение определений типов атомов GAFF2 (`atomtypes_*.itp`) и позиционных ограничений (`posre_*.itp`) из промежуточных файлов ACPYPE.
    *   Коррекция имени молекулы и пути к файлу позиционных ограничений внутри `.itp` файла лиганда.
4.  **Сборка системы:**
    *   Создание уникальной директории для симуляции (`output/ПЕПТИД+ЛИГАНД_Sim`).
    *   Копирование необходимых файлов пептида и лиганда в директорию симуляции.
    *   Объединение координат пептида и лиганда (`insert-molecules`).
    *   Создание симуляционного бокса (`editconf`).
    *   **Автоматическая генерация файла `system.top`:**
        *   Включение необходимых файлов силового поля (`amber14sb.ff`), типов атомов GAFF2 (`gaff2_atomtypes.itp`), лиганда (`ligand.itp`), воды (`tip3p.itp`), ионов (`ions.itp`).
        *   Включение определения `moleculetype` пептида (с его PosRe).
        *   Добавление секций `[ system ]` и `[ molecules ]` (с пептидом и лигандом).
    *   Сольватация системы водой (`solvate`).
    *   Добавление ионов для нейтрализации и достижения нужной концентрации (`grompp`, `genion`).
5.  **Запуск МД симуляций (выполняется отдельным скриптом):**
    *   Перемещение собранной системы в директорию `simulations/`.
    *   Создание индексного файла (`make_ndx`) с группами `Protein_Ligand` и `SOL_Ion`.
    *   Копирование и адаптация MDP файлов (для `minim`, `nvt`, `npt`, `prod`) с правильными группами температуры/давления.
    *   Последовательный запуск этапов:
        *   Минимизация энергии (`grompp`, `mdrun`).
        *   NVT уравновешивание (`grompp`, `mdrun` с PosRes).
        *   NPT уравновешивание (`grompp`, `mdrun` с PosRes).
        *   Продуктивная МД (`grompp`, `mdrun` без PosRes).
6.  **Анализ (выполняется отдельным скриптом):**
    *   Генерация карт контактов между остатками пептида и лигандом на основе траектории продуктивной МД с использованием Python (`MDAnalysis`, `matplotlib`).

## Требования

*   **Conda:** Установленный Miniconda или Anaconda.
*   **Системные зависимости:**
    *   Установленный **NVIDIA CUDA Toolkit** (версия, совместимая с GROMACS, например, 11.8 или новее).
    *   Соответствующие **драйверы NVIDIA**.
    *   Утилиты `wget` и `tar` (обычно присутствуют в Linux/WSL).
    *   Утилита `bc` (обычно присутствует).
*   **Git:** Для клонирования репозитория.

## Установка и запуск

1.  **Клонировать репозиторий:**
    ```bash
    git clone https://github.com/VadimG2/antibiotic_contacs.git
    cd antibiotic_contacs
    ```

2.  **Создать и активировать Conda окружение:**
    ```bash
    conda env create -f env.yaml
    conda activate gromacs_md_env # Используйте имя окружения, указанное в env.yaml
    ```
    *Файл `env.yaml` содержит список всех необходимых пакетов Conda (GROMACS с CUDA, AmberTools, ACPYPE, Open Babel, Python библиотеки и т.д.).*

3.  **Установить силовое поле AMBER14SB:**
    *Пакет GROMACS из Conda Forge не всегда включает все силовые поля. Этот шаг скачивает и устанавливает AMBER14SB в созданное окружение.*
    ```bash
    wget -O ~/miniconda3/envs/gromacs_md_env/share/gromacs/top/amber14sb.ff.tar.gz https://ftp.gromacs.org/contrib/forcefields/amber14sb.ff.tar.gz
    tar -xvf ~/miniconda3/envs/gromacs_md_env/share/gromacs/top/amber14sb.ff.tar.gz -C ~/miniconda3/envs/gromacs_md_env/share/gromacs/top/
    rm ~/miniconda3/envs/gromacs_md_env/share/gromacs/top/amber14sb.ff.tar.gz
    ```

4.  **Подготовить входные файлы:**
    *   Поместите PDB файлы ваших **пептидов** в директорию `input/protein/`.
    *   Поместите PDB файлы ваших **лигандов** (антибиотиков) в директорию `input/ligand/`.

5.  **Очистить предыдущие результаты (если нужно перезапустить всё):**
    ```bash
    rm -rf simulations/* output/*
    # Будьте осторожны, эта команда удалит все предыдущие результаты!
    ```

6  **Запустить основной скрипт пайплайна:**
    ```bash
    bash query.sh # Или как называется ваш основной скрипт
    ```

## Структура директорий

*   **`input/protein/`**: Исходные PDB файлы пептидов.
*   **`input/ligand/`**: Исходные PDB файлы лигандов.
*   **`output/`**: Директория для промежуточных и финальных файлов.
    *   `output/ИМЯ_ПЕПТИДА/`: Файлы подготовленного пептида (`.gro`, `.top`, `_optimized.gro`, `posre_*.itp`).
    *   `output/ИМЯ_ЛИГАНДА/`: Файлы подготовленного лиганда (`.itp`, `_GMX.gro`, `atomtypes_*.itp`, `posre_*.itp`).
    *   `output/ИМЯ_ПЕПТИДА+ИМЯ_ЛИГАНДА_Sim/`: Папка с собранной системой перед запуском симуляции (если скрипт сборки не перемещает ее сразу).
*   **`simulations/`**: Директория для запуска МД симуляций.
    *   `simulations/ИМЯ_ПЕПТИДА+ИМЯ_ЛИГАНДА/`: Рабочая папка для конкретной симуляции, содержит `.tpr`, `.log`, `.xtc`, `.edr`, `index.ndx` и т.д.
*   **`scripts/`** (Рекомендуется): Папка для хранения bash-скриптов пайплайна (`prepare_peptides.sh`, `prepare_ligands.sh`, `assemble_systems.sh`, `run_simulations.sh`, `analyze_contacts.py`).
*   **`mdp_files/`** (Рекомендуется): Папка для хранения шаблонных `.mdp` файлов.
*   **`env.yaml`**: Файл для создания окружения Conda.
*   **`README.md`**: Этот файл.

**ВАЖНО! Чтобы все файлы назывались в формате peptX и ligX, например, pept1 и lig1 (но не обязательно номера, например pept_something и lig_something)**

## Запуск симуляций (Шаг 6)

После того как скрипт сборки системы (`step5_...`) успешно отработает и создаст папки `*_Sim` (и они будут перемещены в `simulations/`), необходимо запустить скрипт симуляций (например, `run_simulation_pipeline.sh v16`, который мы обсуждали).

Перед запуском `run_simulation_pipeline.sh` убедитесь, что:

1.  Папки с собранными системами находятся в директории `simulations/`.
2.  Шаблонные MDP файлы (`step6.*.mdp`) находятся в корневой директории проекта.
3.  Вы находитесь в активированном окружении Conda (`gromacs_md_env`).

```bash
bash run_simulation_pipeline.sh # Запуск минимизации, уравновешивания и продакшена
