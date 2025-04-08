# Повторяем grompp и genion
echo "--> Повторный запуск grompp для ионов..."
gmx grompp -f ions.mdp -c pept1_lig1_bcc_manual_solv.gro -p pept1.top -o pept1_lig1_bcc_manual_ions.tpr -maxwarn 5
if [ $? -ne 0 ]; then echo "Ошибка на этапе grompp для ионов"; exit 1; fi

echo "--> Повторный запуск genion... Выберите группу 'SOL' для замены на ионы, если будет предложено."
echo "SOL" | gmx genion -s pept1_lig1_bcc_manual_ions.tpr -o pept1_lig1_bcc_manual_solv_ions.gro -p pept1.top -pname NA -nname CL -neutral -conc 0.15
if [ $? -ne 0 ]; then echo "Ошибка на этапе genion"; exit 1; fi

# Очистка временных файлов
rm -f ions.mdp \#*# *~

echo "--------------------------------------------------------------------"
echo "Добавление ионов завершено (повторно)."
echo "Готовые файлы для минимизации:"
echo "  Координаты: pept1_lig1_bcc_manual_solv_ions.gro"
echo "  Топология:  pept1.top"
echo "--------------------------------------------------------------------"
