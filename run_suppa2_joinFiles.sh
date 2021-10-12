mkdir -p suppa2_diff

echo tpm
python ~/lab_bin/suppa.py joinFiles -f tpm -i $(echo *18C*quant/tpm.txt) -o suppa2_diff/tpm_18C
python ~/lab_bin/suppa.py joinFiles -f tpm -i $(echo *25C*quant/tpm.txt) -o suppa2_diff/tpm_25C
python ~/lab_bin/suppa.py joinFiles -f tpm -i $(echo *30C*quant/tpm.txt) -o suppa2_diff/tpm_30C


for loc in A3 A5 AF AL MX RI SE isoform; do 
  echo $loc
  python ~/lab_bin/suppa.py joinFiles -f psi -i $(echo *30C*quant/*${loc}.psi) -o suppa2_diff/${loc}_30C; 
  python ~/lab_bin/suppa.py joinFiles -f psi -i $(echo *25C*quant/*${loc}.psi) -o suppa2_diff/${loc}_25C; 
  python ~/lab_bin/suppa.py joinFiles -f psi -i $(echo *18C*quant/*${loc}.psi) -o suppa2_diff/${loc}_18C; 
done
