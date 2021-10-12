diff_dir=$1
if [ -z $diff_dir ]; then diff_dir=$(pwd); fi

for loc in isoform A3 A5 AF AL MX RI SE; do
  echo $loc
  dpsi=${diff_dir}/diffSplice_result_${loc}.dpsi
  psivec=${diff_dir}/diffSplice_result_${loc}.psivec
  sed -i 's/nan/1.0/g' ${dpsi}
  sed -i 's/nan/0.0/g' ${psivec}
  ls -lha $dpsi $psivec
  python ~/lab_bin/suppa.py clusterEvents --dpsi ${dpsi} --psivec ${psivec} --sig-threshold 0.1 --eps 0.05 --min-pts 20 --groups 1-3,4-6,7-9 -o clusterEvents_${loc}
done
