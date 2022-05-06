#/users/j/r/jrboyd/lab_shared/indexes/DM6/SUPPA2/suppa2.dm6.ensGene.gtf..ioi
#-rw-r--r-- 1 sfrietze pi-sfrietze 393K Jan  5 15:17 /users/j/r/jrboyd/lab_shared/indexes/DM6/SUPPA2/suppa2.dm6.ensGene.gtf._A3_strict.ioe
#-rw-r--r-- 1 sfrietze pi-sfrietze 398K Jan  5 15:17 /users/j/r/jrboyd/lab_shared/indexes/DM6/SUPPA2/suppa2.dm6.ensGene.gtf._A5_strict.ioe
#-rw-r--r-- 1 sfrietze pi-sfrietze 593K Jan  5 15:17 /users/j/r/jrboyd/lab_shared/indexes/DM6/SUPPA2/suppa2.dm6.ensGene.gtf._AF_strict.ioe
#-rw-r--r-- 1 sfrietze pi-sfrietze 157K Jan  5 15:17 /users/j/r/jrboyd/lab_shared/indexes/DM6/SUPPA2/suppa2.dm6.ensGene.gtf._AL_strict.ioe
#-rw-r--r-- 1 sfrietze pi-sfrietze 392K Jan  5 15:17 /users/j/r/jrboyd/lab_shared/indexes/DM6/SUPPA2/suppa2.dm6.ensGene.gtf._MX_strict.ioe
#-rw-r--r-- 1 sfrietze pi-sfrietze 231K Jan  5 15:17 /users/j/r/jrboyd/lab_shared/indexes/DM6/SUPPA2/suppa2.dm6.ensGene.gtf._RI_strict.ioe
#-rw-r--r-- 1 sfrietze pi-sfrietze 422K Jan  5 15:17 /users/j/r/jrboyd/lab_shared/indexes/DM6/SUPPA2/suppa2.dm6.ensGene.gtf._SE_strict.ioe

gtf_base=~/lab_shared/indexes/DM6/SUPPA2/suppa2.dm6.ensGene.gtf

diff_dir=$1
if [ -z $diff_dir ]; then diff_dir=$(pwd); fi

for loc in isoform A3 A5 AF AL MX RI SE; do
  ref=${gtf_base}._${loc}_strict.ioe
  if [ ${loc} = isoform ]; then 
    ref=${gtf_base}..ioi
  fi
  if [ ! -f $ref ]; then echo $ref reference file not found. quit.; exit 1; fi
  tpm_files="$diff_dir/tpm_18C.tpm $diff_dir/tpm_25C.tpm $diff_dir/tpm_30C.tpm"
  psi_files="$diff_dir/${loc}_18C.psi $diff_dir/${loc}_25C.psi $diff_dir/${loc}_30C.psi"
  ls -lha $tpm_files $psi_files $ref 
  
  python ~/lab_bin/suppa.py diffSplice -p $psi_files -e $tpm_files -i ${ref} -o diffSplice_result_${loc} -m empirical -a 1000 -l .05 -c
  mv diffSplice_result_${loc}.dpsi $diff_dir/diffSplice_result_${loc}.dpsi
  mv diffSplice_result_${loc}.psivec $diff_dir/diffSplice_result_${loc}.psivec
done
