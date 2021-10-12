


SPLICE_EVENTS = list(
  "SkippingExon" = "SE",
  "Alternative5Prime" = "A5",
  "Alternative3Prime" = "A3",
  "MutuallyExclusiveExon" = "MX",
  "RetainedIntron" = "RI",
  "AlternativeFirstExon" = "AF",
  "AlternativeLastExon" = "AL",
  "Isoform" = 'isoform'
)

SPLICE_EVENTS.DECODE = unlist(SPLICE_EVENTS)
SPLICE_EVENTS.REVERSE = names(SPLICE_EVENTS)
names(SPLICE_EVENTS.REVERSE) = unlist(SPLICE_EVENTS)

SUPPA_PATH = "/slipstream/home/joeboyd/anaconda2/envs/suppa2_env/bin/suppa.py"


#' Title
#'
#' @param input_files
#' @param output_name
#' @param run_TPM
#' @param PSI_todo
#' @param bam_suffix
#' @param output_location
#'
#' @return
#' @export
#'
#' @examples
#' wd = "/slipstream/home/dbgap/data/alignment_RNA-Seq/"
#' bam_files = find_bam_files(wd)
#' input_files = bam_files[1:5]
#' suppa_joinFiles(input_files, "test")
suppa_joinFiles = function(input_files, output_name, run_TPM = TRUE, PSI_todo = unlist(SPLICE_EVENTS), bam_suffix = ".Aligned.+", output_location = "suppa2_diff"){
  dir.create(output_location, showWarnings = FALSE, recursive = TRUE)
  if(all(grepl("bam$", input_files))){
    message("Finding suppa files for bam input...")
    input_files = sub(bam_suffix, ".salmon_quant", input_files)
  }else{
    stop("NYI: currently only bam_files from VACC RNAseq pipeline are accepted.")
  }

  if(run_TPM){
    tpm_files = paste(file.path(input_files, "tpm.txt"), collapse = " ")
    cmd = paste0(SUPPA_PATH, " joinFiles -f tpm -i ", tpm_files, " -o ", file.path(output_location, paste0(output_name, ".tpm")))
    system(paste0("bash inst/extdata/env_eval.sh; bash inst/extdata/wrap_python.sh", " ", cmd))
  }
  for(psi in PSI_todo){
    psi_files = dir(input_files, pattern = paste0(psi, ".psi$"), full.names = TRUE)
    cmd = paste0(SUPPA_PATH, " joinFiles -f psi -i ", paste(psi_files, collapse = " "), " -o ", file.path(output_location, paste0(output_name, ".", psi)))
    system(paste0("bash inst/extdata/env_eval.sh; bash inst/extdata/wrap_python.sh", " ", cmd))


    cmd = paste0(SUPPA_PATH, " joinFiles -f psi -i ", paste(psi_files, collapse = " "), " -o ", file.path(output_location, paste0(output_name, ".", psi)))

  }



}

# mkdir -p suppa2_diff
#
# echo tpm
# python ~/lab_bin/suppa.py joinFiles -f tpm -i $(echo *18C*quant/tpm.txt) -o suppa2_diff/tpm_18C
# python ~/lab_bin/suppa.py joinFiles -f tpm -i $(echo *25C*quant/tpm.txt) -o suppa2_diff/tpm_25C
# python ~/lab_bin/suppa.py joinFiles -f tpm -i $(echo *30C*quant/tpm.txt) -o suppa2_diff/tpm_30C
#
#
# for loc in A3 A5 AF AL MX RI SE isoform; do
# echo $loc
# python ~/lab_bin/suppa.py joinFiles -f psi -i $(echo *30C*quant/*${loc}.psi) -o suppa2_diff/${loc}_30C;
# python ~/lab_bin/suppa.py joinFiles -f psi -i $(echo *25C*quant/*${loc}.psi) -o suppa2_diff/${loc}_25C;
# python ~/lab_bin/suppa.py joinFiles -f psi -i $(echo *18C*quant/*${loc}.psi) -o suppa2_diff/${loc}_18C;
# done
