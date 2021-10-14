#' suppa_joinFiles
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
#'
#' suppa_joinFiles(bam_files[1:5], "test1")
#' suppa_joinFiles(bam_files[6:10], "test2")
suppa_joinFiles = function(input_files,
                           output_name,
                           run_TPM = TRUE,
                           PSI_todo = unlist(SPLICE_EVENTS),
                           bam_suffix = ".Aligned.+",
                           output_location = "suppa2_diff"){
  dir.create(output_location, showWarnings = FALSE, recursive = TRUE)
  if(all(grepl("bam$", input_files))){
    message("Finding suppa files for bam input...")
    input_files = sub(bam_suffix, ".salmon_quant", input_files)
  }else{
    stop("NYI: currently only bam_files from VACC RNAseq pipeline are accepted.")
  }
  message("Running joinFiles...")
  if(run_TPM){
    tpm_files = paste(file.path(input_files, "tpm.txt"), collapse = " ")
    out_file = file.path(output_location, paste0(output_name, ".tpm.tpm"))
    if(!file.exists(out_file)){
      cmd = paste0(SUPPA_PATH, " joinFiles -f tpm -i ", tpm_files, " -o ", sub(".tpm$", "", out_file))
      system(paste0("bash inst/extdata/wrap_python.sh", " ", cmd))
    }else{
      message("Skipping ", out_file, " already exists.")
    }
  }
  for(psi in PSI_todo){
    psi_files = dir(input_files, pattern = paste0(psi, ".psi$"), full.names = TRUE)
    out_file = file.path(output_location, paste0(output_name, ".", psi, ".psi"))
    if(!file.exists(out_file)){
      cmd = paste0(SUPPA_PATH, " joinFiles -f psi -i ", paste(psi_files, collapse = " "), " -o ", sub(".psi$", "", out_file))
      system(paste0("bash inst/extdata/wrap_python.sh", " ", cmd))
    }else{
      message("Skipping ", out_file, " already exists.")
    }
  }
}


#' Title
#'
#' @param wd
#' @param ref_location
#' @param PSI_todo
#'
#' @return
#' @export
#'
#' @examples
#'
#' wd = "suppa2_diff"
#' tpm_files = dir(wd, pattern = ".tpm$", full.names = TRUE)
#' PSI_todo = unlist(SPLICE_EVENTS)
#' psi = PSI_todo[[1]]
#' psi_files = dir(wd, pattern = paste0(psi, ".psi$"), full.names = TRUE)
#' ref_location = "~/indexes/HG38/SUPPA2"
#' suppa_diffSplice(wd, ref_location)
suppa_diffSplice = function(wd,
                            ref_location = "~/indexes/MM10/SUPPA2",
                            PSI_todo = unlist(SPLICE_EVENTS),
                            output_location = wd){
  for(psi in PSI_todo){
    tpm_files = dir(wd, pattern = ".tpm$", full.names = TRUE)
    psi_files = dir(wd, pattern = paste0(psi, ".psi$"), full.names = TRUE)
    if(psi != "isoform"){
      ref_file = dir(ref_location, pattern = paste0(psi, ".+ioe"), full.names = TRUE)
    }else{
      ref_file = dir(ref_location, pattern = ".+ioi", full.names = TRUE)
    }
    out_root = file.path(output_location, paste0("diffSplice_result_", psi))
    if(!(file.exists(paste0(out_root, ".dpsi")) & file.exists(paste0(out_root, ".psivec")))){
      cmd = paste0(SUPPA_PATH, " diffSplice",
                   " -p ", paste(psi_files, collapse = " "),
                   " -e ", paste(tpm_files, collapse = " "),
                   " -i ", ref_file,
                   " -o ", out_root,
                   " -m empirical -a 1000 -l .05 -c")
      system(paste0("bash inst/extdata/wrap_python.sh", " ", cmd))
    }else{
      message("Skipping ", out_root, " already exists.")
    }

  }
}

#' Title
#'
#' @param wd
#' @param PSI_todo
#' @param output_location
#'
#' @return
#' @export
#'
#' @examples
#' wd = "suppa2_diff"
#' PSI_todo = unlist(SPLICE_EVENTS),
#' output_location = wd
#' psi = PSI_todo[1]
suppa_clusterEvents = function(wd,
                               PSI_todo = unlist(SPLICE_EVENTS),
                               output_location = wd){
  # echo $loc
  # dpsi=${diff_dir}/diffSplice_result_${loc}.dpsi
  # psivec=${diff_dir}/diffSplice_result_${loc}.psivec
  # sed -i 's/nan/1.0/g' ${dpsi}

  a = read.table("suppa2_diff/diffSplice_result_A3.dpsi")
  b = read.table("suppa2_diff/diffSplice_result_A3.psivec")
  head(a)
  library(ggplot2)
  library(data.table)
  a$is_nan = is.nan(a$test1.test2_dPSI)
  is_nan_ids = rownames(a[a$is_nan,])
  good_ids = setdiff(rownames(a), is_nan_ids)
  length(is_nan_ids)
  nrow(a)
  a[a$is_nan,]$test1.test2_dPSI = 0
  ggplot(a[good_ids,], aes(x = test1.test2_dPSI, y = test1.test2_p.val, color = is_nan)) +
    geom_point()

  head(b)
  b[is_nan_ids,]
  ggplot(b[is_nan_ids,], aes(x = test1_3, y = test1_4)) +
    geom_point(alpha = .03)

  ggplot(b[good_ids,], aes(x = test1_3, y = test1_4)) +
    geom_point(alpha = .03)
  bdt= as.data.table(b)
  nbins = 50
  bdt[, xbin := round(test1_3 *nbins)]
  bdt[, ybin := round(test1_4 *nbins)]

  bdt = bdt[, .N, .(xbin, ybin)]
  bdt[order(N, decreasing = TRUE)][1:10]
  bdt[is.nan(xbin)][order(N, decreasing = TRUE)][1:10]
  bdt[is.nan(ybin)][order(N, decreasing = TRUE)][1:10]



  # sed -i 's/nan/0.0/g' ${psivec}
  # ls -lha $dpsi $psivec
  # python ~/lab_bin/suppa.py clusterEvents --dpsi ${dpsi} --psivec ${psivec} --sig-threshold 0.1 --eps 0.05 --min-pts 20 --groups 1-3,4-6,7-9 -o clusterEvents_${loc}
  for(psi in PSI_todo){
    dpsi_file = file.path(wd, paste0("diffSplice_result_", psi, ".dpsi"))
    psivec_file = file.path(wd, paste0("diffSplice_result_", psi, ".psivec"))

  out_root = file.path(output_location, paste0("clusterEvents_result_", psi))
  cmd = paste0(SUPPA_PATH, " clusterEvents",
               " --dpsi ", dpsi_file,
               " --psivec ", psivec_file,
               " --sig-threshold 0.1 --eps 0.05 --min-pts 20 ",
               " --groups ", )
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
