get_wrap_script = function(){
  system.file("extdata/wrap_python.sh", package = "ssvSplicing")
}

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
    message("tpm", "...")
    tpm_files = paste(file.path(input_files, "tpm.txt"), collapse = " ")
    out_file = file.path(output_location, paste0(output_name, ".tpm.tpm"))
    if(!file.exists(out_file)){
      cmd = paste0(SUPPA_PATH, " joinFiles -f tpm -i ", tpm_files, " -o ", sub(".tpm$", "", out_file))
      system(paste0("bash ", get_wrap_script(), " ", cmd))
    }else{
      message("Skipping ", out_file, " already exists.")
    }
  }
  for(psi in PSI_todo){
    message(psi, "...")
    psi_files = dir(input_files, pattern = paste0(psi, ".psi$"), full.names = TRUE)
    out_file = file.path(output_location, paste0(output_name, ".", psi, ".psi"))
    if(!file.exists(out_file)){
      cmd = paste0(SUPPA_PATH, " joinFiles -f psi -i ", paste(psi_files, collapse = " "), " -o ", sub(".psi$", "", out_file))
      system(paste0("bash ", get_wrap_script(), " ", cmd))
    }else{
      message("Skipping joinFiles for ", out_file, ", output already exists.")
    }
  }
}

get_tpm_files = function(wd){
  dir(wd, pattern = ".tpm$", full.names = TRUE)
}

get_psi_files = function(wd, psi){
  dir(wd, pattern = paste0(psi, ".psi$"), full.names = TRUE)
}

#' suppa_diffSplice
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
#' ref_location = "~/indexes/HG38/SUPPA2"
#' suppa_diffSplice(wd, ref_location)
suppa_diffSplice = function(wd,
                            ref_location = "~/indexes/MM10/SUPPA2",
                            PSI_todo = unlist(SPLICE_EVENTS),
                            output_location = wd){
  for(psi in PSI_todo){
    message(psi, "...")
    tpm_files = get_tpm_files(wd)
    psi_files = get_psi_files(wd, psi)
    if(psi != "isoform"){
      ref_file = dir(ref_location, pattern = paste0(psi, ".+ioe"), full.names = TRUE)
    }else{
      ref_file = dir(ref_location, pattern = ".+ioi", full.names = TRUE)
    }
    out_root = file.path(output_location, paste0("diffSplice_result_", psi))
    if(!(file.exists(paste0(out_root, ".psivec")))){
      cmd = paste0(SUPPA_PATH, " diffSplice",
                   " -p ", paste(psi_files, collapse = " "),
                   " -e ", paste(tpm_files, collapse = " "),
                   " -i ", ref_file,
                   " -o ", out_root,
                   " -m empirical -a 1000 -l .05 -c")
      system(paste0("bash ", get_wrap_script(), " ", cmd))
    }else{
      message("Skipping diffSplice for ", out_root, ", output already exists.")
    }

  }
}

#' suppa_clusterEvents
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
#' suppa_clusterEvents(wd)
suppa_clusterEvents = function(wd,
                               PSI_todo = unlist(SPLICE_EVENTS),
                               output_location = wd){
  for(psi in PSI_todo){
    message(psi, "...")
    psivec_file = file.path(wd, paste0("diffSplice_result_", psi, ".psivec"))
    stopifnot(file.exists(psivec_file))
    psivec_df.raw = read.table(psivec_file)


    all_dpsi_files = dir(wd, pattern = paste0("diffSplice_result_", psi, "\\.dpsi"), full.names = TRUE)
    stopifnot(all(file.exists(all_dpsi_files)))
    for(dpsi_file in all_dpsi_files){
      dpsi_df.raw = read.table(dpsi_file, header = TRUE, row.names = 1)
      pair_str = sub("_dPSI", "", colnames(dpsi_df.raw)[1])
      message("  ", pair_str)
      pair_a = strsplit(pair_str, "\\.")[[1]][1]
      pair_b = strsplit(pair_str, "\\.")[[1]][2]

      psivec_df.dpsi_matched = psivec_df.raw[, grepl(pair_a, colnames(psivec_df.raw)) | grepl(pair_b, colnames(psivec_df.raw))]

      good_ids = rownames(dpsi_df.raw[!is.nan(dpsi_df.raw[[1]]),])
      psivec_df.dpsi_matched[good_ids,]

      dpsi_file.no_nan = paste0(sub(".dpsi", "", dpsi_file), ".", pair_a, ".", pair_b, ".no_nan.dpsi")
      psivec_file.no_nan = paste0(sub(".psivec", "", psivec_file), ".", pair_a, ".", pair_b, ".no_nan.psivec")


      write.table(dpsi_df.raw[good_ids,], dpsi_file.no_nan, sep = "\t", quote = FALSE)
      write.table(psivec_df.dpsi_matched[good_ids,], psivec_file.no_nan, sep = "\t", quote = FALSE)

      grps = sapply(strsplit(colnames(psivec_df.dpsi_matched), "_"), function(x){
        paste(x[-length(x)], collapse = "_")
      })

      ranges = sapply(unique(grps), function(x)range(which(x == grps)))
      groups_str = paste(paste(ranges[1,], ranges[2,], sep = "-"), collapse = ",")
      out_root = file.path(normalizePath(output_location), paste0("clusterEvents_result_", psi, ".", pair_a, ".", pair_b))
      if(!file.exists(paste0(out_root, ".clustvec"))){
        owd = getwd()
        setwd(dirname(dpsi_file.no_nan))
        cmd = paste0(SUPPA_PATH, " clusterEvents",
                     " --dpsi ", basename(dpsi_file.no_nan),
                     " --psivec ", basename(psivec_file.no_nan),
                     " --sig-threshold 0.1 --eps 0.05 --min-pts 20 ",
                     " --groups ", groups_str,
                     " -o ", out_root)
        message(cmd)
        system(paste0("bash ", get_wrap_script(), " ", cmd))
        setwd(owd)
      }else{
        message("Skipping clusterEvents for ", out_root, ", output already exists.")
      }
    }
  }
}
