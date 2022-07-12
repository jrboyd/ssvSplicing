

strand_cols = c("-" = "blue", "+" = "red", "*" = "gray")
strand_cols = c("-" = "red", "+" = "black")
# #gray bars in plots
# representative_transcript = "ENST00000331340.8"
# #plots some extra green circles, implemented to look at anti sense gene that looked interesting
# additional_transcript = "ENST00000644669.1"

#' filter_readable
#'
#' @param f
#'
#' @return
#' @export
#'
#' @examples
filter_readable= function(f){
  f = f[file.access(f, mode = 4) == 0]
  f
}

#' load_gtf
#'
#' @param gtf_file
#' @param gene_of_interest
#'
#' @return
#' @export
#' @importFrom rtracklayer import.gff
#'
#' @examples
#' gtf_file = "/slipstream/home/joeboyd/indexes/DM6/GTF/dm6.ensGene.gtf"
#' goi = "FBgn0029137"
#' ref_info = load_gtf(gtf_file, goi)
#' #or if you already have exon annotatation loaded, this is much faster
#' ref_info.faster = load_gtf(ref_info$ex_gr, goi)
load_gtf = function(gtf_file, gene_of_interest, expansion = 1.3){
  if(is(gtf_file, "GRanges")){
    ex_gr = gtf_file
  }else{
    stopifnot(file.exists(gtf_file))
    ex_gr = rtracklayer::import.gff(gtf_file, feature.type = "exon", format = "gtf")
  }

  ik_gr = subset(ex_gr, gene_name == gene_of_interest)
  ik_dt = as.data.table(ik_gr)
  tx_widths = ik_dt[, .(total_width = sum(end - start)), .(transcript_id)]
  tx_widths = tx_widths[order(total_width)]
  ik_dt$transcript_id = factor(ik_dt$transcript_id, levels = tx_widths$transcript_id)

  view_gr = range(ik_gr)
  view_gr = resize(view_gr, expansion*width(view_gr), fix = "center")
  return(list(ex_gr = ex_gr, view_gr = view_gr, goi_exon_dt = ik_dt))
}



#' load_splicing_from_bam_files
#'
#' @param bam_files
#'
#' @return
#' @export
#'
#' @examples
#' wd = "/slipstream/home/dbgap/data/alignment_RNA-Seq/"
#' sj_files = find_SJ.out.tab_files(wd)[1:10]
#'
#' features = load_gtf("~/gencode.v36.annotation.gtf", gene_of_interest = "IKZF1")
#'
#' load_splicing_from_SJ.out.tab_files(sj_files, features$view_gr)
load_splicing_from_bam_files = function(bam_files, view_gr, trim_extension = ".Aligned.sortedByCoord.out.bam"){
  qdt = data.table(file = bam_files)
  qdt[, sample := sub(trim_extension, "", basename(file))]
  dt_sp.first = ssvFetchBam(qdt, view_gr, splice_strategy = "splice_count", return_data.table = TRUE, flag = Rsamtools::scanBamFlag(isFirstMateRead = TRUE))
  dt_sp.second = ssvFetchBam(qdt, view_gr, splice_strategy = "splice_count", return_data.table = TRUE, flag = Rsamtools::scanBamFlag(isFirstMateRead = FALSE))
  flip_strand = c("+" = "-", "-" = "+")
  dt_sp.first$strand = flip_strand[dt_sp.first$strand]
  dt_sp = rbind(dt_sp.first, dt_sp.second)
  dt_sp = dt_sp[, .(N = sum(N)), .(which_label, seqnames, start, end, strand, sample)]
  dt_sp[, id := paste(seqnames, start, end, strand)]
  dt_sp
}

load_sj = function(f, view_gr, file_ext = ".SJ.out.tab"){
  cn_sj = c("seqnames", "start", "end", "strand_indicator", "intron_motif", "annotation", "number_unique", "number_multi", "max_overhang")
  dt_sj = fread(f, col.names = cn_sj)

  decode_strand = c("*", "+", "-")
  dt_sj[, strand := decode_strand[strand_indicator+1]]
  dt_sj$strand_indicator = NULL

  decode_motif = c("non-canon",  "GT/AG", "CT/AC", "GC/AG", "CT/GC", "AT/AC", "GT/AT")
  dt_sj[, intron_motif := decode_motif[intron_motif+1]]

  decode_anno = c("unannotated", "annotated")
  dt_sj[, annotation := decode_anno[annotation+1]]

  gr_sj = GRanges(dt_sj)
  dt_sj = suppressWarnings({
    as.data.table(subsetByOverlaps(gr_sj, view_gr, ignore.strand = TRUE))
  })
  dt_sj$sample = sub(file_ext, "", basename(f))
  dt_sj
}

#' find_SJ.out.tab_files
#'
#' @param wd
#'
#' @return
#' @export
#'
#' @examples
#' wd = "/slipstream/home/dbgap/data/alignment_RNA-Seq/"
#' sj_files = find_SJ.out.tab_files(wd)[1:10]
#' sj_files
find_SJ.out.tab_files = function(wd){
  stopifnot(dir.exists(wd))
  sj_files = dir(wd, pattern = "SJ.out.tab$", full.names = TRUE)
  sj_files
}

#' find_bam_files
#'
#' @param wd
#'
#' @return
#' @export
#'
#' @examples
#' wd = "/slipstream/home/dbgap/data/alignment_RNA-Seq/"
#' bam_files = find_bam_files(wd)[1:10]
#' bam_files
find_bam_files = function(wd){
  stopifnot(dir.exists(wd))
  bam_files = dir(wd, pattern = "bam$", full.names = TRUE)
  bam_files
}


#' load_splicing_from_SJ.out.tab_files
#'
#' @param sj_files
#' @param view_gr
#'
#' @return
#' @export
#'
#' @examples
#' wd = "/slipstream/home/dbgap/data/alignment_RNA-Seq/"
#' sj_files = find_SJ.out.tab_files(wd)[1:10]
#'
#' features = load_gtf("~/gencode.v36.annotation.gtf", gene_of_interest = "IKZF1")
#'
#' load_splicing_from_SJ.out.tab_files(sj_files, features$view_gr)
load_splicing_from_SJ.out.tab_files = function(sj_files, view_gr){
  res = pbmcapply::pbmclapply(sj_files, load_sj, view_gr = view_gr)
  dt_sj = rbindlist(res)
  dt_sj[, id := paste(seqnames, start, end, strand)]
  dt_sj
}
