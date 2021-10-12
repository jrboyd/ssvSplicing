

strand_cols = c("-" = "blue", "+" = "red", "*" = "gray")
strand_cols = c("-" = "red", "+" = "black")
# #gray bars in plots
# representative_transcript = "ENST00000331340.8"
# #plots some extra green circles, implemented to look at anti sense gene that looked interesting
# additional_transcript = "ENST00000644669.1"

#' Title
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

#' Title
#'
#' @param gtf_file
#' @param gene_of_interest
#'
#' @return
#' @export
#' @importFrom rtracklayer import.gff
#'
#' @examples
load_gtf = function(gtf_file, gene_of_interest){
  ex_gr = rtracklayer::import.gff("~/gencode.v35.annotation.gtf.gz", feature.type = "exon", format = "gtf")
  ik_gr = subset(ex_gr, gene_name == gene_of_interest)
  ik_dt = as.data.table(ik_gr)
  tx_widths = ik_dt[, .(total_width = sum(end - start)), .(transcript_id)]
  tx_widths = tx_widths[order(total_width)]
  ik_dt$transcript_id = factor(ik_dt$transcript_id, levels = tx_widths$transcript_id)

  view_gr = range(ik_gr)
  view_gr = resize(view_gr, 1.3*width(view_gr), fix = "center")
  return(list(ex_gr = ex_gr, view_gr = view_gr, goi_exon_dt = ik_dt))
}



#' Title
#'
#' @param bam_files
#'
#' @return
#' @export
#'
#' @examples
load_splicing_from_bam_files = function(bam_files, view_gr, trim_extension = ".Aligned.sortedByCoord.out.bam"){
  qdt = data.table(file = bam_files)
  qdt[, sample := sub(trim_extension, "", basename(file))]
  dt_sp.first = ssvFetchBam(qdt, view_gr, splice_strategy = "splice_count", return_data.table = TRUE, flag = Rsamtools::scanBamFlag(isFirstMateRead = TRUE))
  dt_sp.second = ssvFetchBam(qdt, view_gr, splice_strategy = "splice_count", return_data.table = TRUE, flag = Rsamtools::scanBamFlag(isFirstMateRead = FALSE))
  flip_strand = c("+" = "-", "-" = "+")
  dt_sp.first$strand = flip_strand[dt_sp.first$strand]
  dt_sp = rbind(dt_sp.first, dt_sp.second)
  dt_sp = dt_sp[, .(N = sum(N)), .(which_label, seqnames, start, end, strand, sample)]
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
  dt_sj = as.data.table(subsetByOverlaps(gr_sj, view_gr, ignore.strand = TRUE))
  dt_sj$sample = sub(file_ext, "", basename(f))
  dt_sj
}

#' Title
#'
#' @param wd
#'
#' @return
#' @export
#'
#' @examples
find_SJ.out.tab_files = function(wd){
  stopifnot(dir.exists(wd))
  sj_files = dir(wd, pattern = "SJ.out.tab$", full.names = TRUE)
  sj_files
}

#' Title
#'
#' @param wd
#'
#' @return
#' @export
#'
#' @examples
find_bam_files = function(wd){
  stopifnot(dir.exists(wd))
  bam_files = dir(wd, pattern = "bam$", full.names = TRUE)
  bam_files
}


#' Title
#'
#' @param sj_files
#' @param view_gr
#'
#' @return
#' @export
#'
#' @examples
load_splicing_from_SJ.out.tab_files = function(sj_files, view_gr){
  # sj_files = dir("monaco_and_CCLE/", pattern = ".SJ.out.tab$", full.names = TRUE)
  # sj_files = filter_readable(sj_files)
  dt_sj = rbindlist(pbmcapply::pbmclapply(sj_files, load_sj, view_gr = view_gr))
  dt_sj
}

#
# qdt = data.table(file = bam_files)
# qdt[, sample := sub("_rep1.Aligned.sortedByCoord.out.bam", "", basename(file))]
#
#
# f = sj_files[1]
#
#
#
#
# options(mc.cores = 20)
#
# sj_file = paste0("ccle_ik_splice_js.", length(sj_files),".csv")
# length(sj_files)
# sj_lev = sub(".SJ.out.tab", "", basename(sj_files))
#
# if(file.exists(sj_file)){
#   dt_sj = fread(sj_file)
# }else{
#   dt_sj = rbindlist(pbmcapply::pbmclapply(sj_files, load_sj, view_gr = view_gr))
#   fwrite(dt_sj, sj_file)
# }
#
# length(unique(dt_sj$sample))

# ssvRecipes::ssvFetchBamPE.RNA
# dt_sp.first = ssvFetchBam(qdt, view_gr, splice_strategy = "splice_count", return_data.table = TRUE, flag = Rsamtools::scanBamFlag(isFirstMateRead = TRUE))
# dt_sp.second = ssvFetchBam(qdt, view_gr, splice_strategy = "splice_count", return_data.table = TRUE, flag = Rsamtools::scanBamFlag(isFirstMateRead = FALSE))
# flip_strand = c("+" = "-", "-" = "+")
# dt_sp.first$strand = flip_strand[dt_sp.first$strand]
# dt_sp = rbind(dt_sp.first, dt_sp.second)
# dt_sp = dt_sp[, .(N = sum(N)), .(which_label, seqnames, start, end, strand, sample)]

# dt = merge(by = c("seqnames", "start", "end", "strand", "sample"), all = TRUE,
#            dt_sj[, .(seqnames, start, end, strand, sample, number_unique, number_multi)],
#            dt_sp[, .(seqnames, start, end, strand, sample, N)]
# )
#
# dt = copy(dt_sj)
#
# dt[is.na(number_unique), number_unique := 0]
# dt[is.na(number_multi), number_multi := 0]
#
#
#
#
# dt[, number_total := number_unique + number_multi]
# plot(log2(cbind(dt$number_total, dt$number_unique)+16), col = rgb(0,0,0,.5), pch = 16, xlab = "log2 SJ.out.tab total", ylab = "log2 SJ.out.tab unique")
# dt[number_total > 20]
#
# dt = dt[strand %in% c("+", '-')]
#
# dtw = melt(dt[, .(seqnames, start, end, strand, sample, number_unique, number_total)], measure.vars = c("number_unique", "number_total"))
#
# dtw.filled = dcast(dtw, seqnames+start+end+strand+variable~sample, value.var = "value", fill = 0)
# dtw = melt(dtw, id.vars = c("seqnames", "start", "end", "strand", "variable", "sample"))
# dtw$variable.1 = NULL
#
# dtw[, seq_sample := sample]
# dtw[, sample := sub("_rep1", "", seq_sample)]
# dtw = dtw[, .(value = sum(value)), .(seqnames, start, end, strand, sample, variable)]
#
# # saveRDS(dtw, "dtw.RDS")
#
# anno_dt = as.data.table(ex_basic)[, .(start, end)]
# anno_dt$m = ""
# anno_rng = dtw[, .(ymin = 0, ymax = max(value)), .(sample)]
# anno_rng$m = ""
# anno_dt = merge(anno_dt, anno_rng, by = "m", allow.cartesian = TRUE)
#
#
#
# plot_splice_view_full = function(sel_sample, name){
#   ggplot() +
#     geom_rect(data = anno_dt[sample %in% sel_sample], aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax), color = "gray", fill = "gray", size = .2) +
#     geom_arch(GRanges(dtw[value > 10 & variable == "number_total" & sample %in% sel_sample]), aes(height = value, color = strand)) +
#     facet_wrap(sample~., scales = "free_y", ncol = 3) +
#     scale_color_manual(values = strand_cols) +
#     scale_x_continuous(labels = function(x)paste(x/1e3)) +
#     labs(x = "kb") +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
#     labs(title = name) +
#     coord_cartesian(xlim = c(start(view_gr), end(view_gr)))
# }
#
# # donor_samples = unique(dtw[grepl("donor", sample)]$sample)
# # xenograft_samples = unique(dtw[!grepl("donor", sample)]$sample)
#
# # plot_splice_view_full(donor_samples, "donor samples")
# # plot_splice_view_full(xenograft_samples, "xenograft samples")
#
# plot_splice_view_full(unique(anno_dt$sample)[1:9], "monaco and CCLE samples")
#
# plot_splice_heatmap = function(sel_sample, pdf_name, name, min_test_value = 20, win_size = 1, nclust = 4){
#   if(!grepl(".pdf$", pdf_name)){
#     stop("file must be .pdf")
#   }
#   csv_name = sub(".pdf$", ".csv", pdf_name)
#
#   dtw[, id := paste(seqnames, start, end, strand)]
#   dtw[, fraction := value / sum(value), .(sample)]
#
#   sel_id = dtw[, .(test_val = max(value)), .(id)][order(test_val)][test_val >= min_test_value]$id
#
#   sel_dtw = dtw[variable == "number_total" & id %in% sel_id]
#   sel_dtw.filled = dcast(sel_dtw, seqnames+start+end+strand+variable+id~sample, value.var = "value", fill = 0)
#   sel_dtw = melt(sel_dtw.filled, id.vars = c("seqnames", "start", "end", "strand", "variable", "id"), variable.name = "sample")
#   sel_dtw[, fraction := value / sum(value), .(sample)]
#   sel_dtw$id = factor(sel_dtw$id, levels = sel_id)
#
#   sel_dtw = sel_dtw[sample %in% sel_sample]
#
#   set.seed(0)
#   id_lev = levels(ssvSignalClustering(sel_dtw, row_ = "id", column_ = "sample", fill_ = "fraction", facet_ = "", dcast_fill = 0, max_cols = Inf, nclust = 5)$id)
#
#   set.seed(1)
#   sp_clust_dt = ssvSignalClustering(sel_dtw,
#                                     nclust = nclust,
#                                     row_ = "sample", column_ = "id",
#                                     fill_ = "fraction", facet_ = "",
#                                     dcast_fill = 0,
#                                     max_cols = Inf,
#
#                                     max_rows = Inf)
#   sp_clust_dt$id = factor(sp_clust_dt$id, levels = id_lev)
#
#   n_samples = length(unique(dtw$sample))
#   toplot_id = as.character(sel_dtw[, sum(fraction > 0) >= .3*n_samples, id][V1 == TRUE]$id)
#   p_sp_heat = ssvSignalHeatmap(sp_clust_dt[id %in% toplot_id], row_ = "sample", column_ = "id", fill_ = "fraction", facet_ = "", dcast_fill = 0, max_cols = Inf, nclust = 3)
#
#   sp_agg_dt = sp_clust_dt[, .(fraction = mean(fraction), value = mean(value)), .(seqnames, start, end, strand, id, cluster_id)]
#
#   anno_dt_clust = as.data.table(ex_basic)[, .(start, end)]
#   anno_dt_clust$m = ""
#   anno_rng = sp_agg_dt[, .(ymin = 0, ymax = max(fraction)), .(cluster_id)]
#   anno_rng$m = ""
#   anno_dt_clust = merge(anno_dt_clust, anno_rng, by = "m", allow.cartesian = TRUE)
#
#   p_sp_arch = ggplot() +
#     geom_rect(data = anno_dt_clust, aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax), color = "gray", fill = "gray", size = .2) +
#     geom_arch(GRanges(sp_agg_dt), aes(height = fraction, color = strand)) +
#     facet_grid(cluster_id~., scales = "free_y") +
#     scale_color_manual(values = strand_cols) +
#     labs(y = "fraction") +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust= 1)) +
#     coord_cartesian(xlim = c(start(view_gr), end(view_gr)))
#   # coord_cartesian(xlim = c(50375000, 50400000)) +
#   # annotate("point", x = c(ex_add$start, ex_add$end), y = 0, col = "green") +
#   # annotate("segment", x = ex_add$start, xend = ex_add$end, y = 0, yend = 0, col = "green")
#
#   pg_arch = cowplot::plot_grid(p_sp_heat + labs(title = "name"), p_sp_arch)
#
#   ggsave(pdf_name, pg_arch, width = 8, height = 5)
#
#   sp_assign_dt = unique(sp_clust_dt[, .(cluster_id, sample)])[order(cluster_id)]
#   fwrite(sp_assign_dt, csv_name)
#
#   sel_prof_qdt = qdt
#   prof_gr = reduce(ik_gr)
#   prof_gr = range(prof_gr)
#   prof_gr = resize(prof_gr, 1.2*width(prof_gr), fix = "center")
#
#   prof_dt = ssvFetchBam(sel_prof_qdt, prof_gr, win_size = win_size, return_data.table = TRUE, fragLens = NA)
#   prof_dt.add = ssvFetchBam(sel_prof_qdt, prof_gr, win_size = win_size, return_data.table = TRUE, splice_strategy = "add", fragLens = NA)
#
#   prof_dt = prof_dt[, .(y = sum(y)), .(x, seqnames, start, end, id, strand, sample)]
#   prof_dt.add = prof_dt.add[, .(y = sum(y)), .(x, seqnames, start, end, id, strand, sample)]
#
#   sp_assign_dt = unique(sp_clust_dt[, .(cluster_id, sample)])
#
#   prof_dt = merge(sp_assign_dt, prof_dt, by = "sample")
#   prof_dt$sample = factor(prof_dt$sample, levels = levels(sp_assign_dt$sample))
#
#   prof_dt.add = merge(sp_assign_dt, prof_dt.add, by = "sample")
#   prof_dt.add$sample = factor(prof_dt.add$sample, levels = levels(sp_assign_dt$sample))
#
#
#   agg_prof_dt = prof_dt[, .(y = mean(y)), .(cluster_id, x, id, strand, start, end)]
#   agg_prof_dt.add = prof_dt.add[, .(y = mean(y)), .(cluster_id, x, id, strand, start, end)]
#
#   # ggplot(agg_prof_dt, aes(x = x, y = y, color = strand)) +
#   #   geom_path() +
#   #   facet_grid(cluster_id~id, scales = "free_x")
#
#   agg_prof_dt[, xgen := (start + end)/2]
#   agg_prof_dt.add[, xgen := (start + end)/2]
#
#   anno_dt_pile = as.data.table(ex_basic)[, .(start, end)]
#   anno_dt_pile$m = ""
#   anno_rng = agg_prof_dt.add[, .(ymin = 0, ymax = max(y)), .(cluster_id)]
#   anno_rng$m = ""
#   anno_dt_pile = merge(anno_dt_pile, anno_rng, by = "m", allow.cartesian = TRUE)
#
#
#   p_pileup = ggplot() +
#     geom_rect(data = anno_dt_pile, aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax), fill = "lightgray", color = "lightgray") +
#     geom_ribbon(data = agg_prof_dt.add, aes(x = xgen, ymin = 0, ymax = y, group = paste(id, strand)), alpha = .2, fill = "dodgerblue2") +
#     # geom_path(data = agg_prof_dt, aes(x = xgen, y = y, color = strand, group = paste(id, strand))) +
#     geom_ribbon(data = agg_prof_dt, aes(x = xgen, ymin = 0, ymax = y, group = paste(id, strand)), alpha = 1, fill = "dodgerblue2") +
#     geom_arch(GRanges(sp_agg_dt), aes(height = value), color = "black") +
#     facet_grid(cluster_id~., scales = "free") + theme(panel.background = element_blank()) +
#     labs(title = name) +
#     coord_cartesian(xlim = c(start(view_gr), end(view_gr)))
#
#   ggsave(sub(".pdf", ".pileup.pdf", pdf_name), p_pileup, width = 12, height = 10)
#
#   p_pileup.3b4 = p_pileup +
#     coord_cartesian(xlim = c(50365000, 50380000)) +
#     labs(subtitle = "exon3b and exon4 region")
#
#   ggsave(sub(".pdf", ".pileup_ex3b.pdf", pdf_name), p_pileup.3b4, width = 12, height = 10)
#
#   invisible(list(clust_dt = sp_clust_dt, agg_dt = sp_agg_dt, assign_dt = sp_assign_dt))
#
#
# }
#
# # plot_splice_heatmap(donor_samples, "donor_heatmap.pdf", "donors")
# # plot_splice_heatmap(xenograft_samples, "xenograft_heatmap.pdf", "xenografts")
# plot_splice_heatmap(sel_sample = unique(anno_dt$sample),
#                     pdf_name = "monaco_ccle_heatmap.pdf",
#                     name = "Monaco and CCLE",
#                     min_test_value = 10, nclust = 6)
#
