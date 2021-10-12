#' Title
#'
#' @param sel_sample
#' @param pdf_name
#' @param name
#' @param min_test_value
#' @param win_size
#'
#' @return
#' @export
#' @import cowplot
#' @examples
plot_splice_heatmap = function(dtw, sel_sample, pdf_name, name, min_test_value = 20, win_size = 1){
  if(!grepl(".pdf$", pdf_name)){
    stop("file must be .pdf")
  }
  csv_name = sub(".pdf$", ".csv", pdf_name)

  dtw[, id := paste(seqnames, start, end, strand)]
  dtw[, fraction := value / sum(value), .(sample)]

  sel_id = dtw[, .(test_val = max(value)), .(id)][order(test_val)][test_val >= min_test_value]$id

  sel_dtw = dtw[variable == "number_total" & id %in% sel_id]
  sel_dtw.filled = dcast(sel_dtw, seqnames+start+end+strand+variable+id~sample, value.var = "value", fill = 0)
  sel_dtw.filled[1:5, 1:10]
  sel_dtw = melt(sel_dtw.filled, id.vars = c("seqnames", "start", "end", "strand", "variable", "id"), variable.name = "sample")
  sel_dtw[, fraction := value / sum(value), .(sample)]
  sel_dtw$id = factor(sel_dtw$id, levels = sel_id)

  sel_dtw = sel_dtw[sample %in% sel_sample]

  set.seed(0)
  id_lev = levels(ssvSignalClustering(sel_dtw, row_ = "id", column_ = "sample", fill_ = "fraction", facet_ = "", dcast_fill = 0, max_cols = Inf, nclust = 5)$id)

  set.seed(1)
  sp_clust_dt = ssvSignalClustering(sel_dtw,
                                    nclust = 4,
                                    row_ = "sample", column_ = "id",
                                    fill_ = "fraction", facet_ = "",
                                    dcast_fill = 0,
                                    max_cols = Inf,

                                    max_rows = Inf)
  sp_clust_dt$id = factor(sp_clust_dt$id, levels = id_lev)

  n_samples = length(unique(dtw$sample))
  toplot_id = as.character(sel_dtw[, sum(fraction > 0) >= .3*n_samples, id][V1 == TRUE]$id)
  p_sp_heat = ssvSignalHeatmap(sp_clust_dt[id %in% toplot_id], row_ = "sample", column_ = "id", fill_ = "fraction", facet_ = "", dcast_fill = 0, max_cols = Inf, nclust = 3)

  sp_agg_dt = sp_clust_dt[, .(fraction = mean(fraction), value = mean(value)), .(seqnames, start, end, strand, id, cluster_id)]

  anno_dt_clust = as.data.table(ex_basic)[, .(start, end)]
  anno_dt_clust$m = ""
  anno_rng = sp_agg_dt[, .(ymin = 0, ymax = max(fraction)), .(cluster_id)]
  anno_rng$m = ""
  anno_dt_clust = merge(anno_dt_clust, anno_rng, by = "m", allow.cartesian = TRUE)

  p_sp_arch = ggplot() +
    geom_rect(data = anno_dt_clust, aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax), color = "gray", fill = "gray", size = .2) +
    geom_arch(GRanges(sp_agg_dt), aes(height = fraction, color = strand)) +
    facet_grid(cluster_id~., scales = "free_y") +
    scale_color_manual(values = strand_cols) +
    labs(y = "fraction") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust= 1))
  # coord_cartesian(xlim = c(50375000, 50400000)) +
  # annotate("point", x = c(ex_add$start, ex_add$end), y = 0, col = "green") +
  # annotate("segment", x = ex_add$start, xend = ex_add$end, y = 0, yend = 0, col = "green")

  pg_arch = cowplot::plot_grid(p_sp_heat + labs(title = "name"), p_sp_arch)

  ggsave(pdf_name, pg_arch, width = 8, height = 5)

  sp_assign_dt = unique(sp_clust_dt[, .(cluster_id, sample)])[order(cluster_id)]
  fwrite(sp_assign_dt, csv_name)

  sp_assign_dt = unique(sp_clust_dt[, .(cluster_id, sample)])

  invisible(list(clust_dt = sp_clust_dt, agg_dt = sp_agg_dt, assign_dt = sp_assign_dt))
}

plot_splice_pileup = function(sel_prof_qdt, sp_agg_dt, ik_gr, sp_assign_dt, return_data = FALSE){
  prof_gr = reduce(ik_gr)
  prof_gr = range(prof_gr)
  prof_gr = resize(prof_gr, 1.2*width(prof_gr), fix = "center")

  prof_dt = ssvFetchBam(sel_prof_qdt, prof_gr, win_size = win_size, return_data.table = TRUE, fragLens = NA)
  prof_dt.add = ssvFetchBam(sel_prof_qdt, prof_gr, win_size = win_size, return_data.table = TRUE, splice_strategy = "add", fragLens = NA)

  prof_dt = prof_dt[, .(y = sum(y)), .(x, seqnames, start, end, id, strand, sample)]
  prof_dt.add = prof_dt.add[, .(y = sum(y)), .(x, seqnames, start, end, id, strand, sample)]



  prof_dt = merge(sp_assign_dt, prof_dt, by = "sample")
  prof_dt$sample = factor(prof_dt$sample, levels = levels(sp_assign_dt$sample))

  prof_dt.add = merge(sp_assign_dt, prof_dt.add, by = "sample")
  prof_dt.add$sample = factor(prof_dt.add$sample, levels = levels(sp_assign_dt$sample))


  agg_prof_dt = prof_dt[, .(y = mean(y)), .(cluster_id, x, id, strand, start, end)]
  agg_prof_dt.add = prof_dt.add[, .(y = mean(y)), .(cluster_id, x, id, strand, start, end)]

  # ggplot(agg_prof_dt, aes(x = x, y = y, color = strand)) +
  #   geom_path() +
  #   facet_grid(cluster_id~id, scales = "free_x")

  agg_prof_dt[, xgen := (start + end)/2]
  agg_prof_dt.add[, xgen := (start + end)/2]

  anno_dt_pile = as.data.table(ex_basic)[, .(start, end)]
  anno_dt_pile$m = ""
  anno_rng = agg_prof_dt.add[, .(ymin = 0, ymax = max(y)), .(cluster_id)]
  anno_rng$m = ""
  anno_dt_pile = merge(anno_dt_pile, anno_rng, by = "m", allow.cartesian = TRUE)

  if(return_data){
    return(list(
      annotation = anno_dt_pile,
      pileup = agg_prof_dt,
      splicing = sp_agg_dt
    ))
  }

  p_pileup = ggplot() +
    geom_rect(data = anno_dt_pile, aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax), fill = "lightgray", color = "lightgray") +
    geom_ribbon(data = agg_prof_dt.add, aes(x = xgen, ymin = 0, ymax = y, group = paste(id, strand)), alpha = .2, fill = "dodgerblue2") +
    # geom_path(data = agg_prof_dt, aes(x = xgen, y = y, color = strand, group = paste(id, strand))) +
    geom_ribbon(data = agg_prof_dt, aes(x = xgen, ymin = 0, ymax = y, group = paste(id, strand)), alpha = 1, fill = "dodgerblue2") +
    geom_arch(GRanges(sp_agg_dt), aes(height = value), color = "black") +
    facet_grid(cluster_id~., scales = "free") + theme(panel.background = element_blank()) +
    labs(title = name)

  # ggsave(sub(".pdf", ".pileup.pdf", pdf_name), p_pileup, width = 12, height = 10)

  # p_pileup.3b4 = p_pileup +
  #   coord_cartesian(xlim = c(50365000, 50380000)) +
  #   labs(subtitle = "exon3b and exon4 region")

  # ggsave(sub(".pdf", ".pileup_ex3b.pdf", pdf_name), p_pileup.3b4, width = 12, height = 10)
  return(p_pileup)
}
