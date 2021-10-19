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

plot_splice_pileup = function(sel_prof_qdt, sp_agg_dt, ik_gr, sp_assign_dt, return_data = FALSE, pool_by = NULL, pool_FUN = sum, sample_sep = "\n"){
  prof_gr = reduce(ik_gr)
  prof_gr = range(prof_gr)
  prof_gr = resize(prof_gr, 1.2*width(prof_gr), fix = "center")

  if(!is.null(pool_by)){
    stopifnot(pool_by %in% colnames(sel_prof_qdt))
  }

  prof_dt = ssvFetchBam(sel_prof_qdt, prof_gr, win_size = win_size, return_data.table = TRUE, fragLens = NA)
  prof_dt.add = ssvFetchBam(sel_prof_qdt, prof_gr, win_size = win_size, return_data.table = TRUE, splice_strategy = "add", fragLens = NA)

  prof_dt = prof_dt[, .(y = sum(y)), .(x, seqnames, start, end, id, strand, sample)]
  prof_dt.add = prof_dt.add[, .(y = sum(y)), .(x, seqnames, start, end, id, strand, sample)]

  browser()
  if(!is.null(pool_by)){
    new_sample_dt = unique(sel_prof_qdt[, c(pool_by), with = FALSE])
    new_sample_dt$sample = apply(new_sample_dt, 1, function(x){paste(x, collapse = sample_sep)})

    prof_dt = prof_dt[, .(y = pool_FUN(y)), c("x", "seqnames", "start", "end", "id", "strand", pool_by)]
    prof_dt.add = prof_dt.add[, .(y = pool_FUN(y)), c("x", "seqnames", "start", "end", "id", "strand", pool_by)]

    prof_dt = merge(prof_dt, new_sample_dt, by = pool_by)
    prof_dt.add = merge(prof_dt.add, new_sample_dt, by = pool_by)
  }

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

#' Title
#'
#' @param sj_dt
#' @param bam_files
#' @param view_gr
#' @param max_span
#'
#' @return
#' @export
#'
#' @examples
#' options(mc.cores = 20)
#' set.seed(0)
#'
#' # wd = "/slipstream/home/joeboyd/R/SF_Ikaros_splicing/data"
#' wd = "/slipstream/home/dbgap/data/alignment_RNA-Seq/"
#' sj_files = find_SJ.out.tab_files(wd)[1:10]
#' bam_files = find_bam_files(wd)[1:10]
#' features = load_gtf("~/gencode.v36.annotation.gtf", gene_of_interest = "IKZF1")
#'
#' sp_sj_dt = load_splicing_from_SJ.out.tab_files(sj_files, view_gr = features$view_gr)
#' # sp_bam_dt = load_splicing_from_bam_files(bam_files, view_gr = features$view_gr)
#' sp_sj_dt.sel = sp_sj_dt[number_unique > 400]
#'
#' plots = plot_sj_pileups(sp_sj_dt.sel, bam_files[1:4], features$view_gr)
#' plots$assembled
#'
#' bam_qdt = data.table(file = bam_files[1:10], a = LETTERS[rep(c(1, 2), each = 5)], b = LETTERS[rep(c(3, 4), 5)])
#' plot_sj_pileups(sp_sj_dt.sel, bam_qdt, features$view_gr, pool_by = c("a", "b"))
plot_sj_pileups = function(sj_dt, bam_files, view_gr, max_span = 200, pool_by = NULL, pool_FUN = sum, sample_sep = "\n"){
  gr_starts = GRanges(unique(sj_dt[, .(seqnames, start, end = start, strand)]))
  gr_starts$id = paste("start", seq_along(gr_starts), sep = "_")
  names(gr_starts) = start(gr_starts)

  gr_ends = GRanges(unique(sj_dt[, .(seqnames, start = end, end, strand)]))
  gr_ends$id = paste("end", seq_along(gr_ends), sep = "_")
  names(gr_ends) = end(gr_ends)

  sj_dt$start_id = gr_starts[as.character(sj_dt$start)]$id
  sj_dt$end_id = gr_ends[as.character(sj_dt$end)]$id

  if(is.character(bam_files)){
    bam_qdt = data.table(file= bam_files)
    bam_qdt[, sample := sub(".Aligned.sortedByCoord.out.bam", "", basename(file))]
    bam_qdt[, sample_split := gsub("_", "\n", sample)]
  }else if(is.data.frame(bam_files)){
    bam_qdt = bam_files
  }else{
    stop("bam_files must be character or data.frame")
  }

  if(!is.null(pool_by)){
    stopifnot(pool_by %in% colnames(bam_qdt))
  }


  message("Retrieving pileups at splice starts...")
  prof_at_starts = ssvRecipes::ssvFetchBamPE.RNA(bam_qdt, resize(gr_starts, max_span*2, fix = "center"),
                                                 target_strand = "both",
                                                 return_data.table = TRUE,
                                                 win_size = 1,
                                                 n_region_splits = 20)

  message("Retrieving pileups at splice ends...")
  prof_at_ends = ssvRecipes::ssvFetchBamPE.RNA(bam_qdt, resize(gr_ends, max_span*2, fix = "center"),
                                               target_strand = "both",
                                               return_data.table = TRUE,
                                               win_size = 1,
                                               n_region_splits = 20)

  if(!is.null(pool_by)){
    new_sample_dt = unique(bam_qdt[, c(pool_by), with = FALSE])
    new_sample_dt$sample = apply(new_sample_dt, 1, function(x){paste(x, collapse = sample_sep)})

    prof_at_starts = prof_at_starts[, .(y = pool_FUN(y)), c("x", "seqnames", "start", "end", "id", "strand", pool_by)]
    prof_at_ends = prof_at_ends[, .(y = pool_FUN(y)), c("x", "seqnames", "start", "end", "id", "strand", pool_by)]

    prof_at_starts = merge(prof_at_starts, new_sample_dt, by = pool_by)
    prof_at_ends = merge(prof_at_ends, new_sample_dt, by = pool_by)
  }

  plot_stuff = function(prof_at){
    prof_at[, facet := paste(sample, strand)]

    clust_dt = ssvSignalClustering(prof_at, facet_ = "facet", dcast_fill = 0, max_cols = Inf, max_rows = Inf)

    pg_heat = cowplot::plot_grid(ncol = 1,
                                 ssvSignalHeatmap(clust_dt[strand == "+"], fill_limits = c(0, 1e3), facet_ = "sample_split", max_cols = Inf) +
                                   labs(title = "(+) strand") +
                                   scale_x_continuous(breaks = c(-max_span, 0, max_span)) +
                                   labs(fill = "reads", x = "bp", y = "splice junctions"),
                                 ssvSignalHeatmap(clust_dt[strand == "-"], fill_limits = c(0, 1e3), facet_ = "sample_split", max_cols = Inf)  +
                                   labs(title = "(-) strand") +
                                   scale_x_continuous(breaks = c(-max_span, 0, max_span)) +
                                   labs(fill = "reads", x = "bp", y = "splice junctions")
    )

    agg_dt = clust_dt[, .(y = mean(y)), .(x, cluster_id, sample, sample_split, strand)]
    p_side = ggplot(agg_dt, aes(x = x, y = y, color = strand, group = strand)) +
      geom_path() +
      facet_grid(cluster_id~sample_split, scales= "free_y") +
      geom_vline(xintercept = 0, color = "blue") +
      scale_color_manual(values = c("-" = "red", "+" = "blue")) +
      scale_x_continuous(breaks = c(-max_span, 0, max_span)) +
      theme(axis.text.x = element_text(angle = 0, vjust = .5, hjust = 1)) +
      labs(y = "mean reads", x= "bp") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
    list(heatmap = pg_heat, side = p_side, assign_dt = unique(clust_dt[, .(id, cluster_id)]))
  }

  set.seed(0)
  plots_at_starts = plot_stuff(prof_at_starts)
  plots_at_starts$side =  plots_at_starts$side + labs(title = "splice starts")
  set.seed(0)
  plots_at_ends = plot_stuff(prof_at_ends)
  plots_at_ends$side =  plots_at_ends$side + labs(title = "splice starts")

  plots_at_starts$assign_dt[cluster_id %in% c(2, 3)]
  plots_at_ends$assign_dt[cluster_id %in% c(2, 3, 4, 5)]

  pg_side_heatmap = cowplot::plot_grid(ncol = 1,
                                       plots_at_starts$heatmap,
                                       plots_at_ends$heatmap
  )


  pg_side_plots = cowplot::plot_grid(ncol = 1,
                                     plots_at_starts$side,
                                     plots_at_ends$side +
                                       labs(title = "splice ends")
  )

  pg_assembled = cowplot::plot_grid(nrow = 1,
                                    pg_side_heatmap,
                                    pg_side_plots)
  invisible(list(
    plots = list(
      starts_heatmap = plots_at_starts$heatmap,
      ends_heatmap = plots_at_ends$heatmap,
      starts_sideplot = plots_at_starts$side,
      ends_sideplot = plots_at_ends$side

    ),
    side_plots = pg_side_plots,
    heatmaps = pg_side_heatmap,
    assembled = pg_assembled
  ))
}
