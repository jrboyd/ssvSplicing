#' Title
#'
#' @param sel_sample
#' @param min_test_value
#' @param win_size
#'
#' @return
#' @export
#' @import cowplot
#' @examples
#' wd = "/slipstream/home/dbgap/data/alignment_RNA-Seq/"
#' sj_files = find_SJ.out.tab_files(wd)[1:10]
#'
#' features = load_gtf("~/gencode.v36.annotation.gtf", gene_of_interest = "IKZF1")
#'
#' sp_sj_dt = load_splicing_from_SJ.out.tab_files(sj_files, view_gr = features$view_gr)
#' sp_sj_dt.sel = sp_sj_dt[number_unique > 5]
#'
#' plots = plot_splice_heatmap(sp_sj_dt.sel)
#'
#' plots$heatmap
#'
plot_splice_heatmap = function(splice_dt,
                               nclust = 4,
                               min_test_value = 0,
                               min_fraction_present = 0,
                               sel_strand = c("+", "-", "*"),
                               win_size = 1){
  splice_dt = copy(splice_dt)
  sel_id = unique(splice_dt[number_unique > min_test_value]$id)
  # splice_dt =   splice_dt[number_unique > min_test_value]
  splice_dt = splice_dt[id %in% sel_id]
  splice_dt = splice_dt[strand %in% sel_strand]

  dtw = melt(splice_dt[, .(seqnames, start, end, strand, sample, number_unique)], measure.vars = c("number_unique"))
  dtw.filled = dcast(dtw,
                     seqnames+start+end+strand+variable~sample,
                     value.var = "value",
                     fill = 0)
  dtw = melt(dtw.filled, id.vars = c("seqnames", "start", "end", "strand", "variable"), variable.name = "sample", value.name = "value")

  dtw[, id := paste(seqnames, start, end, strand)]
  dtw[, fraction := value / sum(value), .(sample)]

  sel_id = dtw[, .(test_val = max(value)), .(id)][order(test_val)][test_val >= min_test_value]$id

  sel_dtw = dtw[id %in% sel_id]
  sel_dtw[, fraction := value / max(sum(value),1), .(sample)]
  sel_dtw$id = factor(sel_dtw$id, levels = sel_id)


  set.seed(0)
  id_lev = levels(ssvSignalClustering(sel_dtw, row_ = "id", column_ = "sample", fill_ = "fraction", facet_ = "", dcast_fill = 0, max_cols = Inf, nclust = 5)$id)

  set.seed(1)
  sp_clust_dt = ssvSignalClustering(sel_dtw,
                                    nclust = nclust,
                                    row_ = "sample", column_ = "id",
                                    fill_ = "fraction", facet_ = "",
                                    dcast_fill = 0,
                                    max_cols = Inf,

                                    max_rows = Inf)
  set.seed(NULL)
  sp_clust_dt$id = factor(sp_clust_dt$id, levels = id_lev)

  n_samples = length(unique(dtw$sample))
  toplot_id = as.character(sel_dtw[, sum(fraction > 0) >= min_fraction_present*n_samples, id][V1 == TRUE]$id)
  p_sp_heat = ssvSignalHeatmap(sp_clust_dt[id %in% toplot_id], row_ = "sample", column_ = "id", fill_ = "fraction", facet_ = "", dcast_fill = 0, max_cols = Inf, nclust = 3)

  splice_dt.agg = sp_clust_dt[, .(fraction = mean(fraction), value = mean(value)), .(seqnames, start, end, strand, id, cluster_id)]

  sp_assign_dt = unique(sp_clust_dt[, .(cluster_id, sample)])[order(cluster_id)]

  invisible(list(heatmap = p_sp_heat, clust_dt = sp_clust_dt, agg_dt = splice_dt.agg, assign_dt = sp_assign_dt))
}

#' plot_view_pileup_and_splicing
#'
#' @param splice_dt
#' @param return_data
#' @param pool_by
#' @param pool_FUN
#' @param sample_sep
#' @param bam_files
#' @param view_gr
#' @param plot_title
#' @param splice_height_var
#' @param anno_dt
#' @param exons_to_show
#' @param win_size
#' @param prof_dt
#' @param prof_dt.add
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
#' sp_sj_dt.sel = sp_sj_dt[number_unique > 5]
#'
#' plots = plot_view_pileup_and_splicing(sp_sj_dt.sel, bam_files[1:4], features$view_gr, exons_to_show = features$goi_exon_dt[transcript_id == "ENST00000331340.8"])
#'
#' bam_qdt = data.table(file = bam_files)
#' bam_qdt[, sample := sub(".Aligned.sorted.+", "", basename(file))]
#'
#' anno_dt = unique(sp_sj_dt.sel[, .(sample)])
#' anno_dt$group = rep(c("a", "b"), each = 5)
#' anno_dt$group2 = rep(c("c", "d"), 5)
#' debug(plot_view_pileup_and_splicing)
#'
#' plots = plot_view_pileup_and_splicing(sp_sj_dt.sel, bam_files[1:4], features$view_gr,
#'   exons_to_show = features$goi_exon_dt[transcript_id == "ENST00000331340.8"])
#'
#' plots.pooled = plot_view_pileup_and_splicing(sp_sj_dt.sel, bam_qdt[1:4,], features$view_gr,
#'   anno_dt = anno_dt,
#'   pool_by = c("group", "group2"),
#'   exons_to_show = features$goi_exon_dt[transcript_id == "ENST00000331340.8"])
#'
#' cowplot::plot_grid(plots, plots.pooled)
plot_view_pileup_and_splicing = function(splice_dt,
                                         bam_files,
                                         view_gr,
                                         plot_title = "",
                                         splice_height_var = "number_unique",
                                         return_data = FALSE,
                                         anno_dt = NULL,
                                         pool_by = NULL,
                                         pool_FUN = sum,
                                         sample_sep = "\n",
                                         exons_to_show = NULL,
                                         win_size = 50,
                                         prof_dt = NULL,
                                         prof_dt.add = NULL){
  if(is.character(bam_files)){
    bam_qdt = data.table(file= bam_files)
    bam_qdt[, sample := sub(".Aligned.sortedByCoord.out.bam", "", basename(file))]
  }else if(is.data.frame(bam_files)){
    bam_qdt = bam_files
  }else{
    stop("bam_files must be character or data.frame")
  }

  #try to sync sample names
  if(!all(splice_dt$sample %in% bam_qdt$sample)){
    splice_samples = unique(splice_dt$sample)
    k = sapply(splice_samples, function(x)which(grepl(x, bam_qdt$sample)))
    if(any(lengths(k) > 1)){stop("cannot map splice_dt sample names uniquely to bam_qdt sample names (basename of bam_files).")}
    splice_samples = splice_samples[lengths(k) > 0]
    k = unlist(k[lengths(k) > 0])
    splice_dt = splice_dt[sample %in% splice_samples]

    if(nrow(splice_dt) > 0){
      s2s = bam_qdt$sample[k]
      names(s2s) = names(k)
      splice_dt$sample = s2s[splice_dt$sample]
    }else{
      warning("Could not sync sample names between splice_dt and bam_files. Splicing arches will not be plotted.")
    }


  }

  if(!is.null(pool_by)){
    if(is.null(anno_dt)){
      stop("Need anno_dt to provided sample metadata.")
    }
    stopifnot(pool_by %in% colnames(anno_dt))
    if(!all(bam_qdt$sample %in% anno_dt$sample)){
      stop(paste0(
        "Some sample do not have metadata in anno_dt:\n ",
        paste(setdiff(bam_qdt$sample, anno_dt$sample), collapse = "\n ")
      ))
    }
  }

  if(is.null(prof_dt)){
    message("Retrieving direct read pileups...")
    prof_dt.raw = ssvFetchBam(bam_qdt, view_gr, win_size = win_size, return_data.table = TRUE, fragLens = NA)
  }else{
    message("Using provided direct read pileups.")
    prof_dt.raw = prof_dt
  }
  if(is.null(prof_dt.add)){
    message("Retrieving split read pileups...")
    prof_dt.add.raw = ssvFetchBam(bam_qdt, view_gr, win_size = win_size, return_data.table = TRUE, splice_strategy = "add", fragLens = NA)
  }else{
    message("Using provided split read pileups.")
    prof_dt.add.raw = prof_dt.add
  }

  prof_dt = prof_dt.raw[, .(y = sum(y)), .(x, seqnames, start, end, id, strand, sample)]
  prof_dt.add = prof_dt.add.raw[, .(y = sum(y)), .(x, seqnames, start, end, id, strand, sample)]

  if(!is.null(pool_by)){
    new_sample_dt = unique(anno_dt[, c(pool_by), with = FALSE])
    new_sample_dt$sample = apply(new_sample_dt, 1, function(x){paste(x, collapse = sample_sep)})

    prep_attributes = function(.dt, .anno_dt, .new_sample_dt, y_var = "y"){
      for(ko in setdiff(intersect(colnames(.dt), colnames(.anno_dt)), "sample")){
        .dt[[ko]] = NULL
      }
      .dt = merge(.dt, .anno_dt, by = "sample")
      .pool_by = intersect(colnames(.dt), c("x", "seqnames", "start", "end", "id", "strand", pool_by))
      .dt = .dt[, .(y_ = pool_FUN(get(y_var))), .pool_by]
      setnames(.dt, "y_", y_var)
      merge(.dt, .new_sample_dt, by = pool_by)
    }
    prof_dt = prep_attributes(prof_dt, anno_dt, new_sample_dt)
    prof_dt.add = prep_attributes(prof_dt.add, anno_dt, new_sample_dt)
    splice_dt = prep_attributes(splice_dt, anno_dt, new_sample_dt, y_var = splice_height_var)
  }

  agg_prof_dt = prof_dt#prof_dt[, .(y = mean(y)), .(cluster_id, x, id, strand, start, end)]
  agg_prof_dt.add = prof_dt.add#prof_dt.add[, .(y = mean(y)), .(cluster_id, x, id, strand, start, end)]

  # ggplot(agg_prof_dt, aes(x = x, y = y, color = strand)) +
  #   geom_path() +
  #   facet_grid(cluster_id~id, scales = "free_x")

  agg_prof_dt[, xgen := (start + end)/2]
  agg_prof_dt.add[, xgen := (start + end)/2]

  if(is.null(exons_to_show)){
    anno_dt_pile = NULL
  }else{
    anno_dt_pile = as.data.table(exons_to_show)[, .(start, end)]
    anno_dt_pile$m = ""
    anno_rng = agg_prof_dt.add[, .(ymin = 0, ymax = max(y)), .(sample)]
    anno_rng$m = ""
    anno_dt_pile = merge(anno_dt_pile, anno_rng, by = "m", allow.cartesian = TRUE)
  }


  if(return_data){
    return(list(
      annotation = anno_dt_pile,
      pileup = agg_prof_dt,
      splicing = splice_dt
    ))
  }

  p_pileup = ggplot()

  if(!is.null(anno_dt_pile)){
    p_pileup = p_pileup +
      geom_rect(data = anno_dt_pile,
                aes(xmin = start,
                    xmax = end,
                    ymin = ymin,
                    ymax = ymax),
                fill = "lightgray",
                color = "lightgray")
  }
  p_pileup = p_pileup +
    geom_ribbon(data = agg_prof_dt.add,
                aes(x = xgen,
                    ymin = 0,
                    ymax = y,
                    group = paste(id, strand)),
                alpha = .2,
                fill = "dodgerblue2") +
    geom_ribbon(data = agg_prof_dt,
                aes(x = xgen,
                    ymin = 0,
                    ymax = y,
                    group = paste(id, strand)),
                alpha = 1,
                fill = "dodgerblue2") +
    geom_arch(GRanges(splice_dt),
              aes_string(height = splice_height_var),
              color = "black") +
    facet_grid(sample~., scales = "free") +
    theme(panel.background = element_blank()) +
    labs(title = plot_title)

  return(list(plot = p_pileup,
              prof_dt = prof_dt.raw,
              prof_dt.add = prof_dt.add.raw))
}

#' plot_sj_pileups
#'
#' @param sj_dt
#' @param bam_files
#' @param view_gr
#' @param max_span
#' @param pool_by
#' @param pool_FUN
#' @param sample_sep
#' @param prof_at_starts
#' @param prof_at_ends
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
#' sp_sj_dt.sel = sp_sj_dt[number_unique > 100]
#'
#' plots = plot_sj_pileups(sp_sj_dt.sel, bam_files[1:4], features$view_gr)
#' plots$assembled
#'
#' bam_qdt = data.table(file = bam_files[1:10], a = LETTERS[rep(c(1, 2), each = 5)], b = LETTERS[rep(c(3, 4), 5)])
#' plots.pooled = plot_sj_pileups(sp_sj_dt.sel, bam_qdt, features$view_gr, pool_by = c("a", "b"))
plot_sj_pileups = function(sj_dt,
                           bam_files,
                           view_gr,
                           max_span = 200,
                           pool_by = NULL,
                           pool_FUN = sum,
                           sample_sep = "\n",
                           prof_at_starts = NULL,
                           prof_at_ends = NULL,
                           rel_heights.assembly = c(2, 1)){
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
  }else if(is.data.frame(bam_files)){
    bam_qdt = bam_files
  }else{
    stop("bam_files must be character or data.frame")
  }
  if(!is.null(pool_by)){
    stopifnot(pool_by %in% colnames(bam_qdt))
  }
  if(is.null(prof_at_starts)){
    message("Retrieving pileups at splice starts...")
    prof_at_starts.raw = ssvRecipes::ssvFetchBamPE.RNA(bam_qdt, resize(gr_starts, max_span*2, fix = "center"),
                                                       target_strand = "both",
                                                       return_data.table = TRUE,
                                                       win_size = 1,
                                                       n_region_splits = 20)
  }else{
    message("Using provided prof_at_starts.")
    prof_at_starts.raw = prof_at_starts
  }

  if(is.null(prof_at_ends)){
    message("Retrieving pileups at splice ends...")
    prof_at_ends.raw = ssvRecipes::ssvFetchBamPE.RNA(bam_qdt, resize(gr_ends, max_span*2, fix = "center"),
                                                     target_strand = "both",
                                                     return_data.table = TRUE,
                                                     win_size = 1,
                                                     n_region_splits = 20)
  }else{
    message("Using provided prof_at_ends.")
    prof_at_ends.raw = prof_at_ends
  }

  if(!is.null(pool_by)){
    new_sample_dt = unique(bam_qdt[, c(pool_by), with = FALSE])
    new_sample_dt$sample = apply(new_sample_dt, 1, function(x){paste(x, collapse = sample_sep)})

    prof_at_starts = prof_at_starts.raw[, .(y = pool_FUN(y)), c("x", "seqnames", "start", "end", "id", "strand", pool_by)]
    prof_at_ends = prof_at_ends.raw[, .(y = pool_FUN(y)), c("x", "seqnames", "start", "end", "id", "strand", pool_by)]

    prof_at_starts = merge(prof_at_starts, new_sample_dt, by = pool_by)
    prof_at_ends = merge(prof_at_ends, new_sample_dt, by = pool_by)
  }

  plot_stuff = function(prof_at){
    prof_at[, facet := paste(sample, strand)]

    clust_dt = ssvSignalClustering(prof_at, facet_ = "facet", dcast_fill = 0, max_cols = Inf, max_rows = Inf)

    pg_heat = cowplot::plot_grid(ncol = 1,
                                 ssvSignalHeatmap(clust_dt[strand == "+"], fill_limits = c(0, 1e3), facet_ = "sample", max_cols = Inf) +
                                   labs(title = "(+) strand") +
                                   scale_x_continuous(breaks = c(-max_span, 0, max_span)) +
                                   labs(fill = "reads", x = "bp", y = "splice junctions"),
                                 ssvSignalHeatmap(clust_dt[strand == "-"], fill_limits = c(0, 1e3), facet_ = "sample", max_cols = Inf)  +
                                   labs(title = "(-) strand") +
                                   scale_x_continuous(breaks = c(-max_span, 0, max_span)) +
                                   labs(fill = "reads", x = "bp", y = "splice junctions")
    )

    agg_dt = clust_dt[, .(y = mean(y)), .(x, cluster_id, sample, sample, strand)]
    p_side = ggplot(agg_dt, aes(x = x, y = y, color = strand, group = strand)) +
      geom_path() +
      facet_grid(cluster_id~sample, scales= "free_y") +
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

  pg_side_heatmap = cowplot::plot_grid(nrow = 1,
                                       cowplot::plot_grid(ncol = 1, rel_heights = c(1, 10),
                                                          ggplot() + theme_void() + cowplot::draw_text("Splice Site Starts"),
                                                          plots_at_starts$heatmap),
                                       cowplot::plot_grid(ncol = 1, rel_heights = c(1, 10),
                                                          ggplot() + theme_void() + cowplot::draw_text("Splice Site Ends"),
                                                          plots_at_ends$heatmap)
  )


  pg_side_plots = cowplot::plot_grid(nrow = 1,
                                     plots_at_starts$side,
                                     plots_at_ends$side +
                                       labs(title = "splice ends")
  )

  pg_assembled = cowplot::plot_grid(ncol = 1, rel_heights = rel_heights.assembly,
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
    assembled = pg_assembled,
    assign_dt = list(starts = plots_at_starts$assign_dt, ends = plots_at_ends$assign_dt),
    prof_at_starts = prof_at_starts.raw,
    prof_at_ends = prof_at_ends.raw
  ))
}

