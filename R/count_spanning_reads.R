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
#' sj_files = find_SJ.out.tab_files(wd)
#' bam_files = find_bam_files(wd)
#' features = load_gtf("~/gencode.v36.annotation.gtf", gene_of_interest = "IKZF1")
#'
#' sp_sj_dt = load_splicing_from_SJ.out.tab_files(sj_files, view_gr = features$view_gr)
#' # sp_bam_dt = load_splicing_from_bam_files(bam_files, view_gr = features$view_gr)
#' sp_sj_dt.sel = sp_sj_dt[number_unique > 400]
#'
#' plots = plot_sj_pileups(sp_sj_dt.sel, bam_files[1:4], features$view_gr)
#' plots$assembled
plot_sj_pileups = function(sj_dt, bam_files, view_gr, max_span = 200){
  gr_starts = GRanges(unique(sj_dt[, .(seqnames, start, end = start, strand)]))
  strand(gr_starts) %>% table
  gr_starts$id = paste("start", seq_along(gr_starts), sep = "_")
  names(gr_starts) = start(gr_starts)

  gr_ends = GRanges(unique(sj_dt[, .(seqnames, start = end, end, strand)]))
  strand(gr_ends) %>% table
  gr_ends$id = paste("end", seq_along(gr_ends), sep = "_")
  names(gr_ends) = end(gr_ends)

  sj_dt$start_id = gr_starts[as.character(sj_dt$start)]$id
  sj_dt$end_id = gr_ends[as.character(sj_dt$end)]$id

  bam_qdt = data.table(file= bam_files)
  bam_qdt[, sample := sub(".Aligned.sortedByCoord.out.bam", "", basename(file))]
  bam_qdt[, sample_split := gsub("_", "\n", sample)]

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
#' sj_files = find_SJ.out.tab_files(wd)
#' bam_files = find_bam_files(wd)
#' features = load_gtf("~/gencode.v36.annotation.gtf", gene_of_interest = "IKZF1")
#'
#' sp_sj_dt = load_splicing_from_SJ.out.tab_files(sj_files, view_gr = features$view_gr)
#' # sp_bam_dt = load_splicing_from_bam_files(bam_files, view_gr = features$view_gr)
#' sp_sj_dt.sel = sp_sj_dt[number_unique > 400]
#'
#' count_sj_spanning_pairs(sp_sj_dt.sel, bam_files[1:4], features$view_gr)
count_sj_spanning_pairs = function(sj_dt, bam_files, view_gr, max_span = 200){
  gr_starts = GRanges(unique(sj_dt[, .(seqnames, start, end = start, strand)]))
  strand(gr_starts) %>% table
  gr_starts$id = paste("start", seq_along(gr_starts), sep = "_")
  names(gr_starts) = start(gr_starts)

  gr_ends = GRanges(unique(sj_dt[, .(seqnames, start = end, end, strand)]))
  strand(gr_ends) %>% table
  gr_ends$id = paste("end", seq_along(gr_ends), sep = "_")
  names(gr_ends) = end(gr_ends)

  sj_dt$start_id = gr_starts[as.character(sj_dt$start)]$id
  sj_dt$end_id = gr_ends[as.character(sj_dt$end)]$id

  bam_qdt = data.table(file= bam_files)
  bam_qdt[, sample := sub(".Aligned.sortedByCoord.out.bam", "", basename(file))]
  bam_qdt[, sample_split := gsub("_", "\n", sample)]



  names(gr_starts) = NULL
  names(gr_ends) = NULL
  # function(bam_files, )
  message("Retrieving reads at splice starts...")
  reads_at_starts = ssvFetchBam(bam_files, resize(gr_starts, max_span*2, fix = "center"), return_unprocessed = TRUE)
  message("Retrieving reads at splice ends...")
  reads_at_ends = ssvFetchBam(bam_files, resize(gr_ends, max_span*2, fix = "center"), return_unprocessed = TRUE)

  message("Filtering read pairs...")

  nrow(reads_at_starts)/1e6
  nrow(reads_at_ends)/1e6

  common = intersect(reads_at_starts$qname, reads_at_ends$qname)
  length(common)/1e6

  reads_at_starts = reads_at_starts[qname %in% common]
  reads_at_ends = reads_at_ends[qname %in% common]

  sel_qname = reads_at_starts$qname[1]
  a = reads_at_starts[qname %in% sel_qname]
  b = reads_at_ends[qname %in% sel_qname]

  i = 1

  sj_dt.i = sj_dt[i,]
  sj_dt.i$end - sj_dt.i$start
  i_start_id = sj_dt[i,]$start_id
  i_end_id = sj_dt[i,]$end_id

  i_starts_dt = reads_at_starts[id == i_start_id]
  i_ends_dt = reads_at_ends[id == i_end_id]

  sam_flags = as.list(Rsamtools::FLAG_BITNAMES)
  names(sam_flags) = sam_flags

  start_is_first = test_flag(sam_flags$isFirstMateRead, i_starts_dt$flag)
  start_is_first.qname = i_starts_dt$qname[start_is_first]
  end_is_second = test_flag(sam_flags$isSecondMateRead, i_ends_dt$flag)
  end_is_second.qname = i_ends_dt$qname[end_is_second]

  start_is_second = test_flag(sam_flags$isSecondMateRead, i_starts_dt$flag)
  start_is_second.qname = i_starts_dt$qname[start_is_second]
  end_is_first = test_flag(sam_flags$isFirstMateRead, i_ends_dt$flag)
  end_is_first.qname = i_ends_dt$qname[end_is_first]

  good_prime = intersect(start_is_first.qname, end_is_second.qname)
  good_second = intersect(start_is_second.qname, end_is_first.qname)

  i_starts_dt.paired = unique(rbind(
    i_starts_dt[start_is_first & qname %in% good_prime],
    i_starts_dt[start_is_second & qname %in% good_second]
  ))

  i_ends_dt.paired = unique(rbind(
    i_ends_dt[end_is_second & qname %in% good_prime],
    i_ends_dt[end_is_first & qname %in% good_second]
  ))

  nrow(i_starts_dt.paired) / nrow(i_starts_dt)
  nrow(i_ends_dt.paired) / nrow(i_ends_dt)

  paired_qname = intersect(i_starts_dt.paired$qname, i_ends_dt.paired$qname)

  message("Detecting spanning pairs...")

  pbmcapply::pbmclapply(paired_qname, function(t_qname){
    gr_s = expandCigar(i_starts_dt.paired[qname == t_qname])
    gr_e = expandCigar(i_ends_dt.paired[qname == t_qname])
    outside_splice_site =
      end(range(gr_s)) < sj_dt.i$start &
      start(range(gr_e)) > sj_dt.i$end

    outside_splice_site
  })

  data.table(qname = paired_qname)
  i_starts_dt.paired[, which_label := qname]
  cig_dt = expandCigar(i_starts_dt.paired, return_data.table = TRUE)[]
  cig_dt$cigar_type %>% table
  cig_dt_starts = cig_dt[, .(start_min = min(start), end_max = max(end)) , .(which_label)]

  i_ends_dt.paired[, which_label := qname]
  cig_dt = expandCigar(i_ends_dt.paired, return_data.table = TRUE)[]
  cig_dt$cigar_type %>% table
  cig_dt_ends = cig_dt[, .(start_min = min(start), end_max = max(end)) , .(which_label)]

  qname_spans_splice = intersect(
    cig_dt_starts[end_max < sj_dt.i$start]$which_label,
    cig_dt_ends[start_min > sj_dt.i$end]$which_label
  )

  length(qname_spans_splice)

  read_sample_dt = unique(i_starts_dt.paired[, .(qname, sample)])
  counts_sample_dt = read_sample_dt[qname %in% qname_spans_splice][, .N, .(sample)]

  zero_bams = setdiff(bam_files, counts_sample_dt$sample)
  if(length(zero_bams) > 0){
    counts_sample_dt = rbind(
      counts_sample_dt,
      data.table(sample = zero_bams, N = 0)
    )
  }
  counts_sample_dt
}

decode_flags = function(flags_to_test){
  lapply(flags_to_test, decode_flags.single)
}

decode_flags.single = function(flag_to_test){
  Rsamtools::FLAG_BITNAMES[which(intToBits(flag_to_test) > 0)]
}

test_flag = function(flag_to_match, flags_to_test){
  i_to_match = which(Rsamtools::FLAG_BITNAMES == flag_to_match)
  #convert to bits
  bits_to_test = matrix(intToBits(flags_to_test), ncol = 32, byrow = TRUE)
  bits_to_test[, i_to_match] > 0
}

