#https://github.com/comprna/SUPPA#differential-splicing-analysis-for-transcripts-and-local-events
#suppa 2 requires a tpm and psi file per condition
library(data.table)
library(ggplot2)
library(magrittr)
library(GenomicRanges)

### controlling parameters
#location of ReadsPerGene.out.tab files
wd = "/slipstream/home/dbgap/data/alignment_RNA-Seq/"
message(paste(dir(dir(wd, pattern = "salmon", full.names = TRUE)[1]), collapse = "\n"))

# sj_files = find_SJ.out.tab_files(wd)
# bam_files = find_bam_files(wd)
#
# suppa_diff = file.path(wd, "suppa2_diff")

features = load_gtf("~/gencode.v36.annotation.gtf", gene_of_interest = "IKZF1")
dir(suppa_diff, pattern = "dpsi$", full.names = TRUE)
dir(suppa_diff, pattern = "clustvec$", full.names = TRUE)

ref_gr = rtracklayer::import.gff(gtf)

goi_gn = "FBgn0283521" #lola
# goi_gn = "FBgn0027525" #just some other one


# goi_tx = subset(ref_gr, gene_name == goi_gn & type == "transcript")
table(ref_gr$type)
tx_gr = subset(ref_gr, type == "transcript")
names(tx_gr) = tx_gr$transcript_id
gn_gr = unlist(reduce(split(tx_gr, tx_gr$gene_name)))

goi_tx = subset(tx_gr, gene_name == goi_gn)$transcript_id
gn_gr[goi_gn]

ex_gr = subset(ref_gr, type == "exon")

tx_dt = as.data.table(tx_gr[goi_tx])
tx_dt = tx_dt[order(start)]
tx_dt = tx_dt[order(start)]

qgr = range(gn_gr[goi_gn])
qgr = resize(qgr, 1.2*width(qgr), fix = "center")

ex_dt = as.data.table(subset(ex_gr, type == "exon" & transcript_id %in% goi_tx))
ex_dt$transcript_id = factor(ex_dt$transcript_id, levels = unique(ex_dt$transcript_id))
ex_dt[, ymin := as.numeric(transcript_id)]
ex_dt[, ymax := as.numeric(transcript_id) + 1]

lev = levels(ex_dt$transcript_id)

p_exons = ggplot(ex_dt, aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = strand)) +
  geom_rect() +
  scale_y_continuous(breaks = seq_along(lev)+.5, labels = lev) +
  coord_cartesian(xlim = c(start(qgr), end(qgr)))


#isoform A3 A5 AF AL MX RI SE
all_res = list()

td = "isoform"
ds_dt = fread(paste0("/slipstream/home/joeboyd/pipeline_outputs/RNAseq/NSFR_out/NSFR_out_jrb/suppa2_diff/diffSplice_result_", td, ".dpsi"))
ds_dt[, c("gene_name", "transcript_id") := tstrsplit(V1, ";")]
cutoff = .01
ds_dt.sig = ds_dt[`isoform_18C-isoform_25C_p-val` < cutoff | `isoform_18C-isoform_30C_p-val` < cutoff | `isoform_25C-isoform_30C_p-val` < cutoff]

ds_dt[grepl(goi_gn, V1)]
length(unique(ds_dt$gene_name))
length(unique(ds_dt$transcript_id))

dt1 = fread("/slipstream/home/joeboyd/pipeline_outputs/RNAseq/NSFR_out/NSFR_out_jrb/suppa2_diff/tpm_18C.tpm")
dt2 = fread("/slipstream/home/joeboyd/pipeline_outputs/RNAseq/NSFR_out/NSFR_out_jrb/suppa2_diff/tpm_25C.tpm")
dt3 = fread("/slipstream/home/joeboyd/pipeline_outputs/RNAseq/NSFR_out/NSFR_out_jrb/suppa2_diff/tpm_30C.tpm")
dt = merge(merge(dt1, dt2, by = "V1"), dt3, by = "V1")
dt = melt(dt)
dt[, c("temp", "rep") := tstrsplit(variable, "[_\\.]", keep = c(3, 4))]
dt[, lg_val := log2(value + 1)]
dt[, z := (lg_val - mean(lg_val)) / min(sd(lg_val), 1, na.rm = TRUE), V1 ]
ko = dt[is.nan(z)]$V1
dt = dt[!V1 %in% ko]
dt[z > 3, z := 3]
dt[z < -3, z := -3]

ssvSignalHeatmap(dt[V1 %in% ds_dt.sig$transcript_id], row_ = "V1", column_ = "rep", facet_ = "temp", fill_ = "value")
ssvSignalHeatmap(dt[V1 %in% ds_dt.sig$transcript_id], row_ = "V1", column_ = "rep", facet_ = "temp", fill_ = "z")

for(td in c("isoform", "A3", "A5", "AF", "AL", "MX", "RI", "SE")){
  suppressWarnings({
  ds_dt = fread(paste0("/slipstream/home/joeboyd/pipeline_outputs/RNAseq/NSFR_out/NSFR_out_jrb/suppa2_diff/diffSplice_result_", td, ".dpsi"))
  cl_dt = fread(paste0("/slipstream/home/joeboyd/pipeline_outputs/RNAseq/NSFR_out/NSFR_out_jrb/suppa2_diff/clusterEvents_", td, ".clustvec"))
  })
  res = ds_dt[grepl(goi_gn, V1)]

  # res = lapply(goi_tx, function(tx_id){
  #   ds_dt[grepl(tx_id, V1)]
  # }) %>% rbindlist
  all_res[[td]] = res
}

all_res

tpm_files = dir("/slipstream/home/joeboyd/pipeline_outputs/RNAseq/NSFR_out/NSFR_out_jrb/suppa2_diff", pattern = "tpm", full.names = TRUE)

bam_files = dir(wd, pattern = "sort.+bam$", full.names = TRUE)
cfg_dt = data.table(file = bam_files)
cfg_dt[, c("exp", "temp", "rep") := tstrsplit(basename(file), "[_\\.]", keep = 1:3)]
cfg_dt[, sample := paste(temp, rep, sep = "\n")]

cfg_dt$mapped_reads = sapply(cfg_dt$file, function(f){
  temp = Rsamtools::idxstatsBam(f)
  sum(temp[,3])
})

bam_dt = ssvRecipes::ssvFetchBamPE.RNA(cfg_dt, qgr, win_size = 50, return_data.table = TRUE)
bam_dt[, x := (start + end) / 2]

bam_dt[, rpm := y / mapped_reads * 1e6]

ggplot(bam_dt[strand == '-'], aes(x = x, y = y, color = temp)) +
  geom_path() +
  facet_grid(rep~temp+strand)

p_prof = ggplot(bam_dt[strand == "-"], aes(x = x, y = rpm, color = temp)) +
  annotate("rect", xmin = ex_dt$start, xmax = ex_dt$end, ymin = 0, ymax = 50, fill = "gray90") +
  geom_path() +
  facet_grid(rep~temp) +
  theme(panel.background = element_blank(), panel.grid = element_blank()) +
  labs(title = "Lola exons in gray")

p_prof
ggsave(paste0("NSFR_", goi_gn, "_profiles.pdf"), p_prof, width = 10, height = 8)
