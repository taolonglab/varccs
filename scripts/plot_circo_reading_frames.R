
cat('* Loading R libraries...')
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Biostrings))
cat(' OK.\n')

#-------------------- Load files -----------------------------------------------

min_read_coverage_pct_def = 90
min_sec_size_def = 10
max_reads_per_sector_def = 200
reads_per_sector_step_def = 10
t1_norm = TRUE # Normalize reads to the same length (easier to look at comp)
y_lab_size = 0.5


option_list = list(
  make_option(
    c('-ic', '--input-contigs'),
    help = 'Input file with reading frame contigs. This is the output from
            reading_frame_analysis script (*_reading_frame_contigs.txt).'
  ),
  make_option(
    c('-ie', '--input-exons'),
    help = 'Input file with exons aligned to reads (*_exon_table_hmpfix.txt).'
  ),
  make_option(
    c('-if', '--input-fasta'),
    help = 'Input FASTA file with fixed homopolymer errors (*_hmpfix.fasta).'
  ),
  make_option(
    c('-o', '--output-circoplot'),
    help = 'Output filename for the circo plot (PDF).'
  ),
  make_option(
    c('--min-read-cov-pct'),
    default = min_read_coverage_pct_def
  ),
  make_option(
    c('--first-exon'),
    # default = 1,
    default = 'auto',
    help = 'Integer denoting first exon (to be used for tandem analysis).
    Default: "auto" (smallest exon_no from "input-exons" table).'
  ),
  make_option(
    c('--last-exon'),
    default = 'auto',
    help = 'Integer denoting last exon (to be used for tandem analysis).
    Default: "auto" (largest exon_no from "input-exons" table).'
  ),
  make_option(
    c('--min-sector-size'),
    default = min_sec_size_def,
    help = 'Minimum size (angular width) of a sector on a circo plot.'
  ),
  make_option(
    c('--max-reads-per-sector'),
    default = max_reads_per_sector_def,
    help = 'Maximum number of reads to plot per sector.'
  ),
  make_option(
    c('--reads-per-sector-step'),
    default = reads_per_sector_step_def,
    help = 'For sectors exceeding "max-reads-per-sector" reads, plot every
    "reads-per-sector-step"th reads. This is used to reduce plotting too
    many reads (lines) per sector.'
  ),
  make_option(
    c('--plot-overlap-offset'),
    default = 0.5,
    help = 'Offset to overlap reads on the plot so that we dont have whitespace
    between them.'
  ),
  make_option(
    c('-v', '--verbose'),
    default = TRUE,
    help = 'Print progress bars during execution.'
  )
)

opt = parse_args(OptionParser(option_list=option_list))

# Parse input arguments
input_contigs = opt$`input-contigs`
input_exons = opt$`input-exons`
input_reads_fa = opt$`input-fasta`
output_circoplot = opt$`output-circoplot`

min_read_cov = opt$`min-read-cov-pct`/100
first_exon = opt$`first-exon`
last_exon = opt$`last-exon`
min_sec_size = opt$`min-sector-size`
max_reads_per_sector = opt$`max-reads-per-sector`
reads_per_sector_step = opt$`reads-per-sector-step`
plot_overlap_offset = opt$`plot-overlap-offset`
verbose = opt$verbose


# Load files
# reads.dt contains full read lengths calculated from the FASTA file after
# homopolymer error filtering./
cat('* Loading input files')
blast.dt = fread(input_exons)
cat('.')
contigs.dt = fread(input_contigs)
cat('.')
reads.fa = readDNAStringSet(input_reads_fa)
cat('.')
reads.dt = data.table(
   read_id=sapply(names(reads.fa), function (s) strsplit(s, ';')[[1]][1]), 
   seq=as.character(reads.fa))
reads.dt[, seq_len := sapply(seq, nchar)]
reads.dt[, seq := NULL]
rm(reads.fa)
cat('.')
cat(' OK.\n')

# Generate read combinations
if (first_exon == 'auto') {
  first_exon = blast.dt[, min(exon_no)]
  cat('? autodetect first exon: ', first_exon, '\n')

}
if (last_exon == 'auto') {
  last_exon = blast.dt[, max(exon_no)]
  cat('? autodetect last exon: ', last_exon, '\n')
}

cat('* Summarizing reads for plotting')
blast.dt[exon_comb_short_rc %like% paste0(last_exon, ','), 
  tandem_midpoints := paste(
    .SD[exon_id %in% c(paste0('exon_', sprintf('%02d', last_exon)), 
                       paste0('exon_', sprintf('%02d', first_exon), '_rc')
                      ), read_end_pos_align], 
    collapse=','),
  by=read_id]
cat('.') # 1
reads_combs.dt = blast.dt[, .(total_aln_len=sum(alignment_length),
                              tandem_mids=unique(tandem_midpoints)), 
                          by=.(exon_comb_short_rc, read_id, read_count)]
cat('.') # 2
reads_combs.dt = merge(
   reads_combs.dt, 
   reads.dt[, .(read_id, seq_len)], 
   all.x=T, by='read_id'
)
cat('.') # 3
reads_combs.dt[, exon_coverage := total_aln_len/seq_len]
reads_combs.dt = reads_combs.dt[exon_coverage >= min_read_cov]
reads_combs.dt[, exon_comb_avg_len := mean(seq_len), by=exon_comb_short_rc]

reads_combs.dt = reads_combs.dt[order(-exon_comb_avg_len, -seq_len)]
reads_combs.dt[, x := 1:.N]
cat('.') # 4




sector.dt = reads_combs.dt[, .(n=.N, avg_len=unique(exon_comb_avg_len),
                               counts=sum(read_count)), by=exon_comb_short_rc]
sector.dt[, x2 := cumsum(n)]
sector.dt[, x1 := c(1, x2[1:(length(x2)-1)])]

sector.dt[, x1_n := x1]
sector.dt[, x2_n := x2]

# Adjust all sector lengths such that none is smaller than min_sec_size
for (i in 1:sector.dt[, .N]) {
   row = sector.dt[i]
   if (i > 1) sector.dt[i, x1_n := sector.dt[i-1, x2_n] + 1]
   if (row[, x2-x1+1] < min_sec_size) {
      sector.dt[i, x2_n := x1_n + min_sec_size - 1]
   }
   else {
      sector.dt[i, x2_n := x1_n + (x2 - x1)]
   }
}

# Adjust for inter-sector gaps of width "1"
sector.dt[, x2_n := c(x2_n + (1:.N))]
sector.dt[, x1_n := c(x1_n[1], x1_n[2:.N] + (1:(.N-1)))]

sector.dt[, log10_counts := log10(counts)]
sector.dt[, exon_comb_short_rc := factor(exon_comb_short_rc,
   levels = exon_comb_short_rc)]
cat('.') # 5

# Now, for sectors with very many reads show only a sample to 
# make the plot faster.
reads_combs.dt[
  exon_comb_short_rc %in% sector.dt[n>max_reads_per_sector, exon_comb_short_rc],
  remove := rep(c(F, rep(T, reads_per_sector_step - 1)), .N), 
  by=exon_comb_short_rc
]
reads_combs.dt[is.na(remove), remove := F]
reads_combs.dt[remove==F, xf := 1:.N]
cat('.') # 6

# Merge reads_combs with contig information
contigs2.dt = merge(
  contigs.dt, 
  reads_combs.dt[remove==F, .(read_id, x=xf)],
  by='read_id'
)
contigs2.dt[, start_abs := start + read_start_pos_align - 1]
contigs2.dt[, end_abs := end + read_start_pos_align - 1]
cat('.') # 7
# After we have reads_combs with remove==F, calculate sector widths
sector.dt[, n_s := reads_combs.dt[remove==F, .(n_s=.N), 
                                  by=exon_comb_short_rc][, n_s]]
sector.dt[, w := (x2_n-x1_n)/n_s]
cat('.') # 8
cat(' OK.\n')

#
#
# ----------------------- Circo plot -----------------------------
#
#

cat('* Generating circo plot\n')
n_combs = sector.dt[, .N]
# Generate colors for each exon
exon_col_int2 = colorRampPalette(brewer.pal(n=8, name='RdYlBu'))
exon_colors = exon_col_int2(last_exon - first_exon + 1)
names(exon_colors) = first_exon:last_exon
exon_colors = as.list(exon_colors)
cairo_pdf(output_circoplot, width=5, height=5)

circos.clear()
circos.par(cell.padding = c(0, 0, 0, 0), 
           # gap.degree=c(rep(1, sector.dt[, .N-1]), 1), 
           gap.degree = 0.5,
           points.overflow.warning=F,
           start.degree=90)

# exon_colors = shift(exon_colors, 8)

# If we plot each sector of equal width use this
# sector width

circos.genomicInitialize(
   # Sector width proportional to the number of reads in that sector
   sector.dt[, .(exon_comb_short_rc, x1_n, x2_n)],
   plotType=''
)

# ----------------------------- Read length Track -------------------------------
# Show read lengths
tid_len = 1 # Circo index for track with sequence lengths
len_grid_ticks = reads_combs.dt[, pretty(range(exon_comb_avg_len))]
len_grid_tick_len = 2.5

len_grid.dt = sector.dt[, .(x1=c(x1_n, x2_n-len_grid_tick_len), 
                            x2=c(x1_n+len_grid_tick_len, x2_n)), 
                        by=.(exon_comb_short_rc)]
len_grid.dt = len_grid.dt[, .(y=len_grid_ticks), 
                          by=.(exon_comb_short_rc, x1, x2)]

circos.genomicTrackPlotRegion(
   len_grid.dt,
   ylim = c(min(len_grid_ticks), max(len_grid_ticks)),
   # ylim = c(0, 2252)/(c(1, 2252)[t1_norm+1]),
   bg.border = 'black',
   track.height = 0.1,
   panel.fun = function (region, value, ...) {
      circos.genomicLines(region, value, type='segment', ...)
   }
)

x0 = 1
x_off = 0
prev_comb = reads_combs.dt[1, exon_comb_short_rc]

cat('  - Plotting length track\n')
cat('      tick marks: ', len_grid_ticks, '\n')
for (i in 1:reads_combs.dt[remove==F, .N]) {
   # i iterates over reads
   row = reads_combs.dt[remove==F][i]
   seq_len = reads.dt[read_id==row[, read_id], seq_len]
   lwid = sector.dt[exon_comb_short_rc == row[, exon_comb_short_rc], w]
   if (row[, exon_comb_short_rc] == prev_comb & i > 1) {
     x_off = plot_overlap_offset
   } else {
     x0 = sector.dt[exon_comb_short_rc == row[, exon_comb_short_rc], x1_n]
     x_off = 0
   }
   x1 = x0 + lwid
   circos.rect(
      xleft = x0 - x_off, xright = x1,
      ybottom = 0, ytop = row[, seq_len],
      col = adjustcolor('black', alpha.f=1.0), border=NA, 
      # baseline='top',
      track.index = tid_len, sector.index = row[, exon_comb_short_rc]
   )
   prev_comb = row[, exon_comb_short_rc]
   if ((i %% 10 == 0 | i == reads_combs.dt[remove==F, .N]) & verbose == T) {
     cat(paste0('\r      processed read ', i, ' out of ',
         reads_combs.dt[remove==F, .N], '.'))
   }
   x0 = x1
}
cat('\n')

# legend(x='topleft', title='', legend=sector.dt[, exon_comb_short_rc], 
#        y.intersp=0.8, cex=0.7, bty='n',
#        xjust=0, yjust=1)

# --------------Exon combinations and reading frames Tracks --------------------
# To do: remove black border as a background and manually draw the black outline
#        over the exon colors
# Add exon color legend.

cat('  - Plotting exon combination track\n')
tid_combs = 2
t1_ticks = seq(0, 2000, by=500)
t1_other_seq_col = 'black'
t2_other_seq_col = 'white'

# orf_colors = list('0'='black', '1'='#edf8b1', '2'='')
if (t1_norm) {
   t1_ticks = t1_ticks/t1_ticks[length(t1_ticks)]
}

circos.genomicTrackPlotRegion(
   ylim = c(0, 2252)/(c(1, 2252)[t1_norm+1]),
   bg.border = 'black',
   # bg.border = NA,
   # bg.col = exon_colors[['1']],
   bg.col = NA,
   track.height = 0.4
)

fs_colors = list('1'='#edf8b1', '2'='#7fcdbb', '0'='#2c7fb8')
# plot_combs = reads_combs.dt[, unique(exon_comb_short_rc)][1:21]
prev_comb = reads_combs.dt[remove==F][1, exon_comb_short_rc]
x0 = 1
x_off = 0
x1_v = rep(NA, reads_combs.dt[remove==F, .N])

for (i in 1:reads_combs.dt[remove==F, .N]) {
   # i iterates over reads
   row = reads_combs.dt[remove==F][i]
   lwid = sector.dt[exon_comb_short_rc == row[, exon_comb_short_rc], w]
   seq_len = reads.dt[read_id==row[, read_id], seq_len]
   if (row[, exon_comb_short_rc] == prev_comb & i > 1) {
     x_off = plot_overlap_offset
   } else {
     x0 = sector.dt[exon_comb_short_rc == row[, exon_comb_short_rc], x1_n]
     x_off = 0
   }
   coords.dt = blast.dt[read_id==row[, read_id], .(
      exon_no, read_start_pos_align, read_end_pos_align,
      exon_id
   )][order(read_start_pos_align)]
   # If the read is reverse complemented (e.g. skipped PCR primer alignment
   # step, we have to reverse plotting!
   if (coords.dt[1, exon_id] %like% '_rc') {
     rc = TRUE
   } else {
     rc = FALSE
   }
   x1 = x0 + lwid
   x1_v[i] = x1
   y0 = 1
   
   # j iterates over exons within each read
   for (j in 1:coords.dt[, .N]) {

      # Before drawing the line for each exon, check to see if there is a
      # region in the read not aligned to the exon and draw it in gray.
      # fs_cls = fs.dt[read_id==row[, read_id] & exon_id == coords.dt[j, exon_id] &
      #                read_start_pos_align == coords.dt[j, read_start_pos_align],
      #                exon_in_orf_class]
      # if (length(fs_cls) > 0) fs_col = fs_colors[[as.character(fs_cls)]]
      # else fs_col = t2_other_seq_col
      fs_col = t2_other_seq_col

      if (coords.dt[j, read_start_pos_align] > y0) {
         y_coords = c(y0, coords.dt[j, read_start_pos_align])
         if (t1_norm) {
            y_coords = y_coords/seq_len
         }
         if (rc) y_coords = ifelse(t1_norm, 1, seq_len) - y_coords
         circos.rect(
            xleft = x0 - x_off, xright = x1,
            ybottom = y_coords[1], ytop = min(y_coords[2], 1),
            col = t1_other_seq_col, border=NA,
            track.index = tid_combs, sector.index = row[, exon_comb_short_rc]
         )
      }
      y0 = coords.dt[j, read_start_pos_align]
      
      
      # Now draw the known exon in color
      y_coords = c(y0, coords.dt[j, read_end_pos_align + 1])
      if (t1_norm) {
         y_coords = y_coords/seq_len
      }
      if (rc) y_coords = ifelse(t1_norm, 1, seq_len) - y_coords
      circos.rect(
         xleft = x0 - x_off, xright = x1,
         ybottom = y_coords[1], ytop = min(y_coords[2], 1),
         col = exon_colors[[as.character(coords.dt[j, exon_no])]],
         border=NA,
         track.index = tid_combs, sector.index = row[, exon_comb_short_rc]
      )
      y0 = coords.dt[j, read_end_pos_align+1]
   }
   # If there is anything left after last exon, plot that with t1_other_seq_col
   if (coords.dt[j, read_end_pos_align] < seq_len) {
      y_coords = c(y0, seq_len)
      if (t1_norm) y_coords = y_coords/seq_len
      circos.rect(
         xleft = x0 - x_off, xright = x1,
         ybottom = y_coords[1], ytop = min(y_coords[2], 1),
         col = t1_other_seq_col, border=NA,
         track.index=tid_combs, sector.index = row[, exon_comb_short_rc]
      )
   }
   
   prev_comb = row[, exon_comb_short_rc]
   x0 = x1
   if ((i %% 10 == 0 | i == reads_combs.dt[remove==F, .N]) & verbose)  {
     cat(paste0('\r      processed read ', i, ' out of ', 
                reads_combs.dt[remove==F, .N], '.'))
   }
}

# sector.dt[, .(exon_comb_short_rc, x1_n, x2_n)]


reads_combs.dt[remove==F, x1 := x1_v]
if (!('x1' %in% names(contigs2.dt))) {
   contigs2.dt = merge(contigs2.dt, reads_combs.dt[remove==F, .(read_id, x1)],
                       by='read_id')
}
cat('\n')

#----------------------- Open Reading Frame analysis ---------------------------

cat('  - Plotting open reading frames\n')
fs_colors = list('0'='#e31a1c', '1'='#ffffb3', '2'='#6a3d9a', '3'='#a6cee3')
tid_orfs = 3
circos.genomicTrackPlotRegion(
   ylim = c(0, 2252)/(c(1, 2252)[t1_norm+1]),
   bg.border = 'black', 
   # bg.border = NA,
   bg.col = fs_colors[['3']],
   track.height = 0.3
)

k = 1
x_off = 0
prev_comb = contigs2.dt[1, exon_comb_short_rc]
prev_read = contigs2.dt[1, read_id]
read_stop = F

for (i in 1:contigs2.dt[, .N]) {
   row = contigs2.dt[i]
   if (row[, read_id] != prev_read) read_stop = F
   if (read_stop) next
   lwid = sector.dt[exon_comb_short_rc == row[, exon_comb_short_rc], w]
   if (k > 1 & row[, exon_comb_short_rc] == prev_comb) {
     x_off = plot_overlap_offset
   } else { 
     x_off = 0
     k = k + 1 
    }
   x0 = row[, x1 - lwid - x_off]
   x1 = row[, x1]
   y0 = row[, (start_abs - 1)/read_len]
   y1 = row[, min((end_abs + 1)/read_len, 1)]
   if (row[, type] == 2) { # stop codon
      y0 = y0 - 0.03 # Manually overemphasize STOP codon for visibility
      read_stop = T
   }
   circos.rect(
      xleft = x0, xright = x1, ybottom = y0, ytop = y1,
      col = fs_colors[[as.character(row[, type])]], border=NA,
      track.index = tid_orfs, sector.index = row[, exon_comb_short_rc]
   )
   prev_read = row[, read_id]
   prev_comb = row[, exon_comb_short_rc]
   # Progress bar
   if ((i %% 10 == 0 | i == contigs2.dt[, .N]) & verbose) {
     cat(paste0('\r      processed contig ', i, ' out of ',
        contigs2.dt[, .N], '.'))
   } 
}

# Plot black border frame around sectors
for (i in 1:sector.dt[, .N]) {
   circos.rect(
      xleft = sector.dt[i, x1_n], xright = sector.dt[i, x2_n],
      ybottom = 0, ytop = 1, lwd = 0.5,
      track.index = tid_orfs, sector.index = sector.dt[i, exon_comb_short_rc],
      col = NA, border='black'
   )
   circos.rect(
      xleft = sector.dt[i, x1_n], xright = sector.dt[i, x2_n],
      ybottom = 0, ytop = 1, lwd = 0.5,
      track.index = tid_combs, sector.index = sector.dt[i, exon_comb_short_rc],
      col = NA, border='black'
   )
}

cat('\n')
cat('  - Generating legends...')
# Plot a legend on the left for exons
# legend(x=-1.28, y=1.1, title='', legend=names(exon_colors),
#        fill=as.character(exon_colors), y.intersp=0.8, cex=0.7, bty='n',
#        xjust=0, yjust=1, title.adj=0, adj=0)
   
# Plot a legend for reading frame analysis
# legend(x=-1.28, y=-1.0, title='', legend=c('out of frame', 'in frame', 'stop', 'non-translated'),
#        fill=as.character(fs_colors), y.intersp=0.8, cex=0.7, bty='n',
#        xjust=0, yjust=0, title.adj=0, adj=0)
cat(' OK.\n')
cat('  -->', output_circoplot, '\n')

# Generate a summary
# - total number of reads, for each exon combination
# - 
# Here we are ignoring reads with exon coverage <= 0.9 (90%) and exon_start_pos_align
# > 12 (first 4 codons) for reading frame analysis.
cat('* Generating summary...')
sum.dt = contigs.dt[read_id %in% reads_combs.dt[, unique(read_id)],
                    .SD[.N, .(type)], 
                    by=.(exon_comb_short_rc, read_id)]

sum_n.dt = sum.dt[, .N, by=.(exon_comb_short_rc, type)][order(exon_comb_short_rc, N)]

sum_n_d.dt = dcast(sum_n.dt, exon_comb_short_rc ~ paste0('type_', type), value.var = 'N', fill=0)
sum_n_d.dt[, total := type_0 + type_1 + type_2 + type_3]
sum_n_d.dt[, in_frame := type_1 + type_2]
sum_n_d.dt[, in_frame_pct := round(100 * (type_1 + type_2) / total, 1)]
sum_n_d_out.dt = sum_n_d.dt[, .(exon_combination = exon_comb_short_rc, in_frame, total, 
                                in_frame_pct)][order(-in_frame_pct, -total)]

# write.table(sum_n_d_out.dt, paste0('data/', input_base, '_reading_frame_analysis_summary.txt'),
#             sep='\t', quote=F, row.names=F)

# Now count how many of each end either inframe or with stop
# (these we count as 'in frame' reads)
cat(' OK.\n')

