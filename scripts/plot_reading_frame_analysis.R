

# ------------------------ Generate unique read table ------------------------
# For each read we can calculate the total number of variants based on part of
# the read that was aligned to reference exon sequences.

# First, some reads are clipped at beginning (exon_01) or end (exon_18). If the
# exon_start_pos_align > 1 and read_start_pos_align == 1, then just assume the
# first few nts exists like on the exon.
cat('Generating unique read table... ')
blast.dt[, c('read_seq', 'exon_seq') := NULL]

# Filter out reads that may be PCR chimeras.
# blast.dt[, unique(exon_comb_short_rc)]

exon_comb_short_rc_chimeras = c(
   '1-3,1,17-18', '1-3,1-3,17-18', '1-3,16-18,17-18',
   '1-5,16-17,1-3,17-18', '1,1-3,16-18', '1-3,1-2,17-18', '1-3,17-18,17-18',
   '1-3,16-18,2-3,16-18', '1-3,16-18,18,18', '1-3,17-18,18',
   "1-3,17-18,18',18'"
)


# We are left with 68 exon combinations
blast.dt = blast.dt[!(exon_comb_short_rc %in% exon_comb_short_rc_chimeras)]
all_exon_comb = blast.dt[, unique(exon_comb_short_rc)]
# Just single (no tandem) combinations:
# all_exon_comb[!(all_exon_comb %like% '18,1' | all_exon_comb %like% "18'")]

exon_comb_display_manual = c(
   '1-3,17-18', '1-3,16-18', '1-2,17-18'
)

# CLean-up few remaining issues with BLAST alignment
# blast.dt[exon_comb_short_rc=='1-3,16-18,18,18',
#          .SD[alignment_length==max(alignment_length)], by=exon_no]

# blast.dt[read_id=='read_01281' & exon_id=='exon_02' & is.na(read_start_pos_align),
#          read_start_pos_align := 1108]
# blast.dt[read_id=='read_01281' & exon_id=='exon_02' & is.na(read_end_pos_align),
#          read_end_pos_align := 1114]

blast.dt = merge(blast.dt, unique(contigs.dt[, .(read_id, stop_pos)]),
                 by='read_id')
blast.dt = blast.dt[order(read_id, read_start_pos_align)]

# For each exon combination in tandem, label the mid-points on the read as the
# average positions of the two joined tandems (rounded).
# blast.dt[exon_comb_short_rc %like% '18,',
#          tandem_midpoints := paste(.SD[exon_id %in% c('exon_18', 'exon_01_rc'),
#                                        read_end_pos_align],
#                                    collapse=','),
#          by=read_id]


# Ignore tandem repeats for this plot
#blast.dt = blast.dt[!(exon_comb_short_rc %like% '18,1' |
#                      exon_comb_short_rc %like% "18,18'") ]


cat('OK.\n')

# -------------- Reads combinations such that bin width more reflects -------------
#                            the number of unique reads
# "Coarse-grain" sectors with more than this number of reads
cat('Coarse-graining display sectors... ')
max_reads = 300

reads_combs.dt = blast.dt[
   !(exon_comb_short_rc %like% '18,1' | exon_comb_short_rc %like% "18,18'"),
   .(total_aln_len=sum(alignment_length), tandem_mids=unique(tandem_midpoints)),
   by=.(exon_comb_short_rc, read_id, read_count, stop_pos)]
reads_combs.dt = merge(
   reads_combs.dt,
   # reads.dt[, .(read_id, seq_len)], all.x=T, by='read_id'
   unique(blast.dt[, .(read_id, seq_len=read_length)]),
   all.x=T, by='read_id'
)
# reads_combs.dt[, exon_coverage := total_aln_len/seq_len]
# reads_combs.dt = reads_combs.dt[exon_coverage > 0.9]

# Sort first by average length of each sector, then by read seq length
reads_combs.dt[, exon_comb_avg_len := mean(seq_len), by=exon_comb_short_rc]
reads_combs.dt = reads_combs.dt[order(-exon_comb_avg_len, -stop_pos)]
reads_combs.dt[, x := 1:.N]
cat('OK.\n')

# ---------------------- Calculate sector -----------------------------------------
# Now generate sector lengths to initialize circos plot
# from the number of unique reads in each exon combination.
min_sec_size = 50
big_sector_n = 300
read_display_step = 3

cat('Calculating sector coordinates... ')
sector.dt = reads_combs.dt[, .(n=.N, avg_len=unique(exon_comb_avg_len),
                               counts=sum(read_count)), by=exon_comb_short_rc]
sector.dt[, x2 := cumsum(n)]
sector.dt[, x1 := c(1, x2[1:(length(x2)-1)])]

sector.dt[, x2_n := x2]
sector.dt[, x1_n := x1]

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

sector.dt[, log10_counts := log10(counts)]
sector.dt[, exon_comb_short_rc := factor(exon_comb_short_rc,
   levels = exon_comb_short_rc)]

# Now, for sectors with very many reads show only a sample to
# make the plot faster.
reads_combs.dt[exon_comb_short_rc %in% sector.dt[n>big_sector_n, exon_comb_short_rc],
               remove := rep(c(F, rep(T, read_display_step-1)), .N),
               by=exon_comb_short_rc]
reads_combs.dt[is.na(remove), remove := F]
reads_combs.dt[remove==F, xf := 1:.N]

# If we have less than 50 reads per exon_comb_short_rc, then just replicate
# them N times, so that the graph can be seen on the faceted plot.
read_n = 200
small_exon_combs = reads_combs.dt[remove==F, .(read_id, exon_comb_short_rc)][, .N,
   by=exon_comb_short_rc][N < read_n, exon_comb_short_rc]


# Merge reads_combs with contig information
contigs2.dt = merge(contigs.dt, reads_combs.dt[remove==F, .(read_id, x=xf)],
                   by='read_id')
contigs2.dt[, start_abs := start + read_start_pos_align - 1]
contigs2.dt[, end_abs := end + read_start_pos_align - 1]
cat('OK.\n')

# Colors for different exon combinations
# exon_comb_colors = c('black',   '#e41a1c', '#377eb8', '#4daf4a',
#                 '#984ea3', '#a65628', 'pink',    'orange')

exon_comb_colors = list(
   '1-3,17-18' = 'black',   '1-3,16-18' = '#e41a1c',
   '1-2,17-18' = '#377eb8', '1-5,16-18' = '#4daf4a',
   '1-2,18'    = '#984ea3', '1-7,9-18'  = '#a65628',
   '1-3,14-18' = 'pink',    '1,12-18'   = 'orange')

exon_comb_colors = rainbow(sector.dt[, .N])
names(exon_comb_colors) = sector.dt[, exon_comb_short_rc]

exon_colors = list(
    '1'='#a6cee3',  '2'='#1f78b4',  '3'='#b2df8a',  '4'='#33a02c',  '5'='#fb9a99',
    '6'='#d9d9d9',  '7'='#8dd3c7',  '9'='#ccebc5', '10'='#bebada', '11'='#fb8072',
   '12'='#e31a1c', '13'='#fdbf6f', '14'='#ff7f00', '15'='#cab2d6', '16'='#6a3d9a',
   '17'='#ffff99', '18'='#b15928'
)
exon_col_int = colorRampPalette(exon_colors)

shift = function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}


blast_sub.dt = merge(blast.dt,
                     reads_combs.dt[remove==F, .(read_id, exon_comb_avg_len)])
blast_sub.dt[, exon_no_fac := factor(exon_no, levels=sort(unique(exon_no)))]
# Exon present factor: 1 = full present, 2 = partially present
blast_sub.dt[, exon_present := 1]
# Allow ONE NT off at edges due to variants
blast_sub.dt[exon_length > exon_end_pos_align - exon_start_pos_align + 2 &
                !(exon_no %in% c(1, 18)), exon_present := 2]
# blast_sub.dt[, read_length := sum(sapply(read_seq_corrected, nchar)), by=read_id]
#blast_sub.dt = merge(blast_sub.dt, reads.dt[, .(read_id, read_length=seq_len)],
#                     by='read_id')
blast_sub.dt[, c('exon_seq_corrected', 'read_seq_corrected') := NULL]
blast_sub.dt[, total_exons := .N, by=read_id]

blast_sub.dt = blast_sub.dt[, {
   if (.BY[[1]] %in% small_exon_combs) {
      dts = list()
      nreads = .SD[, length(unique(read_id))]
      for (i in 1:floor(read_n/nreads)) {
         dt = copy(.SD)
         dt[, read_id := paste(read_id, '.', sprintf('%04d', i), sep='')]
         dts[[length(dts)+1]] = dt
      }
      rbindlist(dts)
   } else {
      .SD
   }
}, by=exon_comb_short_rc]


# Order reads by length descending
blast_sub.dt = blast_sub.dt[order(-exon_comb_avg_len, -stop_pos, read_id)]
blast_sub.dt[, read_id := factor(read_id,
   levels=blast_sub.dt[, unique(read_id)])]
blast_sub.dt[, exon_comb_short_rc := factor(exon_comb_short_rc,
   levels=blast_sub.dt[, unique(exon_comb_short_rc)])]
blast_sub.dt = blast_sub.dt[!(exon_comb_short_rc %like% '18,1')]




# We draw rectangles for each exon in each read
# If we don't have many read for that exon
exon_no_to_coord = function (n) {
  # Convert exon no n to coordinate
  # Each exon will have length 1.
  nn = ifelse(n > 7, n-1, n)
}

blast_sub.dt[, exon_w := {
  (exon_end_pos_align-exon_start_pos_align+1)/exon_length
}]
blast_sub.dt[, exon_x := {
  # This generates a whole list of ns for all exons/reads
  ns = ifelse(exon_no > 7, exon_no-1, exon_no)
  offset = (exon_start_pos_align-1)/exon_length
  offset + ns - (1-exon_w)*0.5
}]

# Add exon colors
exon_col_int2 = colorRampPalette(brewer.pal(n=8, name='RdYlBu'))
exon_colors = exon_col_int2(17)
# # exon_colors = shift(exon_colors, 8)
names(exon_colors) = c(1:7, 9:18)

# Correct this one due to variant in exon 2
# blast_sub.dt[exon_comb_short_rc=='1-3,17-18' & exon_no==2, exon_present := 1]


#  Plot
square.plot = ggplot(
   #blast_sub.dt[read_id %in% reads_combs.dt[remove==F, read_id] |
   #             gsub('\\.[0-9]+', '', read_id) %in% reads_combs.dt[, read_id]],
   blast_sub.dt,
   aes(x=exon_x, y=read_id, fill=factor(exon_no),
       color=factor(exon_no), width=exon_w, height=1)
   ) +
   geom_tile(show.legend=F) +
   geom_vline(
      xintercept=seq(0.5, 17.5),
      color='gray',
      size=0.1
      #colour='#FFFFFF'
   ) +
   facet_grid(exon_comb_short_rc ~ ., switch='y', space='free_y', scales='free_y') +
   scale_x_continuous(name='Exon', expand=c(0, 0), breaks=1:17,
                    labels=c(1:7, 9:18))+
   scale_y_discrete(name='', expand=c(0, 0)) +
   scale_color_manual(values=unlist(exon_colors)) +
   scale_fill_manual(values=unlist(exon_colors)) +
   theme(
      panel.background = element_rect(fill='white'),
      panel.border = element_rect(colour = "gray", fill=NA, size=0.1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(color='black'),
      panel.spacing = unit(0.0, "cm"),
      # strip.text.y = element_text(angle = 180)
      strip.text=element_blank()
   )
square.plot

ann.plot = ggplot(data=data.table(xmin=15.5, xmax=16.5, ymin=-0.1, ymax=100),
                  mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)) +
   geom_rect(
      inherit.aes = F, data=data.table(xmin=15.5, xmax=16.5, ymin=-0.1, ymax=100),
      mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
      color='red', fill=alpha('red', 0.30)
   ) +
   theme(
      panel.background = element_blank(),
      #panel.border = element_rect(colour = "gray", fill=NA, size=0.1),
      panel.border = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      panel.spacing = unit(0.0, "cm"),
      # strip.text.y = element_text(angle = 180)
      strip.text=element_blank()
   )
ann.plot

# square.plot2 = square.plot + theme_cowplot(line_size = 0)
ann2.plot = ann.plot + theme_void()

ggsave(paste0('plots/', input_base, '_exon_combinations_squares.pdf'),
       square.plot, width = unit(4, 'in'), height = unit(4.5, 'in'))

# Summarize unique and total counts for each exon combination
blast_sum.dt = blast_sub.dt[,
  .(unique_reads = length(unique(gsub('\\.[0-9]+', '', read_id)))),
  by=exon_comb_short_rc]


# Now show in more detail: 1,17-18 and 1-2,18
# square_2.plot = ggplot(
#    blast_sub.dt[exon_comb_short_rc %in% c('1,17-18', '1-2,18') &
#                 ((exon_no <= 4) | (exon_no >= 16))],
#    aes(x=exon_x, y=read_id, fill=factor(exon_no),
#        color=factor(exon_no), width=exon_w, height=1)
# ) +
#    geom_tile(show.legend=F) +
#    geom_vline(
#       xintercept=seq(0.5, 17.5),
#       color='gray',
#       size=0.1
#       #colour='#FFFFFF'
#    ) +
#    facet_grid(exon_comb_short_rc ~ ., space='free', scales='free') +
#    scale_x_continuous(name='Exon', expand=c(0, 0), breaks=c(1:4, 15:17),
#                       labels=c(1:4, 16:18)) +
#    scale_y_discrete(name='', expand=c(0, 0)) +
#    scale_color_manual(values=unlist(exon_colors)) +
#    scale_fill_manual(values=unlist(exon_colors)) +
#    theme(
#       panel.background = element_rect(fill='white'),
#       panel.border = element_rect(colour = "gray", fill=NA, size=0.1),
#       panel.grid.minor = element_blank(),
#       panel.grid.major = element_blank(),
#       axis.title.y = element_blank(),
#       axis.text.y = element_blank(),
#       axis.ticks.y = element_blank(),
#       axis.ticks.x = element_blank(),
#       axis.text.x = element_text(color='black'),
#       panel.spacing = unit(0.0, "cm"),
#       # strip.text.y = element_text(angle = 180)
#       strip.text=element_blank()
#    )
# square_2.plot
#
# square2_2.plot = ggdraw(square_2.plot) + draw_plot(ann2.plot, x=0.836,
# width=0.065, y=0.05,
#                                                height=0.978)
# save_plot(paste0('plots/', input_base, '_exon_combinations_squares_ab_2.pdf'),
# square2_2.plot,
#           base_width = 4, base_height=1.5)
#

# Generat contig coordinates in this weird reference frame
contigs2.dt = merge(contigs2.dt, exon_len.dt[, .(exon_no, exon_length)], by='exon_no')
contigs2.dt[, contig_w := {
   (end-start+1)/exon_length
}]
contigs2.dt[, contig_x := {
   ns = ifelse(exon_no > 7, exon_no-1, exon_no)
   offset = (start+exon_start_pos_align-1)/exon_length
   offset + ns - (1-contig_w)*0.5
}]
contigs2.dt = merge(contigs2.dt,
                    unique(reads_combs.dt[, .(read_id, exon_comb_avg_len)]),
                   by='read_id')



# contigs3.dt = contigs2.dt[, {
#    if (.BY[[1]] %in% small_exon_combs) {
#       dts = list()
#       for (i in 1:read_n) {
#          dt = copy(.SD)
#          dt[, read_id := paste(read_id, '.', i, sep='')]
#          dts[[length(dts)+1]] = dt
#       }
#       rbindlist(dts)
#    } else {
#       .SD
#    }
# }, by=exon_comb_short_rc]

contigs3.dt = contigs2.dt[, {
   if (.BY[[1]] %in% small_exon_combs) {
      dts = list()
      nreads = .SD[, length(unique(read_id))]
      for (i in 1:floor(read_n/nreads)) {
         dt = copy(.SD)
         dt[, read_id := paste(read_id, '.', sprintf('%04d', i), sep='')]
         dts[[length(dts)+1]] = dt
      }
      rbindlist(dts)
   } else {
      .SD
   }
}, by=exon_comb_short_rc]

contigs3.dt[, read_id := factor(read_id,
   levels=blast_sub.dt[, levels(read_id)])]
contigs3.dt[, exon_comb_short_rc := factor(exon_comb_short_rc,
   levels=blast_sub.dt[, levels(exon_comb_short_rc)])]

contigs3.dt = contigs3.dt[order(-exon_comb_avg_len, -read_len, read_id, exon_no,
                                start)]


# Contig types:
# 0 = translated differently
# 1 = translated
# 2 = STOP codon
# 3 = untranslated

# Need to fill the gaps in contigs(2).dt
# In an exon, if STOP is encountered, make sure we have another
# contig from STOP to the end of that exon alignment filled with
# 3s.

# Furthermore, to be able to see exon combinations with very small
# number of reads, copy/paste these reads n times so we can see
# it on display more clearly.

rf_colors = c('0'='#ffffbf', '1'='#1a9850', '2'='#000000', '3'='#fee08b')
grid_color = 'lightgray'
fill_color = 'white'

# Now reading frame analysis in the same fashion
square_rf.plot = ggplot(
   #contigs2.dt[(exon_comb_short_rc=='1-7,9-18')],
   #contigs2.dt[read_id=='read_00554'],
   contigs3.dt,
   aes(x=contig_x, y=read_id, fill=factor(type),
       color=factor(type), width=contig_w, height=1)
   ) +
   geom_tile(show.legend = F) +
   geom_vline(
      xintercept=seq(0.5, 17.5),
      color=grid_color,
      size=0.1
   ) +
   # geom_hline(yintercept=seq(1.5, 17.5), colour='#000000') +
   facet_grid(exon_comb_short_rc ~ ., switch='y', space='free_y', scales='free_y') +
   scale_x_continuous(name='Exon', expand=c(0, 0), breaks=1:17,
                      labels=c(1:7, 9:18))+
   scale_y_discrete(name='', expand=c(0, 0)) +
   scale_color_manual(values=rf_colors) +
   scale_fill_manual(values=rf_colors) +
   theme(
      panel.background = element_rect(fill=fill_color),
      panel.border = element_rect(colour = grid_color, fill=NA,
                                  size=0.1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(color='black'),
      panel.spacing = unit(0.0, "cm"),
      # strip.text.y = element_text(angle = 180)
      strip.text = element_blank()
   )
square_rf.plot

ggsave(paste0('plots/', input_base, '_reading_frames.pdf'),
       square_rf.plot, width = unit(4, 'in'), height = unit(4.5, 'in'))



# # reads with exon_combinations and counts
# reads_c.dt = unique(reads_combs.dt[, .(read_id, exon_comb_short_rc, read_count)])
# reads_c.dt[, exon_comb_count := .N, by=exon_comb_short_rc]
# reads_c.dt = reads_c.dt[!(exon_comb_short_rc %like% '18,1')]
# reads_c.dt = merge(
#    reads_c.dt,
#    unique(blast_sub.dt[, .(read_id, exon_comb_avg_len, read_length)]),
#    by='read_id')
# reads_c.dt[, exon_comb_short_rc := factor(exon_comb_short_rc,
#    levels=blast_sub.dt[, levels(exon_comb_short_rc)])]
#
# blast_sub.dt =  merge(blast_sub.dt, reads_c.dt[, .(read_id, exon_comb_count)],
#                       by='read_id')
#
# # Now also add read abundance for each read
# square_ab.plot = ggplot(
#    reads_c.dt[order(-exon_comb_avg_len, -read_length)],
#    aes(y=exon_comb_count, x=read_id)
#    ) +
#    geom_col(color='black', fill='black') +
#    scale_y_continuous(name='Unique reads', expand=c(0, 0),
#                       breaks=seq(0, 1000, 250)) +
#    scale_x_discrete(name='', expand=c(0, 0)) +
#    facet_grid(exon_comb_short_rc ~ ., switch='y', scales='free_y',
#               space='fixed') +
#    coord_flip() +
#    theme(
#       panel.background = element_rect(fill='white', color='#EBEBEB'),
#       panel.grid.minor.x = element_blank(),
#       panel.grid.major.x = element_blank(),
#       axis.title.y = element_blank(),
#       axis.text.y = element_blank(),
#       axis.ticks.y = element_blank(),
#       # axis.ticks.x = element_blank(),
#       axis.text.x = element_text(color='black'),
#       panel.spacing = unit(0.05, "cm"),
#       strip.text.y = element_text(angle = 180)
#    )
# square_ab.plot
# ggsave(paste0('plots/', input_base, '_exon_combinations_squares_counts.pdf'),
#        square_ab.plot, width = unit(3.5, 'in'), height = unit(3, 'in'))



# Generate a summary
# - total number of reads, for each exon combination
# -
# Here we are ignoring reads with exon coverage <= 0.9 (90%) and
# exon_start_pos_align
# > 12 (first 4 codons) for reading frame analysis.
cat('Generating summary... ')
sum.dt = contigs.dt[read_id %in% reads_combs.dt[, unique(read_id)],
                    .SD[.N, .(type)],
                    by=.(exon_comb_short_rc, read_id)]

sum_n.dt = sum.dt[, .N,
                  by=.(exon_comb_short_rc, type)][order(exon_comb_short_rc, N)]

sum_n_d.dt = dcast(sum_n.dt, exon_comb_short_rc ~ paste0('type_', type),
                   value.var = 'N', fill=0)
# Check for missing reading frame types
if (!('type_0' %in% names(sum_n_d.dt))) {
  sum_n_d.dt[, `type_0` := 0L]
}
if (!('type_1' %in% names(sum_n_d.dt))) {
  sum_n_d.dt[, `type_1` := 0L]
}
if (!('type_2' %in% names(sum_n_d.dt))) {
  sum_n_d.dt[, `type_2` := 0L]
}
if (!('type_3' %in% names(sum_n_d.dt))) {
  sum_n_d.dt[, `type_3` := 0L]
}


sum_n_d.dt[, total := type_0 + type_1 + type_2 + type_3]
sum_n_d.dt[, in_frame := type_1 + type_2]
sum_n_d.dt[, in_frame_pct := round(100 * (type_1 + type_2) / total, 1)]
sum_n_d_out.dt = sum_n_d.dt[, .(exon_combination = exon_comb_short_rc, in_frame,
                                total,
                                in_frame_pct)][order(-in_frame_pct, -total)]


exon_sum.dt = reads_combs.dt[, .SD[1], by=.(exon_comb_short_rc, read_id)]
exon_sum.dt = exon_sum.dt[,
  .(unique_reads = length(read_id),
    total_reads = sum(read_count)),
  by=exon_comb_short_rc
]

write.table(
  exon_sum.dt[, .(`Exon combination`=exon_comb_short_rc,
                  `Unique reads`=unique_reads,
                  `Total reads`=total_reads)],
  paste0('data/', input_base, '_exon_comb_summary.txt'),
  sep='\t', quote=F, row.names=F
)

cat('OK.\n')

#
# # Now count how many of each end either inframe or with stop
# # (these we count as 'in frame' reads)
#
# # Plot a legend on the left for exons
# legend(x=-1.28, y=1.1, title='', legend=names(exon_colors),
#        fill=as.character(exon_colors), y.intersp=0.8, cex=0.7, bty='n',
#        xjust=0, yjust=1, title.adj=0, adj=0)
#
# # Plot a legend for reading frame analysis
# legend(x=-1.28, y=-1.0, title='', legend=c('out of frame', 'in frame', 'stop',
#                                            'non-translated'),
#        fill=as.character(fs_colors), y.intersp=0.8, cex=0.7, bty='n',
#        xjust=0, yjust=0, title.adj=0, adj=0)
#
#

# Now generate annotation for each exon combination to be displayed
# outside of track 1 (lines + labels)

# Check this in illustrator first to see if it even fits.
cat('\n')


# Write these contigs to file
cat('Writing contigs table to file... ')
write.table(contigs.dt, paste0('data/', input_base, '_orf_contigs.txt'),
            sep='\t', quote=F, row.names=F)
write.table(contigs3.dt, paste0('data/', input_base, '_orf_contigs_display.txt'),
            sep='\t', quote=F, row.names=F)
cat('OK.\n')

cat('Writing blast sub-data for plotting... ')
write.table(blast_sub.dt, paste0('data/', input_base, 'blast_sub.txt'),
            sep='\t', quote=F, row.names=F)
cat('OK.\n')

cat('Writing reading frame analysis to file...')
write.table(sum_n_d_out.dt, paste0('data/', input_base,
                                   '_reading_frame_analysis_summary.txt'),
            sep='\t', quote=F, row.names=F)
cat('OK.\n')


# Further analyze only reads which have variants at know AD locations


codon_compare = function (v1, v2) {
   # First vector is always the vector of read seq codons.
   # Second vector are exon seq_codons.
   # Return a list based on matching 1 by 1 codon:
   #  v1 = c('A', 'ATG', 'CTG', 'CCC', 'GG')
   #  v2 = c('G', 'ATG', 'CTG', 'TGA', 'GG')
   # out: c(0, 1, 1, 1, 1, 1, 1, 2, 3)
   # 1 = in frame and equal, 2 = stop, 0 = in frame and not equal,
   # 3 = untranslated (after STOP)
   n = min(length(v1), length(v2))
   # This gives cummulative lengths for each codon
   v1_ind = as.numeric(cumsum(sapply(v1, nchar)))
   # Total number of characters in all codons is the last
   # entry in cummulative vector.
   v1_tot = v1_ind[length(v1_ind)]
   out = rep(0, v1_tot)
   for (i in 1:n) {
      if (i == 1) inds = c(1:v1_ind[1])
      else inds = c((v1_ind[i-1]+1):v1_ind[i])
      if (v1[i] %in% stop_codons) { # STOP codon
         out[inds] = 2
         if (inds[length(inds)] < v1_tot) {
            out[(inds[length(inds)]+1):length(out)] = 3
         }
         return(out)
      } else { # Non-stop
         if (v1[i] == v2[i]) out[inds] = 1
         else out[inds] = 0
      }
   }
   return(out)
}

codon_cmp_fix_stop = function (l) {
   # Input: a list from lapplying codon_compare to multiple
   # sequences.
   # Output. WHen the STOP codon occurs (2), replace every number
   # in reads downstream with 3 (untranslated region.)
   out = list()
   stopped = FALSE
   for (l_i in l) {
      if (stopped) {
         out[[length(out)+1]] = rep(3, length(l_i))
      } else {
         out[[length(out)+1]] = l_i
         if (2 %in% l_i) {
            stopped = TRUE
            next
         }
      }
   }
   return(out)
}

codon_mismatch_pos = function (c1, c2) {
   # Return the mismatch positions between two codons as a list
   c1_v = strsplit(c1, '')[[1]]
   c2_v = strsplit(c2, '')[[1]]
   out = list()
   for (i in 1:length(c1_v)) {
      if (c1_v[i] != c2_v[i]) {
         out[[length(out)+1]] = i
      }
   }
   return(out)
}

str_replace_nth = function (x, n, r) {
   # Replace n-th character in string x by a single character r
   if (n == 1) { # First
      return(paste0(r, substring(x, 2)))
   } else if (n == nchar(x)) { # Last
      return(paste0(substring(x, 1, nchar(x)-1), r))
   } else { # inside
      return(paste0(substring(x, 1, n-1), r, substring(x, n+1)))
   }
}


# Load locations that we need to check for SNVs
# Generate codon to aa reference table from codon.dt
# Do this once we have the reading frame analysis in unique_reads_comb_circo_v2.R.
aa_to_codon.dt = codon.dt[, .(codon = strsplit(codons, ', ')[[1]]), by=aa_1let]

# Calculate codon position mismatch between exon seq and known variant

if (file.exists(input_knownsim)) {
  cat('Loading known and similar AD variants...')
  advars.dt = fread(input_knownsim)
  cat(' OK.\n')
  cat('Processing AD variants...')
  advars.dt = advars.dt[!is.na(codon_new_ref)]
  advars.dt[, codon_mismatch_pos := as.integer(mapply(
     codon_mismatch_pos, codon_new_ref, codon_old))]

  # Generate a new codon
  advars.dt[, codon_new_data := mapply(str_replace_nth, codon_old,
                                       codon_mismatch_pos, mut_new)]

  advars.dt = merge(
     advars.dt,
     aa_to_codon.dt[, .(aa_1let_old=aa_1let, codon_old=codon)], by='codon_old')
  advars.dt = merge(
     advars.dt,
     aa_to_codon.dt[, .(aa_1let_new_ref=aa_1let, codon_new_ref=codon)],
     by='codon_new_ref')
  advars.dt = merge(
     advars.dt,
     aa_to_codon.dt[, .(aa_1let_new_data=aa_1let, codon_new_data=codon)],
     by='codon_new_data')
  advars.dt[, mutation := 'non-syn']
  advars.dt[aa_1let_new_ref == aa_1let_new_data, mutation := 'syn']
  advars.dt[codon_new_data==codon_new_ref, mutation := 'exact']

  # Now check for each read_id and exon_abs_pos
  advars.dt = merge(advars.dt, exon_len.dt[, .(exon_no, exon_start_abs)], by='exon_no')
  advars.dt = merge(advars.dt, unique(blast.dt[, .(read_id, exon_no,
                                                   exon_start_pos_align)]),
                    by=c('read_id', 'exon_no'))
  advars.dt[, exon_rel_pos_trunc := exon_abs_pos - exon_start_abs + 1 -
               exon_start_pos_align + 1]
  advars.dt[, exon_rel_pos := exon_abs_pos - exon_start_abs + 1]
  advars.dt[, codon_new_data_region := -1]
  advars.dt[, codon_new_data_region_len := -1]
  cat(' OK.\n')

  for (i in 1:nrow(advars.dt)) {
     # cat(i)
     codon_region = contigs.dt[
        read_id==advars.dt[i, read_id] & exon_no==advars.dt[i, exon_no] &
        start <= advars.dt[i, exon_rel_pos_trunc] &
        end >= advars.dt[i, exon_rel_pos_trunc],
        type]
     codon_region_len = contigs.dt[
        read_id == advars.dt[i, read_id] & exon_no==advars.dt[i, exon_no] &
        start <= advars.dt[i, exon_rel_pos_trunc] &
        end >= advars.dt[i, exon_rel_pos_trunc],
        end-start+1]

     if (length(codon_region) > 0) {
        advars.dt[i, codon_new_data_region := codon_region]
        advars.dt[i, codon_new_data_region_len := codon_region_len]
     }
  }
  advars.dt[, codon_new_data_orf := 'transl. in different rf']
  advars.dt[codon_new_data_region==1 & codon_new_data_region_len > 3,
            codon_new_data_orf := 'transl. in correct rf']
  advars.dt[codon_new_data_region==0 & codon_new_data_region_len == 3,
            codon_new_data_orf := 'transl. in correct rf']
  advars.dt[codon_new_data_region==3, codon_new_data_orf := 'untranslated']

  # Add more annotations from the original file
  fad.dt = fread(input_knownfad)

  advars.dt = merge(
     advars.dt,
     fad.dt[, .(Mutation_2, `Primary  Papers`, `Mutation Type / Codon Change`)],
     all.x=T, by=c('Primary  Papers', 'Mutation Type / Codon Change'))

  advars.dt[codon_new_data_orf %in% c('transl. in different rf', 'untranslated') &
               mutation != 'exact',
            c('mutation', 'codon_new_data', 'aa_1let_new_data') := c(NA, NA, NA)]
  advars.dt[codon_new_data_orf %in% c('transl. in different rf', 'untranslated'),
            c('codon_new_data', 'aa_1let_new_data') := c(NA, NA)]

  # Add old nucleotides that were mutated
  exons_fa = readDNAStringSet(input_exons_fa)
  exons_fa.dt = data.table(exon_no=1:18, exon_seq=as.character(exons_fa))

  for (i in 1:nrow(advars.dt)) {
     exon = advars.dt[i, exon_no]
     exon_pos = advars.dt[i, exon_rel_pos]
     advars.dt[i, mut_old := exons_fa.dt[
        exon_no==exon,
        substring(exon_seq, exon_pos, exon_pos)]]
  }

  # Add full read sequences
  advars.dt = merge(advars.dt, reads.dt[, .(read_id, seq)], by='read_id')
  advars.dt = advars.dt[order(exon_comb_short_rc, `Mutation Type / Codon Change`,
                              codon_new_data_orf, mutation)]

  # Remove reads
  advars.dt = advars.dt[!is.na(codon_new_data)]
  # Write table
  write.table(advars.dt, paste0('data/', input_base, '_fad_and_similar.txt'),
              sep='\t', quote=F, row.names=F)

  # Generate exon combination unique and total read summary
  uniques.dt = blast.dt[, .SD[1], by=read_id][, .(uniques=.N, totals=sum(read_count)),
                                              by=exon_comb_short_rc]
  # Write this as a reference
  write.table(uniques.dt[, .(`Exon Combination`=exon_comb_short_rc,
                             `Uniques`=uniques, `Total`=totals)][order(-Uniques)],
              paste0('data/', input_base, '_uniques_totals.txt'),
              sep='\t', quote=F, row.names=F)

  # Write text table for Excel
  write.table(
     advars.dt[, .(`Exon Combination`=exon_comb_short_rc,
                   `Translation Type`=codon_new_data_orf,
                   `Mutation Type`=mutation, Type=type, Exon=exon_no,
                   `Position on gencDNA`=exon_abs_pos,
                   `Nucleotide in Exon`=mut_old,
                   `Nucleotide in Ref`=mut_new_ad,
                   `Nucleotide in Data`=mut_new,
                   `Codon in Exon`=codon_old,
                   `Codon in Ref`=codon_new_ref,
                   `Codon in Data`=codon_new_data,
                   `AA in Exon`=aa_1let_old,
                   `AA in Ref`=aa_1let_new_ref,
                   `AA in Data`=aa_1let_new_data,
                   `Clinical Phenotype`, Mutation=Mutation_2,
                   `Primary  Papers`, `Mutation Type / Codon Change`,
                   Read=read_id, `Read Count`=read_count,
                   `Read Sequence`=seq)],
     paste0('data/', input_base, '_fad_and_similar_excel.txt'),
     sep='\t', quote=F, row.names=F
  )

}


cat('Generating SNV summary...')
# Generate a summary for each exon combination
snv_c.dt = unique(
   snv_indel.dt[, .(
      read_id,
      read_count),
      by=exon_comb_short_rc])
snv_counts.dt = snv_c.dt[, .(unique_reads=length(unique(read_id)),
                             total_reads=sum(read_count)), by=exon_comb_short_rc]

# Now find reading frame variants for each exon_comb_short_rc
frame_counts.dt = contigs.dt[exon_no==18 & ((type==2 & end >= 100) |
   (type==1 & end >=97)), .(unique_inframe=length(unique(read_id)),
                            read_id=unique(read_id)), by=exon_comb_short_rc]

# Abeta region here is defined as between 2014 2139 absolute coordinates on the
# reference 1-7,9-18 variant. (It spans parts of Exons 16 and 17.)
ab_coord = c(2014L, 2139L)
blast_ab.dt = merge(
   blast.dt[, .(read_id, exon_no, exon_start_pos_align, exon_end_pos_align,
                exon_comb_short_rc, read_count)],
   exon_len.dt[, .(exon_no, exon_start_abs, exon_end_abs)],
   by='exon_no'
)
blast_ab.dt[, exon_start_pos_align_abs := exon_start_pos_align + exon_start_abs - 1L]
blast_ab.dt[, exon_end_pos_align_abs := exon_end_pos_align + exon_start_abs - 1L]
# First calculate number of nucleotides from exon 16's Ab region we have in each read
# then number of nucleotides from exon 17s Ab region.
# Add them together, then see if this number is >= 50% of ab_coord[2]-ab_coord[1].
blast_ab.dt[, ab_fraction := NULL]
blast_ab.dt[, ab_fraction := {
   ab_region = c(0L, 0L)
   # Do we have exon 16?
   if (16 %in% exon_no) {
      # Does the end of exon 16 alignment reach inside ab region?
      if (.SD[exon_no==16, exon_end_pos_align_abs] >= ab_coord[1]) {
         ab_region[1] = .SD[exon_no==16, exon_end_pos_align_abs] - ab_coord[1] + 1L
      }
      # Does the exon 16 alignment start inside Abeta region? If so, subtract
      # the missing part.
      if (.SD[exon_no==16, exon_start_pos_align_abs] > ab_coord[1]) {
         ab_region[1] = ab_region[1] - (exon_start_pos_align_abs - ab_coord[1])
      }
   }
   if (17 %in% exon_no) {
      # Does the start of exon 17 alignment reach inside ab region?
      if (.SD[exon_no==17, exon_start_pos_align_abs] <= ab_coord[2]) {
         ab_region[2] = ab_coord[2] - .SD[exon_no==17, exon_start_pos_align_abs] + 1L
      }
      # Does the exon 17 alignment end inside Abeta region? If so, subtract
      # the missing part.
      if (.SD[exon_no==17, exon_end_pos_align_abs] < ab_coord[2]) {
         ab_region[2] = ab_region[2] - (ab_coord[2] - exon_end_pos_align_abs)
      }
   }
   sum(ab_region)/(ab_coord[2] - ab_coord[1] + 1L)
}, by=read_id]
# blast_ab.dt = blast_ab.dt[ab_fraction >= 0.5]
blast_ab_count.dt = unique(blast_ab.dt[, .(read_id, read_count, ab_fraction,
                                           exon_comb_short_rc)])
blast_ab_count.dt[, ab := FALSE]
blast_ab_count.dt[ab_fraction >= 0.5, ab := TRUE]
blast_ab_count_sum.dt = blast_ab_count.dt[,
   .(unique_reads_with_ab=length(unique(.SD[ab_fraction >= 0.5, read_id]))),
   by=.(exon_comb_short_rc)]

frame_counts_sum.dt = merge(
   frame_counts.dt,
   blast_ab_count.dt,
   all=T,
   by=c('read_id', 'exon_comb_short_rc')
)
frame_counts_sum.dt[, inframe := TRUE]
frame_counts_sum.dt[is.na(unique_inframe), inframe := FALSE]

# Don't count variants using SNV counts.
counts.dt = merge(
   reads_combs.dt[, .(total_reads = sum(read_count), unique_reads=.N),
                  by=exon_comb_short_rc],
   frame_counts_sum.dt[ab_fraction >= 0.5,
                       .(total_has_ab=sum(read_count),
                         total_has_ab_and_inframe=.SD[inframe==T, sum(read_count)],
                         unique_has_ab=length(unique(read_id)),
                         unique_has_ab_and_inframe=.SD[inframe==T, length(unique(read_id))]),
                       by=exon_comb_short_rc],
   all=T,
   by='exon_comb_short_rc'
)
cat('OK.\n')

cat('Calculating unique exon joins...')
counts.dt = merge(
   counts.dt,
   frame_counts_sum.dt[inframe==T, .(total_inframe = sum(read_count),
                                     unique_inframe = length(unique(read_id))),
                       by=exon_comb_short_rc],
   all=T,
   by='exon_comb_short_rc'
)
for (colname in names(counts.dt)) {
   counts.dt[is.na(get(colname)), colname := 0, with=FALSE]
}

# Now find number of unique junctions
# For each read keep only left and right exons in the intra-exonic join
blast.dt = blast.dt[!(exon_comb_short_rc %like% "-1',17" |
                      exon_comb_short_rc  %like% '-18,17-')]
exon_join.dt = blast.dt[, {
   exon_l = as.integer(gsub('.*([0-9]+),.*', '\\1', .BY[[3]]))
   exon_r = as.integer(gsub('.*,([0-9]+).*', '\\1', .BY[[3]]))
   aln_overlap = .SD[exon_no==exon_l, read_end_pos_align] -
                 .SD[exon_no==exon_r, read_start_pos_align] + 1L
   #cat(.BY[[1]], '\n')
   if (length(aln_overlap) == 0) aln_overlap = as.integer(NA)
   else {
      if (aln_overlap < 0) aln_overlap = 0L
   }
   # Scan aln_overlap + 1 nt on each side for variants
   .(exon_l=exon_l, exon_r=exon_r, aln_overlap=aln_overlap,
     join_seq=paste0(
        str_sub(.SD[exon_no==exon_l, read_seq_corrected], start=-1L-aln_overlap),
        str_sub(.SD[exon_no==exon_r, read_seq_corrected], end=1L)
     )
   )
}, by=.(read_id, read_count, exon_comb_short_rc)]
exon_join.dt = exon_join.dt[!is.na(aln_overlap)]
cat('OK.\n')


# Write exon joins to table
cat('Writing intra-exon joins to table...')
write.table(
   exon_join.dt,
   paste0('data/', input_base, '_intra-exon_joins.txt'),
   sep='\t', quote=F, row.names=F
)
cat('OK.\n')

# Now we can count unique joins
join_count.dt = exon_join.dt[, .(unique_joins=length(unique(join_seq))),
                             by=exon_comb_short_rc]

# Merge with counts.dt
counts.dt = merge(counts.dt, join_count.dt, by='exon_comb_short_rc')

counts.dt = counts.dt[match(contigs3.dt[, levels(exon_comb_short_rc)],
                            exon_comb_short_rc)]


# Write counts of everything to a table
cat('Writing counts table to file... ')
write.table(counts.dt[, .(`Exon combination`=exon_comb_short_rc,
                          `Unique intra-exon junctions`=unique_joins,
                          `Total reads`=total_reads,
                          `Total reads with Ab region`=total_has_ab,
                          `Total reads in-frame`=total_inframe,
                          `Total reads with Ab region and in-frame`=
                             total_has_ab_and_inframe,
                          `Unique reads`=unique_reads,
                          `Unique reads with Ab region`=unique_has_ab,
                          `Unique reads in-frame`=unique_inframe,
                          `Unique reads with Ab region and in-frame`=
                             unique_has_ab_and_inframe
                          )
                      ],
            paste0('data/', input_base, '_counts_per_exoncomb.txt'),
            sep='\t', quote=F, row.names=F)
cat('OK.\n')