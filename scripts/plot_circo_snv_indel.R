# Load required libraries
cat('* Loading R libraries...')
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(RColorBrewer))
cat(' OK.\n')

##----------------------- Parse input arguments ----------------------------------
option_list = list(
  make_option(c('-i', '--input-exon-table'),
              help = 'Output from blast_reads. *_exon_table or _exon_table_hmpfix.'),
  make_option(c('-s', '--input-snvs-indels'),
              help = 'Output from runs_snvs_indels_analysis.R'),
  make_option(c('-l', '--input-exon-lengths'),
              help = 'Input file with exon lengths. Generate in bash then import here.'),
  make_option(c('-f', '--input-fads'), default = ' ',
              help = 'FAD variant matches, output from run_snvs_indels_analysis.R'),
  make_option(c('-o', '--output-circoplot'),
              help = 'Output filename for the circo plot.'),
  make_option(c('--plot-abeta-region'),
              default = 1,
              help = 'Plot red frame / overlay around Abeta region.')
)

opt = parse_args(OptionParser(option_list=option_list))

# Input
blast_fixed = opt$`input-exon-table`
snv_indel_fix = opt$`input-snvs-indels`
exonlen = opt$`input-exon-lengths`
fads = opt$`input-fads`
output_circoplot = opt$`output-circoplot`
plot_abeta = opt$`plot-abeta-region`

##-------------------- Analysis and plotting parameters -------------------------
annotation = c('insertion', 'deletion', 'SNV', '-')
nts = c('A', 'C', 'G', 'T', '-')
count_ivals = c(1,3,10,30,100,300)
no_bins = 5
count_colors = rev(brewer.pal(no_bins, 'RdYlBu'))
link_colors = c(brewer.pal(12, 'Paired'), 'black')
link_colors[11] = 'gray'
# The number of link colors limits how many unique exon combinations we
# can plot. If we have more uniques than colors, some will not be plotted.

##-------------------------- Load files -----------------------------------------
cat('* Loading files')
exonlen.dt = fread(exonlen)
exonlen.dt[, exon_no := as.integer(
   sub('_rc', '', sub('exon_(.*).*', '\\1', exon_id)))]
cat('.')
snv_indel_fix.dt = fread(snv_indel_fix)
cat('.')
blast_fixed.dt = fread(blast_fixed)
if (fads != ' ') {
  var_match.dt = fread(fads)
  cat('.')
}
cat(' OK.\n')

# Calculate absolute exon coordinates
exon_coords.dt = exonlen.dt[order(exon_id)][!(exon_id %like% '_rc')]
exon_coords.dt[, exon_start_abs := c(0, cumsum(exon_length)[1:(.N-1)])+1]
exon_coords.dt[, exon_end_abs := exon_start_abs + exon_length - 1]
exon_coords.dt[, exon_no := as.integer(gsub('_rc', '', gsub('exon_(.*).*','\\1', exon_id)))]

if (fads != ' ') {
  var_matches.dt = var_match.dt[,
    .(count_unique=.N), by=.(exon_no, exon_abs_pos)]
}


#---------------------- Generate separate tables for plotting --------------------

# SNVs
# Let's manually generate a frequency table for each exon and for each exon
# combination.
# First, for each exon, fill missing values?
snv_ft.dt = snv_indel_fix.dt[type=='SNV', .(count_unique=.N,
                                        sum_unique=sum(read_count)),
                         by=c('exon_no', 'exon_pos')]
snv_ft.dt = merge(snv_ft.dt, exonlen.dt[!(exon_id %like% 'rc'),
   .(exon_no, exon_length)], all.x=T, by='exon_no')
snv_ft.dt = merge(snv_ft.dt, exon_coords.dt[,
   .(exon_no, exon_start_abs, exon_end_abs)], all.x=T, by='exon_no')
# For reverse-complemented sequences recalculate position on the + straind
# snv_ft.dt[exon_id %like% 'rc', exon_pos := exon_length - exon_pos]
snv_ft.dt[, exon_pos_abs := exon_start_abs + exon_pos - 1]

# Indels
ins_ft.dt = snv_indel_fix.dt[type=='insertion',
   .(count_unique=.N, sum_unique=sum(read_count)),
                         by=c('exon_no', 'exon_pos')]
ins_ft.dt = merge(ins_ft.dt, exonlen.dt[!(exon_id %like% 'rc'),
   .(exon_no, exon_length)], all.x=T, by='exon_no')
ins_ft.dt = merge(ins_ft.dt, exon_coords.dt[,
   .(exon_no, exon_start_abs, exon_end_abs)], all.x=T, by='exon_no')
ins_ft.dt[, exon_pos_abs := exon_start_abs + exon_pos - 1]

del_ft.dt = snv_indel_fix.dt[, .(count_unique=.N, sum_unique=sum(read_count)),
                         by=c('exon_no', 'exon_pos')]
del_ft.dt = merge(del_ft.dt, exonlen.dt[!(exon_id %like% 'rc'),
   .(exon_no, exon_length)], all.x=T, by='exon_no')
del_ft.dt = merge(del_ft.dt, exon_coords.dt[,
   .(exon_no, exon_start_abs, exon_end_abs)], all.x=T, by='exon_no')
del_ft.dt[, exon_pos_abs := exon_start_abs + exon_pos - 1]

del_nopresuf_ft.dt = snv_indel_fix.dt[type=='deletion',
   .(count_unique=.N, sum_unique=sum(read_count)),
   by=c('exon_no', 'exon_pos')]
del_nopresuf_ft.dt = merge(del_nopresuf_ft.dt,
   exonlen.dt[!(exon_id %like% 'rc'), .(exon_no, exon_length)],
   all.x=T, by='exon_no')
del_nopresuf_ft.dt = merge(del_nopresuf_ft.dt, exon_coords.dt[,
   .(exon_no, exon_start_abs, exon_end_abs)], all.x=T, by='exon_no')
del_nopresuf_ft.dt[, exon_pos_abs := exon_start_abs + exon_pos - 1]


cat('* Combining SNVs and indels into a single table...')
# Recombine snvs, indels and dels into one data table
snv_indel_ft.dt = rbindlist(list(snv_ft.dt, ins_ft.dt, del_ft.dt))
snv_indel_ft.dt[, type := c(rep(annotation[3], nrow(snv_ft.dt)),
                            rep(annotation[1], nrow(ins_ft.dt)),
                            rep(annotation[2], nrow(del_ft.dt)))]

snv_indel_ft_abs.dt = dcast(
   snv_indel_ft.dt,
   exon_no + exon_pos_abs ~ paste('count', type, sep='_'),
   value.var='count_unique',
   fill=0
)
snv_indel_ft_abs.dt = merge(
   snv_indel_ft_abs.dt,
   exon_coords.dt[, .(exon_no, exon_start_abs, exon_end_abs)],
   by='exon_no')

# This table is similar, except that it doesn't have prefix and suffix deletions
# near intra-exon joins.
snv_indel_ft2.dt = rbindlist(list(snv_ft.dt, ins_ft.dt, del_nopresuf_ft.dt))
snv_indel_ft2.dt[, type := c(rep(annotation[3], nrow(snv_ft.dt)),
                            rep(annotation[1], nrow(ins_ft.dt)),
                            rep(annotation[2], nrow(del_nopresuf_ft.dt)))]

snv_indel_ft_abs2.dt = dcast(
   snv_indel_ft2.dt,
   exon_no + exon_pos_abs ~ paste('count', type, sep='_'),
   value.var='count_unique',
   fill=0
)
snv_indel_ft_abs2.dt = merge(
   snv_indel_ft_abs2.dt,
   exon_coords.dt[, .(exon_no, exon_start_abs, exon_end_abs)],
   by='exon_no')


if (!('count_deletion' %in% names(snv_indel_ft_abs2.dt))) {
  snv_indel_ft_abs2.dt[, count_deletion := 0L]
}
if (!('count_insertion' %in% names(snv_indel_ft_abs2.dt))) {
  snv_indel_ft_abs2.dt[, count_insertion := 0L]
}

snv_indel_ft_abs2.dt[, count_total := count_SNV + count_deletion + count_insertion]
snv_indel_ft_abs2.dt[, count_total_log10 := log10(count_total)]
count_range_log10 = snv_indel_ft_abs2.dt[, range(count_total_log10)]
# count_range_log10[2] = count_range_log10[2] + 0.1
# count_ivals_log10 = seq(count_range_log10[1], count_range_log10[2],
# length.out = no_bins)
count_ivals_log10_o = log10(count_ivals)
count_ivals_log10 = list(count_ivals_log10_o[1:(length(count_ivals_log10_o)-1)],
                         count_ivals_log10_o[2:length(count_ivals_log10_o)])
count_ivals_log10 = transpose(count_ivals_log10)

for (i in 1:length(count_ivals_log10)) {
   snv_indel_ft_abs2.dt[count_total_log10 >= count_ivals_log10[[i]][1] &
                       count_total_log10 <= count_ivals_log10[[i]][2],
                       count_color := count_colors[i]]
}
# We can cap the colors at some threshold (3000) and then everything above that
# threshold will have the same color.
snv_indel_ft_abs2.dt[
   count_total_log10 > count_ivals_log10[[length(count_ivals_log10)]][2],
   count_color := count_colors[length(count_ivals_log10)]]

cat(' OK.\n')


cat('* Finding all intra-exon joins...')
# Find all exon joins from the blast_fixed table.
blast_fixed.dt[, exon_no := as.integer(gsub('_rc', '', gsub('exon_(.*).*',
                                                            '\\1', exon_id)))]
blast_fixed.dt = blast_fixed.dt[order(read_id, read_start_pos_align)]

exon_joins.dt = blast_fixed.dt[, {
   # ex = exon_no[1]
   j = 1
   joins = list()
   for (i in 2:.N) {
      if (!(exon_no[i] %in% c(exon_no[i-1]+1, exon_no[i-1]-1))) {
         # This exon is not the prev +- 1. Therefore this is the intra-exon join.
         if (exon_id[i-1] %like% 'rc') {
            exon_l_rc = 1
            exon_l_pos = exon_length[i-1] - exon_end_pos_align[i-1] + 1L
         } else {
            exon_l_rc = 0
            exon_l_pos = exon_end_pos_align[i-1]
         }
         if (exon_id[i] %like% 'rc') {
            exon_r_rc = 1
            exon_r_pos = exon_length[i] - exon_start_pos_align[i] + 1L
         } else {
            exon_r_rc = 0
            exon_r_pos = exon_start_pos_align[i]
         }
         joins[[j]] = list(
            exon_l = exon_no[i-1], exon_r = exon_no[i],
            exon_l_rc = exon_l_rc, exon_r_rc = exon_r_rc,
            exon_l_pos = exon_l_pos, exon_r_pos = exon_r_pos,
            read_count = read_count[1],
            exon_comb_short_rc = exon_comb_short_rc[1])
         j = j + 1
      }
   }
   rbindlist(joins)
}, by=read_id]

# Calculate absolute exon positions. For this we need absolute positions for
# the position of each exon on the full isoform 1-18, which is in exon_coords.dt
exon_coords.dt[, exon_l := as.integer(lapply(exon_id, function (e)
   strsplit(e, '_')[[1]][2]))]
exon_coords.dt[, exon_r := as.integer(lapply(exon_id, function (e)
   strsplit(e, '_')[[1]][2]))]

exon_joins.dt = merge(
   exon_joins.dt,
   exon_coords.dt[,
      .(exon_l = exon_no, exon_start_abs_l = exon_start_abs,
        exon_end_abs_l = exon_end_abs)], by='exon_l')
exon_joins.dt = merge(
   exon_joins.dt,
   exon_coords.dt[,
      .(exon_r = exon_no, exon_start_abs_r = exon_start_abs,
        exon_end_abs_r = exon_end_abs)], by='exon_r')
exon_joins.dt[, exon_l_pos_abs := exon_start_abs_l + exon_l_pos - 1]
exon_joins.dt[, exon_r_pos_abs := exon_start_abs_r + exon_r_pos - 1]

# Now extract only unique joins, since that's what we can plot in the
# circo plot
exon_join_u.dt = exon_joins.dt[,
   .(join_count = sum(read_count), join_uniques = .N),
   by=c('exon_comb_short_rc', 'exon_l', 'exon_r', 'exon_l_pos_abs',
        'exon_r_pos_abs')]

cat(' OK.\n')




cat('* Generating circo plot...')
#--------------------------------- CIRCO PLOT -----------------------------------
#
# Required data tables:
# exon_coords.dt (absolute coordinates for each exon)
# exon_joins.dt  (for each read information about intra-exonic join(s))
# snv_ft.dt      (frequency table of SNVs for each genomic position and exon)
# ins_ft.dt      (same for first letter of insertions)
# del_nopresuf_ft.dt (same for deletions, excluding major missing contig pieces)
# del_presuf_ft.dt   (missing contig pieces of deletions)
# exon_join_u.dt     (unique exon joins)

# Backgrounds for SNVs, deletions and insertions.
# circo_color = alpha(rev(brewer.pal(3, 'Blues')), 0.5)

snv_indel_ft_abs_offset = 0
# Need 8 colors

y_lab_size = 0.5

cairo_pdf(output_circoplot, width=5, height=5)

# ------- prototyping
circo_color = rev(brewer.pal(7, 'Blues')[1:3])
circos.clear()
circos.par("start.degree" = 82, "gap.degree" = c(rep(1, 17), 19),
            points.overflow.warning=F)
circos.genomicInitialize(exon_coords.dt[, .(exon_no, exon_start_abs, exon_end_abs)],
                         plotType = 'labels')
circos.genomicTrackPlotRegion(
   ylim = c(0, 1),
   bg.col = rep(count_colors[1], 18),
   track.height = 0.07
)

# ------------ Plot total changes on each genomic position excluding big dels
for (i in 1:nrow(snv_indel_ft_abs2.dt[])) {
   row = snv_indel_ft_abs2.dt[i]
   if (row[, exon_pos_abs > exon_start_abs + snv_indel_ft_abs_offset &
             exon_pos_abs < exon_end_abs - snv_indel_ft_abs_offset]) {
      # Skip the points at exon boundaries for visibility
      circos.lines(x=c(row[, rep(exon_pos_abs, 2)]), y=c(-0.2, 1.2),
                   col=row[, count_color],
                   track.index=2, sector.index=row[, exon_no])
   }
}

# ------------  Plot SNVs
circos.genomicTrackPlotRegion(
   bg.col = rep(c(circo_color[1]), length.out = 25),
   bg.border = NA, track.height = 0.15,
   ylim = c(0, 3)
)
# for (i in c(0.3, 1.65, 3, 4.35)) {
#    circos.lines(c(-200, 2350), c(i, i), col=alpha('white', 0.5))
# }
for (i in c(0:3)) {
   circos.lines(c(-200, 2350), c(i, i), col=alpha('white', 0.5),
                track.index = 3)
}

# Add annotated SNVs
if (fads != ' ') {
  if (nrow(var_match.dt) > 0) {
     circos.genomicTrackPlotRegion(
        var_matches.dt[, .(exon_no, exon_abs_pos, exon_abs_pos2 = exon_abs_pos,
                           log10(count_unique))],
        track.index=3,
        panel.fun = function (region, value, ...) {
           circos.genomicPoints(region, value, ..., pch=21, cex=0.6,
                                bg=count_colors[length(count_colors)],
                                col=NA
           )
        },
        ylim = c(0, 3)
     )
  }
}

circos.genomicTrackPlotRegion(
   snv_ft.dt[, .(count_unique = log10(as.numeric(sum(count_unique)))),
             by=c('exon_no', 'exon_pos_abs', 'exon_pos_abs')],
   track.index = 3,
   panel.fun = function (region, value, ...) {
      circos.genomicPoints(region, value, ..., pch=16, cex=0.25)
   },
   ylim = c(0, 3)
)

circos.yaxis(at=c(0, 1, 2, 3), sector.index=1, track.index=3, labels.cex=y_lab_size,
             labels = c('1', '10', '100', '1,000'),
             labels.niceFacing = F)



# ------------ Plot insertions
circos.genomicTrackPlotRegion(
   bg.col = rep(c(circo_color[2]), length.out = 18),
   track.height = 0.15,
   ylim = c(0, 2)
)
for (i in c(0:2)) {
   circos.lines(c(0, 2413), c(i, i), col=alpha('white', 0.5),
                track.index=4, sector.index='1')
}

circos.genomicTrackPlotRegion(
   ins_ft.dt[, .(count_unique = log10(as.numeric(sum(count_unique)))),
             by=c('exon_no', 'exon_pos_abs', 'exon_pos_abs')],
   track.index = 4,
   panel.fun = function (region, value, ...) {
      circos.genomicPoints(region, value, ..., pch=16, cex=0.25)
   },
   ylim = c(0, 2)
)

circos.yaxis(at = 0:2, sector.index = 1, track.index = 4, labels.cex=y_lab_size,
             labels = c(1, 10, 100),
             labels.niceFacing = T)

# -------------- Plot deletions

del_range = c(0, 2)

circos.genomicTrackPlotRegion(
   bg.col = rep(c(circo_color[3]), length.out = 18),
   track.height = 0.15,
   ylim = del_range
)
for (i in c(del_range[1]:del_range[2])) {
   circos.lines(c(0, 2413), c(i, i), col=alpha('gray', 0.5),
                track.index=5, sector.index='1')
}
# First the actual deletions
circos.genomicTrackPlotRegion(
   del_nopresuf_ft.dt[, .(count_unique = log10(as.numeric(sum(count_unique)))),
             by=c('exon_no', 'exon_pos_abs', 'exon_pos_abs')],
   track.index = 5,
   panel.fun = function (region, value, ...) {
      circos.genomicPoints(region, value, ..., pch=16, cex=0.25)
   },
   ylim = del_range
)

circos.yaxis(at = c(0:2), sector.index=1, track.index=5,
             labels.cex=y_lab_size,
             labels=c(1, 10, 100),
             labels.niceFacing = T)


# ----------- Add links between intra-exonic junctions
# Display only top two junctions by the total (sum) read count
# exon_join_u_top.dt = exon_join_u.dt[, .SD[order(-join_count)][1:min(2, .N)],
# by=exon_comb_short_rc]
exon_join_u_top.dt = exon_join_u.dt[, .SD[order(-join_count)],
                                    by=exon_comb_short_rc]

exon_join_u_top.dt = exon_join_u_top.dt[!(exon_comb_short_rc %like% '18,[0-9]+')]
exon_join_u_top.dt = exon_join_u_top.dt[!(exon_comb_short_rc %like% ',1[,|-]')]
exon_comb_display = exon_join_u_top.dt[, unique(exon_comb_short_rc)]

for (i in 1:length(link_colors)) {
   circos.genomicLink(
      exon_join_u_top.dt[exon_comb_short_rc == exon_comb_display[i],
         .(exon_l, exon_l_pos_abs, exon_l_pos_abs)],
      exon_join_u_top.dt[exon_comb_short_rc == exon_comb_display[i],
         .(exon_r, exon_r_pos_abs, exon_r_pos_abs)],
                      col=alpha(link_colors[i], 0.9),
                      lwd=1.8)
}

# ------------- Add the annotation for APP beta sequence
if (plot_abeta == 1) {
  aa_prot_coords = c(672, 713)
  abeta_aa_seq = 'DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA'
  e16 = paste0('GTTCTGGGTTGACAAATATCAAGACGGAGGAGATCTCTGAAGTGAAGATGGATGCAGAATT',
               'CCGACATGACTCAGGATATGAAGTTCATCATCAAAAATTG')
  e17 = paste0('GTGTTCTTTGCAGAAGATGTGGGTTCAAACAAAGGTGCAATCATTGGACTCATGGTGGGCG',
               'GTGTTGTCATAGCGACAGTGATCGTCATCACCTTGGTGATGCTGAAGAAGAAACAGTACAC',
               'ATCCATTCATCATGGTGTGGTGGAG')
  abeta_coord = c(aa_prot_coords[1]*3 - 1 - exon_coords.dt[exon_no==16,
                                                           exon_start_abs],
                  aa_prot_coords[2]*3 + 1 - exon_coords.dt[exon_no==17,
                                                           exon_start_abs])
  abeta_dna_seq = paste0(
     str_sub(e16, aa_prot_coords[1]*3 - 1 - exon_coords.dt[exon_no==16,
                                                           exon_start_abs]),
     str_sub(e17, 1, aa_prot_coords[2]*3 + 1 - exon_coords.dt[exon_no==17,
                                                              exon_start_abs]))
  abeta_abs_coord = c(aa_prot_coords[1]*3 - 2, aa_prot_coords[2]*3)
  # draw.sector(start.degree = get.cell.meta.data("cell.start.degree",
  # sector.index='16') +
  #               get.cell.meta.data("cell.top.radius", track.index=2))
  circos.rect(xleft=abeta_abs_coord[1], xright=abeta_abs_coord[2],
              # ybottom=-0.55,
              ybottom = -12,
              ytop=1.55, sector.index='16',
              track.index=2,
              col=alpha('red', 0.30),
              #col=alpha('#ffffbf', 0.4),
              border='red',
              lwd=0.5, lty=1)
  
  # Add dashed line for abeta region inside SNV and indel tracks
  # for (i in 1:2) {
  #    circos.lines(track.index=3, x=c(abeta_abs_coord[i], abeta_abs_coord[i]),
  #                 y=c(-0.5, 4), lty=5, lwd=0.7, sector.index='16')
  #    circos.lines(track.index=4, x=c(abeta_abs_coord[i], abeta_abs_coord[i]),
  #                 y=c(-0.3, 2.7), lty=5, lwd=0.7, sector.index='16')
  #    circos.lines(track.index=5, x=c(abeta_abs_coord[i], abeta_abs_coord[i]),
  #                 y=c(-0.3, 4.3), lty=5, lwd=0.7, sector.index='16')
  # }


  circos.text(x=round(mean(abeta_abs_coord)), y=1,
              labels=expression(paste('A', beta)),
              col = 'red',
              sector.index='17', track.index=1)
  
}

# ------------ Legend for total counts
# Generate intervals from a vector
make_intervals = function (v, min_sym='', max_sym='\u2265') {
   n = length(v)
   ivs = transpose(list(v[1:(n-1)], v[2:n]))
   ivs_str = c()
   if (min_sym != '') {
      ivs_str = c(paste0('\u2264 ', ivs[[1]][1]))
   }
   ivs_str = c(ivs_str, sapply(ivs, function (iv) paste0('[', iv[1], ' - ', iv[2], ']')))
   if (max_sym != '') {
      ivs_str = c(ivs_str, paste0('\u2265 ', ivs[[length(ivs)]][2]))
   }
   return(ivs_str)
}

count_ivals_str = make_intervals(c(0, count_ivals[2:(length(count_ivals)-1)]))

legend(
   # x=1.08, y=-1.1,
   x=-1, y=1.1,
   inset=0, title="total count", yjust=1, xjust=0,
   legend=count_ivals_str,
   fill=count_colors, horiz=F, cex=0.7, y.intersp=0.8, bty='n',
   title.adj=0, adj=0)

# ------------ Legend for exon combinations
# legend(x=0.75, y=1.15, title='', legend=exon_comb_display[1:3],
#        fill=alpha(link_colors[1:4], 0.9), y.intersp=0.8,
#        cex=0.7, bty='n',
#        xjust=1, yjust=1, title.adj=0, adj=0)
# legend(x=1.11, y=1.15, title='', legend=exon_comb_display[4:length(exon_comb_display)],
#        fill=alpha(link_colors[4:length(link_colors)], 0.9), y.intersp=0.8,
#        cex=0.7, bty='n',
#        xjust=1, yjust=1, title.adj=0, adj=0)
legend(x=1.08, y=1.08, title='', legend=exon_comb_display,
       fill=alpha(link_colors, 0.9), y.intersp=0.8,
       cex=0.5, bty='n',
       xjust=1, yjust=1, title.adj=0, adj=0)


# Add labels for inner tracks
circ_font_size = 0.75
circos.text(x=exon_coords.dt[exon_id=='exon_09', round((exon_start_abs+exon_end_abs)/2)],
            y=log10(30), labels='SNV', facing='downward', cex=circ_font_size,
            track.index=3, sector.index=9)
circos.text(x=exon_coords.dt[exon_id=='exon_09', round((exon_start_abs+exon_end_abs)/2)],
            y=log10(10), labels='ins', facing='downward', cex=circ_font_size,
            track.index=4, sector.index=9)
circos.text(x=exon_coords.dt[exon_id=='exon_09', round((exon_start_abs+exon_end_abs)/2)],
            y=log10(10), labels='del', facing='downward', cex=circ_font_size,
            track.index=5, sector.index=9)

# End prototyping circo plot
cat(' OK.\n')
