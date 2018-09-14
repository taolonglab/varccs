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

option_list = list(
  make_option(c('-i', '--input-snvs-indels'),
              help = 'Output from runs_snvs_indels_analysis.R'),
  
)

# Input
blast_fixed = opt$`input-blast`
snv_indel_fix = opt$`input-snv-indel`

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
count_ivals = c(1,3,10,30,100,300)
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


cat('Finding all intra-exon joins...')
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
