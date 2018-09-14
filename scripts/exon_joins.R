#
# Detect and calculate insertions / deletions between exon junctions
#


rm(list=ls(all=T))
cat('Exon-exon join analysis\n')
options(warn = -1)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(optparse))
# source('aux_functions.R')

# Argument defaults
# input_blast_file_def = 'data/job-149_filtered_qs85_hpf_qs30_rep2_unique_blast_results_ws25_go0_gx2_exon_table.txt'
# input_blast_file_def = 'data/job-7871_filtered_qs85_hpf_qs30_rep2_unique_blast_results_ws25_go0_gx2_exon_table.txt'
# exon_orfs_file_def = 'data/app_exons_cds_clean_orfs.txt'
# exon_join_extract_length = 30 # Get 30 nts so 15 frxom each side of the intra-exon join

input_exon_table_def = 'examples/analysis/ex1_exon_table_hmpfix.txt'
# output_exon_table_summary_def = 'exampl'

# Parse user arguments 
option_list = list(
    make_option(c('-i', '--input-blast'),
                help = 'Input blast table file. This is a big blast.dt table 
                        from blast_results_plots.R.'),
    make_option(c('-o', '--output'),
                help = "Filename for the output table of exon-exon joins."),
    make_option(c('--use-homopolymer-fixed-seqs'),
                default = 1,
                help = 'Use exon_seq_corrected and read_seq_corrected columns
                        instead of exon_seq and read_seq in the exon table?'),
    make_option(c('--exon-comb-count-threshold'),
                default = 10,
                help = 'Minimum count threshold for reporting exon combinations.')
    # make_option(c('-o', '--output'),
    #             default = output_folder_def,
    #             help = 'Output folder for plots, e.g. plots/'),
    # make_option(c('-d', '--output_data'),
    #             default = output_data_def,
    #             help = 'Output folder for data, e.g. data/')
)
opt = parse_args(OptionParser(option_list=option_list))

# Load the alignment table
input_blast_file = opt$`input-blast`
output_exonjoins = opt$output
use_hmpfix = opt$`use-homopolymer-fixed-seqs`
exon_comb_th = opt$`exon-comb-count-threshold`
#output_folder = opt$output
#output_data = opt$output_data

cat('* Loading and processing data...')
# Load and generate neccessary tables
blast.dt = fread(input_blast_file, blank.lines.skip = T)
if (use_hmpfix == 1) {
   if ('read_seq_corrected' %in% colnames(blast.dt) &
       'exon_seq_corrected' %in% colnames(blast.dt)) {
      blast.dt[, read_seq := read_seq_corrected]
      blast.dt[, exon_seq := exon_seq_corrected]
      blast.dt[, c('read_seq_corrected', 'exon_seq_corrected') := NULL]
      cat(' OK.\n')
      cat('  - using homopolymer corrected sequences\n')
   } else {
      cat(' OK.\n')
      cat('  ? homopolymer corrected columns not found. Using uncorrected sequences.\n')
   }
}


exon_comb_grp.dt = blast.dt[, list(
   read_id = read_id,
   read_count = read_count,
   read_coverage = read_coverage,
   exon_comb_coverage = exon_comb_coverage
   ), by=exon_comb_short_rc] # changed exon_comb_grp -> exon_comb_short_rc

exon_comb_grp.dt = unique(exon_comb_grp.dt)
# Make a new folder for each exon combination


start.time <- Sys.time()
cat('* Finding and organizing all exon-exon joins... ')
# For each read generate pairs of neighboring exons 
# and note which ones are exon-exon joins (non-consecutive exons)
exon_joins.dt = blast.dt[, {
   c(.SD[1:(.N-1), list(exon_left=exon_no, 
                        exon_left_len=exon_length, 
                        exon_left_start_pos_align=exon_start_pos_align,
                        exon_left_end_pos_align=exon_end_pos_align,
                        exon_left_read_end_pos_align=read_end_pos_align)],
        .SD[2:.N, list(exon_right=exon_no, 
                       exon_right_len=exon_length, 
                       exon_right_start_pos_align=exon_start_pos_align,
                       exon_right_end_pos_align=exon_end_pos_align,
                       read_count=read_count, 
                       exon_comb_short_rc=exon_comb_short_rc,
                       exon_right_read_start_pos_align=read_start_pos_align)])
}, by=read_id]
exon_joins.dt[, exon_comb_type := ifelse(exon_right == exon_left + 1, 1, 0)]

# For each exon-exon join figure out how much is missing (if there is anything missing)
# from left and right exons.

# For left and right exons, find if there are missing nts in the junction
exon_joins.dt[, exon_left_miss := exon_left_len - exon_left_end_pos_align]
exon_joins.dt[, exon_right_miss := exon_right_start_pos_align - 1]
# Exon-exon join strings
exon_joins.dt[, exon_join := paste(exon_left, exon_right, sep='-') ]
#exon_joins.dt[, exon_comb_short := as.character(lapply(exon_comb_order, function (x) comb_text_short(as.integer(strsplit(x, ',')[[1]]))))]
#exon_joins.dt[, exon_comb_short := factor(exon_comb_short, levels=rev(unique(exon_comb_short)))]
exon_joins.dt[, exon_order := 1:.N, by=read_id]
# Calculate the gap between exons. Negative gap means they overlap.
exon_joins.dt[, gap := exon_right_read_start_pos_align - exon_left_read_end_pos_align - 1]
exon_joins.dt[, exon_exon_gap_corrected := ifelse(gap<=0, 0, gap)]

cat('OK.\n')

# cat('* calculating exon-exon join statistics...')
# #
# # For each exon-exon join (in each exon combination isoform) show the histogram
# # of exon_left_end_pos_align and exon_right_pos_align
# #
# exon_join_hist.dt = melt(exon_joins.dt, id.vars=c('exon_join', 'exon_comb_short_rc', 'read_count'), 
#                          measure.vars = c('exon_left_miss', 'exon_right_miss'), variable.name = 'exon_position_side',
#                          value.name = 'exon_missing_nts')
# exon_join_hist_tab.dt = exon_join_hist.dt[, 
#    list(exon_missing_nt_count = sum(read_count)), 
#    by=c('exon_comb_short_rc', 'exon_join', 'exon_position_side', 'exon_missing_nts')]
# 
# 
# # Similarly, for each join show the gaps
# exon_join_gap_tab.dt = exon_joins.dt[, list(exon_exon_gap_count = sum(read_count)), 
#    by=c('exon_comb_short_rc', 'exon_join', 'exon_exon_gap_corrected')]
# 
# cat(' OK.\n')

cat('* saving exon-joins table...')
# Extract only non-consecutive joins (ignore 1,18 and 18,1 for a moment)
write.table(
#   exon_joins.dt[, .(exon_join, read_id, exon_left_read_end_pos_align, 
#                     exon_right_read_start_pos_align,exon_comb_short_rc, read_count)],
   exon_joins.dt,
   output_exonjoins,
   sep='\t', row.names=F, quote=F
)
cat(' OK.\n')
cat('  ->', output_exonjoins, '\n')

