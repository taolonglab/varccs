# Load BLAST results from aligning exons to PacBio reads
#
# Input. 1. BLAST results table
#        2. read length table
#        3. exon length table
# CHeck why is read coverage weird (>1 ?) after the latest run.

cat('* BLAST results analysis\n')
# cat('------------------\n')

# Load required libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(optparse))

# setwd('/Users/igor/cloud/research/app')
options(warn = -1)

# Argument defaults
input_filename_def = 'data/ming1-nonAD-nested-20pass-_blast_results_ws25_go0_gx2.txt'
readlen_filename_def = 'data/ming1-nonAD-nested-20pass-_length_table.txt'

output_folder_def = 'plots/'
exonlen_filename_def = 'data/app_exons_cds_clean_rc_length_table.txt'

# Select only exon combinations with the read count greater than this
exon_comb_dir_th = 0

# Parse user arguments
option_list = list(
    make_option(c('-i', '--input'),
                help = 'Input file from blast results.'),
    make_option(c('-o', '--output'),
                help = 'Output exon alignment table file.'),
    make_option(c('-s', '--output_sum'),
                help = 'Output exon summary file.'),
    make_option(c('-r', '--read_lengths'),
                help = 'File with read lengths: read_id\tread_length'),
    make_option(c('-e', '--exon_lengths'),
                help = 'File with exon lengths: exon_id\texon_length')
)
opt = parse_args(OptionParser(option_list=option_list))

# Load the alignment table
input_filename = opt$input
readlen_filename = opt$read_lengths
exonlen_filename = opt$exon_lengths

output_exon_table = opt$output
output_exon_summary = opt$output_sum

# input_base = gsub('\\..*', '', basename(input_filename))
# 
# # Append slash if not part of folder name
# if (str_sub(output_folder, -1) != '/') {
#    output_folder = paste0(output_folder, '/')
# }

cat('  - loading and processing data...')
blast.dt = fread(input_filename, blank.lines.skip=T)
blast.dt[, exon_no := lapply(exon_id, function (t) as.integer(strsplit(t, '_')[[1]][2]))]
blast.dt[, exon_no := as.integer(exon_no)]
blast.dt[, exon_orientation := 0]
blast.dt[exon_id %like% 'rc', exon_orientation := 1]
blast.dt[, read_direction := '+']
blast.dt[exon_id %like% '_rc', read_direction := '-']
# Extract read_id and read_count from the read_id column
blast.dt[, c('read_id', 'read_count') := tstrsplit(read_id, ';')]
blast.dt[, read_count := as.integer(gsub('count_', '', read_count))]

# Filter out alignments that are literally just primers
# blast.dt = blast.dt[nchar(read_seq) == nchar(exon_seq) & nchar(exon)]

# Load read and exon lengths, then join tables
read_len.dt = fread(readlen_filename)
blast.dt = merge(blast.dt, read_len.dt, by='read_id')
exon_len.dt = fread(exonlen_filename)
blast.dt = merge(blast.dt, exon_len.dt, by='exon_id')
# blast.dt[, c('evalue', 'bitscore') := NULL]

# Calculate aligned exon and read coverage
# For each exon, calculate it's % aligned coverage in the read (how much of the exon is in each read)
blast.dt[, exon_coverage := round((exon_end_pos_align-exon_start_pos_align+1)/exon_length, 2)]
# For each read calculate how much of it is covered by exons
blast.dt[, read_coverage := round(sum(read_end_pos_align-read_start_pos_align+1-no_mismatches-no_gap_openings)/read_length, 2), by=read_id]

# First find the combination in the order in which it occurs
blast.dt = blast.dt[order(read_id, read_start_pos_align)]
blast.dt[, exon_comb_order := paste(exon_no, collapse=','), by=read_id]
blast.dt[, exon_comb_coverage := paste(exon_coverage, collapse=','), by=read_id]

# Filter sequences
# First, remove the sequences that do not start with either exon 1 or exon 18
blast.dt = blast.dt[, if ( (exon_no[1] == 1 & exon_no[.N] == 18) | (exon_no[1] == 18 & exon_no[.N] == 1) ) .SD, by=read_id]

cat('OK.\n')
cat('  - generating various plots...')


# Show histogram of exon numbers per read
exon_no.dt = blast.dt[, .N, by=read_id]
exon_no_sum.dt = exon_no.dt[, list(count=.N), by=N][order(N)]

# Display exon combinations by abundance
read_cov.dt = blast.dt[, list(read_coverage = unique(read_coverage)), by=read_id]

cat('OK.\n')
cat('  - calculating exon statistics...')

# Now group the data table according to unique exon combinations. We will treat the same
# combinations from either 5' or 3' direction. I.e. 1,2,3,4,15 same as 15,4,3,2,1
# Use ' to denote exons from the opposite strand
blast.dt[, exon_no_str := paste(exon_no, ifelse(exon_orientation == 0, '', '\''), sep='')]
blast.dt[, exon_comb_grp := paste(sort(c(paste(exon_no, collapse=','), paste(rev(exon_no), collapse=','))), collapse=';'), by=read_id]
blast.dt[, exon_comb_grp_rc := paste(sort(c(paste(exon_no_str, collapse=','), paste(rev(exon_no_str), collapse=','))), collapse=';'), by=read_id]

# This function reduces some long exon combinations like
# 1,2,3,4,5,6,16,17,18 to 1-6,16-18
#
# Input is a stringsplitted vector of integers

comb_text_short = function (v) {
  x = v[1]
  s = paste0('', v[1])
  d = '' # direction of the current interval
  for (i in 2:length(v)) {
    if ((x+1) == v[i] & d != '-') {
       d = '+'
       if (str_sub(s, -1) != '-') {
         s = paste0(s, '-')
       }
    }
    else if ((x-1) == v[i] & d != '+') {
       d = '-'
       if (str_sub(s, -1) != '-') {
         s = paste0(s, '-')
       }
    }
    else {
       # Update last interval if its open
       d = ''
       if (str_sub(s, -1) == '-') {
         s = paste0(s, v[i-1])
       }
       # Add current number
       s = paste0(s, ',', v[i])
    }
    x = v[i]
  }
  # Check if we have any unclosed intervals
  if (str_sub(s, -1) == '-') {
     s = paste0(s, v[length(v)])
  }
  return(s)
}

comb_text_short_wrc = function (v) {
  #
  # If all sequences are either with exons that are in
  # the same direction, just output them normally.
  if (all(grepl("'", v))) {
     return(comb_text_short(as.integer(gsub('\'', '', v))))
  }
  if (all(!grepl("'", v))) {
     return(comb_text_short(as.integer(v)))
  }
  x = as.integer(gsub('\'', '', v[1]))
  s = paste0('', v[1])
  d = '' # direction of the current interval
  for (i in 2:length(v)) {
    vi = as.integer(gsub('\'', '', v[i]))
    if ((x+1) == vi & d != '-') {
       d = '+'
       if (str_sub(s, -1) != '-') {
         s = paste0(s, '-')
       }
    }
    else if ((x-1) == vi & d != '+') {
       d = '-'
       if (str_sub(s, -1) != '-') {
         s = paste0(s, '-')
       }
    }
    else {
       # Update last interval if its open
       d = ''
       if (str_sub(s, -1) == '-') {
         s = paste0(s, v[i-1])
       }
       # Add current number
       s = paste0(s, ',', v[i])
    }
    x = vi
  }
  # Check if we have any unclosed intervals
  if (str_sub(s, -1) == '-') {
     s = paste0(s, v[length(v)])
  }
  return(s)
}

exon_comb_revcomp = function (x) {
   # Input: string such as
   # Output: reverse complement it, if it starts with a number
   #         followed by prime symbol ', i.e.
   #         1-3,17-18,18'-17',1',1-3,17-18
   o = gsub('([0-9]+)-', '\\1#-', x)
   o = gsub('([0-9]+),', '\\1#,', o)
   o = gsub('([0-9]+)\'', '\\1', o)
   o = gsub('([0-9]+)#', '\\1\'', o)
   return(o)
}

#
# Reorder exon_comb_short to make sure we group together reads from
# sense and anti-sense strands
# blast.dt[, exon_comb_short := as.character(lapply(exon_comb_order, function (x) comb_text_short(as.integer(strsplit(x, ',')[[1]]))))]
blast.dt[, exon_comb_short := as.character(lapply(exon_comb_grp, function (x) comb_text_short(as.integer(strsplit(strsplit(x, ';')[[1]][1], ',')[[1]]))))]
blast.dt[, exon_comb_short_rc := as.character(lapply(exon_comb_grp_rc, function (x) comb_text_short_wrc((strsplit(strsplit(x, ';')[[1]][1], ',')[[1]]))))]
blast.dt[exon_comb_short_rc %like% '^[0-9]+\'', exon_comb_short_rc := exon_comb_revcomp(exon_comb_short_rc)]

exon_comb_grp.dt = blast.dt[, list(
   read_id = read_id,
   read_count = read_count,
   read_coverage = read_coverage,
   exon_comb_coverage = exon_comb_coverage,
   exon_comb_order = exon_comb_order,
   exon_comb_short_rc = exon_comb_short_rc
   ), by=exon_comb_grp]

exon_comb_grp.dt = unique(exon_comb_grp.dt)
exon_comb_grp.dt[, exon_comb_dir := tstrsplit(exon_comb_grp, ';')[[1]]]
exon_comb_grp.dt[, exon_comb_grp := NULL]

# Count how many times each exon combination occurs
# Change dir to undir if we want to count different strands as the same
exon_comb_grp.dt[, read_coverage_good := ifelse(read_coverage > 0.95, read_count, 0.0)]

exon_comb_sum.dt = exon_comb_grp.dt[, list(
   exon_comb_dir_count = sum(read_count),
   read_coverage_good = sum(read_coverage_good),
   exon_comb_short_rc = unique(exon_comb_short_rc)
   ), by=exon_comb_dir][order(-exon_comb_dir_count)]
# Generate short(er) names, which will help display some long exon combinations
exon_comb_sum.dt[, exon_comb_short := as.character(lapply(exon_comb_dir, function (x) comb_text_short(as.integer(strsplit(x, ',')[[1]]))))]
# exon_comb_sum.dt[, exon_comb_short := as.character(lapply(exon_comb_order, function (x) comb_text_short(as.integer(strsplit(x, ',')[[1]]))))]
exon_comb_sum.dt[, exon_comb_short := factor(exon_comb_short, levels=rev(unique(exon_comb_short)))]

# Generate a smaller data where we take the cut on number of exon combinations
exon_comb_sum.subdt = exon_comb_sum.dt[exon_comb_dir_count > exon_comb_dir_th, .(exon_comb_short, exon_comb_dir_count, read_coverage_good)]
exon_comb_sum.subdt[, exon_comb_short := factor(exon_comb_short, levels=rev(unique(exon_comb_short)))]
exon_comb_sum.subdt[, read_coverage_bad := exon_comb_dir_count - read_coverage_good]
exon_comb_sum_melt.subdt = melt(
   exon_comb_sum.subdt, id.vars=c('exon_comb_short'),
   measure.vars=c('read_coverage_good', 'read_coverage_bad'),
   variable.name='read_coverage_type',
   value.name='read_counts'
)
exon_comb_sum_melt.subdt = exon_comb_sum_melt.subdt[order(exon_comb_short, -read_counts)]

# cat('OK.\n')

# Plot read coverage for unique reads
exon_comb_grp.dt[, exon_comb_short := as.character(lapply(exon_comb_dir, function (x) comb_text_short(as.integer(strsplit(x, ',')[[1]]))))]
exon_comb_grp.dt[, exon_comb_short := factor(exon_comb_short, levels=rev(unique(exon_comb_sum.dt$exon_comb_short)))]

# For the most abundant exon combinations (let's say first 6), find the junctions between
# consecutive exon groups, i.e. 1-3,16-18 get both exons on each side of the comma symbol.
# We can extract the borderline

#exon_combinations = exon_comb_sum.dt[exon_comb_dir_count > 10 & !(exon_comb_short %like% "18,"), exon_comb_short]
# exon_combinations = exon_comb_sum.dt[exon_comb_dir_count > 0, exon_comb_short]
#exon_combinations = exon_combinations[1:min(10, length(exon_combinations))]
#exon_comb_ids = sub('.*(\\d+),(\\d+).*', '\\1,\\2', exon_combinations)

cat('OK.\n')

# Write exon combinations to file too
cat('  - writing exon combinations...')
write.table(exon_comb_sum.dt, output_exon_summary, sep='\t', quote=F, row.names=F)
cat(' OK.\n')
cat('    ->', output_exon_summary, '\n')

cat('  - writing final table...')
write.table(blast.dt, output_exon_table, sep='\t', quote=F, row.names=F)
cat(' OK.\n')
# cat('    ->', output_exon_table, '\n')
# cat('  OK.\n')
