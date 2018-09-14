#
# Analyze and make plots for SNVs and indels based on the exon table
# with fixed bad homopolymer runs.
#
# options(datatable.print.nrows = 10)
# Load exon table with fixed homopolymer sequences

rm(list=ls(all=TRUE))

# Load required libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(RColorBrewer))

setwd('~/cloud/research/app')
options(warn = -1)

# Argument defaults
sample_id = 'job43-2AD'
input_filename_def = paste0('data/', sample_id, '/', sample_id, '_filtered_qs85_hpf_qs30_rep2_unique_primerblasted_ws13_unique_blast_results_ws25_go0_gx2_exon_table_fixed.txt')
input_exonsum_def = paste0('data/', sample_id, '/', sample_id, '_filtered_qs85_hpf_qs30_rep2_unique_primerblasted_ws13_unique_blast_results_ws25_go0_gx2_exon_combinations_summary.txt')
input_fasta_def = paste0('reads/', sample_id, '/', sample_id, '_filtered_qs85_hpf_qs30_rep2_unique_primerblasted_ws13_unique.fa')
output_folder_def = paste0('plots/', sample_id, '/')
output_data_def = paste0('data/', sample_id, '/')


# input_filename_def = paste0('data/ad-31-150/ad-31-150_filtered_qs85_hpf_qs30_rep2_unique_blast_results_ws25_go0_gx2_exon_table_fixed.txt')
# input_exonsum_def = paste0('data/ad-31-150/ad-31-150_filtered_qs85_hpf_qs30_rep2_unique_blast_results_ws25_go0_gx2_exon_combinations_summary.txt')
# input_fasta_def = 'reads/ad-31-150/ad-31-150_filtered_qs85_hpf_qs30_rep2_unique.fa'


input_exonlen_def = 'data/app_exons_cds_clean_rc_length_table.txt'


# For ad-31-150 (AD)
exon_comb_display = c(
   '1-3,17-18', '1-3,16-18', '1-2,17-18', '1-5,16-18', '1-2,18',
   '1-7,9-18', '1-3,14-18', '1,12-18', '1-2,16-18', '1,17-18',
   '1,14-18', '1-2,18')

# Combined AD and ND to show the same colors
exon_comb_display = c(
   '1,17-18', '1-2,18', '1-7,9-18', '1-3,17-18', '1-2,16-18', '1-6,16-18',
   '1-6,9-18', '1,14-18', '1-3,18', '1-3,16-18',  '1-5,16-18',
   '1,12-18', '1-3,14-18'
)

# Link colors for CIRCO plot
link_colors = c('black',   '#e41a1c', '#377eb8', '#4daf4a',
                '#984ea3', '#a65628', 'pink', 'orange', 'gray',
                'brown')
link_colors = c(brewer.pal(12, 'Paired'), 'black')
link_colors[11] = 'gray'




# Annotation for insertions, deletions and SNVs in this order
# '-' means both read and exon have '-' so it should be removed.
# (in principle should not happen at this stage.)
annotation = c('insertion', 'deletion', 'SNV', '-')
nts = c('A', 'C', 'G', 'T', '-')

# Parse user arguments
option_list = list(
    make_option(c('-i', '--input'),
       default = input_filename_def,
       help = 'Input file from fasta_homopolymer_fix.R (*_exon_table_hmpfix.txt'),
    make_option(c('-e', '--input_exon'),
       default = input_exonlen_def,
       help = 'Input file with exon lengths. We need this for plotting.'),
    make_option(c('-s', '--input_sum'),
       default = input_exonsum_def,
       help = 'Input file with the counts for each exon combination.'),
    make_option(c('-f', '--input_fa'),
       default = input_fasta_def,
       help = 'Input corrected FASTA file. Will need this for manual exon comb correction.'
    ),
    make_option(c('-o', '--output'),
       default = output_folder_def,
       help = 'Output folder for SNV and indel plots, e.g. plots/.'),
    make_option(c('-d', '--output_data'),
       default = output_data_def,
       help = 'Output folder for computed tables, i.e. data/.')
)
opt = parse_args(OptionParser(option_list=option_list))

# Load the alignment table
input_filename = opt$input
input_exonlen = opt$input_exon
input_exonsum = opt$input_sum
input_codontab = opt$input_codontable
input_fasta = opt$input_fa
output_folder = opt$output
output_data_folder = opt$output_data

input_base = gsub('\\..*', '', basename(input_filename))
input_base = gsub('_blast_results.*', '', input_base)

# Append slash if not part of folder name
if (str_sub(output_folder, -1) != '/') {
   output_folder = paste0(output_folder, '/')
}

cat('Loading FASTA sequence file...')
fa = readDNAStringSet(input_fasta)
fa.dt = data.table(
   read_id=sapply(strsplit(names(fa), ';', fixed=T), `[`, 1),
   seq=as.character(fa)
)
cat(' OK.\n')

cat('Loading exon lengths...')
exonlen.dt = fread(input_exonlen)
exonlen.dt[, exon_no := as.integer(gsub('_rc', '',
                                        gsub('exon_(.*).*', '\\1', exon_id)))]
cat(' OK.\n')

cat('Loading exon table fixed file...')
blast_fixed.dt = fread(input_filename)
blast_fixed.dt = merge(blast_fixed.dt, exonlen.dt[, .(exon_id, exon_length)],
                       by='exon_id')
cat(' OK.\n')

# cat('Loading exon summary file...')
# exonsum.dt = fread(input_exonsum)
# exon_comb_out = c('18,1', '17,1', '3,1-', '1,1-', '3,1,1')
# exon_comb_display = exonsum.dt[!(exon_comb_short_rc %like%
# paste(exon_comb_out, collapse='|')), exon_comb_short_rc]
# 
# cat(' OK.\n')

#
# ---------------------- Manually fix some missing alignments --------------------
# First remove joins: 1-2,16-18,  1-5,17-18,   1-2,17-18
# Then update these:
# 1,17-18 --> 1-2,17-18: add perfectly joined 2 with length 7 everywhere
# 1-3,15-18 --> 1-3,14-18: add perfectly joined 14 with 15 with length 20
#
cat('Manually fixing missing BLAST alignments...')
# 1,17-18 to 1-2,17-18
exon_comb_noisy = c('1-2,16-18', '1-5,17-18', '1-2,17-18', "18'-16',2'-1'",
                    "18'-17',5'-1'", "18'-17',2'-1'", ',18,18,', '18,18$',
                    '1,2,1')
exon_comb_noisy = c(
   #'1-5,17-18',
                    '1-2,17-18', "18'-16',2'-1'",
                    "18'-17',5'-1'", "18'-17',2'-1'", ',18,18,', '18,18$',
                    '1,2,1')
for (i in 1:length(exon_comb_noisy)) {
   #if (!(exon_comb_noisy[i] %in% exon_comb_display)) {
      blast_fixed.dt = blast_fixed.dt[!(exon_comb_short_rc %like% exon_comb_noisy[i])]
   #}
}

adjacent_exon_ids = function (e1, e2, l) {
   # find the indices of l where e1 occurs just before e2
   # e1 = integer, exon 1 to search
   # e2 = integer, exon 2 to search
   # l = vector to search in (usually column from data table)
   for (i in 2:length(l)) {
      if (l[i-1] == e1 & l[i] == e2) return(c(i-1, i))
   }
   return(c(NA, NA))
}

# Sort blast file by read_start_pos_align
blast_fixed.dt = blast_fixed.dt[order(read_id, read_start_pos_align)]
blast_fixed.dt[exon_comb_short_rc=='1-2,17-18', exon_comb_short_rc := {
   if (.SD[exon_id %like% 'exon_02', alignment_length] > 137) {
      '1-2,17-18x'
   }
}, by=read_id]

# blast_fixed.dt[exon_comb_short_rc=='1-3,17-18', exon_comb_short_rc := {
#    if (.SD[exon_id %like% 'exon_02', alignment_length] > 137) {
#       '1-2,17-18x'
#    }
# }, by=read_id]


# This shows up in job-149 non-AD data that needs fixing.
if (nrow(blast_fixed.dt[exon_comb_short_rc %like% '1-2,17-18x' & exon_id %like% 'exon_02' &
                        alignment_length > 137]) > 0) {
   blast_fixed.dt = merge(
      blast_fixed.dt,
      unique(blast_fixed.dt[exon_comb_short_rc %like% '1-2,17-18x', {
        # Find coordinates for the adjacent exon 1 and exon 17 if we have multiple
        # exons 1 and 17 in the subdata datable.
        if (.SD[exon_id %like% 'exon_02', alignment_length] > 137) {
           read = .BY[[1]]
           dt = .SD[exon_id %like% 'exon_02']
           dt2 = .SD[exon_id %like% 'exon_17']
           if (.SD[exon_id=='exon_17', read_start_pos_align] >
               .SD[exon_id=='exon_02', read_end_pos_align+1]) {
              rc = F
              ids = adjacent_exon_ids('exon_02', 'exon_17', .SD[, exon_id])
              if (all(is.na(ids)==T)) {
                 ids = adjacent_exon_ids('exon_17_rc', 'exon_02_rc', .SD[, exon_id])
                 rc = T
              }
              if (!rc) {
                 list(
                    exon_id='exon_03',
                    read_start_pos_align=unique(.SD[ids[1], read_end_pos_align+1]),
                    read_end_pos_align=unique(.SD[ids[2], read_start_pos_align-1]),
                    exon_start_pos_align=1L,
                    exon_end_pos_align=23L,
                    alignment_length=23,
                    read_seq=substr(fa.dt[read_id==read, seq],
                                    unique(.SD[ids[1], read_end_pos_align+1]),
                                    unique(.SD[ids[2], read_start_pos_align-1])),
                    exon_seq='GTCTACCCTGAACTGCAGATCAC',
                    exon_comb_short_rc=gsub('1-2,17-18x', '1-3,17-18', unique(exon_comb_short_rc)),
                    read_count=unique(read_count),
                    exon_seq_corrected='GTCTACCCTGAACTGCAGATCAC',
                    read_seq_corrected=substr(fa.dt[read_id==read, seq],
                                              unique(.SD[ids[1], read_end_pos_align+1]),
                                              unique(.SD[ids[2], read_start_pos_align-1])),
                    exon_length=exonlen.dt[exon_id=='exon_03', exon_length]
                 )
              }
              else {
                 list(
                    # This might not be bug free. If we decide to skip primer blasting
                    # this might need to be fixed.
                    exon_id='exon_03_rc',
                    read_start_pos_align=unique(.SD[ids[1], read_end_pos_align+1]),
                    read_end_pos_align=unique(.SD[ids[2], read_start_pos_align-1]),
                    exon_start_pos_align = exonlen.dt[exon_id=='exon_03', exon_length-22L],
                    exon_end_pos_align = exonlen.dt[exon_id=='exon_03', exon_length],
                    alignment_length=23,
                    read_seq=substr(fa.dt[read_id==read, seq],
                                    unique(.SD[ids[1], read_end_pos_align+1]),
                                    unique(.SD[ids[2], read_start_pos_align-1])),
                    exon_seq='GTGATCTGCAGTTCAGGGTAGAC',
                    exon_comb_short_rc=gsub('1-2,17-18x', '1-3,17-18', unique(exon_comb_short_rc)),
                    read_count=unique(read_count),
                    exon_seq_corrected='GTGATCTGCAGTTCAGGGTAGAC',
                    read_seq_corrected=substr(fa.dt[read_id==read, seq],
                                              unique(.SD[ids[1], read_end_pos_align+1]),
                                              unique(.SD[ids[2], read_start_pos_align-1])),
                    exon_length=exonlen.dt[exon_id=='exon_03', exon_length]
                 )
              }
           }

        }
     }, by=read_id
     ]),
     all=T,
     by=names(blast_fixed.dt)
   )

   blast_fixed.dt[exon_comb_short_rc %like% '1-2,17-18x',
                  exon_comb_short_rc := '1-3,17-18']
}

# Sort blast file by read_start_pos_align
blast_fixed.dt = blast_fixed.dt[order(read_id, read_start_pos_align)]


if (nrow(blast_fixed.dt[exon_comb_short_rc %like% "18'-17',1'"]) > 0) {
   blast_fixed.dt = merge(
     blast_fixed.dt,
     unique(blast_fixed.dt[exon_comb_short_rc %like% "18'-17',1'", {
      # Find coordinates for the adjacent exon 1 and exon 17 if we have multiple
      # exons 1 and 17 in the subdata datable.
      ids = adjacent_exon_ids('exon_17_rc', 'exon_01_rc', .SD[, exon_id])
      list(
      exon_id='exon_02_rc',
      read_start_pos_align=unique(.SD[ids[1], read_end_pos_align-2]),
      read_end_pos_align=unique(.SD[ids[2], read_start_pos_align-1]),
      exon_start_pos_align=1,
      exon_end_pos_align=7,
      alignment_length=7,
      read_seq='TGGGTAC', exon_seq='TGGGTAC',
      exon_comb_short_rc=gsub("18'-17',1'", "18'-17',2'-1'",
                              unique(exon_comb_short_rc)),
      read_count=unique(read_count),
      exon_seq_corrected='TGGGTAC', read_seq_corrected='TGGGTAC',
      exon_length=exonlen.dt[exon_id=='exon_02_rc', exon_length]
      )}, by=read_id]),
     all=T, by=names(blast_fixed.dt))

   # blast_fixed.dt[exon_comb_short_rc %like% , exon_comb_short_rc := '1-2,17-18']
   blast_fixed.dt[exon_comb_short_rc %like% "18'-17',1'",
                  exon_comb_short_rc := gsub("18'-17',1'", "18'-17',2'-1'",
                                             exon_comb_short_rc)]
}

# Sort blast file by read_start_pos_align
blast_fixed.dt = blast_fixed.dt[order(read_id, read_start_pos_align)]


if (nrow(blast_fixed.dt[exon_comb_short_rc %like% '1-3,15-18']) > 0) {
   blast_fixed.dt = merge(
     blast_fixed.dt,
     unique(blast_fixed.dt[exon_comb_short_rc %like% '1-3,15-18', {
      rc = F
      ids = adjacent_exon_ids('exon_03', 'exon_15', .SD[, exon_id])
      if (all(is.na(ids)==T)) {
         ids = adjacent_exon_ids('exon_15_rc', 'exon_03_rc', .SD[, exon_id])
         rc = T
      }
      if (!rc) {
         list(
            exon_id='exon_14',
            read_start_pos_align=unique(.SD[ids[1], read_end_pos_align-6]),
            read_end_pos_align=unique(.SD[ids[2], read_start_pos_align-1]),
            exon_start_pos_align=exonlen.dt[exon_id=='exon_14', exon_length-19L],
            exon_end_pos_align=exonlen.dt[exon_id=='exon_14', exon_length],
            alignment_length=20,
            read_seq='AGCCAACACAGAAAACGAAG', exon_seq='AGCCAACACAGAAAACGAAG',
            exon_comb_short_rc=gsub('1-3,15-18', '1-3,14-18',
                                                   exon_comb_short_rc),
            read_count=unique(read_count),
            exon_seq_corrected='AGCCAACACAGAAAACGAAG',
            read_seq_corrected='AGCCAACACAGAAAACGAAG',
            exon_length=exonlen.dt[exon_id=='exon_14', exon_length]
         )
      }
      else {
         list(
            exon_id='exon_14_rc',
            read_start_pos_align=unique(.SD[ids[1], read_end_pos_align+1]),
            read_end_pos_align=unique(.SD[ids[2], read_start_pos_align+6]),
            # exon_start_pos_align=exonlen.dt[exon_id=='exon_14', exon_length-19L],
            # exon_end_pos_align=exonlen.dt[exon_id=='exon_14', exon_length],
            exon_start_pos_align = 1L,
            exon_end_pos_align = 20L,
            alignment_length=20,
            read_seq='CTTCGTTTTCTGTGTTGGCT', exon_seq='CTTCGTTTTCTGTGTTGGCT',
            exon_comb_short_rc=gsub('1-3,15-18', '1-3,14-18',
                                                   exon_comb_short_rc),
            read_count=unique(read_count),
            exon_seq_corrected='CTTCGTTTTCTGTGTTGGCT',
            read_seq_corrected='CTTCGTTTTCTGTGTTGGCT',
            exon_length=exonlen.dt[exon_id=='exon_14', exon_length]
         )
      }
   }, by=read_id
   ]), all=T, by=names(blast_fixed.dt))
   blast_fixed.dt[exon_comb_short_rc %like% '1-3,15-18',
                  exon_comb_short_rc := gsub('1-3,15-18', '1-3,14-18',
                                             exon_comb_short_rc)]
}

# Sort blast file by read_start_pos_align
blast_fixed.dt = blast_fixed.dt[order(read_id, read_start_pos_align)]

if (nrow(blast_fixed.dt[exon_comb_short_rc %like% "18'-15',3'-1'"]) > 0) {
   blast_fixed.dt = merge(
     blast_fixed.dt,
     unique(blast_fixed.dt[exon_comb_short_rc %like% "18'-15',3'-1'", {
      ids = adjacent_exon_ids('exon_15_rc', 'exon_03_rc', .SD[, exon_id])
      list(
      exon_id='exon_14_rc',
      read_start_pos_align=unique(.SD[ids[1], read_end_pos_align+1]),
      read_end_pos_align=unique(.SD[ids[2], read_start_pos_align+6]),
      # exon_start_pos_align=exonlen.dt[exon_id=='exon_14', exon_length-19L],
      # exon_end_pos_align=exonlen.dt[exon_id=='exon_14', exon_length],
      exon_start_pos_align = 1L,
      exon_end_pos_align = 20L,
      alignment_length = 20L,
      read_seq='CTTCGTTTTCTGTGTTGGCT', exon_seq='CTTCGTTTTCTGTGTTGGCT',
      exon_comb_short_rc=gsub("18'-15',3'-1'", "18'-14',3'-1'",
                                             exon_comb_short_rc),
      read_count=unique(read_count),
      exon_seq_corrected='CTTCGTTTTCTGTGTTGGCT',
      read_seq_corrected='CTTCGTTTTCTGTGTTGGCT',
      exon_length=exonlen.dt[exon_id=='exon_14', exon_length]
      )}, by=read_id
   ]), all=T, by=names(blast_fixed.dt))
   blast_fixed.dt[exon_comb_short_rc %like% "18'-15',3'-1'",
                  exon_comb_short_rc := gsub("18'-15',3'-1'", "18'-14',3'-1'",
                                             exon_comb_short_rc)]
}

# Sort blast file by read_start_pos_align
blast_fixed.dt = blast_fixed.dt[order(read_id, read_start_pos_align)]

if (nrow(blast_fixed.dt[exon_comb_short_rc %like% '1,18']) > 0) {
   blast_fixed.dt = merge(
     blast_fixed.dt,
     unique(blast_fixed.dt[exon_comb_short_rc %like% '1,18', {
      if (.SD[1, exon_id] == 'exon_01') {
        rc = F
        ids = c(1, 2)
      } else {
        rc = T
        ids = c(2, 1)
      }
      # ids = adjacent_exon_ids('exon_01', 'exon_18', .SD[, exon_id])
      # if (all(is.na(ids)==T)) {
      #    ids = adjacent_exon_ids('exon_18_rc', 'exon_01_rc', .SD[, exon_id])
      #    rc = T
      # }
      if (!rc) {
         list(
            exon_id='exon_17',
            read_start_pos_align=unique(.SD[ids[1], read_end_pos_align+1]),
            read_end_pos_align=unique(.SD[ids[2], read_start_pos_align-1]),
            exon_start_pos_align=exonlen.dt[exon_id=='exon_17', exon_length-13L],
            exon_end_pos_align=exonlen.dt[exon_id=='exon_17', exon_length],
            alignment_length=14,
            read_seq='ATGGTGTGGTGGAG', exon_seq='ATGGTGTGGTGGAG',
            exon_comb_short_rc=gsub('1,18', '1,17-18',
                                                   exon_comb_short_rc),
            read_count=unique(read_count),
            exon_seq_corrected='ATGGTGTGGTGGAG',
            read_seq_corrected='ATGGTGTGGTGGAG',
            exon_length=exonlen.dt[exon_id=='exon_17', exon_length]
         )
      }
      else {
         list(
            exon_id='exon_17_rc',
            read_start_pos_align=unique(.SD[ids[1], read_end_pos_align-1]),
            read_end_pos_align=unique(.SD[ids[2], read_start_pos_align+1]),
            # exon_start_pos_align=exonlen.dt[exon_id=='exon_17', exon_length-13],
            # exon_end_pos_align=exonlen.dt[exon_id=='exon_17', exon_length],
            exon_start_pos_align = 1L,
            exon_end_pos_align = 14L,
            alignment_length=14,
            read_seq='CTCCACCACACCAT', exon_seq='CTCCACCACACCAT',
            exon_comb_short_rc=gsub('1,18', '1,17-18',
                                                   exon_comb_short_rc),
            read_count=unique(read_count),
            exon_seq_corrected='CTCCACCACACCAT',
            read_seq_corrected='CTCCACCACACCAT',
            exon_length=exonlen.dt[exon_id=='exon_17', exon_length]
         )
      }
}, by=read_id
   ]), all=T, by=names(blast_fixed.dt))
   blast_fixed.dt[exon_comb_short_rc %like% '1,18',
                  exon_comb_short_rc := gsub('1,18', '1,17-18',
                                             exon_comb_short_rc)]
}


# blast_fixed.dt = merge(blast_fixed.dt,
#   blast_fixed.dt[exon_comb_short_rc=='1,17-18',
#    list(
#    exon_id='exon_02',
#    read_start_pos_align=unique(.SD[exon_id=='exon_01', read_end_pos_align+1]),
#    read_end_pos_align=unique(.SD[exon_id=='exon_17', read_start_pos_align+2]),
#    exon_start_pos_align=1, exon_end_pos_align=7, alignment_length=7,
#    read_seq='GTACCCA', exon_seq='GTACCCA', exon_comb_short_rc='1-2,17-18',
#    read_count=unique(read_count),
#    exon_seq_corrected='GTACCCA', read_seq_corrected='GTACCCA',
#    exon_length=exonlen.dt[exon_id=='exon_02', exon_length]
#    ), by=read_id
# ], all=T, by=names(blast_fixed.dt))
# blast_fixed.dt[exon_comb_short_rc=='1,17-18', exon_comb_short_rc := '1-2,17-18']
#

if (nrow(blast_fixed.dt[exon_comb_short_rc %like% '1-3,15-18']) > 0) {
   blast_fixed.dt = merge(blast_fixed.dt,
     blast_fixed.dt[exon_comb_short_rc=='1-3,15-18',
      list(
      exon_id='exon_14',
      read_start_pos_align=unique(.SD[exon_id=='exon_03', read_end_pos_align-6]),
      read_end_pos_align=unique(.SD[exon_id=='exon_15', read_start_pos_align-1]),
      exon_start_pos_align=exonlen.dt[exon_id=='exon_14', exon_length-20],
      exon_end_pos_align=exonlen.dt[exon_id=='exon_14', exon_length],
      alignment_length=20,
      read_seq='AGCCAACACAGAAAACGAAG', exon_seq='AGCCAACACAGAAAACGAAG',
      exon_comb_short_rc='1-3,14-18',
      read_count=unique(read_count),
      exon_seq_corrected='AGCCAACACAGAAAACGAAG',
      read_seq_corrected='AGCCAACACAGAAAACGAAG',
      exon_length=exonlen.dt[exon_id=='exon_14', exon_length]
      ), by=read_id
   ], all=T, by=names(blast_fixed.dt))
   blast_fixed.dt[exon_comb_short_rc=='1-3,15-18', exon_comb_short_rc := '1-3,14-18']
}

# Re-sort
blast_fixed.dt = blast_fixed.dt[nchar(read_seq_corrected) == nchar(exon_seq_corrected)]
blast_fixed.dt = blast_fixed.dt[order(read_id, read_start_pos_align)]

cat(' OK.\n')

# Now compare exon_seq_corrected and read_seq_corrected for each entry.
# exon_seq corrected is just exon_seq with insertions removed
# because the insertion on exon was made due to incorrect homopolymer run
# in the read.
#
#

cat('Searching for SNVs and indels...')
seq_diff_search = function (s1, s2) {
   # return a data table of differences between s1 and s2
   # we use s1 as exon, s2 as read for snv/indel annotation
   # pos = position of the snv or indel in the aligned string
   #       to obtain position on the exon, we need to add this
   #       to exon_start_pos_align
   # type = one of types in annotation (SNV, indel, ...)
   # nt = if type is indel or SNV write the new nt
   # hmp = if its indel, check the previous letter in the other string, i.e.
   #       exon: ATG--CCTAAAGT (s1)
   #       read: ATGGGCCTAA-CT (s2)
   #                ^     ^
   #               ins    del
   #       for ins, check read (s2) at i-1 or i+1 position
   #                vs at i position
   #       for del, check exon (s1) at i-1 or i+1 position
   #                vs at i position
   #       if they are the same then indel is a possible
   #       homopolymer error and set it to 1
   l = list()
   j = 1  # index for the list l, increment it each time you add to l manually.
   k = 0  # numer of skipped indices due to insertions in the read.
   if (nchar(s1) != nchar(s2)) {
      stop("Lengths of input sequences are not identical.")
   }
   s1.v = strsplit(s1, '')[[1]]
   s2.v = strsplit(s2, '')[[1]]
   n = length(s1.v)
   for (i in 1:n) {
      if (s1.v[i] != s2.v[i]) {
         if (s1.v[i] == '-') { # Insertion (in exon)
            # If it is still the same insertion as from the previous nt, skip it.
            k = k + 1
            if (s1.v[i-1] == '-') {
               next
            }
            hmp = ifelse(s2.v[i]==s2.v[i-1] | s2.v[i]==s2.v[i+1], 1, 0)
            l[[j]] = list(i-k, annotation[1], s1.v[i], s2.v[i], hmp)
            j = j + 1
            next
         }
         if (s2.v[i] == '-') { # Deletion
            hmp = ifelse(s1.v[i]==s1.v[i-1] | s1.v[i]==s1.v[i+1], 1, 0)
            l[[j]] = list(i-k, annotation[2], s1.v[i], '', hmp)
            j = j + 1
            next
         }
         # ----------------- SNV ------------------------------
         l[[j]] = list(i-k, annotation[3], s1.v[i], s2.v[i], 0)
         j = j + 1
      } else { # This should not occur, but catch it anyway.
         if (s1.v[i] == '-' & s2.v[i] == '-') {
            l[[j]] = list(i-k, annotation[4], s1.v[i], '-', 0)
            j = j + 1
         }
      }
   }
   dt = rbindlist(l)
   #dt = l
   if (length(dt) > 0) {
      names(dt) = c('pos', 'type', 'old_nt', 'nt', 'hmp')
      return(dt)
   }
}


snv_indel.dt = blast_fixed.dt[, seq_diff_search(exon_seq_corrected,
                                                read_seq_corrected),
                          by=c('read_id', 'exon_id', 'read_count',
                               'exon_start_pos_align',
                               'exon_comb_short_rc', 'read_start_pos_align')]
snv_indel.dt[, exon_pos := pos + exon_start_pos_align - 1]

# Find and extract big deletions which are typically near exon joins.
if (nrow(blast_fixed.dt[exon_end_pos_align < exon_length]) > 0) {
  dels_suf.dt = blast_fixed.dt[exon_end_pos_align < exon_length,
                       .(exon_pos = seq(exon_end_pos_align+1, exon_length),
                         type='deletion', nt='', hmp=0),
                       by=c('read_id', 'exon_id', 'read_start_pos_align',
                            'read_count', 'exon_start_pos_align',
                            'exon_comb_short_rc')]
  # dels_suf.dt[, c('read_start_pos_align') := NULL]
  dels_suf.dt[, pos := exon_pos]
}

if (nrow(blast_fixed.dt[exon_start_pos_align > 1]) > 0) {
  dels_pre.dt = blast_fixed.dt[exon_start_pos_align > 1,
     .(exon_pos = seq(1, exon_start_pos_align), type='deletion', nt='', hmp=0),
     by=c('read_id', 'exon_id', 'read_start_pos_align', 'read_count',
          'exon_start_pos_align', 'exon_comb_short_rc')]
  # dels_pre.dt[, c('read_start_pos_align') := NULL]
  dels_pre.dt[, pos := exon_pos]
}
cat(' OK.\n')

# Generate an empty table first.
# Exclude exon 8?
# exon_coords.dt = exonlen.dt[order(exon_id)][!(exon_id %like% '_rc'
# | exon_id %like% 'exon_08')]
exon_coords.dt = exonlen.dt[order(exon_id)][!(exon_id %like% '_rc')]
exon_coords.dt[, exon_start_abs := c(0, cumsum(exon_length)[1:(.N-1)])+1]
exon_coords.dt[, exon_end_abs := exon_start_abs + exon_length - 1]
exon_coords.dt[, exon_no := as.integer(gsub('_rc', '', gsub('exon_(.*).*',
                                                            '\\1', exon_id)))]
# exon_rc_coords.dt = exonlen.dt[order(-exon_id)][exon_id %like%
# '_rc' & !(exon_id %like% 'exon_08')]
# exon_rc_coords.dt[, exon_start_abs := c(0, cumsum(exon_length)[1:(.N-1)])+1]
setkey(exon_coords.dt, exon_no)

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



# Once we are at this, fix the snv_indel.dt table as well and export it
# fixed.
snv_indel_fix.dt = copy(snv_indel.dt)
snv_indel_fix.dt[, exon_no := as.integer(gsub('_rc', '',
   gsub('exon_([0-9]+).*', '\\1', exon_id)))]
snv_indel_fix.dt = merge(
   snv_indel_fix.dt,
   exonlen.dt[!(exon_id %like% 'rc'), .(exon_no, exon_length)],
   all.x=T, by='exon_no'
)
snv_indel_fix.dt = merge(
   snv_indel_fix.dt,
   exon_coords.dt[, .(exon_no, exon_start_abs, exon_end_abs)],
   all.x=T, by='exon_no'
)
snv_indel_fix.dt[exon_id %like% 'rc', exon_pos := exon_length - exon_pos + 1]
snv_indel_fix.dt[, exon_abs_pos := exon_start_abs + exon_pos - 1]
snv_indel_fix.dt[, mut_new := nt]
snv_indel_fix.dt[exon_id %like% 'rc', mut_new := chartr('ACGT', 'TGCA', mut_new)]
snv_indel_fix.dt[exon_id %like% 'rc', old_nt := chartr('ACGT', 'TGCA', old_nt)]
snv_indel_fix.dt[, nt := NULL]


# Let's manually generate a frequency table for each exon and for each exon
# combination.
# First, for each exon, fill missing values?
snv_ft.dt = snv_indel_fix.dt[type=='SNV', .(count_unique=.N,
                                        sum_unique=sum(read_count)),
                         by=c('exon_no', 'exon_pos')]
#snv_ft.dt[, exon_no := as.integer(gsub('_rc', '', gsub('exon_(.*).*', '\\1',
#                                                       exon_id)))]
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
#ins_ft.dt[, exon_no := as.integer(gsub('_rc', '', gsub('exon_(.*).*',
#                                                       '\\1', exon_id)))]
ins_ft.dt = merge(ins_ft.dt, exonlen.dt[!(exon_id %like% 'rc'),
   .(exon_no, exon_length)], all.x=T, by='exon_no')
ins_ft.dt = merge(ins_ft.dt, exon_coords.dt[,
   .(exon_no, exon_start_abs, exon_end_abs)], all.x=T, by='exon_no')
# ins_ft.dt[exon_id %like% 'rc', exon_pos := exon_length - exon_pos]
ins_ft.dt[, exon_pos_abs := exon_start_abs + exon_pos - 1]

# Add missing chunks of each exon as deletions!
# dels_all.dt = rbindlist(list(snv_indel.dt[type=='deletion'],
# dels_pre.dt, dels_suf.dt))
# dels_all.dt = merge(snv_indel_fix.dt[type=='deletion'], dels_pre.dt, all=T, by=intersect(names(snv_indel_fix.dt),
#                                                                                         names(dels_pre.dt)))
# dels_all.dt = merge(dels_all.dt, dels_suf.dt, all=T)

del_ft.dt = snv_indel_fix.dt[, .(count_unique=.N, sum_unique=sum(read_count)),
                         by=c('exon_no', 'exon_pos')]
#del_ft.dt[, exon_no := as.integer(gsub('_rc', '', gsub('exon_(.*).*', '\\1',
#                                                       exon_id)))]
del_ft.dt = merge(del_ft.dt, exonlen.dt[!(exon_id %like% 'rc'),
   .(exon_no, exon_length)], all.x=T, by='exon_no')
del_ft.dt = merge(del_ft.dt, exon_coords.dt[,
   .(exon_no, exon_start_abs, exon_end_abs)], all.x=T, by='exon_no')
# del_ft.dt[exon_id %like% 'rc', exon_pos := exon_length - exon_pos]
del_ft.dt[, exon_pos_abs := exon_start_abs + exon_pos - 1]

del_nopresuf_ft.dt = snv_indel_fix.dt[type=='deletion',
   .(count_unique=.N, sum_unique=sum(read_count)),
   by=c('exon_no', 'exon_pos')]
#del_nopresuf_ft.dt[, exon_no := as.integer(gsub('_rc', '',
#   gsub('exon_(.*).*', '\\1', exon_id)))]
del_nopresuf_ft.dt = merge(del_nopresuf_ft.dt,
   exonlen.dt[!(exon_id %like% 'rc'), .(exon_no, exon_length)],
   all.x=T, by='exon_no')
del_nopresuf_ft.dt = merge(del_nopresuf_ft.dt, exon_coords.dt[,
   .(exon_no, exon_start_abs, exon_end_abs)], all.x=T, by='exon_no')
# del_nopresuf_ft.dt[exon_id %like% 'rc', exon_pos := exon_length - exon_pos]
del_nopresuf_ft.dt[, exon_pos_abs := exon_start_abs + exon_pos - 1]
#


ceil = function (x) {
   # My ceiling function
   # Calculate number of significant digits
   dig = floor(log10(x))
   return(ceiling(x/10^dig)*10^dig)
}

flr = function (x) {
   # My flooring function
   # Calculate number of significant digits
   dig = floor(log10(x))
   return(floor(x/10^dig)*10^dig)
}

# Generate a picture of the cDNA 1-7,9-18, then on it show SNVs indels
ymax = 1
exons_blocks.dt = exon_coords.dt[,
   .(x = c(exon_start_abs, exon_start_abs, exon_end_abs, exon_end_abs),
   y = c(0, ymax, ymax, 0)), by=exon_id]
# Generate text and coordinates to put in the blocks
exons_text.dt = exons_blocks.dt[, .(x=mean(x), y=0.7,
   text=as.integer(gsub('exon_', '', unique(exon_id)))), by=exon_id]
plot_text_size = 20
plot_pt_size = 0.5
exon_colors = rainbow(length(exons_text.dt[, text]), end=0.8)
names(exon_colors) = exons_text.dt[, text]
yminplot = min(snv_ft.dt[, min(count_unique)],
               del_ft.dt[, min(count_unique)],
               ins_ft.dt[, min(count_unique)])
ymaxplot = max(snv_ft.dt[, max(count_unique)],
               del_ft.dt[, max(count_unique)],
               ins_ft.dt[, max(count_unique)])

# Plot all indels and SNVs on one plot for unique reads
# ft.plot = ggplot(
#    # snv_ft.dt[exon_comb_short_rc %in% exon_comb_display],
#    snv_ft.dt,
#    aes(x=exon_pos_abs, y=count_unique)) +
#    annotation_logticks(sides='l') +
#    geom_polygon(data=exons_blocks.dt,
#                 mapping=aes(x=x, y=y*ceil(ymaxplot), group=exon_id, fill=exon_id),
#                 alpha=0.2) +
#    geom_polygon(data=exons_blocks.dt,
#                 mapping=aes(x=x, y=y, group=exon_id, fill=exon_id)) +
#    geom_vline(data=exon_coords.dt[, .(exon_midpt = exon_end_abs[1:(.N-1)])],
#               mapping=aes(xintercept=exon_midpt),
#               color='white') +
#    # Add vertical white lines between exons for legibility
#    # geom_text(data=exons_text.dt, mapping=aes(x=x, y=y, label=text)) +
#    geom_point(data=ins_ft.dt,
#               # color=rgb( 39/255, 164/255, 255/255),
#               mapping=aes(color='insertions'),
#               size=plot_pt_size) + # Plot insertions
#    geom_point(data=del_ft.dt,
#               #color=rgb(247/255, 118/255, 109/255),
#               mapping=aes(color='deletions'),
#               size=plot_pt_size) + # Plot deletions
#    geom_point(mapping=aes(color='SNVs'), size=plot_pt_size) +
#    # geom_point(mapping=aes(color=exon_comb_short_rc)) + # Plots SNVs
#    scale_x_continuous(name='exon', expand=c(0, 0), breaks=exons_text.dt[, x],
#                       labels=exons_text.dt[, text]) +
#    scale_y_continuous(name='unique counts', trans='log10', breaks=10^(0:3),
#                       #limits=c(flr(yminplot), ceil(ymaxplot)),
#                       expand=c(0, 0)) +
#    # scale_color_manual(values=rainbow(length(exon_comb_display), end=0.8)) +
#    # scale_shape_discrete(solid=F) +
#    scale_fill_manual(values=rep(c('gray', 'white'),
#                                 length.out = length(exons_text.dt[, text])),
#                      guide=F) +
#    scale_color_manual(name='', values=c('insertions'='red',
#                                         'deletions'='blue',
#                                         'SNVs'='black')) +
#    #$ scale_fill_discrete(guide=F) +
#    coord_cartesian(ylim=c(yminplot-0.3, ceil(ymaxplot))) +
#    theme(
#          text = element_text(size=plot_text_size),
#          axis.text = element_text(color='black'),
#          # axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
#          axis.ticks.x = element_blank(),
#          panel.background = element_rect(fill=NA, color='black'),
#          panel.ontop = T,
#          # legend.position = 'none',
#          panel.grid.minor = element_blank(),
#          panel.grid.major = element_blank()
#    )
# ft.plot
# ggsave(paste0(output_folder, input_base, '_snv_indel_map_log10.pdf'),
#        ft.plot, width=10, height=3)
# 
cat('OK.\n')

# Save this snv_indel.dt table
cat('Saving SNV/indel table...')
write.table(snv_indel_fix.dt, paste0(gsub('exon_table_fixed.txt', '', input_filename),
   'snv_indel_table.txt'), sep='\t', quote=F, row.names=F)
cat(' OK.\n')

# Show length distributions of insertions and deletions
cat('Generating plots for length distributions...')
ins_len.c = Reduce(c, lapply(str_locate_all(
   blast_fixed.dt[exon_seq_corrected %like% '-', exon_seq_corrected], '-+'),
                             function (x) return(x[,2]-x[,1]+1)))

if (!is.null(ins_len.c)) {
   ins_len_ft.dt = data.table(table(ins_len.c))
   ins_len_ft.dt[, length := as.integer(ins_len.c)]
   ins_len_hist.plot = ggplot(ins_len_ft.dt, aes(x=length, y=N)) +
      geom_bar(stat='identity') +
      scale_x_continuous(name='length (nt)', breaks=1:max(ins_len_ft.dt[, length])) +
      scale_y_continuous(name='count')
   ggsave(paste0(output_folder, input_base, '_insertions_length_hist.pdf'),
          ins_len_hist.plot, width=5, height=4)
}

# Now find length distribution for deletions.
# Find contiguous deletions, i.e. from vector v
# v = (1  2  3  6  9 13 14 17 25 28 29 30)
# Output: 3 1 1 2 1 1 3 as a list
find_lengths = function (v) {
   if (length(v) == 1) return (list(1))
   lens = list()
   run_len = 1 # Keep track of run length
   j = 1 # List lens index
   for (i in 1:(length(v)-1)) {
      current = v[i]
      if (current == v[i+1]-1) {
         run_len = run_len + 1
         next
      } else {
         lens[[j]] = run_len
         j = j + 1
         run_len = 1
      }
   }
   lens[[j]] = run_len
   return(lens)
}

# Find all deletion lengths
if (nrow(snv_indel.dt[type=='deletion']) > 0) {
  del_len.dt = snv_indel.dt[type=='deletion', .(del_len = as.integer(find_lengths(sort(exon_pos))),
                               exon_comb_short_rc = unique(exon_comb_short_rc)),
                           by=c('read_id', 'exon_id')]

  # Count them into a frequency table manually
  del_len_ft.dt = del_len.dt[, .(count = .N), by=c('del_len', 'exon_comb_short_rc')]
  # Set missing values to zero
  del_len_ft2.dt = rbindlist(list(
     del_len_ft.dt,
     del_len_ft.dt[, .(del_len = setdiff(seq(1,max(del_len)), del_len), count = 0),
                   by=exon_comb_short_rc]
     ))
  del_len_ft2.dt = del_len_ft2.dt[order(del_len)]

  # Plot them
  del_len_ft.plot = ggplot(del_len_ft.dt, aes(x=del_len, y=count)) +
     geom_bar(stat='identity') +
     scale_x_continuous(name='length (nt)') +
     scale_y_continuous(name='count')
  ggsave(paste0(output_folder, input_base, '_deletions_length_hist.pdf'),
         del_len_ft.plot, width=5, height=4)
}


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

cat('Combining SNVs and indels into a single table...')
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



no_bins = 5
# count_colors_all = brewer.pal(no_bins+1, 'YlGnBu')
# count_colors = count_colors_all[2:length(count_colors_all)]
# count_colors = brewer.pal(no_bins, 'Greys')

# Colors for total counts
#count_colors = brewer.pal(8)
colfunc = colorRampPalette(c('#2D71AB', 'white'))
# count_colors = rev(colfunc(5))
count_colors = rev(brewer.pal(no_bins, 'RdYlBu'))

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


# ----------------------- Load known AD variants ----------------------------
ad_var.dt = fread('data/app_ad_known_variants_3.txt')

# Load FASTA file with exon sequences
library(Biostrings)
exons_fa = readDNAStringSet('fa/app_exons_cds_clean.fa')
exon_names = sapply(names(exons_fa), function (x) strsplit(x, ' ')[[1]][1], USE.NAMES=F)
exons_seq.dt = data.table(exon_id=exon_names, seq=as.character(exons_fa))

# Amino acid to nucleotide
# Amino acid 201,
ad_var.dt[, exon_abs_pos := aa_pos*3 - 2 + mut_offset]
ad_var.dt[, exon_id := exon_coords.dt[exon_abs_pos >= exon_start_abs &
                            exon_abs_pos <= exon_end_abs, exon_id], by=1:nrow(ad_var.dt)]
ad_var.dt[, exon_pos := exon_abs_pos -
             exon_coords.dt[exon_abs_pos >= exon_start_abs &
                            exon_abs_pos <= exon_end_abs, exon_start_abs] + 1,
          by=1:nrow(ad_var.dt)]

test = c()
for (i in 1:(nrow(ad_var.dt))) {
   row = ad_var.dt[i]
   test = c(test, str_sub(exons_seq.dt[exon_id==row$exon_id]$seq,
                          start = row$exon_pos - row$mut_offset,
                          end = row$exon_pos - row$mut_offset + 2))
}

# Find variant matches by inner joining alzforum list with our results
# var_matches.dt = merge( ad_var.dt,
#                         snv_indel.dt[, .(
#                            count_unique = .N,
#                            sum_unique = sum(read_count)
#                            ),
#                            by=.(exon_pos, exon_id, nt, type, hmp)],
#                         by=c('exon_id', 'exon_pos')
#                       )
# # Extract only positions where the actual nt is matched and mut_type
# var_matches.dt = var_matches.dt[(mut_new == nt | (mut_new == '-' & nt == '')) & mut_type == type]
# var_matches.dt[, exon_no := as.integer(strsplit(exon_id, '_')[[1]][2])]

# Find reads with these variants
# var_snv_indel.dt = rbindlist(lapply(1:nrow(var_matches.dt), function (i) {
#    dt = var_matches.dt[i]
#    snv_indel_fix.dt[exon_pos==dt$exon_pos &
#                     exon_no==dt$exon_no &
#                     nt==dt$nt &
#                     type==dt$type]
# }))

ad_var.dt[, exon_no := as.integer(gsub('exon_([0-9]+).*', '\\1', exon_id))]
var_match.dt = merge(
   ad_var.dt[, .(exon_no, exon_abs_pos, type=mut_type, mut_new,
                      `Clinical Phenotype`, `Primary  Papers`)],
   snv_indel_fix.dt[, .(exon_no, read_id, read_count, exon_comb_short_rc, type,
                        exon_abs_pos, mut_new)],
   by=c('exon_no', 'exon_abs_pos', 'type', 'mut_new'))
var_match.dt = var_match.dt[`Clinical Phenotype` %like% 'Alzheimer']

# Check for other variants at these positions
var_match_similar.dt = merge(
   ad_var.dt[, .(exon_no, exon_abs_pos, type=mut_type, mut_new_ad=gsub('-', '', mut_new),
                      `Clinical Phenotype`, `Primary  Papers`, `Mutation Type / Codon Change`)],
   snv_indel_fix.dt[, .(exon_no, read_id, read_count, exon_comb_short_rc, type,
                        exon_abs_pos, mut_new)],
   by=c('exon_no', 'exon_abs_pos', 'type'))
var_match_similar.dt = var_match_similar.dt[`Clinical Phenotype` %like% 'Alzheimer']

# Extract old/new codons from AD reference
var_match_similar.dt[, codon_old := gsub('.*([A|C|G|T]{3}) to.*', '\\1', `Mutation Type / Codon Change`)]
var_match_similar.dt[, codon_new_ref := gsub('.*to ([A|C|G|T]{3}).*', '\\1', `Mutation Type / Codon Change`)]

# Fix missing old/new codons
var_match_similar.dt[nchar(codon_new_ref) > 3, codon_new_ref := NA]
var_match_similar.dt[nchar(codon_old) > 3, codon_old := NA]
var_match_similar.dt = var_match_similar.dt[, {
   if (.N > 1) {
      .SD[mut_new_ad==mut_new]
   } else {
      .SD
   }
}, by=.(read_id, exon_abs_pos)]

# Generate codon to aa reference table from codon.dt
# Do this once we have the reading frame analysis in unique_reads_comb_circo_v2.R.
# aa_to_codon.dt = codon.dt[, .(codon = strsplit(codons, ', ')[[1]]), by=aa_1let]

# Manually check the new


# Check
# var_match_similar.dt[, codon_new_data := codon_old]


if (nrow(var_match.dt) > 0) {
   cat('Saving identified known variants...')
   write.table(var_match.dt, paste0(gsub('exon_table_fixed.txt', '', input_filename),
      'known_ad_variants.txt'), sep='\t', quote=F, row.names=F)
   cat(' OK.\n')

   var_matches.dt = var_match.dt[,
      .(count_unique=.N), by=.(exon_no, exon_abs_pos)]
}

if (nrow(var_match_similar.dt) > 0) {
   cat('Saving identified known and similar variants...')
   write.table(var_match_similar.dt, paste0(gsub('exon_table_fixed.txt', '', input_filename),
      'known_similar_ad_variants.txt'), sep='\t', quote=F, row.names=F)
   cat(' OK.\n')
}


cat('Generating circo plot...')
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

cairo_pdf(paste0(output_folder, input_base, '_circos_plot.pdf'), width=5, height=5)


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

# dels_big.dt = del_presuf_ft.dt[exon_pos_abs > exon_start_abs + (snv_indel_ft_abs_offset + 2) &
#                     exon_pos_abs < exon_end_abs - (snv_indel_ft_abs_offset + 2),
#                     .(count_unique = log10(as.numeric(sum(count_unique)))),
#                     by=c('exon_no', 'exon_pos_abs', 'exon_pos_abs')]

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

# Plot big deletions (commented after talking to Tao)
# for (i in 1:length(count_ivals_log10)) {
#    iv = count_ivals_log10[[i]]
#    if (i == length(count_ivals_log10)) iv[2] = 1e6
#    dt = dels_big.dt[count_unique >= iv[1] & count_unique < iv[2]]
#    for (j in dt[, unique(exon_no)]) {
#       circos.points(x=dt[exon_no==j,exon_pos_abs], y=rep(-0.08, dt[exon_no==j, .N]),
#                     col=count_colors[i],
#                     cex=0.15, track.index=2, sector.index=j)
#    }
# }

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

# Then the "deletions" which represent missing chunks of these exons
# (commented out after talking to Tao)
# circos.genomicTrackPlotRegion(
#    del_presuf_ft.dt[exon_pos_abs > exon_start_abs + (snv_indel_ft_abs_offset + 2) &
#                     exon_pos_abs < exon_end_abs - (snv_indel_ft_abs_offset + 2),
#                     .(count_unique = log10(as.numeric(sum(count_unique)))),
#                     by=c('exon_no', 'exon_pos_abs', 'exon_pos_abs')],
#    track.index = 5,
#    panel.fun = function (region, value, ...) {
#       circos.genomicPoints(region, value, ..., pch=16, cex=0.25,
#                            col=count_colors[1])
#    },
#    ylim = c(0, 4)
# )
# circos.yaxis(at = c(0, 2, 4), sector.index = 1, track.index = 5,
#              labels.cex=y_lab_size,
#              labels = c(1, 100, '10,000'),
#              labels.niceFacing = T)

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
# exon_comb_display = exon_join_u_top.dt[, unique(exon_comb_short_rc)]

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
dev.off()
cat(' OK.\n')


# Write output tables
cat('Writing output tables...')
write.table(blast_fixed.dt, paste0(output_data_folder, input_base,
   '_blast_results_fixed.txt'), sep='\t', quote=F, row.names=F)
write.table(exon_coords.dt, file.path(output_data_folder, paste0(input_base, '_exon_coordinates.txt')),
   sep='\t', quote=F, row.names=F)
write.table(var_match.dt, file.path(output_data_folder, paste0(input_base, '_variant_matches.txt')),
            sep='\t', quote=F, row.names=F)
# snv_indel.dt[, exon_no := sapply(exon_id,
#    function (e) as.numeric(strsplit(e, '_')[[1]][2]))]
# write.table(snv_indel.dt, paste0(output_data_folder, input_base,
#    'all_snv_indels.txt', sep='\t', quote=F, row.names=F))
cat(' OK.\n')

# Generate all unique sequences with exon joins?
tmp = blast_fixed.dt[, {
   if (exon_no[1] == 18) {
      as.character(reverseComplement(DNAString(paste(read_seq_corrected, collapse=''))))
   } else {
      paste0(read_seq_corrected, sep='')
   }
}, by=read_id]

