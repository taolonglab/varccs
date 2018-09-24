#
# Fix homopolymer runs with exon sequences where the quality is poor. This is only used
# to find SNVs and count unique sequences.
#

# rm(list=ls(all=TRUE))

# Load required libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Biostrings))

options(warn = -1)

# Argument defaults
# input_filename_def = 'data/ad-31-150_filtered_qs85_hpf_qs30_rep2_unique_blast_results_ws25_go0_gx2_exon_table.txt'
# input_fasta_def = 'reads/ad-31-150_filtered_qs85_hpf_qs30_rep2_unique.fa'
# output_folder_def = 'reads/'

# Parse user arguments
option_list = list(
    make_option(c('-i', '--input'),
                help = 'Input file from blast results (output from scripts/blast_reads_vs_exons)'),
    make_option(c('-f', '--input-fasta'),
                help = 'Input FASTA file to fix homopolymer runs with low quality.'),
    make_option(c('-t', '--output-table'),
                help = 'Output table filename.'),
    make_option(c('-o', '--output-fasta'),
                hel = 'Output fasta filename.')
)
opt = parse_args(OptionParser(option_list=option_list))

# Load the alignment table
input_filename = opt$input
input_fasta = opt$`input-fasta`
output_table = opt$`output-table`
output_fasta = opt$`output-fasta`
# output_folder = opt$output
# readlen_filename = opt$read_lengths
# exonlen_filename = opt$exon_lengths

# input_base = gsub('\\..*', '', basename(input_filename))

# Append slash if not part of folder name
# if (str_sub(output_folder, -1) != '/') {
#    output_folder = paste0(output_folder, '/')
# }
cat('* Fixing low quality homopolymer runs\n')
cat('  - loading BLAST file...')
blast.dt = fread(input_filename)
cat(' OK.\n')

cat('  - loading FASTA file...')
fa = readDNAStringSet(input_fasta)
fa_meta = names(fa)
fa_readids = gsub('^(.*);.*', '\\1', fa_meta)
fa_counts = as.integer(gsub('.*;count_([0-9]+).*', '\\1', fa_meta))
fa.dt = data.table(read_id=as.character(fa_readids), count=fa_counts, seq=as.character(fa))
setkey(fa.dt, read_id)
cat(' OK.\n')

# Process FASTA file, generate meta-data

if (nrow(blast.dt[read_seq %like% 'N']) > 0) {
   # Generate a new corrected read based on the same criteria
   hmp_fix.dt = blast.dt[read_seq %like% 'N',
            .(read_id, exon_id, exon_no, read_start_pos_align, read_end_pos_align,
              exon_start_pos_align, exon_end_pos_align, alignment_length,
              read_seq, exon_seq, exon_comb_short_rc, read_count,
              read_coverage, exon_length, exon_comb_coverage,
              exon_comb_order,
              read_seq_corrected = {
      x_reg = gregexpr('(-+[ACGT]N|[ACGT]N-+|[ACGT]-+N)', read_seq)
      # List where each element is another list containing lengths of each match
      x_len = lapply(x_reg, attr, 'match.length')
      # List where each element is a list containing coordinates of first and last
      # letter for a regex match.
   
      x_int = mapply(function (r, l) list(transpose(list(r, r+l-1))), x_reg, x_len)
      # Total number of sequences to check = .N
      n = .N
      #print(length(read_seq))
      seqs = list()
      for (i in 1:n) { # i = read sequence index, can use it to get exons too
         x_inti = x_int[[i]]
         # x_inti = x_int
         if (x_inti[[1]][1] == -1) {
            # No regex match, just return the sequence
            seqs[[i]] = read_seq[[i]]
            next
         }
         m = length(x_inti) # number of index in this sequence
         # First extract the prefix that is not changed
         if (x_inti[[1]][1] > 1) {
            seq = str_sub(read_seq[i], 1, x_inti[[1]][1]-1)
         } else {
            seq = ""
         }
         # Then every change just take the exonic sequence instead.
   
         for (j in 1:m) {
            seq = paste0(seq,
               gsub('-', '', str_sub(exon_seq[i], x_inti[[j]][1], x_inti[[j]][2]))
            )
            # Sequence from the read after this match until the next match with Ns
            if (j < m) {
               seq = paste0(seq, str_sub(read_seq[i], x_inti[[j]][2]+1, x_inti[[j+1]][1]-1))
            }
         }
         # Finally finish with the read_seq suffix that is unchanged
         if (x_inti[[m]][2] < nchar(read_seq[i])) {
            seq = paste0(seq, str_sub(read_seq[i], x_inti[[m]][2]+1))
         }
         seqs[[i]] = seq
      }
      seqs
   },
         exon_seq_corrected = {
      x_reg = gregexpr('(-+[ACGT]N|[ACGT]N-+|[ACGT]-+N)', read_seq)
      # List where each element is another list containing lengths of each match
      x_len = lapply(x_reg, attr, 'match.length')
      # List where each element is a list containing coordinates of first and last
      # letter for a regex match.
      x_int = mapply(function (r, l) list(transpose(list(r, r+l-1))), x_reg, x_len)
      # Total number of sequences to check = .N
      n = .N
      #print(length(read_seq))
      seqs = list()
      for (i in 1:n) { # i = read sequence index, can use it to get exons too
         x_inti = x_int[[i]]
         if (x_inti[[1]][1] == -1) {
            # No regex match, just return the sequence
            seqs[[i]] = exon_seq[[i]]
            next
         }
         m = length(x_inti) # number of indexes in this sequence
         # First extract the prefix that is not changed
         if (x_inti[[1]][1] > 1) {
            seq = str_sub(exon_seq[i], 1, x_inti[[1]][1]-1)
         } else {
            seq = ""
         }
         # Then every change just take the exonic sequence instead
         for (j in 1:m) {
            # If the exonic sequence has '-' that means read sequence had an insertion.
            # Then remove '-' from exon
            seq = paste0(seq,
               gsub('-', '', str_sub(exon_seq[i], x_inti[[j]][1], x_inti[[j]][2]))
            )
            # Sequence from the read after this match until the next match with Ns
            if (j < m) {
               seq = paste0(seq, str_sub(exon_seq[i], x_inti[[j]][2]+1, x_inti[[j+1]][1]-1))
            }
         }
         # Finally finish with the read_seq suffix that is unchanged
         if (x_inti[[m]][2] < nchar(exon_seq[i])) {
            seq = paste0(seq, str_sub(exon_seq[i], x_inti[[m]][2]+1))
         }
         seqs[[i]] = seq
      }
      seqs
         })]
   
   # Generate a new corrected read based on the same criteria
   hmp_fix.dt = hmp_fix.dt[,
            .(read_id, exon_id, exon_no, read_start_pos_align, read_end_pos_align,
              exon_start_pos_align, exon_end_pos_align, alignment_length,
              read_seq, exon_seq, exon_comb_short_rc, read_count,
              read_coverage, exon_length, exon_comb_coverage,
              exon_comb_order,
              read_seq_corrected = {
      x_reg = gregexpr('[ACGT][ACGT]N[ACGT]', read_seq_corrected)
      # List where each element is another list containing lengths of each match
      x_len = lapply(x_reg, attr, 'match.length')
      # List where each element is a list containing coordinates of first and last
      # letter for a regex match.
      x_int = mapply(function (r, l) list(transpose(list(r, r+l-1))), x_reg, x_len)
      # Total number of sequences to check = .N
      n = .N
      #print(length(read_seq))
      seqs = list()
      for (i in 1:n) { # i = read sequence index, can use it to get exons too
         x_inti = x_int[[i]]
         # x_inti = x_int
         if (x_inti[[1]][1] == -1) {
            # No regex match, just return the sequence
            seqs[[i]] = read_seq_corrected[[i]]
            next
         }
         m = length(x_inti) # number of index in this sequence
         # First extract the prefix that is not changed
         if (x_inti[[1]][1] > 1) {
            seq = str_sub(read_seq_corrected[i], 1, x_inti[[1]][1]-1)
         } else {
            seq = ""
         }
         # Then every change just take the exonic sequence instead
         for (j in 1:m) {
            seq = paste0(seq, gsub('-', '', str_sub(exon_seq[i], x_inti[[j]][1], x_inti[[j]][2])))
            # Sequence from the read after this match until the next match with Ns
            if (j < m) {
               seq = paste0(seq, str_sub(read_seq_corrected[i], x_inti[[j]][2]+1, x_inti[[j+1]][1]-1))
            }
         }
         # Finally finish with the read_seq suffix that is unchanged
         if (x_inti[[m]][2] < nchar(read_seq_corrected[i])) {
            seq = paste0(seq, str_sub(read_seq_corrected[i], x_inti[[m]][2]+1))
         }
         seqs[[i]] = seq
      }
      seqs
   },
      exon_seq_corrected = {
      x_reg = gregexpr('[ACGT][ACGT]N[ACGT]', read_seq_corrected)
      # List where each element is another list containing lengths of each match
      x_len = lapply(x_reg, attr, 'match.length')
      # List where each element is a list containing coordinates of first and last
      # letter for a regex match.
      x_int = mapply(function (r, l) list(transpose(list(r, r+l-1))), x_reg, x_len)
      # Total number of sequences to check = .N
      n = .N
      #print(length(read_seq))
      seqs = list()
      for (i in 1:n) { # i = read sequence index, can use it to get exons too
         x_inti = x_int[[i]]
         # x_inti = x_int
         if (x_inti[[1]][1] == -1) {
            # No regex match, just return the sequence
            seqs[[i]] = exon_seq_corrected[[i]]
            next
         }
         m = length(x_inti) # number of index in this sequence
         # First extract the prefix that is not changed
         if (x_inti[[1]][1] > 1) {
            seq = str_sub(exon_seq_corrected[i], 1, x_inti[[1]][1]-1)
         } else {
            seq = ""
         }
         # Then every change just take the exonic sequence instead
         for (j in 1:m) {
            seq = paste0(seq, gsub('-', '', str_sub(exon_seq[i], x_inti[[j]][1], x_inti[[j]][2])))
            # Sequence from the read after this match until the next match with Ns
            if (j < m) {
               seq = paste0(seq, str_sub(exon_seq_corrected[i], x_inti[[j]][2]+1, x_inti[[j+1]][1]-1))
            }
         }
         # Finally finish with the read_seq suffix that is unchanged
         if (x_inti[[m]][2] < nchar(exon_seq_corrected[i])) {
            seq = paste0(seq, str_sub(exon_seq_corrected[i], x_inti[[m]][2]+1))
         }
         seqs[[i]] = seq
      }
      seqs
      })]
   
   hmp_fix.dt[, read_seq_corrected := as.character(read_seq_corrected)]
   hmp_fix.dt[, exon_seq_corrected := as.character(exon_seq_corrected)]
   
} else {
   hmp_fix.dt = blast.dt[, .(read_id, exon_id, read_start_pos_align, read_end_pos_align,
              exon_start_pos_align, exon_end_pos_align, alignment_length,
              read_seq, exon_seq, exon_comb_short_rc, read_count,
              read_coverage, exon_length, exon_comb_coverage,
              exon_comb_order)]
   hmp_fix.dt[, read_seq_corrected := read_seq]
   hmp_fix.dt[, exon_seq_corrected := exon_seq]
}

cat('  - correcting reads...')
cat(' OK.\n')
setkey(hmp_fix.dt, read_id)
# print(head(hmp_fix.dt))

fa1.dt = hmp_fix.dt[fa.dt]
# fa1.dt[is.na(read_seq_corrected), read_seq_corrected := seq]
# print(head(fa1.dt))

fa2.dt = unique(fa1.dt[
   ,
   # !is.na(read_seq_corrected),
   .(seq, count, exon_comb_short_rc, seq_c = {
    if (is.na(read_seq_corrected)) {
       seq
    } else {
       s = ''
       if (read_start_pos_align[1] > 1) {
          s = paste0(s, str_sub(seq, 1, read_start_pos_align[1]-1))
       }
       for (i in 1:(.N-1)) {
          s = paste0(s, read_seq_corrected[i])
          s = paste0(s, str_sub(seq, read_end_pos_align[i]+1, read_start_pos_align[i+1]-1))
       }
       s = paste0(s, read_seq_corrected[.N])
       if (read_end_pos_align[.N] < nchar(seq)) {
          s = paste0(s, str_sub(seq, read_end_pos_align[.N]))
       }
       s
    }
}), by=read_id])

cat('  - saving output FASTA...')
fa_out = DNAStringSet(fa2.dt[, seq_c])
names(fa_out) = fa2.dt[, paste(read_id, ';count_', count, sep='')]
writeXStringSet(fa_out, output_fasta)
cat(' OK.\n')
cat('    ->', output_fasta, '\n')

# Just save the final table here and do the SNV / indel annotation
# in a separate script.
cat('  - saving output homopolymer fixed table...')
# write.table(hmp_fix.dt, paste0(gsub('.txt', '', input_filename),
#                                '_hmpfix.txt'),
#             sep='\t', quote=F, row.names=F)

blast_fixed.dt = blast.dt[!(read_seq %like% 'N'),
         .(read_id, exon_id, exon_no, read_start_pos_align, read_end_pos_align,
           exon_start_pos_align, exon_end_pos_align, alignment_length,
           read_seq, exon_seq, exon_comb_short_rc, read_count,
           read_coverage, exon_length, exon_comb_coverage,
           exon_comb_order)]
blast_fixed.dt[, exon_seq_corrected := exon_seq]
blast_fixed.dt[, read_seq_corrected := read_seq]

blast_fixed.dt = merge(blast_fixed.dt, hmp_fix.dt, by=names(hmp_fix.dt), all=T)

write.table(blast_fixed.dt, output_table,
            sep='\t', quote=F, row.names=F)
cat(' OK.\n')
cat('    ->', output_table, '\n')
