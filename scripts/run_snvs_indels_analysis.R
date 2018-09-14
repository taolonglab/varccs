# Analyze and make plots for SNVs and indels based on the exon table
# with fixed bad homopolymer runs.
#
# options(datatable.print.nrows = 10)
# Load exon table with fixed homopolymer sequences

cat('* Initializing R libraries...')
# Load required libraries
suppressPackageStartupMessages(library(data.table))
# suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(optparse))
# suppressPackageStartupMessages(library(Biostrings))
# options(warn = -1)


# Annotation for insertions, deletions and SNVs in this order
# '-' means both read and exon have '-' so it should be removed.
# (in principle should not happen at this stage.)
annotation = c('insertion', 'deletion', 'SNV', '-')
nts = c('A', 'C', 'G', 'T', '-')

#------------------------------- Parse user arguments -----------------------------

option_list = list(
    make_option(c('-i', '--input-exon-table'),
       help = 'Input file from fasta_homopolymer_fix.R (*_exon_table_hmpfix.txt'),
    make_option(c('-l', '--input-exon-lengths'),
       help = 'Input file with exon lengths. Generate in bash then import here.'),
    make_option(c('-f', '--input-fads'),
                # action = 'store_true',
                default = ' ',
                help = 'Manually prepared file of Familiar AD mutations (FADs),
                based on the Alzforum online table.'),
    make_option(c('-os', '--output-snvs-indels'),
                help = 'Output table of SNVs and indels.'),
    make_option(c('-of', '--output-fads'),
                default = ' ',
                help = 'Output file for Familial AD mutations (if any).')
)
opt = parse_args(OptionParser(option_list=option_list))

# Input / output filenames
input_exonlen = opt$`input-exon-lengths`
input_exontab = opt$`input-exon-table`
input_fads = opt$`input-fads`

output_snvs = opt$`output-snvs-indels`
output_fads = opt$`output-fads`
cat(' OK.\n')

#------------------------------- Load files -----------------------------------
cat('* Loading input files')
exonlen.dt = fread(input_exonlen)
exonlen.dt[, exon_no := as.integer(
   sub('_rc', '', sub('exon_(.*).*', '\\1', exon_id)))]
cat('.')
blast_fixed.dt = fread(input_exontab)
cat('.')
# blast_fixed.dt = merge(blast_fixed.dt, exonlen.dt[, .(exon_id, exon_length)],
#                        by='exon_id')
# cat('.')
if (input_fads != ' ') {
   cat(input_fads)
   if (output_fads != ' ') {
      ad_var.dt = fread(input_fads)
      cat('.')
   } else {
      stop('! ERROR: Output filename for any found FADs is not given.')
   }
}
# exons_fa = readDNAStringSet('fa/app_exons_cds_clean.fa')
# cat('.')
cat(' OK.\n')


# Load FASTA file with exon sequences


# Sort blast file by read_start_pos_align
blast_fixed.dt = blast_fixed.dt[order(read_id, read_start_pos_align)]


# Now compare exon_seq_corrected and read_seq_corrected for each entry.
# exon_seq corrected is just exon_seq with insertions removed
# because the insertion on exon was made due to incorrect homopolymer run
# in the read.
#
#

cat('* Searching for SNVs and indels...')
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
# cat(' OK.\n')

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
cat(' OK.\n')

# Save this snv_indel.dt table
cat('* Saving SNV/indel table...')
write.table(snv_indel_fix.dt, output_snvs, sep='\t', quote=F, row.names=F)
cat(' OK.\n')
cat('  ->', output_snvs, '\n')


# ----------------------- Load known AD variants ----------------------------
# exon_names = sapply(names(exons_fa), function (x) strsplit(x, ' ')[[1]][1], 
#                     USE.NAMES=F)
# exons_seq.dt = data.table(exon_id=exon_names, seq=as.character(exons_fa))

# Amino acid to nucleotide
# Amino acid 201,
if (input_fads != ' ') {
   cat('* Analyzing data for Familial AD mutations...')
   ad_var.dt[, exon_abs_pos := aa_pos*3 - 2 + mut_offset]
   ad_var.dt[, exon_id := exon_coords.dt[exon_abs_pos >= exon_start_abs &
                               exon_abs_pos <= exon_end_abs, exon_id], 
             by=1:nrow(ad_var.dt)]
   ad_var.dt[, exon_pos := exon_abs_pos -
                exon_coords.dt[exon_abs_pos >= exon_start_abs &
                               exon_abs_pos <= exon_end_abs, exon_start_abs] + 1,
             by=1:nrow(ad_var.dt)]
   
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
      ad_var.dt[, .(exon_no, exon_abs_pos, type=mut_type, 
                    mut_new_ad=gsub('-', '', mut_new),
                    `Clinical Phenotype`, `Primary  Papers`, 
                    `Mutation Type / Codon Change`)],
      snv_indel_fix.dt[, .(exon_no, read_id, read_count, exon_comb_short_rc, type,
                           exon_abs_pos, mut_new)],
      by=c('exon_no', 'exon_abs_pos', 'type'))
   var_match_similar.dt = var_match_similar.dt[`Clinical Phenotype` %like% 'Alzheimer']
   
   # Extract old/new codons from AD reference
   var_match_similar.dt[, codon_old := gsub('.*([A|C|G|T]{3}) to.*', '\\1', 
                                            `Mutation Type / Codon Change`)]
   var_match_similar.dt[, codon_new_ref := gsub('.*to ([A|C|G|T]{3}).*', '\\1', 
                                                `Mutation Type / Codon Change`)]
   
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
   cat(' OK.\n')
   
   if (nrow(var_match.dt) > 0) {
      cat('* Saving identified known variants...')
      write.table(var_match.dt, output_fads, sep='\t', quote=F, row.names=F)
      cat(' OK.\n')
      cat('  ->', output_fads, '\n')
   
      # var_matches.dt = var_match.dt[,
      #    .(count_unique=.N), by=.(exon_no, exon_abs_pos)]
   }
   
   # This outputs mutations occuring at same positions as FADs, but
   # not being identical (so 'similar'):
   #
   # if (nrow(var_match_similar.dt) > 0) {
   #    cat('* Saving identified known and similar variants...')
   #    write.table(var_match_similar.dt, paste0(gsub('exon_table_fixed.txt', '', 
   #                                                  input_filename),
   #       'known_similar_ad_variants.txt'), sep='\t', quote=F, row.names=F)
   #    cat(' OK.\n')
   # }
  
}

