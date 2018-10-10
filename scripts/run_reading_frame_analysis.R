#
# -------------------------- Load required libraries ------------------------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(optparse))

#------------------------ Analysis parameters -------------------------------
stop_codons = c('TAG', 'TAA', 'TGA')
codon.dt = data.table(
   aa_3let = c('Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His',
               'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp',
               'Tyr', 'Val', 'STOP'),
   aa_1let = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
               'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'STOP'),
   codons = c('GCT, GCC, GCA, GCG', 'CGT, CGC, CGA, CGG, AGA, AGG',
              'AAT, AAC', 'GAT, GAC', 'TGT, TGC', 'CAA, CAG', 'GAA, GAG',
              'GGT, GGC, GGA, GGG', 'CAT, CAC', 'ATT, ATC, ATA',
              'TTA, TTG, CTT, CTC, CTA, CTG', 'AAA, AAG', 'ATG', 'TTT, TTC',
              'CCT, CCC, CCA, CCG', 'TCT, TCC, TCA, TCG, AGT, AGC',
              'ACT, ACC, ACA, ACG', 'TGG', 'TAT, TAC', 'GTT, GTC, GTA, GTG',
              'TAA, TGA, TAG')
)

# Output contig types:
#  1 = in frame
#  0 = in frame, but different codon
#  2 = stop codon
#  3 = untranslated region after stop codon

# ------------- Command line arguments --------------------------------------

# Parse user arguments
option_list = list(
    make_option(c('-i', '--input-exon-table'),
       help = 'Input blast results file with fixed low qual homopolymer runs.'),
    make_option(c('-n', '--input-exon-orfs'),
       help = 'Input file with ORF information for each exon.'),

    make_option(c('-o', '--output-contigs'),
       help = 'Output contig file.')
)
opt = parse_args(OptionParser(option_list=option_list))

# Load the alignment table
input_filename = opt$`input-exon-table`
input_exon_orfs = opt$`input-exon-orfs`

output_contigfile = opt$`output-contigs`

# --------------------------- Load files -------------------------------------
cat('* Loading fixed BLAST results...')
blast.dt = fread(input_filename)
cat(' OK.\n')


# --------------------------- Reading frame analysis ---------------------------

seq_to_codons = function (x, off=0) {
   # Convert sequence x to codons (3-mers) using offset 'off', e.g.
   # x = ATGCTGCCCGG
   # A, TGC, TGC, CCG, G
   if (off == 0) {
      codon_pos = seq(1, nchar(x), 3)
   } else {
      codon_pos = c(1, seq(1+off, nchar(x), 3))
   }
   codons = substring(x, codon_pos[1:(length(codon_pos)-1)],
                      codon_pos[2:length(codon_pos)]-1)
   codon_pos_end = c(codon_pos[length(codon_pos)],
                     nchar(x))[(codon_pos[length(codon_pos)] < nchar(x))+1]
   return(c(
      codons, 
      substring(x, codon_pos[length(codon_pos)], codon_pos_end)
   ))
}

cat('* Loading exon ORF information...')
exon_orf.dt = fread(input_exon_orfs)
# exon_orf.dt[, seq := as.character(exon.fa)]
exon_first = exon_orf.dt[, min(exon_no)]
exon_last = exon_orf.dt[, max(exon_no)]
cat('OK.\n')

cat('* Filtering out possible PCR chimeras...')
blast.dt = blast.dt[!(exon_comb_short_rc %like% "-2,[12][,'-]")]
blast.dt = blast.dt[!(exon_comb_short_rc %like% paste0(exon_last, ',', 
                                                       exon_last, ','))]
blast.dt = blast.dt[!(exon_comb_short_rc %like% paste0(exon_last, ',',
                                                       exon_first+1))]
blast.dt = blast.dt[!(exon_comb_short_rc %like% paste0("[0-9]+,",
                                                       exon_first, 
                                                       "[,'-]"))]
cat('OK.\n')


# Do this on a subset of blast.dt, for each read separately
# First do only reads for which the first entry is exon_01.
cat('* Calculating positions for reverse-complemented reads...')
blast.dt = merge(blast.dt, 
                 exon_orf.dt[, .(exon_no, orf_position)], by='exon_no')
blast.dt = blast.dt[order(read_id, read_start_pos_align)]
blast.dt[, read_length := .SD[.N, read_end_pos_align] - 
                          .SD[1, read_start_pos_align] + 1L, 
         by=read_id]

# Now recalculate the read coordinates for the reverse complemented reads
# rc_reads = blast.dt[, if (.SD[1, exon_id] == 'exon_18_rc') .SD[1,],
#                     by=read_id][, read_id]
# blast.dt[!(read_id %in% rc_reads), exon_seq_corrected := exon_seq]
# blast.dt[read_id %in% rc_reads,
#          c('read_start_pos_align', 'read_end_pos_align', 
#            'read_seq_corrected',
#            'exon_start_pos_align', 'exon_end_pos_align',
#            'exon_seq_corrected') := list(
#             read_length - read_end_pos_align + 1L,
#             read_length - read_start_pos_align + 1L,
#             reverse(chartr('ACGT', 'TGCA', read_seq_corrected)),
#             exon_length - exon_end_pos_align + 1L,
#             exon_length - exon_start_pos_align + 1L,
#             reverse(chartr('ACGT', 'TGCA', exon_seq))
#             ), by=read_id]
blast.dt = blast.dt[order(read_id, read_start_pos_align)]
cat('OK.\n')

# Remove tandem repeats from analysis?
#
# Keep tandem because we will need them for some of the known FADs!
#blast.dt = blast.dt[!(exon_comb_short_rc %like% '18,1') &
#                    !(exon_comb_short_rc %like% '18,18\''), ]

orf_offset = function (i) {
   # This function finds the orf offset for exon 1, without using any
   # modulo 3 division, where i is exon_start_pos_align. This means
   # if i > imax, it will return NA and the read will not be processed.
   # Codon 1
   if (i %in% seq(1,40,3)) {
      return(0L)
   } else if (i %in% seq(2,41,3)) {
      return(2L)
   } else if (i %in% seq(3,42,3)) {
      return(1L)
   } else {
      return(NA)
   }
}


codon_compare = function (v1, v2) {
   # First vector is always the vector of read seq codons.
   # Second vector are exon seq_codons.
   # Return a list based on matching 1 by 1 codon:
   #  v1 = c('A', 'ATG', 'CTG', 'CCC', 'GG')
   #  v2 = c('G', 'ATG', 'CTG', 'TGA', 'GG')
   # out: c(0, 1, 1, 1, 1, 1, 1, 2, 0)
   # 1 = in frame, 2 = stop, 0 = not equal
   n = min(length(v1), length(v2))
   # This gives cummulative lengths for each codon
   v1_ind = as.integer(cumsum(sapply(v1, nchar)))
   # Total number of characters in all codons is the last
   # entry in cummulative vector.
   v1_tot = v1_ind[length(v1_ind)]
   out = rep(0, v1_tot)
   for (i in 1:n) {
      if (i == 1) inds = c(1:v1_ind[1])
      else inds = c((v1_ind[i-1]+1):v1_ind[i])
      if (v1[i] %in% stop_codons) { # STOP codon
         out[inds] = 2L
         if (inds[length(inds)] < v1_tot) {
            out[(inds[length(inds)]+1):length(out)] = 3L
         }
         return(out)
      } else { # Non-stop
         if (v1[i] == v2[i]) out[inds] = 1L
         else out[inds] = 0L
      }
   }
   return(out)
}

vector_contigs = function (v) {
   # For a vector like c(0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 0, 0)
   # find all contiguous runs of numbers. Return the result as a list
   # list(list(1, 1, 0), list(2, 7, 1), list(8, 10, 2), list(11, 12, 0))
   if (length(v) == 1) return(list(c(1L, 1L, as.integer(v[1]))))
   prev = v[0]
   out = list(c(1L, 1L, as.integer(v[1])))
   j = 1L # Contig index
   for (i in 2:length(v)) {
      if (v[i] != out[[j]][3]) {
         out[[j+1]] = c(i, i, as.integer(v[i]))
         out[[j]][2] = i-1L
         j = j + 1L
      }
   }
   if (out[[length(out)]][2] != length(v)) out[[length(out)]][2] = length(v)
   return(out)
}


read_offset = function (x) {
   if (x == 1) return(2L)
   else if (x == 0) return(0L)
   else if (x == 2) return(1L)
   else return(NA)
}

# Reading frame check
#
#              exon 1           exon 2
# exon:    ATG CTG CCC GGT   T TTC CGG AGT
# read: XX XTG CTG CCC GGA AGT TTC CGG AGT GTA AAG ATT
#           ^----------------^- read_start_pos_align
#
# Offset for read to find codons is obtained from exon 1.
#

# blast.dt[, read_orf_offset := read_start_pos_align, by=read_id]

cat('* Generating table of reading frame contiguous segments...')
# test_reads = blast.dt[, unique(read_id)][1:10]
blast_no_reads = as.numeric(gsub('read_([0-9]+)', '\\1',
                                 blast.dt[order(-read_id)][1, read_id]))

blast.dt = blast.dt[, if (all(read_end_pos_align > read_start_pos_align)) 
   .SD, by=read_id]

contigs.dt = blast.dt[, {
   # Establish the first position of the ORF for the read based on exon 1
   #
   #                v--- exon_start_pos_align = 2, exon_orf_position == 1
   # exon:         ATG CTG CCC GGT
   # read: XXX XXX XTG CTG CCC GGA ATG
   #                TG CTG CCC GGA ATG
   #                ^--- read_start_pos_align = 8
   #                offset from rspa = 2
   #                (espa-1) %% 3 = 1
   #                (eop-1) %% 3 = 0
   #                ((espa-eop) %% 3) = 1
   #
   #                 v--- exon_start_pos_align = 2, exon_orf_position = 2
   # exon:         A TGC TGC CCG GT
   # read: X XXX XXX TGC TGC CCG GAA TG
   #                 ^--- read_start_pos_align = 8
   #                 offset from rspa = 0
   #                 (espa-1)%%3 = 1,  (eop-1)%%3 = 1
   #                 (espa-eop)%%3 = 0
   #
   #                v--- exon_start_pos_align = 2, exon_orf_position = 3
   # exon:         AT GCT GCC CGG T
   # read: XX XXX XXT GCT GCC CGG AAT G
   #                ^--- read_start_pos_align = 8
   #                offset from rspa = 1
   #                (espa-1)%%3 = 1, (eop-1)%%3 = 2
   #                (espa-eop)%%3 = 2

   # Loop over each exon
   # If exon 1 is not aligned within the first 3 nts, just skip that read
   # altogether. This means we cannot use exon_id in the by part of data 
   # frame
   # clause.
   dt = copy(.SD)
   # dt = blast.dt[read_id=='read_00183']
   read_orf_start_off = orf_offset(dt[1, exon_start_pos_align])
   # If exon 01 is aligned too far from start codon, just skip the
   # entire read. If not, run reading frame analysis.

   if (dt[1, exon_no] == 1 & !is.na(read_orf_start_off)) {
      # For intra-exonic joins where exons overlap in alignment, adjust
      # read_start_pos_align for this purpose.
      for (i in 2:dt[, .N]) {
         read_pos_align_diff = dt[i, read_start_pos_align] -
            dt[i-1, read_end_pos_align]
         if (read_pos_align_diff <= 0) {
            dt[i, read_start_pos_align := read_start_pos_align -
                  read_pos_align_diff + 1L]
            dt[i, exon_start_pos_align := exon_start_pos_align -
                  read_pos_align_diff + 1L]
            dt[i, read_seq_corrected := str_sub(read_seq_corrected, 1L -
                                                   read_pos_align_diff + 1L)]
            dt[i, exon_seq_corrected := str_sub(exon_seq_corrected, 1L -
                                                   read_pos_align_diff + 1L)]
         }
      }

      # Remove any indels from the alignment
      read_seqs = lapply(dt[, read_seq_corrected], 
                         function (r) gsub('-', '', r))
      exon_seqs = lapply(dt[, exon_seq_corrected], 
                         function (r) gsub('-', '', r))

      # Calculate codon shift for exon 1.
      codon_shift = dt[1, as.integer((exon_start_pos_align - 
                                         orf_position) %% 3)]
      read_offset_1 = read_offset(codon_shift)
      dt[, read_orf_offset := sapply(
         as.integer((read_start_pos_align - read_start_pos_align[1] +
                        codon_shift) %% 3), read_offset)]

      # Calculate reading frame offset for each aligned exon
      dt[, exon_orf_offset := sapply(
         as.integer((exon_start_pos_align - orf_position) %% 3), read_offset)]

      # Generate codons in the correct reading frame
      # We use offset from exon 1 to set the reading frame for the entire read.
      reads_codons = mapply(function (x_s, x_off) seq_to_codons(x_s, x_off),
                            read_seqs, dt[, read_orf_offset])
      exons_codons = mapply(function (x_s, x_off) seq_to_codons(x_s, x_off),
                            exon_seqs, dt[, exon_orf_offset])

      codon_cmp = mapply(codon_compare, reads_codons, exons_codons)
      # names(codon_cmp) = dt[, exon_no]
      contigs_list = lapply(codon_cmp, function (l)
         as.data.table(transpose(vector_contigs(l))))
      for (i in 1:length(contigs_list)) {
         contigs_list[[i]][, V4 := dt[i, exon_no]]
         contigs_list[[i]][, V5 := dt[i, read_start_pos_align]]
         contigs_list[[i]][, V6 := dt[i, exon_start_pos_align]]
         contigs_list[[i]][, V7 := dt[i, exon_end_pos_align]]
         # contigs_list[[i]][, V5 := dt[i, read_end_pos_align]]
      }
      contigs = rbindlist(contigs_list)
      names(contigs) = c('start', 'end', 'type', 'exon_no',
                         'read_start_pos_align',
                         'exon_start_pos_align',
                         'exon_end_pos_align'
                         )
      # codon_cmp = codon_compare(read, exon)
      # contigs = vector_contigs(codon_cmp)
      contigs
   }
   #  2 = stop codon in the read (this needs to be tested first)
   # If it is not stop then:
   #  0 = they're not the same but could be translated
   #  1 = they're exactly the same
   # Anything before ATG on exon 1 should be like '2' ie black
}, by=.(read_id, exon_comb_short_rc)]


# Now fix the intervals after the stop codon for each read
contigs.dt[, contig_id := 1:.N]
contig_id_update = contigs.dt[, {
   first_stop = .SD[, which(type[1:(.N-1)] == 2)[1]]
   if (!is.na(first_stop)) .SD[(first_stop+1):.N, .(contig_id = contig_id)]
}, by=.(read_id)][, contig_id]
contigs.dt[contig_id %in% contig_id_update, type := 3]


contigs.dt = merge(contigs.dt, 
                   unique(blast.dt[, .(read_id, read_len=read_length)]),
                   by='read_id')

contigs.dt[, start_abs := start + read_start_pos_align - 1L]
contigs.dt[, end_abs := end + read_start_pos_align - 1L]

# Find position of stop codon within the read and record it. Later
# we will use this position to sort reads for readability.
contigs.dt[, stop_pos := read_len]
contigs.dt[, stop_pos := {
   .SD[type==2, end_abs-3L] - .SD[1, start_abs+1L]
}, by=read_id]

cat('OK.\n')

cat('* Savint output table...')
write.table(contigs.dt, output_contigfile, sep='\t', quote=F,
            row.names=F)
cat(' OK.\n')
cat('  ->', output_contigfile, '\n')





