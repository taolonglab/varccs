library(stringr)


comb_short_to_long = function (s) {
   # Take 1-3,17-18
   # and output 1,2,3,17-18
   # Also be able to parse 1-3,17-18,18'-17',3'-1',1-3,17-18
   # to: 1,2,3,17,18,18,17,3,2,1,1,2,3,17,18
   runs = strsplit(as.character(s), ',')[[1]]
   unlist(lapply(runs, function (r) {
      xs = strsplit(gsub("'", "", r), '-', fixed=T)[[1]]
      if (length(xs) == 1) {
         xs[1]
      } else {
         seq(xs[1], xs[2])
      }
   }))
}

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

comb_pairs_from_str = function (s) {
   # Take string s which has comma delimited numbers
   # and return all consecutive pairs of these numbers
   # as matrix, with columns:
   # 1. first exon in pair
   # 2. second exon in pair
   # 3. exon-exon interface type:
   #        0: two consecutive, i.e. 1,2 or 2,1
   #        1: exon-exon join i.e. 3,16
   x_v = as.integer(strsplit(s, ',')[[1]])
   x = x_v[1]
   d = ''
   ss = paste0('', x_v[1])
   pairs = matrix(nrow=length(x_v)-1, ncol=3)
   for (i in 2:length(x_v)) {
      pairs[i-1,1] = as.integer(x_v[i-1])
      pairs[i-1,2] = as.integer(x_v[i])
      if ((x+1) == x_v[i] & d != '-') {
         # Two consecutively increasing numbers
         d = '+'
         if (str_sub(ss, -1) != '-') {
            ss = paste0(ss, '-')
         }
         pairs[i-1,3] = 0
      }
      else if ((x-1) == x_v[i] & d != '+') {
         # Two consecutively decreasing numbers
         d = '-'
         if (str_sub(ss, -1) != '-') {
            ss = paste0(ss, '-')
         }
         pairs[i-1,3] = 0
      }
      else {
         # Two non-consecutive numbers. Close 
         # any previously open interval
         d = ''
         if (str_sub(ss, -1) == '-') {
            ss = paste0(ss, x_v[i-1])
         }
         # Start a new possible interval
         # (the - sign will be added only if
         # the next number is consecutively incr/decr)
         ss = paste0(ss, ',', x_v[i])
         pairs[i-1,3] = 1
      }
      x = x_v[i]
   }
   
   return(pairs)
}

exon_joins_from_short_str = function (s) {
   # From shortened exon combination string s
   # such as 1-3,17-18,1-3,17-18, extract pairs of
   # all exon joins: 3,17  18,1  3,17
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



tandem_combinations = function (s) {
   # Input: strings e.g. 1-3,17-18,18',2'-1',1-2,17-18
   s_joins_0 = strsplit(s, ',')[[1]]
   # Fill out the runs
   s_joins = c()
   for (i in 1:length(s_joins_0)) {
      s_run = strsplit(s_joins_0[i], '-')[[1]]
      if (length(s_run) == 1) {
         s_joins = c(s_joins, s_run)
      } else {
         if (s_run %like% "'") {
            s_runs_i = as.integer(gsub("'", '', s_run))
            s_runs_str = sapply(
               (s_runs_i[1]:s_runs_i[2]), 
               function (j) paste0(as.character(j), "'"))
            s_joins = c(s_joins, s_runs_str)
         } else {
            s_joins = c(s_joins, as.character(c(s_run[1]:s_run[2])))
         }
      }
   }
   # 1-3    17-18    18'   2'-1'   1-2  17-18
   comb = F
   combs = list()
   current_comb = c()
   j = 1
   for (i in 1:length(s_joins)) {
      e_no = as.integer(gsub("'", "", s_joins[i]))
      if (e_no == 1) {
         if (!comb) {
            # Combination starts with 1
            comb = T
            current_comb = c('1')
         } else {
            # Combination finished with 1 (this is rc)
            comb = F
            combs[[j]] = c(current_comb, "1'")
            j = j + 1
            current_comb = c()
         }
      } else if (e_no == 18) {
         if (!comb) {
            comb = T
            # Combination starts with 18 (this is rev comp)
            current_comb = c("18'")
         } else {
            # Combination finished
            combs[[j]] = c(current_comb, '18')
            j = j + 1
            current_comb = c()
            comb = F
         }
      } else {
         current_comb = c(current_comb, s_joins[i])
      }
      
   }
   return(combs)
}

exon_comb_c_to_str = function (v) {
   # Input exon combination as a vector of strings
   # Output string with , and -
   v_prev = as.integer(gsub("'", '', v[1]))
   out = v[1]
   run = F
   for (i in 2:length(v)) {
      v_i = as.integer(gsub("'", '', v[i]))
      if (v_i == (v_prev+1) | v_i == (v_prev-1)) {
         if (!run) {
            out = paste0(out, '-')
         }
         run = T
      } else {
         if (run) {
            out = paste0(out, v_prev, ',', v[i])
         } else {
            out = paste0(out, ',', v[i])
         }
         run = F
      }
      v_prev = v_i
   }
   # Do we have an unfinished run?
   if (run) {
      out = paste0(out, v_prev)
   }
   return(out)
}

