#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Filter homopolymers that start with the base that has poor
quality score and replace them with run_letter (default N).



Created on Wed Apr 19 11:34:52 2017

@author: igor
"""

import sys
import re
import os.path


def main(input_fq, output_st, output_fq, rep=2, run_letter='N', qs_letter='~', qs_threshold=30):
    
    # Define and compile regular expressions for
    # detection of repeated bases.    
    regex = re.compile(r"A{"+str(rep)+",}|C{"+str(rep)+",}|G{"+str(rep)+",}|T{"+str(rep)+",}")
    
    with open(input_fq, 'r') as fq, open(output_fq, 'w') as fo, open(output_st, 'w') as st:
        
        for i, fq_line in enumerate(fq):
            if i % 4 == 0:
                # Meta-data line
                meta = fq_line.strip()
            elif i % 4 == 1:
                # Sequence line
                seq = fq_line.strip()
            elif i % 4 == 3:
                # Quality score line
                s, q = '', ''
                s_i = 0
                no_matches = 0
                for m in regex.finditer(seq):
                    qs = ord(fq_line[m.start()]) - 33
                    if qs < qs_threshold:
                        s += seq[s_i:m.start()+1]+run_letter
                        q += fq_line[s_i:m.start()+1]+qs_letter
                        # s_i = m.end()
                    else:
                        s += seq[s_i:m.end()]
                        q += fq_line[s_i:m.end()]
                    s_i = m.end()
                    no_matches += 1
                    st.write(str(qs)+'\t'+seq[m.start()]+'\n')
                # Add everything as is after the last match
                s += seq[s_i:]
                q += fq_line[s_i:]

                if no_matches == 0:
                    # No matches found
                    fo.write(meta+'\n')
                    fo.write(seq+'\n')
                    fo.write('+\n')
                    fo.write(fq_line)
                else:
                    # Write to output
                    fo.write(meta+'\n')
                    fo.write(s+'\n')
                    fo.write('+\n')
                    fo.write(q)    
                    
                
    


if __name__ == '__main__':
    
    # verbose=False
    # if verbose:
    #     print('FASTQ homopolymer run filter.')    
    # input_fq = 'reads/job-7871_filtered_qs85.fastq'
    # input_fq = 'test.fastq'
    input_fq = sys.argv[1] # Input FASTQ file
    input_qs = sys.argv[2] # Quality score threshold
    input_rl = sys.argv[3] # Minimum number of repeats
    # if verbose:
    #     print('Input file: '+input_fq)
    #     print('Minimum quality score: '+str(input_qs))
    #     print('Minimum homopolymer run length: '+str(input_rl))
    
    # output_fq = os.path.splitext(input_fq)[0]+'_hpf_qs'+str(input_qs)+'_rep'+str(input_rl)+'.fastq'
    # output_st = os.path.splitext(input_fq)[0]+'_hpf_qs'+str(input_qs)+'_rep'+str(input_rl)+'_homopolymer_qs.txt'
    output_fq = sys.argv[4]
    output_st = sys.argv[5]


    main(input_fq, output_st, output_fq, rep=int(input_rl), qs_threshold=int(input_qs))
    # if verbose:
    #     print('Done.')