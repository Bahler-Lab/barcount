import sys
import argparse
from warnings import warn
import pandas as ps
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid

def levenshtein_distance(s1, s2):
    """
    Python version of Levenshtein distance for compatability. The third-party
    library is much faster and recommended. Taken from recipe at:
    http://code.activestate.com/recipes/576874-levenshtein-distance/
    """
    l1 = len(s1)
    l2 = len(s2)

    matrix = [range(l1 + 1)] * (l2 + 1)
    for zz in xrange(l2 + 1):
        matrix[zz] = range(zz, zz + l1 + 1)
    for zz in range(0, l2):
        for sz in range(0, l1):
            if s1[sz] == s2[zz]:
                matrix[zz + 1][sz + 1] = min(matrix[zz + 1][sz] + 1,
                                    matrix[zz][sz + 1] + 1, matrix[zz][sz])
            else:
                matrix[zz + 1][sz + 1] = min(matrix[zz + 1][sz] + 1,
                                    matrix[zz][sz + 1] + 1, matrix[zz][sz] + 1)
    return matrix[l2][l1]

try:
    import Levenshtein
    levenshtein_distance = Levenshtein.distance
except ImportError:
    warn("python-Levenshtein package not found. The package is " +
                  "recommended as it makes barcode splitting much faster. " +
                  "Using native Python Levenshtein distance function instead.")

def hamming2(s1, s2):
    """Calculate the Hamming distance between two bit strings"""
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

    
def run_analysis(barcode_table, out, fastq, flanking_left, flanking_right, min_read_length=None, max_read_length=None, rmdup=False, max_distance_flanks=0, max_distance_barcode=1, verbose=False, debug=False, seqfile_format='fastq', save_extracted_barcodes=False, umiA_start=None, umiA_end=None, umiB_start=None, umiB_end=None):

    #Read the table with gene-barcode matchings
    barcode_data = ps.read_csv(barcode_table, header=None, names=['gene',], index_col=1)
    barcode_data.loc['unassigned'] = np.nan

    if len(barcode_data.index) != len(list(set(barcode_data.index))):
        raise IOError('Duplicated entry in barcode table or multiple barcodes map to the same gene.')
        
    if len(barcode_data.columns) != 1:
        raise IOError('Error reading barcode table. Make sure the file is comma-separated and has exactly two columns: gene id, barcode')
        
    if verbose : 'Read barcode table containing %i barcode-gene pairs'%len(barcode_data.index)
    barcodes = barcode_data.index.tolist()
    
    #Setup output stats
    stats = ps.Series()
    stats['fastq_entries'] = 0
    stats['failed_findBarcode'] = 0
    stats['unmatched'] = 0
    stats['sucessfully_matched'] = 0

    force_read_length = True if (min_read_length or max_read_length) else False
    if force_read_length:
        stats['failed_minReadLength'] = 0
        stats['failed_maxReadLength'] = 0

    barcode_data['count_postFilters'] = 0
    if rmdup:
        barcode_data['failed_rmdup'] = 0
    
    if umiA_start:
        barcode_data['failed_sameUMI'] = 0    
        
    #Set up rmdup cache
    rmdup_cache = []
    umi_cache = {bc : [] for bc in barcode_data.index}
    
    #Matching quality
    matchingqs = {bc : [] for bc in barcode_data.index}
    barcodelengths = {bc : [] for bc in barcode_data.index}

    #Open files
    read_inf = open(fastq, 'r')
    unmachted_fastq = open(out+'_unmatched.fastq', 'w')
    no_barcode_fastq = open(out+'_noBarcode.fastq', 'w')
    if save_extracted_barcodes:
            barcodes_fastq = open(out+'_barcodes.fastq', 'w')
            
    if debug:
        debug_file = open(out+'_debug.csv', 'w')
        debug_log = ['fastq_entry',#0
                'fastq_header',#1
                'sequence',#2
                'sequence_length',#3
                'forceReadLengthFilter',#4
                'left_matchingq',#5
                'left_start',#6
                'left_end',#7
                'right_matchingq',#8
                'right_start',#9
                'right_end',#10
                'barcode',#11
                'found_gene_match',#12
                'matched_gene',#13
                'matching_quality',#14
                'checksum',#15
                'rmdupFilter',#16
                'umiA',#17
                'umiB',#18
                'sameUMIFilter',#19
                ]
        
    ###MAIN LOOP###
    seq_count = 0
    for r in SeqIO.parse(read_inf, seqfile_format):
        #write the debug line from last iteration and reset
        if debug: debug_file.write(','.join(map(str, debug_log))+'\n')
        debug_log = ['NA']*20
        
        seq_count += 1
        if verbose:
            if seq_count%5000 == 0:
                print('Processed %i reads'%seq_count)
        debug_log[0] = str(seq_count)
        debug_log[1] = str(r.description)
        debug_log[2] = str(r.seq)
        debug_log[3] = len(str(r.seq))
        
        #Apply readlength filter, move on to next sequence if it fails
        if force_read_length:
            if (min_read_length and (len(r.seq) < min_read_length)):
                stats['failed_minReadLength'] += 1
                debug_log[4] = 'FAIL'
                continue
            elif max_read_length and (len(r.seq) > max_read_length):
                stats['failed_maxReadLength'] += 1
                debug_log[4] = 'FAIL'
                continue
            else:
                debug_log[4] = 'PASS'                


        #Try to find the barcode
        barcode, debug_log[5:11] = find_region(str(r.seq), flanking_left, flanking_right, max_distance=max_distance_flanks)
        debug_log[11] = barcode
        if barcode == 'FAIL':
            stats['failed_findBarcode'] += 1
            SeqIO.write([r,], no_barcode_fastq, 'fastq')
            continue
        if save_extracted_barcodes:
            SeqIO.write([r[debug_log[7]:debug_log[9]],], barcodes_fastq, 'fastq')
            
        #Match barcode to gene
        matched_barcode, matchq = match_barcode(barcode, barcodes, max_distance_barcode)
        if matched_barcode == 'unassigned':
            SeqIO.write([r,], unmachted_fastq, 'fastq')
            stats['unmatched'] += 1
            debug_log[12] = 'FAIL'
            continue
        else:
            gene = barcode_data.loc[matched_barcode,'gene']
            matchingqs[matched_barcode].append(matchq)
            debug_log[12] = 'PASS'
            debug_log[13] = gene
            debug_log[14] = str(matchq)
            stats['sucessfully_matched'] += 1
            barcodelengths[matched_barcode].append(len(barcode))
      
        #If rmdup filter is active, check if we have seen this sequence before
        if rmdup:
            checksum = seguid(r.seq)
            debug_log[15] = checksum
            if checksum in rmdup_cache:
                debug_log[16] = 'FAIL'
                barcode_data.loc[matched_barcode, 'failed_rmdup'] += 1
                continue
            else:
                debug_log[16] = 'PASS'
                rmdup_cache.append(checksum)
                
        #If umi filter is active, check if we have seen this UMI before
        if umiA_start:
            umiA = str(r.seq)[umiA_start:umiA_end]
            debug_log[17] = umiA
            if umiB_start:
                umiB = str(r.seq)[umiB_start:umiB_end]
            else:
                umiB = ''
            debug_log[18] = umiB
            umiA += umiB
            if umiA in umi_cache[matched_barcode]:
                debug_log[19] = 'FAIL'
                barcode_data.loc[matched_barcode, 'failed_sameUMI'] +=1
                continue
            else:
                debug_log[19] = 'PASS'
                umi_cache[matched_barcode].append(umiA)   
        
        #If loop has not been continued yet, count this as a true hit
        barcode_data.loc[matched_barcode, 'count_postFilters'] += 1
            
        
    #Write debug file for final iteration
    if debug: debug_file.write(','.join(map(str, debug_log))+'\n')

    #Close all files
    read_inf.close()
    unmachted_fastq.close()
    no_barcode_fastq.close()
    if save_extracted_barcodes:
        barcodes_fastq.close()
        
    if debug: debug_file.close()

    if verbose: print('Processed fastq file sucessfully. Generating output files.')

    #Export data
    barcode_data['mean_matching_distance'] = [np.mean(matchingqs[k]) if len(matchingqs[k])>0 else np.nan for k in barcode_data.index]
    barcode_data['mean_barcode_lengths'] = [np.mean(barcodelengths[k]) if len(barcodelengths[k])>0 else np.nan for k in barcode_data.index]

    barcode_data = barcode_data.sort_values('count_postFilters', ascending=False)
    barcode_data.to_csv(out+'.csv')

    stats['fastq_entries'] = seq_count
    stats.to_csv(out + '_stats.csv')

    if verbose: print('Done.')
    return barcode_data, stats
    
    
def find_region(rs, fleft, fright, max_distance=0):
    log = ['NA']*6 #left_matchingq, left_start, left_end, right_matchingq, right_start, right_end       
    log[:3] = search_substring(rs, fleft, max_distance)
    log[3:] = search_substring(rs, fright, max_distance)
    
    try:
        return rs[log[2]:log[4]], log
    except Exception:
        return 'FAIL', log
        
        
def search_substring(rss, subrss, max_distance=0):
    '''This function finds a substring in a string using a custom, greedy search algorith which returns if a match was found without necessarily exploring the entire search space. Barcode reads will in >80% of the cases contain exact matches, so these are looked for first. If the distance is set to 1, substrings with one single base substitution are looked for. Indels are considered last. This means the behaviour for multiple matches is slightly inconsistent: If there are multiple exact matches, the first one is returned. If there are multiple matches with one SNP, the last one is returned. If there are multiple matches with one indel, the first one is returned.'''

    assert max_distance in [0,1]   

    lowest = [100, None, None]#initialise: matchingq, start, end
    
    #Start with substrings of the same length. If there is an exact match, return immediately
    for startpos in range(len(rss)-len(subrss)-1):
        temp_substring = rss[startpos:startpos+len(subrss)]
        temp_distance = hamming2(temp_substring, subrss)#Hamming distance is much quicker than levenshtein, and in 90% of cases, this is all that's needed
        
        if temp_distance < lowest[0]:
            lowest = [temp_distance, startpos, startpos+len(subrss)]
            
        if lowest[0] == 0:
            return lowest

    #If there was a hit with 1 mismatch, return that one. I.e. this algorithm favours SNPs over indels.
    if max_distance == 1:
        if lowest[0] == 1:
            return lowest
    
    #If no match was found yet, look at substrings with 1 length difference. Return immediately when there was a hit
    if max_distance == 1:
        #Start with substrings with one insertion
        for startpos in range(len(rss)-len(subrss)-2):
            temp_substring = rss[startpos:startpos+len(subrss)+1]
            temp_distance = levenshtein_distance(temp_substring, subrss)
            
            if temp_distance < lowest[0]:
                lowest = [temp_distance, startpos, startpos+len(subrss)+1]#Here it is strictly not necessary to keep track of the lowest (because it returns immediately when there is a hit, but it also doesn't do any harm.
               
            if lowest[0] == 1:
                return lowest   
                
        #Next, substrings with one deletion
        for startpos in range(len(rss)-len(subrss)):
            temp_substring = rss[startpos:startpos+len(subrss)-1]
            temp_distance = levenshtein_distance(temp_substring, subrss)
            
            if temp_distance < lowest[0]:
                lowest = [temp_distance, startpos, startpos+len(subrss)-1]
        
            if lowest[0] == 1:
                return lowest

                
    #If no match was found, return 'FAIL'
    return ['FAIL', 'NA', 'NA']
        
def match_barcode(barcode, barcode_list, max_distance):
    """Return the best matching barcode and the matching quality"""
    if barcode in barcode_list:
        return barcode, 0
    
    else:
        # find the Levenshtein distance to each
        distances = [levenshtein_distance(barcode, b) for b in barcode_list]
        best = min(distances)
        # check if there's a tie or the distance is too great:
        if best > max_distance or distances.count(best) > 1:
            return 'unassigned', best
        else:
            # otherwise, return the best one, after caching it for future use
            return barcode_list[distances.index(best)], best
            
            
def print_dict_long(indict):
    outstr = ''
    for k,v in indict.items():
        outstr += str(k) + ': ' + str(v) + '\n'
    return outstr
    
    
def cli():
    #Set up parsing of command line arguments with argparse
    parser = argparse.ArgumentParser(description='Welcome to barcount v0.9. Count barcodes in fastq files and match to genes. Written by stephan.kamrad@crick.ac.uk and maintained at https://github.com/Bahler-Lab/barcount.')
    
    parser.add_argument('--barcode_table', type=str, required=True, help='Path to barcode database. This is a file containing a gene id and its barcode, separated by comma without additional whitespace on each line. There must be no header.')
    
    parser.add_argument('--min_read_length', type=int, help='Minimum read length.')
    parser.add_argument('--max_read_length', type=int, help='Maximum read length.')

    parser.add_argument('--rmdup', default=False, action='store_true', help='If set, consider only the first occurence of the entire read in the fastq file.')
    
    parser.add_argument('--flanking_left', type=str, required=True, help='Flanking sequence of the barcode region on the left.')
    parser.add_argument('--flanking_right', type=str, required=True, help='Flanking sequence of the barcode region on the right.')
    
    parser.add_argument('--max_distance_flanks', type=int, choices=[0,1], help='Number of single-letter permutations (Levenshtein distance) to allow between barcodes. The default is 0. The maximum is 1.', default=0)

    parser.add_argument('--max_distance_barcode', type=int, help='Number of single-letter permutations (Levenshtein distance) to allow between barcodes. The default is 1.', default=1)
    parser.add_argument('--out', type=str, required=True, help='Output files will be saved under this path/name. This includes the count table, qc plots and additional stats (required).')
    
    parser.add_argument('--verbose', default=False, action='store_true', help='If set, use verbose output.')
    parser.add_argument('--debug', default=False, action='store_true', help='If set, make debug file True/False.')

    parser.add_argument('--seqfile_format', default='fastq', help='Format of the sequence file (fastq by default). This will be passed on to Biopython SeqIO.')    
    
    parser.add_argument('--save_extracted_barcodes', default=False, action='store_true', help='If set, save an additional fastq file with the extracted barcodes .')

    parser.add_argument('--fastq', type=str, required=True, help='Path to the fastq file to analyse.')
    
    #UMI options  
    parser.add_argument('--umiA_position', type=str, help='Find the first UMI by position in the read. The position must be in the format start:stop following python-like indexing, which is 0-based and does not include the stop posiion in the slice. Prefilter read file for correct readlength or use with --max_read_length or --min_read_length.')
    parser.add_argument('--umiB_position', type=str, help='Find the second UMI by position in the read. The position must be in the format start:stop following python-like indexing, which is 0-based and does not include the stop posiion in the slice.  Prefilter read file for correct readlength or use with --max_read_length or --min_read_length.')
    
    args = parser.parse_args()
    
    if args.verbose : print('Welcome to barcount v0.92 (20.6.2020). Count barcodes in fastq files and match to genes. Written by stephan.kamrad@crick.ac.uk and maintained at https://github.com/Bahler-Lab/barcount.')

    if args.verbose: print('Command line call: '+' '.join(sys.argv))
        
    if args.umiA_position:
        try:
            umiA_start = int(args.umiA_position.split(':')[0])
            umiA_end = int(args.umiA_position.split(':')[1])
        except Exception:
            raise ValueError('umiA position is not specified correctly. The position must be in the format start:stop following python-like indexing, which is 0-based and does not include the stop posiion in the slice. Argument received: %s'%args.umiA_position)
    else:
        umiA_start = None
        umiA_end = None
            
    if args.umiB_position:
        try:
            umiB_start= int(args.umiB_position.split(':')[0])
            umiB_end = int(args.umiB_position.split(':')[1])
        except Exception:
            raise ValueError('umiB position is not specified correctly. The position must be in the format start:stop following python-like indexing, which is 0-based and does not include the stop posiion in the slice. Argument received: %s'%args.umiB_position)
    else:
        umiB_start = None
        umiB_end = None

    if umiA_start or umiB_start:
        warn('You have set options to finc UMIs by position. This is only recommended if you first make sure that reads are of the correct length and it is your repsonsibility to the the min_red_length and max_read_length accordingly.')
        
    if args.umiB_position and not args.umiA_position:
        raise ValueError('umiB positions are specified but not umiA positions. If your reads contain only one UMI, please consider this umiA.')
        
    if args.verbose: print('Running barcount with the following arguments:\n' + print_dict_long(vars(args)))
        
    #Run analysis
    _,_ = run_analysis(args.barcode_table, args.out, args.fastq, args.flanking_left, args.flanking_right, min_read_length=args.min_read_length, max_read_length=args.max_read_length, rmdup=args.rmdup, max_distance_flanks=args.max_distance_flanks, max_distance_barcode=args.max_distance_barcode, verbose=args.verbose, debug=args.debug, seqfile_format=args.seqfile_format, save_extracted_barcodes=args.save_extracted_barcodes, umiA_start=umiA_start, umiA_end=umiA_end, umiB_start=umiB_start, umiB_end=umiB_end)

