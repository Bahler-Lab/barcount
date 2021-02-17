# barcount
Python script for counting DNA barcodes in next generation sequencing data of pooled library experiments

### Is barcount the right tool for my analysis?

Barcount was developed with three particular goals: 
1. It is a simple all-in-one solution that produces counts directly from fastq files.
2. It provides many options and collects a lot of details during the analysis. In our experiments, we were often dealing with mis-matches in the flanking regions or the barcode. Barcount can be configured to be more or less strict with respect to what is counted will keep track of the number of mis-matches.
3. Barcount implements checking for PCR duplicates based on two unique molecular identifier (UMI) regions which are identified by their flanking regions. Reads with the same barcode and the same UMIs are recorded but not counted towards the final result.

These features also come with drawbacks, so you might want to consider if barcount is the right solution to your analysis:
1. The workflow is relatively complex, e.g. compared with an alternative approach of blasting your fastq file against a database of barcodes. Barcount plays to its strengths only when UMIs are used and the UMIs and barcodes must be identifiable by invariable flanking regions with known sequences. 
2. Barcount is relatively slow, certainly compared to many other NGS data analysis tools, and it can take a few hours to analyse a medium-sized fastq file. It is also not parallelised, so only makes use of a single CPU core. However, the memory footprint is low and multiple fastq files can be analysed in parallel on a standard computer or server using the bash script provided. 


### Getting started
1. Barcount requires Python 3 and a few common packages (including [biopython](https://biopython.org/)), available through the [anaconda distribution](https://www.anaconda.com/distribution/).
2. Download the python script from this repository
3. Open a new terminal and run 'python barcount.py -h' which will show the in-built help page.
4. Make sure your barcode database is in the correct format (two-column csv file, barcode and gene id, no header).

### Barcount manual

```
usage: barcount [-h] --barcode_table BARCODE_TABLE
                [--min_read_length MIN_READ_LENGTH]
                [--max_read_length MAX_READ_LENGTH] [--rmdup] --flanking_left
                FLANKING_LEFT --flanking_right FLANKING_RIGHT
                [--max_distance_flanks {0,1}]
                [--max_distance_barcode MAX_DISTANCE_BARCODE] --out OUT
                [--verbose] [--debug] [--seqfile_format SEQFILE_FORMAT]
                [--save_extracted_barcodes] --fastq FASTQ
                [--umiA_position UMIA_POSITION]
                [--umiB_position UMIB_POSITION]

Welcome to barcount v0.9. Count barcodes in fastq files and match to genes.
Written by stephan.kamrad@crick.ac.uk and maintained at
https://github.com/Bahler-Lab/barcount.

optional arguments:
  -h, --help            show this help message and exit
  --barcode_table BARCODE_TABLE
                        Path to barcode database. This is a file containing a
                        gene id and its barcode, separated by comma without
                        additional whitespace on each line. There must be no
                        header.
  --min_read_length MIN_READ_LENGTH
                        Minimum read length.
  --max_read_length MAX_READ_LENGTH
                        Maximum read length.
  --rmdup               If set, consider only the first occurence of the
                        entire read in the fastq file.
  --flanking_left FLANKING_LEFT
                        Flanking sequence of the barcode region on the left.
  --flanking_right FLANKING_RIGHT
                        Flanking sequence of the barcode region on the right.
  --max_distance_flanks {0,1}
                        Number of single-letter permutations (Levenshtein
                        distance) to allow between barcodes. The default is 0.
                        The maximum is 1.
  --max_distance_barcode MAX_DISTANCE_BARCODE
                        Number of single-letter permutations (Levenshtein
                        distance) to allow between barcodes. The default is 1.
  --out OUT             Output files will be saved under this path/name. This
                        includes the count table, qc plots and additional
                        stats (required).
  --verbose             If set, use verbose output.
  --debug               If set, make debug file True/False.
  --seqfile_format SEQFILE_FORMAT
                        Format of the sequence file (fastq by default). This
                        will be passed on to Biopython SeqIO.
  --save_extracted_barcodes
                        If set, save an additional fastq file with the
                        extracted barcodes .
  --fastq FASTQ         Path to the fastq file to analyse.
  --umiA_position UMIA_POSITION
                        Find the first UMI by position in the read. The
                        position must be in the format start:stop following
                        python-like indexing, which is 0-based and does not
                        include the stop posiion in the slice. Prefilter read
                        file for correct readlength or use with
                        --max_read_length or --min_read_length.
  --umiB_position UMIB_POSITION
                        Find the second UMI by position in the read. The
                        position must be in the format start:stop following
                        python-like indexing, which is 0-based and does not
                        include the stop posiion in the slice. Prefilter read
                        file for correct readlength or use with
                        --max_read_length or --min_read_length.
						
```

### Barcount algorithm

![pyphe logo](https://github.com/Bahler-Lab/barcount/blob/master/documentation/algorithm.png)

Explanation of main loop: The next read in the sequence file is loaded using BioPython SeqIO module. The user has the ability to specify the read file format which will be passed onto the SeqIO’s.parse() function, the default is fastq. Next, barcount checks if the read-length is within the upper and lower limits set by the user (do not necessarily have to be used simultaneously). Next, the barcode is found by finding the constant regions surrounding the barcode (flanks). Users should take care to provide flanks which are long enough to not appear at random. The algorithm first looks for exact matches of the flank in the sequence, then allowing one mismatch, then one insertion and one deletion. If both flanks are found, the sequence between is extracted as the barcode and compared to the database. We check first if the database contains the exact barcode, in which case it is immediately assigned to the gene. If not, the Levenshtein distance of the extracted barcode to all database barcodes is computed. The best match is chosen and the corresponding gene is assigned if the matching distance is within the user-defined cut-off (again, care should be taken to a threshold appropriate for the length of the barcode and considering the distribution of distances within the database). If there is a tie between two database entries, no barcode is assigned. If the rmdup filter is set, a hash of the sequence is computed using the Bio.SeqUtils.CheckSum.seguid() function which is checked against a cache of previously seen hashes. Note that this can increase the runtime of barcount significantly so it might be worth considering alternative tools for this step if needed (example?). Finally, the UMIs are identified by position (please make sure to specify those as described in the readme, using Pythonic indexing). If two UMIs are present (--umiA_position and –umiB_position set), they are concatenated and checked against a cache of previously seen UMIs. The cache is specific to each database entry, which means the maximum signal per gene is at 4^length(UMI) and significant saturation will be seen below that. Please consider this danger when using short UMIs and or very large read files. 

