# Genome Order Optimization using a Genetic Algorithm (GOOGA) (now published: https://doi.org/10.1371/journal.pcbi.1006949)


This project implements a genetic algorithm (GA) designed to take genotype data from a segregating population and use it to identify an optimal scaffold order from a perhaps not so well assembled genome.
It does this by permuting scaffold orders from the genome assembly and fitting a HMM with genotype data to estimate recombination rates and the ultimate likelihood of a particular order.  Then this is fed into the GA to search for a near optimal order.

The code runs in python and was developed for python2.7 and requires the scipy package to be installed.

To run on the test data, first download this repo.  Then you need to extract the test data:

` gunzip data/LVR.cross.tar.gz` 

` tar xf data/LVR.cross.tar` 

This gives several files, including the genotype files (e.g. g.F2.163.txt), which look like:

`head g.F2.163.txt`

`F2.163	1	0	AB`

`F2.163	1	100000	AB`

`F2.163	1	200000	AB`

`F2.163	1	300000	AB`

`F2.163	1	400000	AB`

`F2.163	1	500000	AB`


where the 4 columns are line_name, chromosome, interval, genotype (AA, AB, BB, or NN for missing/unknown)

Also the test data contains a list of all segregating lines you wish to consider for ordering scaffolds, it's called `LVR.f2set.txt`, an estimate of the intitial intra-scaffold recombinantion rates called `LVR.isr.txt`, and a file of error rates called `LVR.er2.txt`.

Finally, the genetic algorithm code is called `parallel.genet.alg.final.py`. It was written to be run on a multiprocessor system and can make use of parallelism.  It takes several flags at runtime.  To get help on these flags simply run:

`python2 parallel.genet.alg.final.py --help`

Usage: parallel.genet.alg.final.py marker_file [options]


Options:

  -h, --help            show this help message and exit

  -e ERROR_FILE         File with precaculated error rates -
                        default=error.rates.txt
                        
  -f F2_FILE            File with F2s to use in analysis -
                        default=test.f2group.txt
                        
  -i INTRASCAFF_RATES_FILE
                        File with precalculated intra-scaffold recombination
                        rates - default=intrascaff_rates.txt
                        
  -c NCPU               Number of CPUs - default=16
  
  -g NGEN               Number of generations to run the GA - default=10^10
  
  -l ELITE              Number of elite in each generation default=3
  
  -t TERMINATION        Number of generations with no improvement before
                        termination - default=1000 
  
  -r                    Use -r flag if your mapping pop is RILs - default is
                        F2 seg. pop


NCPU should be set according to your computer.  The code was developed for a 16 core machine, hence the default setting.    

ELITE designates the number of contigs orders (in a GA they are called individuals) to be carried over to the next generation.  I've had good luck setting this between 2-4.

The GA will automatically select it's population size based on NCPU and ELITE. Specifically, it will be NCPU + ELITE.  So for example, given the default settings, the GA population size will be 19.

You can also control how long the GA runs, using both NGEN and TERMINATION.  NGEN simply sets a max number of generations for the GA, whereas TERMINATION allows you to tell the GA to stop after a specified number of generations where the best scaffold order (individual) has not been improved upon. For example, if NGEN is set to 500 and TERMINATION is set to 200, the GA will stop at 500 generations total or 200 gen without improvement, whichever comes first.  As a default the NGEN is set astronomically high, to rely only on TERMINATION to end the run.

Finally the -r flag is a boolean, which when used switches the HMM backend to accomadate a RIL population (no hets) rather than a standard seg pop, like a F2 (has hets).  If you don't use this flag, the GA defaults to standard F2 seg. pop.

The only mandatory input file is a starting a genomic scaffold order ("marker_file").  There's an example in markers.2.v2.txt.  This file correpsonds to Chromosome 2 in Mimulus guttatus.  The GA uses this order as a starting place.

Then run the code:

`python2 parallel.genet.alg.final.py markers.2.v2.txt -c 8 -l 1 -i LVR.isr.txt -f LVR.f2set.txt -e LVR.er2.txt > some.output.file`

This runs on 8 CPUs with 1 elite on the test files noted above. And at each generation it will output various run statistics, including all current scaffold orders and likelihoods that the algorithm is grinding away on.

At some point it will stop, and you might want to run:

`python2 run.stats.py some.output`

This step requires that you have the matplotlib package installed.  It gives a plot of various run statistics so you can see how it progressed.  It gives some sense of whether you need to tune things or if it was pretty successful.

You'll also want to check out the last generation of the output file.  This will contain the full likelihood of the best individual and the corresponding contig order.
