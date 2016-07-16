This project implements a algorithm to designed to take genotype data from a segregating population and use it to identify an optimal scaffold order from a perhaps not so well assembled genome.
It does this by permuting scaffold orders from the genome assembly and fitting a HMM with genotype data to estimate recombination rates and the ultimate likelihood of a particular order.  Then this is fed into a genetic algorithm that searches for a optimal order.

The code runs in python and was developed for python2.7 and requires the scipy package to be installed.

To run on the test data, first download this repo.  Then you need to extract the test data:

` bunzip2 Lex2.tar.bz2` 

` tar xf Lex2.tar` 

This gives several files, including the genotype files (e.g. `Genotypes.imswc922.txt`), which look like:

`head Genotypes.imswc922.txt`

imswc922	1	0	AB

imswc922	1	100000	AB

imswc922	1	200000	AB

imswc922	1	300000	AB

imswc922	1	400000	AB

imswc922	1	500000	AB

imswc922	1	600000	NN

imswc922	1	700000	AB

imswc922	1	800000	AB

imswc922	1	900000	AB

where the columns are line_name, chromosome, interval, genotype

Also the test data contains a list of all lines you wish to consider for ordering scaffolds, it's called `test.f2group.txt`, and a file of error rates called `error.rates.txt`.

Finally before running you'll need to edit the main genetic algorithm code.  It's called `parallel.genet.alg.optimize.select.contigs.py`. It was written to be run on a multiprocessor system and can make use of parallelism.  The file as several global variables set at the top:

`POP_SIZE = 16  #pop size`

`NGEN = 1000000000 #generations to run`

`MUTATION = [0, 1, 2, 3]  #randomly select one value from this list to determine the number of mutations an indiv. pass on to next gen.`

`ELITE = 3  #number best individuals to save at each generation`

`TERMINATION = 100 #if the most fit line doesn't change for 100 generation, end the run`

`NCPUs = 16 #number of cpus to run on.  Remember that after the 1st generation you will have POP_SIZE - ELITE novel indv`
          
          `#since we save past results, on a machine with 10 CPUs, if ELITE=2, you may want to do a popsize of 12, because that will max out all 10 CPUs after Gen 1 `

NCPUs should be set according to your computer.  I was running on a 16 core machine, hence the settings.  It usually makes sense to set the genetic algorithm population size (POP_SIZE) to be slightly larger than the number of cores.   Specifically it should be set to NCPUs + ELITE.

And ELITE designates the number of contigs orders (in a genet. alg. they are called individuals) to be carried over to the next generation.  I've had good luck setting this between about 2-4.

You can also control how long the genet. alg. runs, using both NGEN and TERMINATION.  NGEN simply sets a max number of generations for the genetic algorithm, whereas TERMINATION allows you to tell the genetic algorithm to stop after a specified number of generations where the best scaffold order (individual) has not be improved upon.
For example, if NGEN is set to 500 and TERMINATION is set to 200, the genetic algorithm will stop at 500 generations total or 200 gen without improvement, whichever comes first.

Finally there is the per indiviudal mutation rate.  It's specified by the list given to MUTATION.  It is currently set to randomly sample between 0 and 3 mutations per indvidiual per generation.  This proved useful for chromosomes with approx. 20 scaffolds.  It may need to be increased or descreased for for different genomes.

The final input file is a starting a genomic scaffold order.  THere's an example in markers.2.LM.txt.  The genetic algorithm uses this as a staring place.

Then run the code:
`python parallel.genet.alg.optimize.select.contigs.py markers.2.LM.txt > some.output.file`

And at each generation it will output various run statistics, including all current scaffold orders and likelihoods that the algorithm is grinding away on.

At some point it will stop, and you might want to run:

`python run.stats.py some.output`

This gives a plot of various run statistics so you can see how it progressed.  It gives some sense of whether you need to tune things or if it was pretty successful.
