POP_SIZE = 12  #pop size
NGEN = 500 #generations to run
MUTATION = [0, 1, 2, 3]  #randomly select one value from this list to determine the number of mutations an indiv. pass on to next gen.
ELITE = 2  #number best individuals to save at each generation
NCPUs = 10 #number of cpus to run on.  Remember that after the 1st generation you will have POP_SIZE - ELITE novel indv
          #since we save past results, on a machine with 10 CPUs, if ELITE=2, you may want to do a popsize of 12, because that will max out all 10 CPUs after Gen 1 

#import statements
from fitness_optimize_select_contigs import fitness ##repackaging of John's likelihood calc. code
from common import sampler, weighted_sampler, get_file, defaultdict
import os, sha, sys, re, random, copy, time, multiprocessing
from cPickle import dump

### FUNCTIONS FOR TRACKING REC. RATES AND RECYCLING PREV. ESTIMATED VALUES GREEDILY##########
def greedy_slices(x):
    out = []
    idx = len(x)
    while idx > 1:
        i = 0
        j = i+idx
        while j <= len(x):
            out.append(tuple(x[i:j]))
            i+=1
            j+=1
        idx-=1
    return out

def update_subset_memo(contig_ord, r_rates):
    def flip_scaff_ord(s):
        return tuple([i*-1 for i in reversed(s)])    
    for i,j in zip(*map(greedy_slices, [contig_ord,r_rates[1:]])):
        if i not in subset_memo: subset_memo[i] = j
        flip_i = flip_scaff_ord(i)
        if flip_i not in subset_memo: subset_memo[flip_i] = tuple(reversed(j))

def fill_in_rates_return_UGaps_and_new_R_rates(contig_ord, subset_memo, intra_scaff = False):
    out = {}
    lkp_pos = {i:idx for idx, i in enumerate(contig_ord)}
    for span in greedy_slices(contig_ord):
        if span in subset_memo:
            for idx, scaff in enumerate(span[:-1]):
                pos = lkp_pos[scaff]
                if pos not in out: out[pos] = subset_memo[span][idx]
    #print 'out', out
    xout = [0.01 for ii in range(len(contig_ord)-1)]
    #need_to_estimate = sorted(set(range(len(contig_ord)-1)).difference(out.keys()))
    pre_estimated =[i+1 for i in sorted( out.keys() )]
    for i in out:
        xout[i] = out[i]
    #print 'xout', xout
    if intra_scaff:
        xout.insert(0, intra_scaff)
        pre_estimated.insert(0,0)
    else:
        xout.insert(0, 0.1)
    UGaps_out, R_rates_out = [pre_estimated,len(contig_ord)], xout
    print 'subset_memo len and memo len'
    print len(subset_memo), len(memo)
    print "UGaps and R_rates"
    print UGaps_out, R_rates_out
    print '#########################################################################'
    return UGaps_out, R_rates_out
#####END##############################################################################

###FUNCTIONS FOR PERFORMING MUTATION AND RECOMBINATION IN GENET. ALGO.################
def swap_mutation(contig_ord, memo, subset_memo, mut_rate = MUTATION):
    #for swap mutation protocol see http://www.theprojectspot.com/tutorial-post/applying-a-genetic-algorithm-to-the-travelling-salesman-problem/5
    a = sampler(mut_rate, 1)[0]
    new = [i for i in contig_ord.chrom_list]
    pos = range(len(new))
    random.shuffle(pos)
    for i in range(a):
        sign = random.choice([-1, 1]) 
        p1 = pos.pop(0)
        p2 = pos.pop(0)
        tmp = new[p1]
        new[p1] = new[p2]*sign
        new[p2] = tmp*sign
    c = ContigOrder(new, chrom_dict, scaff_lookup, memo, subset_memo)
    return c

def recombination(contig_ord1, contig_ord2, memo, subset_memo):
    #use order crossover (see: http://www.theprojectspot.com/tutorial-post/applying-a-genetic-algorithm-to-the-travelling-salesman-problem/5)
    parents = [[i for i in contig_ord1.chrom_list], [i for i in contig_ord2.chrom_list]]
    random.shuffle(parents)
    contig_ord1, contig_ord2 = parents #just scramble which parent is donor, so p2 isn't always
    new = ['-' for i in contig_ord1]
    pos = range(len(new))
    a,b = sorted(sampler(pos, 2))
    sval = set([abs(i) for i in contig_ord2[a:b+1]])
    sidx = set(range(a,b+1))
    for i in range(a, b+1): new[i] = contig_ord2[i]
    picks = [i for i in contig_ord1 if abs(i) not in sval]
    for i in [j for j in range(len(contig_ord1)) if j not in sidx]:
        p = picks.pop(0)
        new[i] = p
    c = ContigOrder(new, chrom_dict, scaff_lookup, memo, subset_memo)
    return c
############################END#######################################################################

#####READ IN DATA FILES AND PREP FOR GA###############################################################
init_file = get_file(sys.argv[1])
chrom_dict = defaultdict(list)
scaff_order = []
for i,j in init_file:
    chrom_dict[i].append(j)
    if i not in scaff_order: scaff_order.append(i)

chrom_list = []
scaff_lookup = {}
idx = 1
for i in scaff_order:
    s = sorted(chrom_dict[i], key = lambda s: int(s))
    if chrom_dict[i] == s:
        chrom_list.append(idx)
    else: chrom_list.append(idx*-1)
    chrom_dict[i] = s
    scaff_lookup[idx] = i
    idx+=1
#print chrom_list
#print [scaff_lookup[abs(i)] for i in chrom_list]

class ContigOrder:
    #made a class to handle scaffold order (oops, called it contig)
    #the class stores one order, and has basic functions to track fitness, randomize order, and output
    def __init__(self, chrom_list, chrom_dict, scaff_lookup, memo, subset_memo):
        self.chrom_list = [i for i in chrom_list]
        self.chrom_dict = {i:[j for j in chrom_dict[i]] for i in chrom_dict}
        self.scaff_lookup = {i:scaff_lookup[i] for i in scaff_lookup}
        self.Fitness = -99999999
        self.tag = sha.sha(str(self.chrom_list)).hexdigest()
        #print self.tag
        self.memo = {i:j for i,j in memo.items()}
        self.subset_memo = {i:j for i,j in subset_memo.items()}
        mygaps, myrates = fill_in_rates_return_UGaps_and_new_R_rates(self.chrom_list, self.subset_memo)
        self.Rates = myrates
        self.UGaps = mygaps
        self.info = {}
        #print 'init'
        #print self.chrom_list
        #print self.Rates
        #print self.UGaps
    #
    def shuffle(self):
        #use to scramble order
        #c = ContigOrder(chrom_list, chrom_dict, scaff_lookup)
        #c.shuffle()
        #now c is randomized
        j = [i for i in self.chrom_list]
        random.shuffle(j)
        self.chrom_list = [i*random.choice([-1, 1]) for i in j]
        self.tag = sha.sha(str(self.chrom_list)).hexdigest()
    #
    def write_file_and_test_fitness(self):
        #this writes a file, and calls JKK's fitness code and points it to the file
        #and retrieves the lnL (aka fitness) and then cleans up the file
        #also uses a cache ("memo") to check if the fitness has already been calculated for a particular order
        #this saves on the expensive compute
        ##6/27/16, further modifications also now allow it to search through precomputed local values and pass in those
        ##in place of full estimates.
        fnm = sha.sha(str(self.chrom_list)).hexdigest()
        self.tag = fnm
#        print 'running', self.chrom_list, self.tag
        if fnm in self.memo:
            print 'hit:', fnm
            self.Fitness = self.memo[fnm][0]
            self.Rates = self.memo[fnm][1]
            #update_subset_memo(self.chrom_list, self.Rates)
        else:  #we haven't tested this one yet
            b = open(fnm, 'w')
            for i in self.chrom_list:
                scaff = self.scaff_lookup[abs(i)]
                output = self.chrom_dict[scaff]
                if i < 0: output = list(reversed(output))
                for j in output: b.write(scaff+'\t'+j+'\n')
            b.close()
#            print self.tag
            myUGaps, my_R_rates = fill_in_rates_return_UGaps_and_new_R_rates(self.chrom_list, self.subset_memo,  intra_scaff = 0.31)
            #print 'run!'
            #print self.chrom_list
            #print self.Rates
            #print self.UGaps
            rates_and_lnLk = fitness(fnm, my_R_rates, myUGaps)
            full_fitness_and_rates = fitness(fnm, my_R_rates, [[], 17])
            info_dump = {'carryoverFitness': rates_and_lnLk[-1], 'fullFitness':full_fitness_and_rates[-1],
                         'carryoverRates': rates_and_lnLk[:-1], 'fullRates':full_fitness_and_rates[:-1], 'ord':self.chrom_list,
                         'num_carryover':len(myUGaps[0]) }
            self.info = info_dump
            #print rates_and_lnLk
            #print subset_memo
            #print memo
            rates = rates_and_lnLk[:-1]
            myfitness = rates_and_lnLk[-1]
            self.Fitness = myfitness
            self.Rates = rates
            os.remove(fnm)
            #memo[fnm] = [myfitness, rates]
            #update_subset_memo(self.chrom_list, self.Rates)
    #
    def output_scaff_order(self):
        #translate and print scaff order
        output = []
        for i in self.chrom_list:
            s = self.scaff_lookup[abs(i)]
            if i < 0: strand = '-'
            else: strand = '+'
            output.append(strand+s)
        return output

#print population[0].chrom_list
#population[0].shuffle()
#print population[0].chrom_list
#population[0].write_file_and_test_fitness()
#print memo
#print population[0].fitness  #it all works!


#print [i.tag for i in population]

def functionalize_write_file_and_test_fitness(x):
    x.write_file_and_test_fitness()
    return x

memo, subset_memo = {}, {}###THE FIRST MEMO TRACKS FULLY PRECOMPUTED CONTIG_ORDERS, WHEREAS THE SECOND "SUBSET_MEMO" TRACKS PARTIAL FRAGMENTS
carryover_info = {}
if __name__ == '__main__':
    c = ContigOrder(chrom_list, chrom_dict, scaff_lookup, memo, subset_memo) 
    population = [copy.copy(c) for i in xrange(POP_SIZE)]  #setting up the popualtion
    for i in range(1, POP_SIZE): #shuffle all but the first one, leave that at whatever is in the file (prob. v2 genome order)
        population[i].shuffle()  #randomize
    print '\nrunning population size='+str(POP_SIZE)+ ', for '+str(NGEN)+' generations\nwith mutation per generation='+str(MUTATION)+', and saving the best '+str(ELITE)+' individuals from each generation', '\non N='+str(NCPUs)+' CPUs' 
    start  = time.time()
    P = multiprocessing.Pool(NCPUs)
    for gen in xrange(NGEN):
        population = P.map(functionalize_write_file_and_test_fitness, population)
        #for i in range(POP_SIZE):
        #    print "generation="+str(gen+1)+';', 'individual='+str(i+1)+';', ' elapsed time (sec):', time.time()-start
        #    population[i].write_file_and_test_fitness()
        print "generation="+str(gen+1)+';', 'elapsed time (sec):', time.time()-start
        population.sort(key = lambda s: s.Fitness*-1)
        for i in population:
            update_subset_memo(i.chrom_list, i.Rates)
            memo[i.tag] = [i.Fitness, i.Rates]
            loc_info = {i:j for i,j in i.info.items()}
            TAG = i.tag
            loc_info['gen'] = gen
            carryover_info[TAG] = loc_info
        dump(carryover_info, open('intra-scaff.carryover.info.pickle', 'w'))
        weights = list(reversed(range(1, len(population)+1)))
        weight_dict = {i:weights[i] for i in range(len(population))} #using ranked based selection (see :http://www.obitko.com/tutorials/genetic-algorithms/selection.php)
        new_population = population[:ELITE]#save the elites
        for i in range(ELITE, POP_SIZE):
            c1 = weighted_sampler(weight_dict)
            c2 = c1
            while c2 == c1:
                c2 = weighted_sampler(weight_dict)
            #c1 and c2 are the 2 indivduals to be recombined
            cnew = swap_mutation(population[c1], memo, subset_memo)# first mutate c1.  could do c2 too, but I didn't for now.  
            cnew = recombination(cnew, population[c2], memo, subset_memo)
            new_population.append(cnew)#add this onto the next generation
        seen = []
        tmp_pop = []
        for i in new_population: #occasionally identical orders get in, this removes duplicates and fills in with a new random order
            if i.tag not in seen:
                seen.append(i.tag)
                tmp_pop.append(i)
            else:
                cnew  = ContigOrder(chrom_list, chrom_dict, scaff_lookup, memo, subset_memo)
                cnew.shuffle()
                tmp_pop.append(cnew)
        new_population = tmp_pop #lazy, but make the next generation from this temp, de-duplicated table
        print "generation="+str(gen+1), 'results:'
        for i in range(len(population)):
            print 'individual='+str(i+1), 'fitness='+str(population[i].Fitness)
            print '\t\t\tscaffold order=', ' '.join(population[i].output_scaff_order())
            print '\t\t\traw list order=', population[i].chrom_list
        population = new_population
        for i in population:
            i.memo = {i:j for i,j in memo.items()}
            i.subset_memo = {i:j for i,j in subset_memo.items()}
