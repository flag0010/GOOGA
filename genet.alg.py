POP_SIZE = 10
NGEN = 1000000000
MUTATION = [0, 1, 2, 3]
ELITE = 2
print 'running population size='+str(POP_SIZE), ', for '+str(NGEN)+' generations, with mutation per generation='+str(MUTATION)+', and saving the best '+str(ELITE)+' individuals at from each generation'
from fitness import fitness
from common import sampler, weighted_sampler, get_file, defaultdict
import os, sha, sys, re, random, copy, time

def swap_mutation(contig_ord, mut_rate = MUTATION):
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
    c = ContigOrder(new, chrom_dict, scaff_lookup)
    return c

def recombination(contig_ord1, contig_ord2):
    #use order crossover (see: http://www.theprojectspot.com/tutorial-post/applying-a-genetic-algorithm-to-the-travelling-salesman-problem/5)
    parents = [[i for i in contig_ord1.chrom_list], [i for i in contig_ord2.chrom_list]]
    random.shuffle(parents)
    contig_ord1, contig_ord2 = parents #just scramble which parent is donor, so p1 isn't always
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
    c = ContigOrder(new, chrom_dict, scaff_lookup)
    return c

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
memo = {}

class ContigOrder:
    def __init__(self, chrom_list, chrom_dict, scaff_lookup):
        self.chrom_list = [i for i in chrom_list]
        self.chrom_dict = {i:[j for j in chrom_dict[i]] for i in chrom_dict}
        self.scaff_lookup = {i:scaff_lookup[i] for i in scaff_lookup}
        self.Fitness = -9999999999999999
        self.tag = sha.sha(str(self.chrom_list)).hexdigest()
    #
    def shuffle(self):
        j = [i for i in self.chrom_list]
        random.shuffle(j)
        self.chrom_list = [i*random.choice([-1, 1]) for i in j]
        self.tag = sha.sha(str(self.chrom_list)).hexdigest()
    #
    def write_file_and_test_fitness(self):
        fnm = sha.sha(str(self.chrom_list)).hexdigest()
        self.tag = fnm
        if fnm in memo: self.Fitness = memo[fnm]
        else:  #havent tested this one yet
            b = open(fnm, 'w')
            for i in self.chrom_list:
                scaff = self.scaff_lookup[abs(i)]
                output = self.chrom_dict[scaff]
                if i < 0: output = list(reversed(output))
                for j in output: b.write(scaff+'\t'+j+'\n')
            b.close()
            myfitness = fitness(fnm)[-1]
            self.Fitness = myfitness
            os.remove(fnm)
            memo[fnm] = myfitness
    #
    def output_scaff_order(self):
        output = []
        for i in self.chrom_list:
            s = scaff_lookup[abs(i)]
            if i < 0: strand = '-'
            else: strand = '+'
            output.append(strand+s)
        return output



c = ContigOrder(chrom_list, chrom_dict, scaff_lookup) 
population = [copy.copy(c) for i in xrange(POP_SIZE)]
#print population[0].chrom_list
#population[0].shuffle()
#print population[0].chrom_list
#population[0].write_file_and_test_fitness()
#print memo
#print population[0].fitness  #it all works!

for i in range(1, POP_SIZE): #shuffle all but the first one, leave that at whatever is in the file (prob. v2 genome order)
    population[i].shuffle()  #randomize
#print [i.tag for i in population]

start  = time.time()
for gen in xrange(NGEN):
    for i in range(POP_SIZE):
        print "generation="+str(gen+1)+';', 'individual='+str(i+1)+';', ' elapsed time (sec):', time.time()-start
        population[i].write_file_and_test_fitness()
    population.sort(key = lambda s: s.Fitness*-1)
    weights = list(reversed(range(1, len(population)+1)))
    weight_dict = {i:weights[i] for i in range(len(population))} #use ranked based selection (see :http://www.obitko.com/tutorials/genetic-algorithms/selection.php)
    new_population = population[:ELITE]
    for i in range(ELITE, POP_SIZE):
        c1 = weighted_sampler(weight_dict)
        c2 = c1
        while c2 == c1:
            c2 = weighted_sampler(weight_dict)
        cnew = swap_mutation(population[c1])
        cnew = recombination(cnew, population[c2])
        new_population.append(cnew)
    seen = []
    tmp_pop = []
    for i in new_population:
        if i.tag not in seen:
            seen.append(i.tag)
            tmp_pop.append(i)
        else:
            cnew  = ContigOrder(chrom_list, chrom_dict, scaff_lookup)
            cnew.shuffle()
            tmp_pop.append(cnew)
    new_population = tmp_pop
    print "generation="+str(gen+1), 'results:'
    for i in range(len(population)):
        print 'individual='+str(i+1), 'fitness='+str(population[i].Fitness), 'order=', ' '.join(population[i].output_scaff_order())
    population = new_population

#found in 11 gen
#11 0 10 -416.58822347 [-1, 2, -3, 4, 6, 10, -7, 8, -9, 5, 11, 12, -13, 14, 15, 16, 17]
#about 2.2 hrs in