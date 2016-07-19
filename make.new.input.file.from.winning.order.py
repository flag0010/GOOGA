from common import *
import sys, json

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

class ContigOrder:
    #made a class to handle scaffold order (oops, called it contig)
    #the class stores one order, and has basic functions to track fitness, randomize order, and output
    def __init__(self, chrom_list, chrom_dict, scaff_lookup):
        self.chrom_list = [i for i in chrom_list]
        self.chrom_dict = {i:[j for j in chrom_dict[i]] for i in chrom_dict}
        self.scaff_lookup = {i:scaff_lookup[i] for i in scaff_lookup}
    #
    def write_file_and_test_fitness(self):
        #this writes a file, and calls JKK's fitness code and points it to the file
        #and retrieves the lnL (aka fitness) and then cleans up the file
        #also uses a cache ("memo") to check if the fitness has already been calculated for a particular order
        #this saves on the expensive compute
        ##6/27/16, further modifications also now allow it to search through precomputed local values and pass in those
        ##in place of full estimates.
        if 0:1 
        else: 
            for i in self.chrom_list:
                scaff = self.scaff_lookup[abs(i)]
                output = self.chrom_dict[scaff]
                if i < 0: output = list(reversed(output))
                for j in output: print scaff+'\t'+j

#print chrom_list
#print chrom_dict
#print scaff_lookup

q = ContigOrder(chrom_list, chrom_dict, scaff_lookup)

new_ord =  json.loads(sys.argv[2])

q = ContigOrder(new_ord, chrom_dict, scaff_lookup)

q.write_file_and_test_fitness()
