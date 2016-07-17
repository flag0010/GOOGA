from matplotlib import pyplot as plt
from common import get_file
import sys
from collections import deque, defaultdict
import itertools

def moving_average(iterable, n=3):
    it = iter(iterable)
    d = deque(itertools.islice(it, n-1))
    d.appendleft(0)
    s = sum(d)
    for elem in it:
        s += elem - d.popleft()
        d.append(elem)
        yield s / float(n)

def yield_N(x, n):
    idx, xidx = 0, 0
    k = {}
    k[idx] = []
    while xidx < len(x):
        while len(k[idx]) < n and xidx < len(x):
            k[idx].append(tuple(x[xidx]))
            xidx+=1
        idx+=1
        k[idx] = []
    for i in k.keys():
        if not k[i]: del(k[i])
    return k

x = get_file(sys.argv[1])
#for i in x: print i
etime = [i for i in x if 'elapsed' in i]
gen = [int(i[0].replace('generation=', '').replace(';', '')) for i in etime]
sec_cum = [float(i[-1]) for i in etime]
last = 0
sec_per = []
for i in sec_cum:
    sec_per.append(i-last)
    last = i

indv_times = [i for i in x if "run" in i and 'time' in i and '(sec)=' in i]
qx = [j for j in x if j]
indvs = [i for i in qx if 'individual=' in i[0]]

g = defaultdict(list)

for i,j in zip(indvs, indv_times):
    g[i[0]].append(float(j[-1])) 

worst_pos = []
worst_per_gen = []
for i in range(len(g[g.keys()[0]])):
    bad = 0
    for indv in g:
        if g[indv][i] > bad: bad = g[indv][i]
    worst_pos.append( i + 1 )
    worst_per_gen.append( bad )

mv_ave = list(moving_average(sec_per, n=10))
plt.subplot(3,1,1)
plt.plot(gen, sec_per)

plt.ylabel('per generation runtime (sec)')
plt.plot(range(1, len(mv_ave)+1), mv_ave, color='r')

plt.plot(worst_pos, worst_per_gen, color='g')


elite_line = [i for i in x if 'saving' in i][0]
elite_num = int(elite_line[elite_line.index('best')+1])


indv_to_grab = ['individual='+str(i) for i in range(1,elite_num+1)]
#print indv_to_grab
hit_lines, elite_change = [], []
for idx, i in enumerate(x):
    for j in indv_to_grab:
        if j in i: hit_lines.append(idx+1)

hit_gen = yield_N([x[i] for i in hit_lines], elite_num)
curr_elite = set()
new_elite = []
for gen in sorted(hit_gen):
    n = len(set(hit_gen[gen]) - curr_elite)
    new_elite.append(n)
    curr_elite = set(hit_gen[gen])

plt.subplot(3,1,2)
#print new_elite
gen = [int(i[0].replace('generation=', '').replace(';', '')) for i in etime]
plt.plot(gen, new_elite)
plt.ylabel('new elite lines added per gen')
#plt.subplot(4,1,3)
#plt.plot(gen, sec_cum)
#plt.xlabel('generation')
#plt.ylabel('cumulative runtime (sec)')


best = [float(i[1].replace('fitness=', '')) for i in x if 'individual=1' in i]
plt.subplot(3,1,3)
plt.plot(gen, best)
#print best
plt.ylabel('fitness of best order (lnLk)')
plt.show()
print max(best)
