from matplotlib import pyplot as plt
import sys
from common import *
from pprint import pprint
import maps
final_maps, lendict = maps.final_maps, maps.lendict

def plot_ord(scaff_ord, lendict, ypos):
    lenlist = [lendict[i] for i in map(lambda s: s[1:], scaff_ord)]
    colorlist = []
    for i in scaff_ord:
        if '-' in i: colorlist.append('r')
        else: colorlist.append('g')
    pos = 0
    out = {}
    for i in range(len(scaff_ord)):
        scaff, lenx, col = scaff_ord[i], lenlist[i], colorlist[i]
        plt.plot([pos, pos+lenx], [ypos, ypos], color=col, linewidth=2, alpha=.74)
        plt.plot([pos+lenx, pos+lenx], [ypos+1, ypos-1], color='k', alpha=.74)
        midpoint = (lenx * .5) + pos
        if (i+1) % 2:
            plt.text(midpoint-1e5,ypos-2, s=scaff, rotation=90, size=7)
        else:
            plt.text(midpoint-1e5, ypos+2, s = scaff, rotation=90, size=7)
        out[scaff] = midpoint
        pos+=lenx
    plt.plot([0, 0], [ypos+1, ypos-1], color='k', alpha=.74)
    return out

def connections(scaff_ord1, scaff_ord2, mids1, mids2, y1, y2):
    #pprint(scaff_ord1)
    #pprint(scaff_ord2)
    #pprint(mids1)
    #pprint(mids2)
#    def amend_mids(x, newx):
#        for i in x.keys():
#            print i, x[i]
#            newx[i] = x[i]
#            if '-' in x:
#                j = i.replace('-', "+")
#                newx[j] = x[i]
#            else:
#                j = i.replace('+', "-")
#                newx[j] = x[i]
#            print j, newx[i], newx[j]
#        return(newx)
#    print len(mids1), len(mids2)
#    mids1x = amend_mids(mids1, {})
#    mids2x = amend_mids(mids2, {})
    def flip(x):
        if '-' in x: return x.replace('-', '+')
        elif '+' in x: return x.replace('+', '-')
    #print len(mids1), len(mids2)
    #print 
    for i in scaff_ord1:
        plotIT = 0
        if i in mids1:
            m1 = mids1[i]
            plotIT+=1
        elif flip(i) in mids1:
            m1 = mids1[flip(i)]
            plotIT+=1
        if i in mids2:
            m2 = mids2[i]
            plotIT+=1
        elif flip(i) in mids2:
            m2 = mids2[flip(i)]
            plotIT+=1
        #print i
        #pprint(mids1)
        #pprint(mids2)
        #print plotIT
        if plotIT > 1: plt.plot([m1, m2], [y1, y2], color='k', alpha=.2)
 
 
pops = 'v2_genome IMPR_rils IM_ase IM_mse LVR SF SWC'.split()   
#pprint(final_maps.keys())

for chrom in '1 2 3 4 5 6 7 8 9 10 11 12 13 14'.split():
    p = pops
    yidx = 100
    mx = 0
    p0_mid = plot_ord(final_maps[p[0]][chrom], lendict, yidx)
    p1_mid = plot_ord(final_maps[p[1]][chrom], lendict, yidx-10)
    p2_mid = plot_ord(final_maps[p[2]][chrom], lendict, yidx-20)
    p3_mid = plot_ord(final_maps[p[3]][chrom], lendict, yidx-30)
    p4_mid = plot_ord(final_maps[p[4]][chrom], lendict, yidx-40)
    p5_mid = plot_ord(final_maps[p[5]][chrom], lendict, yidx-50)
    p6_mid = plot_ord(final_maps[p[6]][chrom], lendict, yidx-60)
    for i in [p0_mid, p1_mid, p2_mid, p3_mid, p4_mid, p5_mid, p6_mid]:
        if max(i.values()) > mx: mx = max(i.values())
    no_sign = lambda s: [i.replace('-','').replace('+','') for i in s]
    complete = lambda s, ovlp: [i for i in s if i.replace('-','').replace('+','') in ovlp] 
    m0, m1, m2, m3, m4, m5, m6 = map(no_sign, [final_maps[p[i]][chrom] for i in range(7)])
    f0, f1, f2, f3, f4, f5, f6 = [final_maps[p[i]][chrom] for i in range(7)]
    ovlp01, ovlp12, ovlp23, ovlp34, ovlp45, ovlp56 = set(m0).intersection(m1), set(m1).intersection(m2), set(m2).intersection(m3), set(m3).intersection(m4), set(m4).intersection(m5), set(m5).intersection(m6)
    m0complete, m1complete = complete(f0, ovlp01), complete(f1, ovlp01)
    print m0, m1, ovlp01, m0complete, m1complete
    connections(m0complete, m1complete, p0_mid, p1_mid, yidx, yidx-10)
    
    m1complete, m2complete = complete(f1, ovlp12), complete(f2, ovlp12)
    connections(m1complete, m2complete, p1_mid, p1_mid, yidx-10, yidx-20)
    
    m2complete, m3complete = complete(f2, ovlp23), complete(f3, ovlp23)
    connections(m2complete, m3complete, p2_mid, p3_mid, yidx-20, yidx-30)
    
    m3complete, m4complete = complete(f3, ovlp34), complete(f4, ovlp34)
    connections(m3complete, m4complete, p3_mid, p4_mid, yidx-30, yidx-40)
    
    m4complete, m5complete = complete(f4, ovlp45), complete(f5, ovlp45)
    connections(m4complete, m5complete, p4_mid, p5_mid, yidx-40, yidx-50)
    
    m5complete, m6complete = complete(f5, ovlp56), complete(f6, ovlp56)
    connections(m5complete, m6complete, p5_mid, p6_mid, yidx-50, yidx-60)
    
    plt.xlim(xmin=-1250000, xmax = mx+1250000)
    plt.title('Chromosome '+ chrom)
    plt.show()

    
for chrom in '1 2 3 4 5 6 7 8 9 10 11 12 13 14'.split():
    yidx = 100
    mx = 0
    for p in pops[1:]:
        v2_mid = plot_ord(final_maps['v2_genome'][chrom], lendict, yidx)
        p_mid = plot_ord(final_maps[p][chrom], lendict, yidx-10)
        #m1, m2 = [i.replace('-','').replace('+','') for i in v2_mid.keys()],  [i.replace('-','').replace('+','') for i in p_mid.keys()]
        #for i in v2_mid.keys():
        #    if i not in p_mid: del(v2_mid[i])
        #for i in p_mid.keys():
        #    if i not in v2_mid: del(p_mid[i])
        #pprint([(i, v2_mid[i]) for i in final_maps['v2_genome'][chrom]])
        #pprint([(i, p_mid[i]) for i in final_maps[p][chrom]])
        #print '*'*40
        if max(v2_mid.values()) > mx: mx =  max(v2_mid.values())
        if max(p_mid.values()) > mx: mx =  max(p_mid.values())
        #print map(len,[v2_mid, p_mid])
        #print map(len, [final_maps['v2_genome'][chrom], final_maps[p][chrom]])
        #print 
        m1, m2 = [i.replace('-','').replace('+','') for i in final_maps['v2_genome'][chrom]],  [i.replace('-','').replace('+','') for i in final_maps[p][chrom]]
        ovlp = set(m1).intersection(m2)
        v2complete = [i for i in final_maps['v2_genome'][chrom] if i.replace('-','').replace('+','') in ovlp]
        pcomplete = [i for i in final_maps[p][chrom] if i.replace('-','').replace('+','') in ovlp]
        #print map(len, [v2complete, pcomplete])
        connections(v2complete, pcomplete, v2_mid, p_mid, yidx, yidx-10)     
        yidx-=24
    plt.xlim(xmin=-1250000, xmax = mx+1250000)
    plt.title('Chromosome '+ chrom)
    plt.show()

# z = 0
# for i in first:
#     #print i, z
#     if z == 1:
#         orig = i[0].replace('scaffold order=', '').split()
#         z+=1
#     if 'individual=1 fitness=' in i[0]:
#         z +=1
# 
# 
# fbest = [idx for idx, i in enumerate(first) if 'individual=1 fitness=' in i[0]]
# sbest = [idx for idx, i in enumerate(second) if 'individual=1 fitness=' in i[0]]
# pprint(sbest)
# 
# origfitf, origfits = first[fbest[0]][0].replace('individual=1 fitness=', ''), second[sbest[0]][0].replace('individual=1 fitness=', '')
# ff, sf = first[fbest[-1]][0].replace('individual=1 fitness=', ''), second[sbest[-1]][0].replace('individual=1 fitness=', '')
# 
# fend =  first[fbest[-1]+1][0].replace('scaffold order=', '').split()
# send =  second[sbest[-1]+1][0].replace('scaffold order=', '').split()
# 
# print orig, origfitf, origfits
# print fend, ff
# print send, sf
# orig_mids = plot_ord(orig, lendict, 20)
# first_mids = plot_ord(fend, lendict, 10)
# second_mids = plot_ord(send, lendict, 0)
# connections(orig, fend, orig_mids, first_mids, 20, 10)
# connections(fend, send, first_mids, second_mids, 10, 0)
# 
