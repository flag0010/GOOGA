# likelihoods based on window calls as input
# all plants for a scaff together in all*txt
# reports likelihood of observed states in forward_backward

#v3.1:: homogenenous recombination within v1scaffs
#v3.1b:: work on exponential scale for rates

import numpy
from scipy import optimize
from scipy.optimize import minimize_scalar
from scipy.special import gammaln
from math import exp,log
import sys

LG=sys.argv[1]

gradval=1e-8
maxrbp=0.1/1000000.0  # recombination rate per bp (morgans / megabase)
prci=10000000.0 # low values mean higher precision default : factr=10000000.0


def calc_v0(r_rates): 
	def scipy_ln_like0(x):
		return -LL(x)

	#recombination rate = max*(10.0**-x)

	bounds = [(0.0,10.0)] # different calculation for intra-scaff rates (per bp)
	for k in range(1,len(r_rates)):
		bounds.append( (0.00000001,0.5) )
	
	best, val, d = optimize.fmin_l_bfgs_b(scipy_ln_like0, r_rates, factr=prci, epsilon=gradval, approx_grad=True, bounds=bounds)
	#print d
	solution = list(best)
	ln_l = -scipy_ln_like0(solution)
	solution.append(ln_l)
	#zbob=ln_like0(parents,famStr,RRL,RAL,AAL,FLnum,1, list(best),matplant)
	return solution



def foward_backward(obs, states, start_p,transition_probability,er):

	alpha=[{} for j in range(len(obs))] # forward:: alpha[j][X] is probability that true genotye is X at marker j (starts at 0)
	lnFactor=0.0

	#print transition_probability
	#print obs
	for y in states:
		alpha[0][y] = start_p[y] * emission_probability(y,obs[0],er)

	for t in xrange(1, len(obs)):
        	for y in states:
			alpha[t][y] = 0.0
			for y0 in states: # y0 is state at t-1

				alpha[t][y] +=alpha[t-1][y0] * transition_probability[t-1][y0][y] * emission_probability(y,obs[t],er)

		normalizer = max(alpha[t]['AA'],alpha[t]['AB'],alpha[t]['BB'])
		#print t, alpha[t]['AA'],alpha[t]['AB'],alpha[t]['BB']
		lnFactor+=log(normalizer)
        	for y in states:
			alpha[t][y] = alpha[t][y]/normalizer

				
	LLobs=lnFactor+log(alpha[len(obs)-1]['AA']+alpha[len(obs)-1]['AB']+alpha[len(obs)-1]['BB'])



	#return alpha,beta,LLobs
	return LLobs

def emission_probability(genotype,calledG,ER): # cc [ AA,AB,BB,NN ] 

	e1 = ER[0] # probability of sequencing error to het
	e2 = ER[1]
	beta=ER[2]

	if calledG == 'NN':
		return 1.0
	elif calledG =='AA':
		if genotype=='AA':
			prob = 1 - e1 - e2 
		elif genotype=='AB':
			prob = beta/2
		elif genotype=='BB':
			prob = e2
	elif calledG =='AB':
		if genotype=='AA' or genotype=='BB':
			prob = e1 
		elif genotype=='AB':
			prob = 1-beta
	elif calledG =='BB':
		if genotype=='AA':
			prob = e2 
		elif genotype=='AB':
			prob = beta/2
		elif genotype=='BB':
			prob = 1-e1-e2

	return prob
	

 
def LL(x):

	# transition probs a global
	transition_probability=[{} for j in xrange(TotalMarkers-1)] # global that is updated within LL(x)
	bsum=0
	for x1 in xrange(TotalMarkers-1): # recom rates
		if distX[x1]==-9:
 			bsum+=1
			r = x[bsum]
		else:
			#dist=abs(Position[plantID][v1s][x1+1]-Position[plantID][v1s][x1])
			r = float(distX[x1])*maxrbp*(10.0**-x[0])

		transition_probability[x1] ={'AA' : {'AA':(1-r)**2.0,'AB':2*r*(1-r),'BB':r**2.0}, 'AB' : {'AA':r*(1-r),'AB':(1-r)**2.0 + r**2.0,'BB':r*(1-r)}, 'BB' : {'AA':r**2.0,'AB':2*r*(1-r),'BB':(1-r)**2.0} }

	#print transition_probability
	
	Total_LL=0.0
	for j in range(len(f2plants)):
		plantID=f2plants[j] # 
		
		ER=[Error_Rates[plantID][0],Error_Rates[plantID][1],Error_Rates[plantID][2]]
		llx=foward_backward(obsA[plantID],states,start_probability,transition_probability,ER)



		Total_LL+=llx
	#print x,Total_LL
	return Total_LL


####################################################################################################
### Main Program

states = ('AA','AB','BB')
start_probability = {'AA':0.25,'AB':0.5,'BB':0.25}


Error_Rates={}
inZ = open("error.rates.txt","rU")
for line_idx, line in enumerate(inZ):
	cols = line.replace('\n', '').split('\t') 
#imswc001 482 0.0608584706859 1e-05 0.0180081186063 -159.05915623
	Error_Rates[cols[0]]=[float(cols[2]),float(cols[3]),float(cols[4])]


Position={}
inY = open("markers."+LG+".LM.txt","rU")
TotalMarkers=0
mxz=[]
Breaks=[]
distX=[]
for line_idx, line in enumerate(inY):
	cols = line.replace('\n', '').split('\t') 
	if line_idx==0:
		cscaff=cols[0]
		lastpos=int(cols[1])
	else:
		if cols[0] != cscaff:
			Breaks.append(line_idx)
			distX.append(-9)
		else:
			distX.append(abs(int(cols[1])-lastpos))
		cscaff=cols[0]
		lastpos=int(cols[1])
	Position[cols[0]+"_"+cols[1]]=line_idx
	mxz.append(cols[0]+"_"+cols[1])
	TotalMarkers+=1

#print distX
interScaff_intervals=len(Breaks)
R_rates=[0.1]  # recombination rate per bp (morgans / megabase) within scaffs
for k in range(interScaff_intervals):
	R_rates.append(0.01) # recomb rates between scaffs

print "No recomb parameters ",len(R_rates)
print "Recomb parameters ",R_rates

obsA={}
f2plants=[]

srx  =open("test.f2group.txt", "rU")
#srx  =open("x30.txt", "rU")
for line_idx, line in enumerate(srx):
	colx = line.replace('\n', '').split('\t') 
	plantID=colx[0]
	f2plants.append(plantID)
	src  =open("Genotypes."+plantID+".txt", "rU")
	for line_idx, line in enumerate(src):
		cols = line.replace('\n', '').split('\t') 
	# isg480	1	400000	AB
		if plantID!=cols[0]:
			print "Whoa"
		if line_idx==0:
			obsA[plantID]=["NN" for j in range(TotalMarkers)]

		try:
			qx=Position[cols[1]+"_"+cols[2]]
			obsA[plantID][qx]=cols[3]
		except KeyError:
			pass


zsol= calc_v0(R_rates) 
out1 = open("MLE.1."+LG+".txt","w")
for j in range(len(zsol)):
	out1.write(str(zsol[j])+'\n')
print zsol



