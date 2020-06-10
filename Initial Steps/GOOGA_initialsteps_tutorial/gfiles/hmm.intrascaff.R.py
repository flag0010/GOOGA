# likelihoods based on window calls as input
# all plants for a scaff together in all*txt
# reports likelihood of observed states in forward_backward

#v2:: bounded optimization of error rates

from scipy import optimize
from scipy.optimize import minimize_scalar
from scipy.special import gammaln
from math import exp,log
import sys

#genotyping error probs
zy=0.00001 # edge for bounds
MaxRR = 0.25


def calc_v0(r_rates): 
	def scipy_ln_like0(x):
		return -LL(x)

	bounds = [ (0.0,MaxRR) for k in range(len(r_rates)) ]
	best, val, d = optimize.fmin_l_bfgs_b(scipy_ln_like0, r_rates, approx_grad=True, bounds=bounds)
	solution = list(best)
	ln_l = -scipy_ln_like0(solution)
	solution.append(ln_l)
	#zbob=ln_like0(parents,famStr,RRL,RAL,AAL,FLnum,1, list(best),matplant)
	return solution


def foward_backward(obs, states, start_p,transition_probability,er):

	alpha=[{} for j in range(len(obs))] # forward:: alpha[j][X] is probability that true genotye is X at marker j (starts at 0)
	lnFactor=0.0

	for y in states:
		alpha[0][y] = start_p[y] * emission_probability(y,obs[0],er)

	for t in xrange(1, len(obs)):
        	for y in states:
			alpha[t][y] = 0.0
			for y0 in states: # y0 is state at ( t-1 )
				alpha[t][y] +=alpha[t-1][y0] * transition_probability[t-1][y0][y] * emission_probability(y,obs[t],er)

		normalizer = max(alpha[t]['AA'],alpha[t]['AB'],alpha[t]['BB']) #print alpha[t]['AA'],alpha[t]['AB'],alpha[t]['BB']
		lnFactor+=log(normalizer)
        	for y in states:
			alpha[t][y] = alpha[t][y]/normalizer

				
	# Likelihood of observed states
	LLobs=lnFactor+log(alpha[len(obs)-1]['AA']+alpha[len(obs)-1]['AB']+alpha[len(obs)-1]['BB'])


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
	transition_probability=[{} for j in xrange(len(x))] # global that is updated within LL(x)
	for x1 in xrange(len(x)): # recom rates
		#dist=abs(Position[plantID][v1s][x1+1]-Position[plantID][v1s][x1])
		r = x[x1] #rbp*float(dist)
		transition_probability[x1] ={'AA' : {'AA':(1-r)**2.0,'AB':2*r*(1-r),'BB':r**2.0}, 'AB' : {'AA':r*(1-r),'AB':(1-r)**2.0 + r**2.0,'BB':r*(1-r)}, 'BB' : {'AA':r**2.0,'AB':2*r*(1-r),'BB':(1-r)**2.0} }

	Total_LL=0.0
	for j in range(len(f2plants)):
		plantID=f2plants[j] # updated for each scaff
		ER=[Error_Rates[plantID][0],Error_Rates[plantID][1],Error_Rates[plantID][2]]
		llx=foward_backward(obsA[plantID],states,start_probability,transition_probability,ER)
		#print plantID,"LL= ",llx
		#print "forward ",fprbs
		#print "backward ",rprbs
		#postProb=[{} for j in range(len(obsA[plantID][v1s]))] # forward:: alpha[j][X] is probability that true genotye is X at marker j (starts at 0)

		#for j in range(len(fprbs)):
		#	denom=0.0
		#	for y in states: 
		#		denom+=(fprbs[j][y]*rprbs[j][y])
		#	for y in states: 
		#		postProb[j][y]=(fprbs[j][y]*rprbs[j][y])/denom

		#print postProb


		Total_LL+=llx
	#print x,Total_LL
	return Total_LL


####################################################################################################
### Main Program

states = ('AA','AB','BB')
start_probability = {'AA':0.25,'AB':0.5,'BB':0.25}


Error_Rates={}
inZ = open(sys.argv[1],"rU") # "nic.f2.error.rates.txt"
for line_idx, line in enumerate(inZ):
	cols = line.replace('\n', '').split('\t') 
#F2_53	683	0.11882	0.00001	0.00630	-266.92
	Error_Rates[cols[0]]=[float(cols[2]),float(cols[3]),float(cols[4])]

#any valid genotype file
inK=open(sys.argv[2],"rU") # "g.F2_412.txt"
v1contig={}
for line_idx, line in enumerate(inK):
	cols = line.replace('\n', '').split('\t') 
	try:
		v1contig[cols[1]].append(cols[2])
	except KeyError:
		v1contig[cols[1]]=[cols[2]]

	# F2_392	1	73254	BB



out1 = open(sys.argv[4],"w") # "MLE.iscaff.nic.f2.txt"

for vcx in v1contig:

	if len(v1contig[vcx])>1: # more than one marker

		Position={}
		TotalMarkers=len(v1contig[vcx])

		for j in range(TotalMarkers): 
			Position[vcx+"_"+v1contig[vcx][j]]=j

		R_rates=[0.01 for j in range(TotalMarkers-1)] # recomb rates to be estimated

		obsA={}
		f2plants=[]

		srx  =open(sys.argv[3], "rU") # sys.argv[] "low.error.f2s.txt"
		for line_idx, line in enumerate(srx):
			colx = line.replace('\n', '').split('\t') 
			# F2_53
			plantID=colx[0]
			f2plants.append(plantID)
			src  =open("g."+plantID+".txt", "rU")
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
		srx.close()

		zsol= calc_v0(R_rates) 

		
		for j in range(len(zsol)):
			out1.write(vcx+'\t'+v1contig[vcx][j]+'\t'+str(zsol[j])+'\n')
		print vcx,zsol



