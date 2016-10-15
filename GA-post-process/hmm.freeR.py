# likelihoods based on window calls as input
# all plants for a scaff together in all*txt
# reports likelihood of observed states in forward_backward

#v2:: bounded optimization of recom rates

from scipy import optimize
from scipy.optimize import minimize_scalar
from scipy.special import gammaln
from math import exp,log
import sys, os, common


GAOutFile=sys.argv[1]
MarkerInit = sys.argv[2]
F2file = sys.argv[3]
Errorfile = sys.argv[4]
#genotyping error probs
zy=0.00001 # edge for bounds

def calc_v0(r_rates): 
	def scipy_ln_like0(x):
		return -LL(x)

	bounds = [ (0.0,0.25) for k in range(len(r_rates)) ]
	best, val, d = optimize.fmin_l_bfgs_b(scipy_ln_like0, r_rates, approx_grad=True, bounds=bounds)
	solution = list(best)
	ln_l = -scipy_ln_like0(solution)
	solution.append(ln_l)
	return solution



def foward_backward(obs, states, start_p,transition_probability,er):

	alpha=[{} for j in range(len(obs))] # forward:: alpha[j][X] is probability that true genotye is X at marker j (starts at 0)
	lnFactor=0.0

	for y in states:
		alpha[0][y] = start_p[y] * emission_probability(y,obs[0],er)

	for t in xrange(1, len(obs)):
        	for y in states:
			alpha[t][y] = 0.0
			for y0 in states: # y0 is state at t-1

				alpha[t][y] +=alpha[t-1][y0] * transition_probability[t-1][y0][y] * emission_probability(y,obs[t],er)

		normalizer = max(alpha[t]['AA'],alpha[t]['AB'],alpha[t]['BB'])
		#print alpha[t]['AA'],alpha[t]['AB'],alpha[t]['BB']
		lnFactor+=log(normalizer)
        	for y in states:
			alpha[t][y] = alpha[t][y]/normalizer

				
	# Likelihood of observed states
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

	print x
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


		Total_LL+=llx
	#print x,Total_LL
	return Total_LL


####################################################################################################
### Main Program

states = ('AA','AB','BB')
start_probability = {'AA':0.25,'AB':0.5,'BB':0.25}


Error_Rates={}
inZ = open(Errorfile,"rU")
for line_idx, line in enumerate(inZ):
	cols = line.replace('\n', '').split() 
#imswc001 482 0.0608584706859 1e-05 0.0180081186063 -159.05915623
	Error_Rates[cols[0]]=[float(cols[2]),float(cols[3]),float(cols[4])]


Qi =  common.get_file(GAOutFile, 'adfadfa')

xx = [idx for idx, i in enumerate(Qi) if 'individual=1' in i[0]]
#print xx
final =  Qi[xx[-1]+2][0].replace('raw list order= ', '')
os.popen('python make.new.input.file.from.winning.order.py '+MarkerInit+' '+'"'+final+'" > better.'+ MarkerInit)

Position={}
inY = open('better.'+ MarkerInit,"rU")
TotalMarkers=0
for line_idx, line in enumerate(inY):
	cols = line.replace('\n', '').split('\t') 
	Position[cols[0]+"_"+cols[1]]=line_idx
	TotalMarkers+=1

R_rates=[0.01 for j in range(TotalMarkers-1)] # recomb rates to be estimated

obsA={}
f2plants=[]

srx  =open(F2file, "rU")
for line_idx, line in enumerate(srx):
	colx = line.replace('\n', '').split('\t') 
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

zsol= calc_v0(R_rates) 

out1 = open("MLE.free."+MarkerOrder,"w")
for j in range(len(zsol)):
	out1.write(str(zsol[j])+'\n')
print zsol




