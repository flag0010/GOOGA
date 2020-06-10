# Given map (and rates), produce genotype posterior probabilities
import sys
from math import exp,log

def foward_backward(obs, states, start_p,transition_probability,er):


	alpha=[{} for j in range(len(obs))] # forward:: alpha[j][X] is probability that true genotye is X at marker j (starts at 0)
	beta= [{} for j in range(len(obs))] # backward:: beta[j][X] is probability that true genotye is X at marker j (starts at 0)

	lnFactor=0.0

	for y in states:
		alpha[0][y] = start_p[y] * emission_probability(y,obs[0],er)
		beta[len(obs)-1][y] = 1.0

	if Cross_type==0: # F2
		for t in xrange(1, len(obs)):
			for y in states:
				alpha[t][y] = 0.0
				for y0 in states: # y0 is state at t-1
					alpha[t][y] +=alpha[t-1][y0] * transition_probability[t-1][y0][y] * emission_probability(y,obs[t],er)

			normalizer = 1.0 # max(alpha[t]['AA'],alpha[t]['AB'],alpha[t]['BB'])
			lnFactor+=log(normalizer)
			for y in states:
				alpha[t][y] = alpha[t][y]/normalizer

		LLobs=lnFactor+log(alpha[len(obs)-1]['AA']+alpha[len(obs)-1]['AB']+alpha[len(obs)-1]['BB'])


		for t in xrange(len(obs)-2,-1,-1):
			for y in states:
				beta[t][y] = 0.0 # y is state at t
				for y0 in states: # y0 is state at t+1
					beta[t][y] +=beta[t+1][y0] * transition_probability[t][y][y0] * emission_probability(y0,obs[t+1],er)



	elif Cross_type==1: # RIL

		for t in xrange(1, len(obs)):
			for y in states:
				alpha[t][y] = 0.0
				for y0 in states: # y0 is state at t-1
					alpha[t][y] +=alpha[t-1][y0] * transition_probability[t-1][y0][y] * emission_probability(y,obs[t],er)


			normalizer = max(alpha[t]['AA'],alpha[t]['BB'])
			lnFactor+=log(normalizer)

			for y in states:

				alpha[t][y] = alpha[t][y]/normalizer

	

		LLobs=lnFactor+log(alpha[len(obs)-1]['AA']+alpha[len(obs)-1]['BB'])	


	return LLobs,alpha,beta




def emission_probability(genotype,calledG,ER): # cc [ AA,AB,BB,NN ] 


	if Cross_type==0:

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

	elif Cross_type==1:
		e2 = ER[0] 

		if calledG == 'NN':
			return 1.0

		elif calledG =='AA':
			if genotype=='AA':
				prob = 1 - e2 

			elif genotype=='BB':
				prob = e2

		elif calledG =='BB':

			if genotype=='AA':
				prob = e2 

			elif genotype=='BB':
				prob = 1-e2

	return prob



######################################## Main program

LG=sys.argv[1]
error_rate_file=sys.argv[2]
plant_list=sys.argv[3] # "all.f2s.txt"
map_file=sys.argv[4] # MLE.v2.txt"

Cross_type=0 # 0 = f2s
TotalMarkers=0

out1 = open(LG+".pp.txt","w")  # relevant map


inZ = open(error_rate_file,"rU")
srx  =open(plant_list, "rU")


Error_Rates={}
for line_idx, line in enumerate(inZ):
	cols = line.replace('\n', '').split('\t') 
	Error_Rates[cols[0]]=[float(cols[2]),float(cols[3]),float(cols[4])]

r_list=[]
Position={}
MINFO=[]
cc=0
inY = open(map_file,"rU")  # relevant map
for line_idx, line in enumerate(inY): 
# 11	224041	0.0224242323879
# 11	634494	0.0527302850596
# map.1.txt	79_0	0.014655172413793
	cols = line.replace('\n', '').split('\t')
	if cols[0]==LG:
		Position[cols[0]+"_"+cols[1]]=cc
		MINFO.append([cols[1]])
		cc+=1
		if float(cols[2])>=0.0: # not a LL
			r_list.append(float(cols[2]))
		TotalMarkers+=1
inY.close()

obsA={}
f2plants=[]

for idx1, line in enumerate(srx):

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



states = ('AA','AB','BB')
start_probability = {'AA':0.25,'AB':0.5,'BB':0.25}

transition_probability=[{} for j in xrange(TotalMarkers-1)] # global that is updated within LL(x)
bsum=0
for x1 in xrange(TotalMarkers-1): # recom rates
	r = r_list[x1]  # This is rate
	transition_probability[x1] ={'AA' : {'AA':(1-r)**2.0,'AB':2*r*(1-r),'BB':r**2.0}, 'AB' : {'AA':r*(1-r),'AB':(1-r)**2.0 + r**2.0,'BB':r*(1-r)}, 'BB' : {'AA':r**2.0,'AB':2*r*(1-r),'BB':(1-r)**2.0} }

T_LL=0.0
for j in range(len(f2plants)):
	plantID=f2plants[j] # 
	ER=[Error_Rates[plantID][0],Error_Rates[plantID][1],Error_Rates[plantID][2]]
	llx,fprbs,rprbs=foward_backward(obsA[plantID],states,start_probability,transition_probability,ER)
	postProb=[{} for j in range(len(fprbs))] # forward:: alpha[j][X] is probability that true genotye is X at marker j (starts at 0)

	for j in range(len(fprbs)):
		denom=0.0
		for y in states: 
			denom+=(fprbs[j][y]*rprbs[j][y])
		for y in states: 
			postProb[j][y]=(fprbs[j][y]*rprbs[j][y])/denom

		out1.write(LG+'\t'+plantID+'\t'+str(MINFO[j][0])+'\t'+str(obsA[plantID][j])+'\t'+str(postProb[j]["AA"])+'\t'+str(postProb[j]["AB"])+'\t'+str(postProb[j]["BB"])+'\n')

	T_LL+=llx



print LG,T_LL



