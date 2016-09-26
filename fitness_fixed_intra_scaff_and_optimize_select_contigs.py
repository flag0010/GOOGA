# likelihoods based on window calls as input
# all plants for a scaff together in all*txt
# reports likelihood of observed states in forward_backward

#v3.1:: homogenenous recombination within v1scaffs
#v3.1b:: work on exponential scale for rates

# fir = fixed intra-scaff rates

import numpy
from scipy import optimize
from scipy.optimize import minimize_scalar
from scipy.special import gammaln
from math import exp,log

gradval=1e-8
prci=10000000.0 # low values mean higher precision default : factr=10000000.0

def fitness(filename, R_rates, UGaps,
			lines_file = "test.f2group.txt",
			error_rates = "error.rates.txt",
			intra_scaff_rates_file="intrascaff_rates.txt"):
	#filename points to the intial contig order.  the test example was the file markers.2.LM.txt
	#R_rates is the total list of all estimated inter-scaffold rec. rates, include those that need to be changed and those that don't
	#need to rename the list of F2s for each population as test.f2.group.txt	
	#UGaps is a list of 2 items, first another list with the gap positions that stay the same (note the position is zero-based indexing)
	#and second a int which gives the total number of recombination parameters (which is 1 minus the number of scaffolds plus 1 for
	#the global intra scaffold parameter, which always comes first)  i.e. N + 1 - 1  simple!  :)
	#for example, if R_rates = [0.1, 0.25, 0.3]
	#and if UGaps = [[1], 3]
	#the program would assume the first (0.1) and the last (0.3) rates are to be estimated, and won't attempt to optimize the middle rate (0.25).
	#and the program will ignore whatever value you put in the first and last position (i.e. 0.1 and 0.3 are irrelevant)
	#print filename, R_rates, UGaps
	
	def calc_v0(r_rates): 
		def scipy_ln_like0(x):
			return -LL(x)
		def final_scipy_ln_like0(x):
			##this function is a fudge to get out the really real -ln_lk without using the version above that has been
			#polluted by mucking around with the global params R_rates and UGaps
			#not exactly pretty, but it works
			return -LL_without_global_tweaking(x)
		bounds = [] # CHANGE HERE
		for k in range(0,len(r_rates)):
			bounds.append( (0.00000001,0.1) )
		#print len(bounds), len(r_rates)
		best, val, d = optimize.fmin_l_bfgs_b(scipy_ln_like0, r_rates, factr=prci, epsilon=gradval, approx_grad=True, bounds=bounds)
		#print d
		solution = list(best)
		ln_l = -final_scipy_ln_like0(solution)
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
			#pass in rates we will be changing and/or those we will leave.  
			#could use a global to tell which rates change and which we leave 
			#once i wrap this whole thing in a big old function, global could come in there
		position, tot_len = UGaps
		xnew = [-999 for i in xrange(tot_len)]
		for r, p in zip(R_rates, position): xnew[p] = r
		for p, r in enumerate(x):
			if p not in position: xnew[p] = r
		#print x
		#print len(x), len(xnew)
		x = xnew
		#print x
		#ok, the rates in x are now patched with the available carry overs
		transition_probability=[{} for j in xrange(TotalMarkers-1)] # global that is updated within LL(x)
		bsum=0
		for x1 in xrange(TotalMarkers-1): # recom rates
			if distX[x1]==-9: # we have hit a new v1 scaff 
				r = x[bsum]  # This is inter-scaff rate for this position: bsum = number of v1 scaffs into LG 
				bsum+=1  # CHANGE HERE
			else:
				r = distX[x1] # CHANGE HERE
	
			transition_probability[x1] ={'AA' : {'AA':(1-r)**2.0,'AB':2*r*(1-r),'BB':r**2.0}, 'AB' : {'AA':r*(1-r),'AB':(1-r)**2.0 + r**2.0,'BB':r*(1-r)}, 'BB' : {'AA':r**2.0,'AB':2*r*(1-r),'BB':(1-r)**2.0} }
	
		Total_LL=0.0
		for j in range(len(f2plants)):
			plantID=f2plants[j] # 		
			ER=[Error_Rates[plantID][0],Error_Rates[plantID][1],Error_Rates[plantID][2]]
			llx=foward_backward(obsA[plantID],states,start_probability,transition_probability,ER)
			Total_LL+=llx
#		Calls_to_LL+=1
		return Total_LL
	
	def LL_without_global_tweaking(x):
		transition_probability=[{} for j in xrange(TotalMarkers-1)] # global that is updated within LL(x)
		bsum=0
		for x1 in xrange(TotalMarkers-1): # recom rates
			if distX[x1]==-9: # we have hit a new v1 scaff 
				r = x[bsum]  # This is inter-scaff rate for this position: bsum = number of v1 scaffs into LG 
				bsum+=1  # CHANGE HERE
			else:
				r = distX[x1] # CHANGE HERE
	
			transition_probability[x1] ={'AA' : {'AA':(1-r)**2.0,'AB':2*r*(1-r),'BB':r**2.0}, 'AB' : {'AA':r*(1-r),'AB':(1-r)**2.0 + r**2.0,'BB':r*(1-r)}, 'BB' : {'AA':r**2.0,'AB':2*r*(1-r),'BB':(1-r)**2.0} }
	
		Total_LL=0.0
		for j in range(len(f2plants)):
			plantID=f2plants[j] # 		
			ER=[Error_Rates[plantID][0],Error_Rates[plantID][1],Error_Rates[plantID][2]]
			llx=foward_backward(obsA[plantID],states,start_probability,transition_probability,ER)
			Total_LL+=llx
#		Calls_to_LL+=1
		return Total_LL
	####################################################################################################
	### Main Program
	
#	Calls_to_LL=0 # global for testing
	
	states = ('AA','AB','BB')
	start_probability = {'AA':0.25,'AB':0.5,'BB':0.25}
	
	Error_Rates={}
	inZ = open(error_rates,"rU")
	for line_idx, line in enumerate(inZ):
		cols = line.replace('\n', '').split(' ') 
		Error_Rates[cols[0]]=[float(cols[2]),float(cols[3]),float(cols[4])]
	
	## NEW STUFF
	forward_rate={}
	reverse_rate={}
	iISR=open(intra_scaff_rates_file,"rU")
	for line_idx, line in enumerate(iISR):
		cols = line.replace('\n', '').split('\t') 
		#print cols
		forward_rate[cols[0]+"_"+cols[1]]=float(cols[2])
		reverse_rate[cols[0]+"_"+cols[1]]=float(cols[3])
	
	Position={}
	inY = open(filename,"rU") 
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
			if cols[0] != cscaff: # this marker starts a new v1 scaffold
				Breaks.append(line_idx)
				distX.append(-9) # this is code for program to estimate a free recombination rate
			else: # this marker is within a v1 scaffold
				if int(cols[1]) > lastpos: # forward orientation of scaff
					distX.append( reverse_rate[cols[0]+"_"+cols[1]] ) # recombination rate 
				else:
					distX.append( forward_rate[cols[0]+"_"+cols[1]] ) 
			cscaff=cols[0]
			lastpos=int(cols[1])
		Position[cols[0]+"_"+cols[1]]=line_idx
		mxz.append(cols[0]+"_"+cols[1])
		TotalMarkers+=1
	## CHANGE HERE
	interScaff_intervals=len(Breaks)
	#R_rates=[0.01 for k in range(interScaff_intervals)]  # Initial recombination rate per bp (morgans / megabase) within scaffs
	
#	print "No recomb parameters ",len(R_rates)
#	print "Recomb parameters ",R_rates
	
	obsA={}
	f2plants=[]
	srx  =open(lines_file, "rU")
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
	
	zsol= calc_v0(R_rates) 
	return zsol
	#print zsol
#	print Calls_to_LL


