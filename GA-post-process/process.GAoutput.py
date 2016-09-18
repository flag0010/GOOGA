# take output from GA ordering programs and make input file for hmm

# v3 pull order from GA output file

import sys

inputfile=sys.argv[1]
in1 = open(inputfile,"rU")
out1 = open("markers."+inputfile,"w")

bestLL=-9999999999
getnext=0
scorder=''
for line_idx, line in enumerate(in1):
	cols = line.replace('\n', '').split('\t') 
	if len(cols)==4 and getnext==1:
# 			scaffold order= +79 +853 +69 +181 +559 +60 -333 -314 +308 -362 +428 -165 +1184 +240 -26b +563 -254b +232 +427 -80a +1044 +671 +83a +187 +146 -214 -82 +115b +432 -376 -74
		if cols[3][:5]=="scaff":
			#print bestLL,cols
			scorder=cols[3].split("= ")[1]

	elif cols[0][:10]=="individual":
# individual=1 fitness=-3561.00638608
		gg = cols[0].split(' ')
		if len(gg)==2:
			gg1=gg[1].split('=')
			if float(gg1[1])>bestLL:
				if bestLL==-9999999999:
					stLL=float(gg1[1])
				bestLL=float(gg1[1])
				getnext=1
				
			else:
				getnext=0



v1s=[]
orientations=[]

cols = scorder.replace('\n', '').split(' ')
#print scorder
totscaffs=len(cols)
for j in range(totscaffs):
	v1s.append(cols[j][1:])
	if cols[j][0]=="+":
		orientations.append('pos')	
	elif cols[j][0]=="-":
		orientations.append('neg')	
	else:
		print "ugh ",cols[j]


v1marks={}
iny = open("Genotypes.F2.002.txt","rU") # any genotype file from this mapping pop
for line_idx, line in enumerate(iny):
	cols = line.replace('\n', '').split('\t') 	
# imswc700	1	0	NN
	try:
		v1marks[cols[1]].append(cols[2])
	except KeyError:
		v1marks[cols[1]]=[cols[2]]

for j in range(totscaffs):
	chosenv1=v1s[j]
	try:
		zz=v1marks[v1s[j]]
		#print v1s[j],v1marks[v1s[j]]
		if orientations[j]=='pos':
			for k in range(len(zz)):
				out1.write(v1s[j]+'\t'+zz[k]+'\n')
				#print v1s[j],zz[k]
		elif orientations[j]=='neg':
			for k in range(len(zz)):
				out1.write(v1s[j]+'\t'+zz[len(zz)-k-1]+'\n')		

	except KeyError:
		print v1s[j],"missing"

print stLL,bestLL 

