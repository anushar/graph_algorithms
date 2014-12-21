
"""
--function validate_over_represented(g,edgeAttr,sig,correct,isDirected):
identifies the statistically significant edges in a graph
#input arguments:
#1 g - a graph object
#2 edgeAttr- the edge property (i.e. transactions)
#3 sig - the univariate significance level
#4 - correct - boolean to toggle multivariate correction
#5 - isDirected - boolean to toggle directed/no directed graph
"""

def validate_over_represented(g,edgeAttr,sig,correct,isDirected):
	
	#get a list of all weights and find the u
	weights = [e[2][edgeAttr]for e in g.edges(data=True)]
	sumWeights = int(np.sum(weights))
	
	#three dictionaries (i.e., maps, hash-maps, key/val pairs) to store output
	pmfHyper = {} #key: vertex pair, value: probability mass function of the hypergeometric distribution 
	pval = {} #key  vertex pair, value: associated p-value
	validatedDict = {} #key  vertex pair, value: associated p-value of validated links
	
	#perform Bonferroni by dividing by the number of test (number of edges)
	if correct == True:
		multivariateSignificanceCorrection = sig/float(g.number_of_edges())#Bonferroni
	else:
		multivariateSignificanceCorrection = sig
	
	#find the probability of each weight for the hypergeometric null model
	
	#for all edges
	for e in g.edges_iter(data=True):
		source = e[0] #get source/target and weight
		target = e[1]
		weight = e[2][edgeAttr]
		
		#if isDirected get out strength of source and in strength of target
		if isDirected:
			sout = g.out_degree(source,weight=edgeAttr)
			sin = g.in_degree(target,weight=edgeAttr)
		#otherwise get total strength (sum of in and out strengths)
		else:
			sout = g.degree(source,weight=edgeAttr)
			sin = g.degree(target,weight=edgeAttr)

		#find the probability mass function at the value of the observed
		#weight, strengths and sumWeights
		pmfHyper[(source,target)] = hypergeom.pmf(weight,sumWeights ,sout, sin, loc=0)
	
	#now find the p-value
	for e in g.edges_iter(data=True):
		
		#as before
		source = e[0]
		target = e[1]
		weight = e[2][edgeAttr]
		
		#as before
		if isDirected:
			sout = g.out_degree(source,weight=edgeAttr)
			sin = g.in_degree(target,weight=edgeAttr)
		else:
			sout = g.degree(source,weight=edgeAttr)
			sin = g.degree(target,weight=edgeAttr)
		
		
		#find the lower and upper limits for the hypergeometric pmf sum
		lowerSumLim = int(weight)
		upperSumLim = int(sin)
		if sout < sin:
			upperSumLim = int(sout)
		
		#initialize pval for the pair
		pval[(source,target)] = 0
		
		#sum the hypergeometric pmf in the limits found to calculate the p-val
		for X in range(lowerSumLim,upperSumLim+1):
			pval[(source,target)] += hypergeom.pmf(X,sumWeights ,sout, sin, loc=0)
	
	#test agains the correction to form the collection of validated links
	for source,target in pval:
		if pval[(source,target)] < multivariateSignificanceCorrection:
			validatedDict[(source,target)] = pval[(source,target)]
	
	#return the links
	return validatedDict.keys()		
	
