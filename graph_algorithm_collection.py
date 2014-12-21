import math
import random
import networkx as nx
import collections
import numpy as np
import itertools
import collections
from scipy.stats import hypergeom


def aggressor_out_strength(g,u):
	
	aggressorOutStrength = 0
	for v in g.successors(u):
		aggressorOutStrength += g[u][v]['aggressor_transactions']
	return aggressorOutStrength

def quoter_out_strength(g,u):
	
	quoterOutStrength = 0
	for v in g.successors(u):
		quoterOutStrength += g[u][v]['quoter_transactions']
	return quoterOutStrength
	
def aggressor_in_strength(g,u):
	
	aggressorInStrength = 0
	for v in g.predecessors(u):
		aggressorInStrength += g[v][u]['aggressor_transactions']
	return aggressorInStrength

def quoter_in_strength(g,u):
	
	quoterInStrength = 0
	for v in g.predecessors(u):
		quoterInStrength += g[v][u]['quoter_transactions']
	return quoterInStrength
		
		

def compose_weighted(g,h,weightOperation='add'):
	"""
	composition is the union of node sets and edge sets
	weights for edges that exist in both graphs are added
	"""
	r = g.__class__()
	r.add_nodes_from(g.nodes())
	r.add_nodes_from(h.nodes())
	
	edges = g.edges(data=True)
	edges.extend(h.edges(data=True))
	
	print edges

	for edge in edges:

		u = edge[0]
		v = edge[1]
		weightDict = edge[2]
		#if edge in r.edges():	
		if v in r[u]:
			print "edge found, updating weights with addition"
			for key,val in weightDict.items():
				r[u][v][key] += val
		else:
			print "edge not found"
			r.add_edge(u,v)
			for key,val in weightDict.items():
				r[u][v][key] = val
				
	return r
	#TO DO update node attributes

def compose_weighted_all(graphs,weightOperation='add'):
	"""
	composition is the union of node sets and edge sets
	weights for edges that exist in both graphs are added
	"""
	#TO DO. Add checks that all graphs are the same class
	r = graphs[0].__class__()
	for g in graphs:
		r.add_nodes_from(g.nodes())
	
	edges = []
	for g in graphs:
		e = g.edges(data=True)
		edges.extend(e)
	
	for edge in edges:

		u = edge[0]
		v = edge[1]
		weightDict = edge[2]
		#if edge in r.edges():	
		if v in r[u]:
			#print "edge found, updating weights with addition"
			for key,val in weightDict.items():
				r[u][v][key] += val
		else:
			#print "edge not found"
			r.add_edge(u,v)
			for key,val in weightDict.items():
				r[u][v][key] = val
				
	return r
			
		
	
	#TO DO update node attributes


def degree_averaged_cc(ccs,degrees):
	
	degreeValues = list(set(degrees))
	k_mean_cc = {k:[] for k in degreeValues}


	for i in range(len(degrees)):
		k = degrees[i]
		k_mean_cc[k].append(ccs[i])

	
	for k in degreeValues:
		k_mean_cc[k] = np.mean(k_mean_cc[k])
	
	for k in sorted(k_mean_cc.keys()):
		print k,k_mean_cc[k]
	return k_mean_cc

def kNN(g,source,target,edgeAttr):
	
	#print "source: ", source, "target", target
	#if type(g) == nx.Graph():
		#return {u:np.sum([g.degree(v,weight=w) for v in g.neighbors(u)])/float(g.degree(u,weight=w)) for u in g}
	#if type(g) == nx.DiGraph():
		
	if source == 'in' and target == 'in':

		knn = {u:0 for u in g if g.in_degree(u) > 0}
		for u in knn.keys():
			for v in g.predecessors(u):
				knn[u] += g.in_degree(v,weight=edgeAttr)
			knn[u] /= float(g.in_degree(u,weight=None))
	
	if source == 'in' and target == 'out':
		knn = {u:0 for u in g if g.in_degree(u) > 0}
		for u in knn.keys():
			for v in g.predecessors(u):
				knn[u] += g.out_degree(v,weight=edgeAttr)
			knn[u] /= float(g.in_degree(u,weight=None))
	
	if source == 'out' and target == 'in':
		knn = {u:0 for u in g if g.out_degree(u) > 0}
		for u in knn.keys():
			for v in g.successors(u):
				knn[u] += g.in_degree(v,weight=edgeAttr)
			knn[u] /= float(g.out_degree(u,weight=None))

	if source == 'out' and target == 'out':
		knn = {u:0 for u in g if g.out_degree(u) > 0}
		for u in knn.keys():
			for v in g.successors(u):
				knn[u] += g.out_degree(v,weight=edgeAttr)
			knn[u] /= float(g.out_degree(u,weight=None))

	return knn
	
	

def knn_k(g,source,target,edgeAttr):
	
	if source == 'in':
		knn_k = {k:[] for k in g.in_degree(weight=edgeAttr).values() if k > 0}
	if source == 'out':
		knn_k = {k:[] for k in g.out_degree(weight=edgeAttr).values() if k > 0}
		
	knn = kNN(g,source,target,edgeAttr)
	
	
	if source == 'in':
		for u,k in g.in_degree(weight=edgeAttr).items():
			if k > 0:
				knn_k[k].append(knn[u])
		
	if source == 'out':
		for u,k in g.out_degree(weight=edgeAttr).items():
			if k > 0:
				knn_k[k].append(knn[u])
	
	for k in knn_k:
		knn_k[k] = np.mean(knn_k[k])

	return knn_k

	
def sum_of_weights(g,edgeAttr):

	sumWeights = 0
	
	for u in g:
		for v in g.successors(u):
			sumWeights += g[u][v][edgeAttr] 
	
	sumWeights = float(sumWeights)
	return sumWeights

def to_undirected_weighted(g):

	f = nx.Graph()
	f.add_nodes_from(g.nodes())
	
	for u,v,edgeDict in g.edges(data=True):
		
		if f.has_edge(u,v):
			for item in edgeDict.items():
			#f[u][v][edgeAttr] += edgeDict[edgeAttr]
				f[u][v][item[0]] += item[1]
		else:
			f.add_edge(u,v,edgeDict)
	return f


def degree_distribution(g):
	
	
	if g.is_directed():
		inDegDist = {k:0 for k in g.in_degree().values()}
		outDegDist = {k:0 for k in g.out_degree().values()}
		for u in g:
			inDegDist[g.in_degree(u)] += 1
			outDegDist[g.out_degree(u)] += 1
		return (inDegDist,outDegDist)
	
	degDist = {k:0 for k in g.degree().values()}
	for u in g:
		degDist[g.degree(u)] += 1
	return degDist
		

def normalize_out_strength(g,weightString):
	
	for u in g:
		sOut = float(g.out_degree(u,weight=weightString))
		for v in g.successors(u):
			g[u][v][weightString] /= sOut
	
def max_weight(g,weightString):
	maxWeight = 0.0
	for u in g:
		for v in g.successors(u):
			if g[u][v][weightString] > maxWeight:
				maxWeight = g[u][v][weightString]
	return maxWeight

def sum_weight(g,weightString):
	sumWeight = 0.0
	for u in g:
		for v in g.successors(u):
			sumWeight += g[u][v][weightString]
	return sumWeight
	
def normalize_weights(g,weightString,norm='max'):
	w = max_weight(g,weightString)
	if norm == 'total':
		w = sum_weight(g,weightString)
	for u in g:
		for v in g.successors(u):
			g[u][v][weightString]  /= float(w)
			
def add_edge_from_transaction(gd,transaction):

		quoter = transaction[8]
		aggressor = transaction[10]
		verb = transaction[12]
		vol = float(transaction[5])
		
		gd.add_node(quoter)
		gd.add_node(aggressor)
	
		if verb == "Sell":
			#if edge not already there add to the graph with one transaction		
			if quoter not in gd[aggressor]:
				gd.add_edge(aggressor,quoter,transactions=1,volume=vol)
			else:
				#otherwise increment the transactions/volume
				gd[aggressor][quoter]['transactions'] += 1
				gd[aggressor][quoter]['volume'] += vol
		
		if verb == "Buy":
			if aggressor not in gd[quoter]:
				gd.add_edge(quoter,aggressor,transactions=1,volume=vol)
			else:
				gd[quoter][aggressor]['transactions'] += 1
				gd[quoter][aggressor]['volume'] += vol

#computes a dictionary of borrow preference index values keyed by edge
def preference_index(g,attributeString):
	
	#dict from edge to borrow(lend) preference index
	
	
	#lpi = {e:0 for e in g.edges()}
	#bpi = {e:0 for e in g.edges()}
	
	#bpi=dict.fromkeys(g.edges(),0)
	#lpi=dict.fromkeys(g.edges(),0)
	
	#for all edges
	for e in g.edges():
	
		#get access to the source(0) and target(1) vertices
		u=e[0]
		v=e[1]
		#borrow preference index
		inStrength = g.in_degree(v,weight=attributeString)
		bpi[e] = g[u][v][attributeString] / (1.0 *inStrength)
		
		
		#lender preference index
		outStrength = g.out_degree(u,weight=attributeString)
		lpi[e] = g[u][v][attributeString] / (1.0 *outStrength)
		
	
	#print bpi
	return  (bpi,lpi)
	
	#for all nodes
	"""
	for u in g:
		#get the node in-strength
		inStrength = g.in_degree(u,weight=attributeString)
		
		for v in g.predecessors(u):
			transactions = g[v][u][attributeString]
	"""

#ratio of transactions on the edge to average number of transactions on link
# lpi_i =  W_ij / S_i / K_i = (W_ij*K_i) / S_i
def alt_preference_index(g,attributeString):
	
	#lpi = {e:0 for e in g.edges()}
	#bpi = {e:0 for e in g.edges()}
	
	bpi=dict.fromkeys(g.edges(),0)
	lpi=dict.fromkeys(g.edges(),0)
	
	#for all edges
	for e in g.edges():
	
		#get access to the source(0) and target(1) vertices
		u=e[0]
		v=e[1]
		
		#lender preference index
		si = g.out_degree(u,weight=attributeString)
		ki = g.out_degree(u)
		meanW = si / (1.0 * ki)
		lpi[e] = 1.0*g[u][v][attributeString] / (meanW)
		#print lpi[e]
		
	return  (bpi,lpi)
	
#computes a dictionary of Y2 keyed by vertex 
def Y2(g,attributeString):
	
	#create the dictionary to hold the output
	#Y2Borrow(Lend) can only be defined for vertices that have in(out) degrees respectively
	y2Lend   = {u:0 for u in g if len(g.successors(u)) > 0}
	y2Borrow = {u:0 for u in g if len(g.predecessors(u)) > 0}
	#print y2
	
	for u in g:
		#get the node strength
		outStrength = g.out_degree(u,weight=attributeString)
		
		for v in g.successors(u):
			transactions = g[u][v][attributeString]
			y2Lend[u] += pow(transactions/(1.0*outStrength),2.0)
		
		inStrength = g.in_degree(u,weight=attributeString)
		
		for v in g.predecessors(u):
			transactions = g[v][u][attributeString]
			y2Borrow[u] += pow(transactions/(1.0*inStrength),2.0)
		
	
	#y2Lend = np.mean(y2Lend.values())
	#y2Borrow = np.mean(y2Borrow.values())
	
	lendVals = [ val for val in y2Lend.values() if val > 0]
	borrowVals = [ val for val in y2Borrow.values() if val > 0]
	
	#return the average of Y2
	return (np.mean(lendVals),np.mean(borrowVals))
	
def multigraph_to_weighted(gdNew,gdm,nodeattr,edgeattr):
	
	edgeList = [e for e in gdm.edges(data=True)]
	
	n=gdm.nodes(data=True)
	gdNew.add_nodes_from(n)

	for item in edgeList:
		#print item
		u=item[0]
		v=item[1]

		#gdm.add_node(u,group=gd.node[u][nodeattr])
		#gdm.add_node(v,group=gd.node[v][nodeattr])

		if v in gdNew[u]:
			#print "yes"
			gdNew[u][v][edgeattr] += 1
		else:
			gdNew.add_edge(u,v,transactions=1)
		
def weighted_to_multigraph(gd,gdm,nodeattr,edgeattr):

	edgeList = [e for e in gd.edges(data=True)]
	
	for item in edgeList:
		source = item[0]
		target = item[1]
		gdm.add_node(source,group=gd.node[source][nodeattr])
		gdm.add_node(target,group=gd.node[target][nodeattr])
		numTransactions = item[2][edgeattr]
		for i in range(numTransactions):
			gdm.add_edge(source,target,transactions=1)


def edge_based_directed_multigraph_double_edge_swap(g,nswap=1,ntries=100):

	
	swapcount = 0	
	while swapcount < nswap:
		
		
		edges = list(g.edges_iter())
		e1 = random.choice(edges)
		u = e1[0]
		v = e1[1]
		e2 = random.choice(edges)
		x = e2[0]
		y = e2[1]
		#dont create self loops
		while u == y or x == v:
			e2 = random.choice(edges)
			print u,v,x,y
		
		print u,v,x,y
		g.add_edge(u,y)
		g.add_edge(x,v)
		g.remove_edge(u,v)
		g.remove_edge(x,y)
		swapcount += 1
		print swapcount

def directed_edge_sampling(g,nswap=1):
	
	edges = list(g.edges_iter())
	freq = {e:0 for e in edges}
	numEdges = float(g.number_of_edges())
	
	

	swapcount = 0	
	while swapcount < nswap:
				
		e = random.choice(edges)
		freq[e] += 1	
		swapcount += 1
	
	for e in edges:
		#pEdge = (g.out_degree(e[0])*g.in_degree(e[1]))/numEdges
		pEdge = (g.out_degree(e[0]))/(numEdges*g.in_degree(e[1]))
		print freq[e]/float(nswap),pEdge
		#print g.out_degree(e[0]),g.in_degree(e[1]),numEdges,pEdge
		#print freq[e]/numEdges,(g.out_degree(e[0]))/numEdges


#selects two random edges u-->v and x-->y
#deletes them and creates u-->y and x-->v
def directed_multigraph_double_edge_swap(g,nswap=1,ntries=100):

	#probability weighted by degree
	n=0
	swapcount=0
	

	keys,degrees=zip(*g.out_degree().items())
	cdf=nx.utils.cumulative_distribution(degrees)

	
	while swapcount < nswap:
	
		
		#select two random vertices proportional to their out-degree

		(ui,xi) = nx.utils.discrete_sequence(2,cdistribution=cdf)
		if ui==xi:
			continue
		u=keys[ui]
		x=keys[xi]

		
		#select two source vertices completely at random
		"""
		u=random.choice(g.nodes())
		x=random.choice(g.nodes())
		
		if u == x:
			continue
		if g.out_degree(u) == 0:
			continue
		if g.out_degree(x) == 0:
			continue
		"""
		#choose target uniformly from neighbours
		#v=random.choice(list(g[u]))
		v=random.choice([t[1] for t in list(g.out_edges_iter(u))])
		#y=random.choice(list(g[x]))
		y=random.choice([t[1] for t in list(g.out_edges_iter(x))])
		
		
		if u==y:
			continue #dont create self loops
		if x==v:
			continue #dont create self loops

		g.add_edge(u,y)
		g.add_edge(x,v)
		g.remove_edge(u,v)
		g.remove_edge(x,y)
		swapcount += 1
		
		

#selects two random edges u-->v and x-->y
#deletes them and creates u-->y and x-->v
def directed_double_edge_swap(g,nswap=1,ntries=100):


	#probability weighted by degree

	swapcount=0
	keys,degrees=zip(*g.out_degree().items())
	
	cdf=nx.utils.cumulative_distribution(degrees)
	

	while swapcount < nswap:
	
		(ui,xi) =nx.utils.discrete_sequence(2,cdistribution=cdf)

		if ui==xi: #same source skip, swap with no effect
			continue

		u=keys[ui]
		x=keys[xi]
		
		#choose target uniformly from neighbours
		#un = [u for g.in_degree(u)]
		#yn = [y for g.in_degree(y)]
		
		v=random.choice(list(g[u]))
		y=random.choice(list(g[x]))
		
		
		if v==y:
			continue #same target skip, swap with no effect
		
		if u==y:
			continue #dont create self loops
		if x==v:
			continue #dont create self loops
		if(y not in g[u])  and (v not in g[x]):#dont create parallel edges
			g.add_edge(u,y)
			g.add_edge(x,v)
			g.remove_edge(u,v)
			g.remove_edge(x,y)
			swapcount += 1
		
		

def double_edge_swap(g,nswap=1):


	#probability weighted by degree
	n=0
	swapcount=0
	keys,degrees=zip(*g.degree().items())
	
	cdf=nx.utils.cumulative_distribution(degrees)
	
	while swapcount < nswap:
	
		(ui,xi) =nx.utils.discrete_sequence(2,cdistribution=cdf)
		if ui==xi: #same source, skip
			continue
		u=keys[ui]
		x=keys[xi]
		
		#choose target uniformly from neighbours
		v=random.choice(list(g[u]))
		y=random.choice(list(g[x]))
		
		
		if v==y:
			continue #same target, skip
		if u==y:
			continue #dont create self loops
		if x==v:
			continue #dont create self loops
		if(y not in g[u])  and (v not in g[x]):#dont create parallel edges
			g.add_edge(u,y)
			g.add_edge(x,v)
			g.remove_edge(u,v)
			g.remove_edge(x,y)
			swapcount += 1
			

#selects two random edges u-->v and x-->y
#deletes them and creates u-->y and x-->v
def directed_double_weight_swap(g,nswap=1,ntries=100,edgeAttribute='transactions'):


	#probability weighted by node strength
	n=0
	swapcount=0
	keys,degrees=zip(*g.out_degree(weight=edgeAttribute).items())
	
	cdf=nx.utils.cumulative_distribution(degrees)
	
	while swapcount < nswap:
	
		#print "swapcount %d" %(swapcount)
	
		(ui,xi) =nx.utils.discrete_sequence(2,cdistribution=cdf)
		
		#if ui==xi:
		#	continue #same source skip
		
		u=keys[ui]
		x=keys[xi]
		
		#choose target uniformly from neighbours
		v=random.choice(list(g[u]))
		y=random.choice(list(g[x]))
		
		
		#if v==y:
		#	continue #same target skip
		if u==y:
			continue #dont create self loops
		if x==v:
			continue #dont create self loops
		if(y not in g[u])  and (v not in g[x]):#dont create parallel edges
			
			if (g[u][v]['transactions']== 1) and (g[x][y]['transactions']== 1):
				"""
				print "case 1"
				print "edge 1"
				print u,v,g[u][v]['transactions']
				print "edge 2"
				
				prevKout_u = g.out_degree(u,weight='transactions')
				prevKout_x = g.out_degree(x,weight='transactions')
				prevKin_v = g.in_degree(v,weight='transactions')
				prevKin_y = g.in_degree(y,weight='transactions')
				
				print g.out_degree(u,weight='transactions')
				print g.out_degree(x,weight='transactions')
				print g.in_degree(v,weight='transactions')
				print g.in_degree(y,weight='transactions')
				print x,y,g[x][y]['transactions']
				"""
				g.add_edge(u,y,transactions=1)
				g.add_edge(x,v,transactions=1)
				g.remove_edge(u,v)
				g.remove_edge(x,y)
				swapcount += 1
				"""
				Kout_u = g.out_degree(u,weight='transactions')
				Kout_x = g.out_degree(x,weight='transactions')
				Kin_v = g.in_degree(v,weight='transactions')
				Kin_y = g.in_degree(y,weight='transactions')
				
				if (prevKout_u != Kout_u) or (prevKout_x != Kout_x) or (prevKin_y != Kin_y) or (prevKin_v != Kin_v):
					input("error")
				else:
					print("ok")
				"""
				
				
			elif (g[u][v]['transactions'] > 1) and (g[x][y]['transactions'] > 1):
				"""
				print "case 2"
				print "edge 1"
				print u,v,g[u][v]['transactions']
				print "edge 2"
				print x,y,g[x][y]['transactions']
				
				prevKout_u = g.out_degree(u,weight='transactions')
				prevKout_x = g.out_degree(x,weight='transactions')
				prevKin_v = g.in_degree(v,weight='transactions')
				prevKin_y = g.in_degree(y,weight='transactions')
				
				print u,g.out_degree(u,weight='transactions')
				print x,g.out_degree(x,weight='transactions')
				print v,g.in_degree(v,weight='transactions')
				print y,g.in_degree(y,weight='transactions')
				"""
				g.add_edge(u,y,transactions=1)
				g.add_edge(x,v,transactions=1)
				g[u][v]['transactions'] -= 1
				g[x][y]['transactions'] -= 1
				swapcount += 1
				"""
				Kout_u = g.out_degree(u,weight='transactions')
				Kout_x = g.out_degree(x,weight='transactions')
				Kin_v = g.in_degree(v,weight='transactions')
				Kin_y = g.in_degree(y,weight='transactions')
				
				if (prevKout_u != Kout_u) or (prevKout_x != Kout_x) or (prevKin_y != Kin_y) or (prevKin_v != Kin_v):
					input("error")
				else:
					print("ok")
				"""

			elif (g[u][v]['transactions'] > 1) and (g[x][y]['transactions'] == 1):
				"""
				print "case 3"
				print "edge 1"
				print u,v,g[u][v]['transactions']
				print "edge 2"
				print x,y,g[x][y]['transactions']

				prevKout_u = g.out_degree(u,weight='transactions')
				prevKout_x = g.out_degree(x,weight='transactions')
				prevKin_v = g.in_degree(v,weight='transactions')
				prevKin_y = g.in_degree(y,weight='transactions')

				print u,g.out_degree(u,weight='transactions')
				print x,g.out_degree(x,weight='transactions')
				print v,g.in_degree(v,weight='transactions')
				print y,g.in_degree(y,weight='transactions')
				"""
				g.add_edge(u,y,transactions=1)
				g.add_edge(x,v,transactions=1)
				g[u][v]['transactions'] -= 1
				g.remove_edge(x,y)
				swapcount += 1
				"""
				Kout_u = g.out_degree(u,weight='transactions')
				Kout_x = g.out_degree(x,weight='transactions')
				Kin_v = g.in_degree(v,weight='transactions')
				Kin_y = g.in_degree(y,weight='transactions')
				
				if (prevKout_u != Kout_u) or (prevKout_x != Kout_x) or (prevKin_y != Kin_y) or (prevKin_v != Kin_v):
					input("error")
				else:
					print("ok")
				"""


			elif (g[u][v]['transactions'] == 1) and (g[x][y]['transactions'] > 1):
				"""
				print "case 4"
				print "edge 1"
				print u,v,g[u][v]['transactions']
				print "edge 2"
				print x,y,g[x][y]['transactions']

				prevKout_u = g.out_degree(u,weight='transactions')
				prevKout_x = g.out_degree(x,weight='transactions')
				prevKin_v = g.in_degree(v,weight='transactions')
				prevKin_y = g.in_degree(y,weight='transactions')
				"""
				g.add_edge(u,y,transactions=1)
				g.add_edge(x,v,transactions=1)
				g[x][y]['transactions'] -= 1
				g.remove_edge(u,v)
				swapcount += 1
				"""
				Kout_u = g.out_degree(u,weight='transactions')
				Kout_x = g.out_degree(x,weight='transactions')
				Kin_v = g.in_degree(v,weight='transactions')
				Kin_y = g.in_degree(y,weight='transactions')
				
				if (prevKout_u != Kout_u) or (prevKout_x != Kout_x) or (prevKin_y != Kin_y) or (prevKin_v != Kin_v):
					input("error")
				else:
					print("ok")
				
				"""
				
#unfinished
def directed_double_edge_triangle_swap(g,nswap=1,ntries=100):


	n=0
	swapcount=0
	keys,degrees=zip(*g.out_degree().items())
	
	cdf=nx.utils.cumulative_distribution(degrees)
	
	while swapcount < nswap:
	
		(ui,xi) =nx.utils.discrete_sequence(2,cdistribution=cdf)
		if ui==xi:
			continue
		u=keys[ui]
		x=keys[xi]
		
		v=random.choice(list(g[u]))
		y=random.choice(list(g[x]))
		
		
		if v==y:
			continue
		if u==y:
			continue
		if x==v:
			continue
		if(y not in g[u])  and (v not in g[x]):
			g.add_edge(u,y)
			g.add_edge(x,v)
			g.remove_edge(u,v)
			g.remove_edge(x,y)
			swapcount += 1
			
		#do the triangle swap
		
		#find the middleman triangles
		middlemen=[]
		for i in g:
			for j in g[i]:
				for k in g[j]:
					if i in g[k]:
						middlemen.append((i,j,k))
						#print i,j,k

		#pick a random triangle
		randomTriangle = random.choice(middlemen)
		i=randomTriangle[0]
		j=randomTriangle[1]
		k=randomTriangle[2]
		
		g.remove_edge(i,j)
		g.add_edge(j,i)
		
		g.remove_edge(j,k)
		g.add_edge(k,j)
		
		g.remove_edge(k,i)
		g.add_edge(i,k)
		
def edge_transaction_correlation(g):

	tCorrel = []

	for u in g:
		for v in g[u]:
			trans = g[u][v]['transactions']
			sourceBorrowingStrength = g.out_degree(u,weight='transactions')
			targetLendingStrength   = g.in_degree(v,weight='transactions')
			tCorrel.append((trans*1.0/sourceBorrowingStrength,trans*1.0/targetLendingStrength))

	#unzip the list of tuples and calculate the correlatio between the lists
	lhs,rhs=zip(*tCorrel)
	correl = np.corrcoef(lhs,rhs)[0][1]
	return correl
	
def edge_density(g,directed):
	
	if directed == True:
		return g.number_of_edges() / (1.0*g.number_of_nodes()*(g.number_of_nodes()-1))
	return g.number_of_edges() / (2.0*g.number_of_nodes()*(g.number_of_nodes()-1))

#calculalates how many of all nodes are reachable in a directed graph
#and returns a value between 0 and 1
def reachability_index(g):
	
	pathLengths = nx.all_pairs_shortest_path_length(g)
	reachableNodesSum = 0
	reachableNodesMean = 0
	
	for u in g:
		#calculate how many other nodes are reachable via a directed path excluding self
		reachableNodesSum += len(pathLengths[u])-1
		
	reachableNodesMean = reachableNodesSum
	#average reachable nodes per node
	reachableNodesMean /= (1.0*g.number_of_nodes())
	
	#average reachable nodes per node in units of system size
	reachableNodesMean /= (1.0*g.number_of_nodes())
	
	#find the fraction of pairs that are reachable
	#reachableNodesSum /= (1.0*g.number_of_nodes()*(g.number_of_nodes()-1))
	#return reachableNodesSum
	
	#print reachableNodesMean
	return reachableNodesMean 

#returns true if two graphs with the same set of nodes have Sin(u,g1) =  Sin(u,g2)
#and Sout(u,g1) =  Sout(u,g2) for all u
def is_strength_conserved(g1,g2):

	strengthConserved = True
	
	for u in g2:
				
		kOut_1 = g1.out_degree(u,weight='transactions')
		kIn_1 = g1.in_degree(u,weight='transactions')
			
		kOutS_1 = g2.out_degree(u,weight='transactions')
		kInS_2 = g2.in_degree(u,weight='transactions')
		
		#if either of the in or out strength are not equal
		#strength conservation fails	
		if kOut_1 != kOutS_1 or kIn_1 != kInS_2:
			strengthConserved = False
		
		return strengthConserved

#return the number of parallel edges in a multigraph
def number_of_parallel_edges(g):
	
	numParallel = 0
	for u in g:
		for v in g[u]:
			numParallel += g.number_of_edges(u,v)-1

	return numParallel

def delta_function(arg1,arg2):
	if arg1 == arg2:
		return 1
	return 0
	


def gnm_multigraph(n,e,edgeattr):
	
	g = nx.MultiDiGraph()
	g.add_nodes_from(range(n))
	
	edgeCount = 0
	while edgeCount < e:
		u = random.choice(g.nodes())
		v = random.choice(g.nodes())
		while v == u:
			v = random.choice(g.nodes())
		g.add_edge(u,v)
		g[u][v][edgeattr] = 1
		edgeCount += 1
	return g
	
def reachable_pairs(g):
	
	nActual  = 0
	n = g.number_of_nodes()
		
	if g.is_directed():
		nMax = n*(n-1)
	else:
		nMax = (n*(n-1))/2

	for u in g:
		for v in g:
			if v != u:
				nActual += nx.has_path(g,u,v)
	return nActual / float(nMax)
	
	
	
def parallel_degree(p,vertexToInt):
	
	kp = {i:0 for i in vertexToInt.values()}
	for i in vertexToInt.values():
		for j in vertexToInt.values():
			if j != i:
				kp[i] += p[j][i]*p[i][j]
	return kp	

def k_parallel(g):
	
	"""
	returns the number of parallel edges for each vertex in a dict
	"""
	
	kp = {u:0 for u in g}
	for u in g:
		for v in g.successors(u):
			if u in g.successors(v):
				kp[u] += 1
	return kp

def num_parallel_edges(g,fraction=True):

	"""
	returns the number of pairs (u,v) such that (u,v) in g.edges() and (v,u) in g.edges()
	"""
	eParallel = 0
	for e in list(g.edges_iter()):
		if e[0] in g.successors(e[1]):
			eParallel += 1
	if fraction:
		return eParallel / (2.0* g.number_of_edges())
	
	return eParallel / 2.0 #correct double counting

		
def directed_clustering_coefficient(g,vertexToInt,weighted=False,weightString=None):

	n = g.number_of_nodes()
	from numpy import linalg as LA
	#adjacency matrix
	#sortedVertices = sorted(g.nodes())
	#A = nx.to_numpy_matrix(g,nodelist=sortedVertices,weight=None)
	A = nx.to_numpy_matrix(g)
	#A = nx.to_numpy_matrix(g,nodelist=vertexToInt.keys(),weight=None)
	
	#A = np.zeros((g.number_of_nodes(),g.number_of_nodes()))
	
	
	AT = A.T
	Asum = AT + A
	Asum3 = LA.matrix_power(Asum,3)
	A2 = LA.matrix_power(A,2)
	A3 = LA.matrix_power(A,3)
	AATA =  A*AT*A
	#AATA =  AT*A
	#AATA =  A*AATA


	ATA2 = AT*A2
	A2AT = A2*AT
	#cVector = [i for i in range(n)]
	#cVector = [i for i in vertexToInt.values()]
	#cVector = np.asmatrix(cVector)
	
	"""
	kin = {i:np.dot(AT[i],cVector.T) for i in range(n)}
	kout = {i:np.dot(A[i],cVector.T)for i in range(n)}
	ktot = {i:kin[i] + kout[i] for i in range(n)}
	kparallel = {i:np.dot(Asum[i],cVector.T)for i in range(n)}
	"""
	"""
	kin = {i:np.dot(AT[i],cVector.T) for i in vertexToInt.values()}
	kout = {i:np.dot(A[i],cVector.T)for i in vertexToInt.values()}
	ktot = {i:kin[i] + kout[i] for i in vertexToInt.values()}
	kparallel = {i:np.dot(Asum[i],cVector.T)for i in vertexToInt.values()}
	"""
	"""
	kin = {u:np.dot(AT[vertexToInt[u]],cVector.T) for u in vertexToInt}
	kout = {u:np.dot(A[vertexToInt[u]],cVector.T)for u in vertexToInt}
	ktot = {u:kin[vertexToInt[u]] + kout[vertexToInt[u]] for u in vertexToInt}
	kparallel = {u:np.dot(Asum[vertexToInt[u]],cVector.T)for u in vertexToInt}
	"""
	
	kin = g.in_degree()
	kout = g.out_degree()
	ktot = {u:kin[u] + kout[u] for u in g}
	kparallel = k_parallel(g)
	
	

	
	ccycle = {u:0 for u in g}
	cmiddle = {u:0 for u in g}
	cin = {u:0 for u in g}
	cout = {u:0 for u in g}
	ctot = {u:0 for u in g}

	
	#----------normalize weights and raise to 1/3-----
	if weighted:
	
		normalize_weights(g,weightString)

		for u in g:
			for v in g.successors(u):
				g[u][v][weightString] **= 1/3.0
	#-------------------------------------------------	
		
	#weight matrix
	"""
	W = nx.to_numpy_matrix(g,nodelist=vertexToInt.keys(),weight=weightString)
	WT = W.T
	W2 = LA.matrix_power(W,2)
	W3 = LA.matrix_power(W,3)
	"""
	
	if weighted:		
		MMTM =  W*WT*W
		MTM2 = WT*W2
		M2MT = W2*WT
		M3 = W3
	
	else:
		A2 = LA.matrix_power(A,2)
		M3 = LA.matrix_power(A,3)
		MMTM =  A*AT*A
		MTM2 = AT*A2
		M2MT = A2*AT
	
	#for i in vertexToInt.values():#range(n):
	threshold = 1.0
	for u in g:
		
		"""
		posCycleTriangles = float((kin[i]*kout[i] - kparallel[i]))
		posMiddleTriangles = posCycleTriangles
		posInTriangles = float((kin[i]*(kin[i]-1)))
		posOutTriangles = float((kout[i]*(kout[i]-1)))
		posDirectedTriangles = float((ktot[i]*(ktot[i]-1)) - 2.0*kparallel[i])
		"""
		posCycleTriangles = float((kin[u]*kout[u] - kparallel[u]))
		posMiddleTriangles = float((kin[u]*kout[u] - kparallel[u]))
		posInTriangles = float((kin[u]*(kin[u]-1)))
		posOutTriangles = float((kout[u]*(kout[u]-1)))
		posDirectedTriangles = float((ktot[u]*(ktot[u]-1)) - 2.0*kparallel[u])
		
		#if u =='GB0012':
		"""
		print vertexToInt[u],A3[vertexToInt[u],vertexToInt[u]],posCycleTriangles
		print vertexToInt[u],AATA[vertexToInt[u],vertexToInt[u]],posMiddleTriangles
		print vertexToInt[u],ATA2[vertexToInt[u],vertexToInt[u]],posInTriangles
		print vertexToInt[u],A2AT[vertexToInt[u],vertexToInt[u]], posOutTriangles
		"""
		
		
		if posCycleTriangles >= threshold:
			ccycle[u] = A3[vertexToInt[u],vertexToInt[u]] / posCycleTriangles
			#assert ccycle[u] <= 1

		if posMiddleTriangles >= threshold:
			cmiddle[u] = AATA[vertexToInt[u],vertexToInt[u]] / posMiddleTriangles
			#assert cmiddle[u] <= 1
			
			
		if posInTriangles >= threshold: 
			cin[u] = ATA2[vertexToInt[u],vertexToInt[u]] / posInTriangles
			assert cin[u] <= 1
			
		if posOutTriangles >= threshold: 
			cout[u] = A2AT[vertexToInt[u],vertexToInt[u]] / posOutTriangles
			assert cout[u] <= 1
		
		"""
		if posDirectedTriangles > 0 :
			#ctot[sortedVertices[i]] = Asum3[i,i] / posDirectedTriangles
			ctot[i] = Asum3[i,i] / posDirectedTriangles
		"""
		
		#print sortedVertices[i],ctot[sortedVertices[i]]
		#input()
	
	#return (np.mean(ccycle.values()),np.mean(cmiddle.values()),np.mean(cin.values()),np.mean(cout.values()))
	#return five dictionaries
	#return (ccycle,cmiddle,cin,cout,ctot)
	return (ccycle,cmiddle,cin,cout)
def directed_weighted_clustering(g,weightString):
	
	
	n = g.number_of_nodes()
	from numpy import linalg as LA
	#adjacency matrix
	A = nx.to_numpy_matrix(g,nodelist=g.nodes(),weight=None)
	A2 = LA.matrix_power(A,2)
	AT = A.T
	Asum = AT + A
	cVector = [i for i in range(n)]
	cVector = np.asmatrix(cVector)
	
	kin = {i:np.dot(AT[i],cVector.T) for i in range(n)}
	kout = {i:np.dot(A[i],cVector.T)for i in range(n)}
	kparallel = {i:np.dot(Asum[i],cVector.T)for i in range(n)}

	#print "kin"
	#print kin
	#weight matrix
	W = nx.to_numpy_matrix(g,nodelist=g.nodes(),weight=weightString)
	WT = W.T
	W2 = LA.matrix_power(W,2)
	W3 = LA.matrix_power(W,3)
			
	WWTW =  W*WT*W
	WTW2 = WT*W2
	W2WT = W2*WT
	
	ccycle = {i:0 for i in range(n)}
	cmiddle = {i:0 for i in range(n)}
	cin = {i:0 for i in range(n)}
	cout = {i:0 for i in range(n)}

	for i in range(n):
			
			if kin[i]*kout[i]  - kparallel[i] > 0:
				ccycle[i] = W3[i,i] / float((kin[i]*kout[i] - kparallel[i]))
				cmiddle[i] = WWTW[i,i] / float((kin[i]*kout[i] - kparallel[i]))
			if kin[i] > 1: 
				cin[i] = WTW2[i,i] / float((kin[i]*(kin[i]-1)))
			if kout[i] > 1: 
				cout[i] = W2WT[i,i] / float((kout[i]*(kout[i]-1))) 
	#print type((np.mean(ccycle.values()),np.mean(cmiddle.values()),np.mean(cin.values()),np.mean(cout.values())))
	#print "here"
	#input()
	#return (np.mean(ccycle.values()),np.mean(cmiddle.values()),np.mean(cin.values()),np.mean(cout.values()))
	return (ccycle,cmiddle,cin,cout)


def directed_weighted_clustering_coefficient(g,vertexToInt,edgeAttr):

	n = g.number_of_nodes()
	from numpy import linalg as LA
	#----------normalize weights and raise to 1/3-----
	
	#norm = 'total'
	#normalize_weights(g,edgeAttr,norm)
	#norm = float(sum_of_weights(g,edgeAttr))
	norm = float(max_weight(g,edgeAttr))
	#norm = 1.0
	#print "norm: ", norm
	
	for u in g:
		for v in g.successors(u):
			#print g[u][v][edgeAttr]
			#input()
			g[u][v][edgeAttr] /= norm
			g[u][v][edgeAttr] **= (1/3.0)
			


	#-------------------------------------------------	
		
	#weight matrix

	W = nx.to_numpy_matrix(g,weight=edgeAttr)
	WT = W.T
	W2 = LA.matrix_power(W,2)
	W3 = LA.matrix_power(W,3)
	WWTW =  W*WT*W
	WTW2 = WT*W2
	W2WT = W2*WT
	W3 = W3
	
	kin = g.in_degree()
	kout = g.out_degree()
	ktot = {u:kin[u] + kout[u] for u in g}
	kparallel = k_parallel(g)
	

	ccycle = {u:0 for u in g}
	cmiddle = {u:0 for u in g}
	cin = {u:0 for u in g}
	cout = {u:0 for u in g}
	ctot = {u:0 for u in g}
	
	
	threshold = 0.0
	for u in g:
		
		posCycleTriangles = float((kin[u]*kout[u] - kparallel[u]))
		posMiddleTriangles = float((kin[u]*kout[u] - kparallel[u]))
		posInTriangles = float((kin[u]*(kin[u]-1)))
		posOutTriangles = float((kout[u]*(kout[u]-1)))
		posDirectedTriangles = float((ktot[u]*(ktot[u]-1)) - 2.0*kparallel[u])
		

		if posCycleTriangles > threshold:
			ccycle[u] = W3[vertexToInt[u],vertexToInt[u]] / posCycleTriangles
			assert ccycle[u] <= 1

		if posMiddleTriangles > threshold:
			cmiddle[u] = WWTW[vertexToInt[u],vertexToInt[u]] / posMiddleTriangles
			assert cmiddle[u] <= 1
			
			
		if posInTriangles > threshold: 
			cin[u] = WTW2[vertexToInt[u],vertexToInt[u]] / posInTriangles
			assert cin[u] <= 1
			
		if posOutTriangles > threshold: 
			cout[u] = W2WT[vertexToInt[u],vertexToInt[u]] / posOutTriangles
			assert cout[u] <= 1
		
		"""
		if posDirectedTriangles > 0 :
			#ctot[sortedVertices[i]] = Asum3[i,i] / posDirectedTriangles
			ctot[i] = Asum3[i,i] / posDirectedTriangles
		"""
	
	#print np.mean(ccycle.values()),np.mean(cmiddle.values()),np.mean(cin.values()),np.mean(cout.values())
	#input()
	return (ccycle,cmiddle,cin,cout)


	
def directed_clustering_coefficient_slow(g,vertexToInt):

	#make adjacency matrix
	p={i:{j:0.0 for j in vertexToInt.values()} for i in vertexToInt.values()}
	for u in g:
		for v in g:
			if v in g.successors(u):
				p[vertexToInt[u]][vertexToInt[v]] = 1.0

	
	kin={}
	kout={}
	for u in g:
		kin[vertexToInt[u]] = g.in_degree(u)
		kout[vertexToInt[u]] = g.out_degree(u)
	
	kp = parallel_degree(p,vertexToInt)


	cCycleAlt = {i:0 for i in vertexToInt.values()}
	cMiddleAlt = {i:0 for i in vertexToInt.values()}
	cInAlt = {i:0 for i in vertexToInt.values()}
	cOutAlt = {i:0 for i in vertexToInt.values()}
	
	for i in vertexToInt.values():
		for j in vertexToInt.values():
				for k in vertexToInt.values():
					if j != i and k != i and k != j:
						cInAlt[i] += (p[j][k]*p[j][i]*p[k][i]) 
						cOutAlt[i] += (p[i][k]*p[i][j]*p[j][k]) 
						cCycleAlt[i] += (p[i][j]*p[j][k]*p[k][i])
						cMiddleAlt[i] += (p[i][k]*p[j][k]*p[j][i])
		"""
		if i == 7:
			print 'bad bank'
			print i,cCycleAlt[i],kin[i]*kout[i] - kp[i]
			print i,cMiddleAlt[i],kin[i]*kout[i] - kp[i]
			print i,cInAlt[i],kin[i]*(kin[i]-1)
			print i,cOutAlt[i],kout[i]*(kout[i]-1)
			input()
		"""

		if (kin[i]*(kin[i]-1)) > 0.0: #should be one really but correcting from fp error, i.e 1.000000001
			cInAlt[i] /= kin[i]*(kin[i]-1)
		else:
			cInAlt[i] = 0
		#if cInAlt[i] > 1:
			#print "error: ",kin[i]
			
		if (kout[i]*(kout[i]-1)) > 0.0:
			cOutAlt[i] /= kout[i]*(kout[i]-1)
		else:
			cOutAlt[i] = 0
			
		if (kin[i]*kout[i] - kp[i]) > 0:
			cCycleAlt[i] /= (kin[i]*kout[i] - kp[i])
			cMiddleAlt[i] /= (kin[i]*kout[i] - kp[i])
		else: 
			cCycleAlt[i] = 0
			cMiddleAlt[i] = 0

		


	
	return (cCycleAlt,cMiddleAlt,cInAlt,cOutAlt)
	
	
def make_vertex_to_int(g):

	vertexToInt={}
	intIndex = 0
	for u in g:
		vertexToInt[u] = intIndex
		intIndex += 1
	return vertexToInt
	
	
def mobility_impact(g):

	kIn = g.in_degree(g)
	kOut = g.out_degree(g)
	kMean = np.mean(kIn.values())
	kMeanSquared = kMean*kMean
	kInMax = max(kIn.values())
	kOutMax = max(kOut.values())
	
	kInSquared = [kIn[u]*kIn[u] for u in kIn]
	kOutSquared = [kOut[u]*kOut[u] for u in kOut]
	
	kInSquaredMean = np.mean(kInSquared)
	kOutSquaredMean = np.mean(kOutSquared)
	
	val = (1/kMean) +(2/kMeanSquared)*(kInMax*kOutSquaredMean + kOutMax*kInSquaredMean)
	return (g.number_of_nodes(),val)

def swap_mobility(g):
		
	from numpy import linalg as LA
	
	A = nx.to_numpy_matrix(g)	
	AT = A.transpose()
	kSum= np.sum(g.degree().values())**2
	kProdInOut_ii = 0
	kProdInOut_ij = 0
	for u in g:
		kProdInOut_ii += g.in_degree(u)*g.out_degree(u)
		for v in g.successors(u):
			kProdInOut_ii += g.out_degree(u)*g.in_degree(v)
			
		
	
	swapMobility = 0.5*np.trace(A*AT*A*AT)  - kProdInOut_ij + np.trace(A*AT*A) \
	+ 0.5*kSum + 0.5*np.trace(LA.matrix_power(A,2)) - kProdInOut_ii
	
	return float(swapMobility)
	
	
def cycle_mobility(g):
	
	from numpy import linalg as LA
	
	A = nx.to_numpy_matrix(g)
	A2 = LA.matrix_power(A,2)
	A3 = LA.matrix_power(A,3)
	
	
	Ap = A
	for i in range(g.number_of_nodes()):
		for j in range(g.number_of_nodes()):
			Ap[i,j] = 0
			Ap[i,j] = A[i,j]*A[j,i]
	
	Ap3 = LA.matrix_power(Ap,3)	
	cycleMobility = 0.5*A3.trace() - np.trace(Ap*A2) + np.trace(Ap*Ap*A2) - 0.33333333333 *np.trace(Ap3)
	
	return float(cycleMobility)


def get_cycle_triangles(g,parallel):
	"""
	represent a cycle triangle as a list of 2-tuples of length three
	return a list of such lists
	v such that v = source(e1) and v = target(e2) or
	v = source(e2) and v = target(e1), i.e. the configurations
	u---->v---->w or w---->v---->u
	"""
	candidateNodes = [u for u in g if g.in_degree(u) > 0 and g.out_degree(u) > 0]
	#cycles = {u:0 for u in candidateNodes}
	#print g.number_of_nodes(),len(candidateNodes)
	triangles = []
	for u in candidateNodes:
		outSet = g.successors(u)
		inSet = g.predecessors(u)
		for v in outSet:
			for w in inSet:
				if g.has_edge(v,w): #cycle triangle found
					#cycles[u] += 1
					if parallel:
						triangles.append([(v,w),(u,v),(w,u)])
					else:
						if (v not in g[w]) and (u not in g[v]) and (w not in g[u]):
							triangles.append([(v,w),(u,v),(w,u)])
	return triangles				
	#return (cycles,triangles)
		
def cycle_edge_reversion(g,parallel=False):
	
	#get the list of cycle triangles, select a random triangle and access its edges

	cycles = get_cycle_triangles(g,parallel)
	cycle = random.choice(cycles)

	for edge in cycle:
		source = edge[0]
		target = edge[1]
		g.remove_edge(source,target)
		g.add_edge(target,source)
	


def filter_graph_by_weight_percentile(g,per,edgeAttr):


	edgesToRemove = []

	weights = [e[2][edgeAttr] for e in g.edges(data=True)]
	percentileCuttoff = np.percentile(weights,per)

	for e in g.edges(data=True):
		if e[2][edgeAttr] < percentileCuttoff:
			edgesToRemove.append(e)
	g.remove_edges_from(edgesToRemove)
	
def write_pajek(g,edgeAttr,gName,isDirected=True):
	
	f = open(gName, 'wt')
	f.write('*Vertices')
	f.write(' %s' %(g.number_of_nodes()))
	f.write('\n')
	#write vertices
	#map vertices to integers
	intToVertex={}
	i = 1
	for u in g:
		f.write('%s'%i)
		f.write(' \"'+u+'\"')
		f.write('\n')
		intToVertex[u] = i
		i += 1



	#write edges
	if isDirected:
		f.write('*edges')
	else:
		f.write('*Arcs')
	f.write(' %s' %(g.number_of_edges()))
	for e in g.edges():
		f.write('\n')
		u = e[0]
		v = e[1]
	#	print u,v
		f.write('%s %s %s'%(intToVertex[u],intToVertex[v],g[u][v][edgeAttr]))
	f.close()

#---------statistical--link--validation
def validate_over_represented(g,edgeAttr,sig,correcting,isDirected):

	weights = [e[2][edgeAttr]for e in g.edges(data=True)]
	sumWeights = int(np.sum(weights))
	
	pmfHyper = {}
	pval = {}
	validatedDict = {}
	
	if correcting == True:
		multivariateSignificanceCorrection = sig/float(g.number_of_edges())#Bonferroni
	else:
		multivariateSignificanceCorrection = sig

	#find the probability of each weight for the hypergeometric null model
	for e in g.edges_iter(data=True):
		source = e[0]
		target = e[1]
		weight = e[2][edgeAttr]
		if isDirected:
			sout = g.out_degree(source,weight=edgeAttr)
			sin = g.in_degree(target,weight=edgeAttr)
		else:
			sout = g.degree(source,weight=edgeAttr)
			sin = g.degree(target,weight=edgeAttr)
		pmfHyper[(source,target)] = hypergeom.pmf(weight,sumWeights ,sout, sin, loc=0)
	
	#now find the p-value
	for e in g.edges_iter(data=True):
		
		source = e[0]
		target = e[1]
		weight = e[2][edgeAttr]
		
		#print source,target,weight
		
		pmfHyper[(source,target)] = hypergeom.pmf(weight,sumWeights ,sout, sin, loc=0)
		if isDirected:
			sout = g.out_degree(source,weight=edgeAttr)
			sin = g.in_degree(target,weight=edgeAttr)
		else:
			sout = g.degree(source,weight=edgeAttr)
			sin = g.degree(target,weight=edgeAttr)
		weight = e[2][edgeAttr]
		lowerSumLim = int(weight)
		upperSumLim = int(sin)
		if sout < sin:
			upperSumLim = int(sout)
		pval[(source,target)] = 0
		
		for X in range(lowerSumLim,upperSumLim+1):
			pval[(source,target)] += hypergeom.pmf(X,sumWeights ,sout, sin, loc=0)
		
	
	#now apply the statistical correction for performing num edges tests
	for source,target in pval:
		if pval[(source,target)] < multivariateSignificanceCorrection:
			validatedDict[(source,target)] = pval[(source,target)]
	return validatedDict.keys()		

#---------statistical--link--validation
def validate_over_represented_fast(g,edgeAttr,sig,correct,isDirected):

	weights = [e[2][edgeAttr] for e in g.edges(data=True)]
	sumWeights = int(np.sum(weights))
	
	pmfHyper = {}
	pval = {}
	validatedDict = {}
	
	if correct == True:
		multivariateSignificanceCorrection = sig/float(g.number_of_edges())#Bonferroni
	else:
		multivariateSignificanceCorrection = sig

	#find the probability of each weight for the hypergeometric null model
	for source,nbrsdict in g.adjacency_iter():
		for target,keydict in nbrsdict.iteritems():
			for key,eattr in keydict.iteritems():
				
				if key == edgeAttr:
					
					if isDirected:
						sout = g.out_degree(source,weight=edgeAttr)
						sin = g.in_degree(target,weight=edgeAttr)
					else:
						sout = g.degree(source,weight=edgeAttr)
						sin = g.degree(target,weight=edgeAttr)
				
					pmfHyper[(source,target)] = hypergeom.pmf(eattr,sumWeights ,sout, sin, loc=0)
					

	#now find the p-value
	for source,nbrsdict in g.adjacency_iter():
		for target,keydict in nbrsdict.iteritems():
			for key,eattr in keydict.iteritems():
			
				if key == edgeAttr:
					
					if isDirected:
						sout = g.out_degree(source,weight=edgeAttr)
						sin = g.in_degree(target,weight=edgeAttr)
					else:
						sout = g.degree(source,weight=edgeAttr)
						sin = g.degree(target,weight=edgeAttr)

					
					lowerSumLim = int(eattr)
					upperSumLim = int(sin)
					if sout < sin:
						upperSumLim = int(sout)
				
					pval[(source,target)] = 0
		
					for X in range(lowerSumLim,upperSumLim+1):
						pval[(source,target)] += hypergeom.pmf(X,sumWeights ,sout, sin, loc=0)
		
	#now apply the statistical correction for performing num edges tests
	for source,target in pval:
		if pval[(source,target)] < multivariateSignificanceCorrection: #reject the null hypothesis if True
			validatedDict[(source,target)] = pval[(source,target)]
	return validatedDict.keys()		

#---------statistical--link--validation
def validate_over_under_represented_fast(g,edgeAttr,sig,correct,isDirected):

	weights = [e[2][edgeAttr] for e in g.edges(data=True)]
	sumWeights = int(np.sum(weights))
	
	pmfHyper = {}
	pvalOver = {}
	pvalUnder = {}
	validatedOver = {}
	validatedUnder = {}
	e = g.number_of_edges()
	nB = len([u for u in g if g.in_degree(u) > 0])
	nL = len([u for u in g if g.out_degree(u) > 0])
	nLB = len([u for u in g if g.out_degree(u) > 0 and g.in_degree(u) > 0])
	T = e + nB*nL - nLB
	
	print "e: " + str(e)
	print "nB: " + str(nB)
	print "nL: " + str(nL)
	print "nLB: " + str(nLB)
	print "T: " + str(T)

	if correct == True:
		multivariateSignificanceCorrection = sig/float(T)#Bonferroni
	else:
		multivariateSignificanceCorrection = sig

	#find the probability of each weight for the hypergeometric null model
	for source,nbrsdict in g.adjacency_iter():
		for target,keydict in nbrsdict.iteritems():
			for key,eattr in keydict.iteritems():
				
				if key == edgeAttr:
					
					if isDirected:
						sout = g.out_degree(source,weight=edgeAttr)
						sin = g.in_degree(target,weight=edgeAttr)
					else:
						sout = g.degree(source,weight=edgeAttr)
						sin = g.degree(target,weight=edgeAttr)
				
					pmfHyper[(source,target)] = hypergeom.pmf(eattr,sumWeights ,sout, sin, loc=0)
					

	#now find the p-value
	for source,nbrsdict in g.adjacency_iter():
		for target,keydict in nbrsdict.iteritems():
			for key,eattr in keydict.iteritems():
			
				if key == edgeAttr:
					
					if isDirected:
						sout = g.out_degree(source,weight=edgeAttr)
						sin = g.in_degree(target,weight=edgeAttr)
					else:
						sout = g.degree(source,weight=edgeAttr)
						sin = g.degree(target,weight=edgeAttr)

					#do the over validation
					lowerSumLim = int(eattr)
					upperSumLim = int(sin)
					if sout < sin:
						upperSumLim = int(sout)
				
					pvalOver[(source,target)] = 0
		
					for X in range(lowerSumLim,upperSumLim+1):
						pvalOver[(source,target)] += hypergeom.pmf(X,sumWeights ,sout, sin, loc=0)
					
					#do the under validation
					lowerSumLim = 0
					upperSumLim = int(eattr)
				
					pvalUnder[(source,target)] = 0
		
					for X in range(lowerSumLim,upperSumLim+1):
						pvalUnder[(source,target)] += hypergeom.pmf(X,sumWeights ,sout, sin, loc=0)
		
	#now apply the statistical correction for performing num edges tests
	for source,target in pvalOver:
		if pvalOver[(source,target)] < multivariateSignificanceCorrection: #reject the null hypothesis if True
			validatedOver[(source,target)] = pvalOver[(source,target)]
		if pvalUnder[(source,target)] < multivariateSignificanceCorrection: #reject the null hypothesis if True
			validatedUnder[(source,target)] = pvalUnder[(source,target)]
	
	#sorted_pvalOver = sorted(pvalOver.iteritems(), key=operator.itemgetter(1))

	return (validatedOver.keys(),validatedUnder.keys())


def write_pajek(g,edgeAttr,gName,isDirected=True):
	
	f = open(gName, 'wt')
	f.write('*Vertices')
	f.write(' %s' %(g.number_of_nodes()))
	f.write('\n')
	#write vertices
	#map vertices to integers
	intToVertex={}
	i = 1
	for u in g:
		f.write('%s'%i)
		f.write(' \"'+u+'\"')
		f.write('\n')
		intToVertex[u] = i
		i += 1
	#write edges
	if isDirected:
		f.write('*edges')
	else:
		f.write('*Arcs')
	f.write(' %s' %(g.number_of_edges()))
	for e in g.edges():
		f.write('\n')
		u = e[0]
		v = e[1]
	#	print u,v
		f.write('%s %s %s'%(intToVertex[u],intToVertex[v],g[u][v][edgeAttr]))
	f.close()
	

def out_degree_distribution(g,weight=None,cumulative=False):

	assert nx.is_directed(g)
	out_degrees = g.out_degree(weight=weight).values()
	n = float(g.number_of_nodes())
	kmax = np.max(out_degrees)+1
	kmin = 0
	dd = {k:out_degrees.count(k)/n for k in range(kmin,kmax)}
	
	if cumulative:
		cum = dict.fromkeys(dd,0)
		for k in cum:
			cum[k] = np.sum([dd[K] for K in range(k)])
		return cum
		
	return dd
	
def digraph_completeness(g):
	assert nx.is_directed(g)
	n = float(g.number_of_nodes())
	return (g.number_of_edges / (n*(n-1.0)))

def digraph_reciprocity(g,completenessCorrected=True):
	
	numParallel = 0
	for u in g:
		for v in g:
			numParallel += nx.has_edge(u,v)*nx.has_edge(v,u)
	rho = numParallel / float(g.number_of_edges())

	if completenessCorrected:
		alpha = digraph_completeness(g)
		return (rho - alpha) / (1 - alpha)
	return rho

