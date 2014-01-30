''' Re-implementation of multiple goal part algorithm from [1] Supplemental 
Listings, with functions for efficiently performing multiple runs,  
validation of assembly plans, and computing the minimum number of stages.

@author: Jonathan Blakes <jonathan.blakes@nottingham.ac.uk>
@copyright (c) 2013-present Jonathan Blakes
@date August 31th 2013
@license: LGPLv3, http://www.gnu.org/licenses/lgpl-3.0.txt

References:

[1]	Densmore, D., Hsiau, T. H.-C., Kittleson, J. T., DeLoache, W., Batten, C. 
	and Anderson, J. C. 2010. Algorithms for automated DNA assembly. 
	Nucleic Acids Research 38, 8 (Mar. 2010), 2607-2616. 
	DOI=http://dx.doi.org/10.1093/nar/gkq165

[2]	Blakes, J., Raz, O., Feige, U., Bacardit, J. Ben-Yehezkel, T., Shapiro, E. 
	and Krasnogor, N. A heuristic for maximizing DNA reuse in synthetic 
	DNA library assembly. ACS Synthetic Biology (submitted 2013).

'''


from collections import Counter
from random import Random
from time import time

hashable = lambda seq: tuple(map(tuple, [(part,) for part in seq]))

import math
def min_stages(goal_part):
	return int(math.ceil(math.log(len(goal_part), 2)))


#function create_asm_graph_sgp( part, hash_mem, [slack], [hash_sharing] )
def create_asm_graph_sgp(part, hash_mem={}):
	''' Supplemental Listing 1: Create Assembly Graph for a Single Goal Part
	
	For use of the optional parameters see function create_asm_graph_sgp2 below.   
	'''
	
	part = tuple(part)
	
	# Memoization Case : memoization hash already has desired part
	if part in hash_mem:		# if ( hash_mem.has_graph(part, [slack]) )
		return hash_mem[part]	# 	return hash_mem.get_graph(part)
	
	# Base Case : part is primitive part
	if len(part) == 1:			# if ( part.length = 1 )
		return asm_graph(part)	# 	return new graph with just the given part
	
	# Recursive Step : iteratively partition part and recurse
	graph_best = None	# graph_best := initialize empty graph
	# for i := 0 to part.length do
	for i in xrange(1, len(part)):	# i must be > 0
		
		# Find best graph for left and right partitions
		subpart_L = part[:i]	# subpart_L := part.subpart ( 0, i )
		subpart_R = part[i:]	# subpart_R := part.subpart ( i, p.length )
		# graph_L := create_asm_graph_sgp( subpart_L, hash_mem, [slack-1], [hash_sharing] )
		graph_L = create_asm_graph_sgp(subpart_L, hash_mem)
		# graph_R := create_asm_graph_sgp( subpart_R, hash_mem, [slack-1], [hash_sharing] )
		graph_R = create_asm_graph_sgp(subpart_R, hash_mem)
		
		# Combine left and right graphs into new graph for intermediate part
		# graph_new := combine_graphs( graph_L, graph_R, [hash_sharing] )
		graph_new = combine_graphs(graph_L, graph_R)
		
		# If cost of new graph is the best so far save the graph
		# graph_best := min_cost( graph_new, graph_best, [slack] )
		graph_best = min_cost(graph_new, graph_best)
		
	# end for
	
	# Add best graph to hash table and return
	hash_mem[part] = graph_best	# hash_mem.insert( part, graph_best )
	return graph_best
# end function


class asm_graph(object):
	
	def __init__(self, part, stages=[]):
		self.part = part
		self._stages = stages
		self.stages = 0
		self.steps = 0
		
		self.sharing = 0
		self.subgraphs = []
		
	def __iter__(self):
		return self.subgraphs.__iter__()

	def reset_cost(self):
		self._steps = self.steps
		self.steps = 0
	
	
# function combine_graphs( graph_L, graph_R ) from Listing 2
def combine_graphs(graph_L, graph_R):
	''' Return graph created from combining two child graphs 
	'''
	# graph_new := merge graph_L and graph_R data structures
	stages = merge_stages(graph_L._stages, graph_R._stages)
	step = (graph_L.part, graph_R.part)
	stages.append([step])
	part = graph_L.part + graph_R.part # == step[0] + step[1]
	graph_new = asm_graph(part, stages)
	
	# Calculate the cost of new graph
	graph_new.stages = max(graph_L.stages, graph_R.stages) + 1
	graph_new.steps = graph_L.steps + graph_R.steps + 1

	# Set subgraphs of graph_new to subgraphs of graph_L and graph_R + itself   
	graph_new.subgraphs.extend(graph_L.subgraphs)
	graph_new.subgraphs.extend(graph_R.subgraphs)
	graph_new.subgraphs.append(graph_new)
	
	return graph_new
# end function


def merge_stages(stages_0, stages_1):
	stages = [set() for _ in range(max(len(stages_0), len(stages_1)))]
	for i, stage in enumerate(stages_0):
		stages[i].update(stage)
	for i, stage in enumerate(stages_1):
		stages[i].update(stage)
	return stages


# function combine_graphs( graph_L, graph_R, hash_sharing ) from Listing 3
def combine_graphs2(graph_L, graph_R, hash_sharing):
	''' Return graph created from combining two child graphs factoring in 
	sharing
	'''
	# Call original combine_graphs function from Listing 2
	graph_new = combine_graphs(graph_L, graph_R)
	
	# Calculate sharing value of new graph
	# graph_new.sharing := graph_L.sharing + graph_R.sharing + hash_sharing.get_sharing( graph_new.part )
	graph_new.sharing = graph_L.sharing + graph_R.sharing + hash_sharing.get(graph_new.part, 0)
	
	return graph_new
# end function 


# function min_cost(graph_new, graph_best) from Listing 2
def min_cost(graph_0, graph_1):
	''' Return graph with minimum cost 
	'''
	if graph_1 is None:
		return graph_0
	
	# Number of stages always take priority
	if graph_0.stages < graph_1.stages: return graph_0
	if graph_1.stages < graph_0.stages: return graph_1

	# If number of stages equal, then graph with less steps is lower cost
	if graph_0.steps < graph_1.steps: return graph_0
	if graph_1.steps < graph_0.steps: return graph_1

	# Graphs have identical cost so arbitrarily choose one
	return graph_0
#end function


# function min_cost(graph_new, graph_best, slack) from Listing 3
def min_cost2(graph_0, graph_1, slack):
	''' Return graph with minimum cost factoring in slack 
	'''
	if graph_1 is None:
		return graph_0
	
	# If either graph has more stages than slack, then call original min_cost function from Listing 2
	if graph_0.stages > slack or graph_1.stages > slack:
		return min_cost(graph_0, graph_1)
	
	# Graph has fewer stages than slack so no need to consider stages
	# Allow sharing to compensate for a graph with more steps
	adjusted_steps_0 = graph_0.steps - graph_0.sharing
	adjusted_steps_1 = graph_1.steps - graph_1.sharing
	if adjusted_steps_0 < adjusted_steps_1: return graph_0
	if adjusted_steps_1 < adjusted_steps_0: return graph_1

	# Graphs have identical cost so arbitrarily choose one
	return graph_0
#end function


#function create_asm_graph_mgp( list_goal_parts, hash_part_library )
def create_asm_graph_mgp(list_goal_parts, hash_part_library={}, extra_slack=0, maxstages=None, hash_sharing=None):
	''' Supplementary Listing 4: Multi-Goal-Part Assembly Algorithm 
	'''
	
	list_goal_parts = map(tuple, list_goal_parts)

	if maxstages is None:	
		# Run single goal part algorithm on each part to find max stages
		maxstages = compute_maxstage(list_goal_parts)
	maxstages += extra_slack

	if hash_sharing is None:
		# Calculate the sharing factors for each possible intermediate part
		hash_sharing = compute_sharing(list_goal_parts)
	
	# Iterate across all goal parts until we have a result graph for each goal part
	list_result_graphs = []
	hash_pinned = {}
	
	while len(list_goal_parts) > 0:
		
		# Reinitialize memoization hash with part library and pinned graphs each with zero cost
		# hash_mem := hash_part_library + hash_pinned
		hash_mem = dict(hash_part_library)
		hash_mem.update(hash_pinned)

		# Call single-goal-part algorithm for each each goal part and determine which graph to pin
		# graph_pinned := initialize empty graph
		graph_pinned = None
		for part in list_goal_parts:
			
			# IR, slack and sharing
			graph_new = create_asm_graph_sgp2(part, hash_mem, slack=maxstages, hash_sharing=hash_sharing)
#			# emulate IR without slack and sharing
#			graph_new = create_asm_graph_sgp2(part, hash_mem, slack=0, hash_sharing={})
			
			# Pin graph with most stages (up to maxstages)
			if graph_pinned is None or (
				graph_new.stages > graph_pinned.stages	# correct as published
#				graph_new.stages < graph_pinned.stages	# prefer shallower graphs  
#				graph_new.steps < graph_pinned.steps	# prefer fewer steps
			):
				graph_pinned = graph_new
				
		# Add pinned graph and graph for each intermediate part to our hash of pinned graphs
		for subgraph in graph_pinned:
			subgraph.reset_cost()
			hash_pinned[subgraph.part] = subgraph
		# end for

		# Remove pinned graph from goal part list and add to result list
		list_goal_parts.remove(graph_pinned.part)
		list_result_graphs.append(graph_pinned)

	# end while
	return list_result_graphs
# end function


def compute_maxstage(list_goal_parts):
#	# Run single goal part algorithm on each part to find max stages
#	maxstages = 0
#	hash_mem = {}
#	for part in list_goal_parts:
#		stages = create_asm_graph_sgp(part, hash_mem).stages
#		if stages > maxstages:
#			maxstages = stages
#	assert max(map(min_stages, list_goal_parts)) == maxstages 
#	return maxstages
	# don't ^ because it can simply be computed:
	return max(map(min_stages, list_goal_parts))

def compute_sharing(list_goal_parts):
	''' according to page 2611, last paragraph of [1]: " 
The sharing factor is the number of redundant times each
sequence of primitive parts appears across all goal parts.
For the goal parts abcde and abcgh, the sharing factor
for all intermediate parts is zero except for the parts ab,
bc and abc. These intermediate parts have a sharing
factor of one, since they appear in both goal parts. ..."
	'''
	counter = Counter()
	for part in list_goal_parts:
		m = len(part) + 1
		for n in xrange(2, m):
			for i in xrange(0, m - n):
				counter[part[i:i+n]] += 1
	hash_sharing = dict((part, count - 1) for part, count in counter.items())
	return hash_sharing


def create_asm_graph_sgp2(part, hash_mem, slack, hash_sharing):
	'''Supplementary Listing 1 (2): Create Assembly Graph for a Single Goal Part 
	'''
	
	part = tuple(part)
	
	# Memoization Case : memoization hash already has desired part
	if part in hash_mem:		# if ( hash_mem.has_graph(part, [slack]) )
		return hash_mem[part]	# 	return hash_mem.get_graph(part)
	
	# Base Case : part is primitive part
	if len(part) == 1:	# if ( part.length = 1 )
		return asm_graph(part)	# return new graph with just the given part
	
	# Recursive Step : iteratively partition part and recurse
	# graph_best = initialize empty graph
	graph_best = None
	# for i := 0 to part.length do
	for i in xrange(1, len(part)):	# i must be > 1
		
		# Find best graph for left and right partitions
		subpart_L = part[:i]	# subpart_L := part.subpart ( 0, i )
		subpart_R = part[i:]	# subpart_R := part.subpart ( i, p.length )
		graph_L = create_asm_graph_sgp2(subpart_L, hash_mem, slack - 1, hash_sharing)
		graph_R = create_asm_graph_sgp2(subpart_R, hash_mem, slack - 1, hash_sharing)
		
		# Combine left and right graphs into new graph for intermediate part
		graph_new = combine_graphs2(graph_L, graph_R, hash_sharing)
		
		# If cost of new graph is the best so far save the graph
		graph_best = min_cost2(graph_new, graph_best, slack)
		
	# end for
	
	# Add best graph to hash table and return
	# hash_mem.insert( part, graph_best )
	hash_mem[part] = graph_best
	return graph_best
# end function


def create_asm_graph_mgp_optimised_for_multiple_runs(list_goal_parts, 
	hash_part_library={}, 
	extra_slack=0,
	runs=1,
	random=Random(),
	shuffle_first=False):

	_start = time()

	list_goal_parts = map(tuple, list_goal_parts) # for efficient compute_sharing

	# unchanging
	maxstages = compute_maxstage(list_goal_parts)
	hash_sharing = compute_sharing(list_goal_parts)

	# to return
	all_list_result_graphs = []
	time_runs = []
	
	for run in xrange(runs):
		
		start = time()
		
		# shuffle input to (possibly) get a different output
		if run > 1 or shuffle_first:
			random.shuffle(list_goal_parts)	# Fisher-Yates O(n)

		list_result_graphs = create_asm_graph_mgp(list_goal_parts, 
			dict(hash_part_library), 
			extra_slack, 
			maxstages, 
			hash_sharing)

		stop = time()
		time_runs.append(stop - start)

		all_list_result_graphs.append(list_result_graphs)

	_stop = time()
	time_total = _stop - _start

	return all_list_result_graphs, time_runs, time_total


def all_list_result_graphs_to_all_stages(all_list_result_graphs):
	# each list_result_graphs is an assembly plan
	# each stages is a list of merged sets of binary concatenations from graphs
	# 	equivalent to output of an Anew run, which can be visualised similarly  
	all_stages = {}
	stages_to_list_result_graphs = {}
	for list_result_graphs in all_list_result_graphs:
		
		# determine quality
		merged = _merge_result_graphs(list_result_graphs)
		quality = (merged.stages, merged.steps)
		num_stages, sum_steps = quality

		# stages are merger of 
		stages = tuple(frozenset(steps) for steps in merged._stages)
		# map stages to list_result_graphs so we can validate plan by assembling
		stages_to_list_result_graphs[stages] = list_result_graphs
		# all_stages is {num_stages, {sum_steps, [stages1, stages2]}} 
		all_stages.setdefault(num_stages, dict()).setdefault(sum_steps, []).append(stages)
		
	return all_stages, stages_to_list_result_graphs


def _merge_result_graphs(list_result_graphs):
	stages = reduce(merge_stages, map(lambda graph: graph._stages, list_result_graphs)) 
	merged = asm_graph(map(lambda graph: graph.part, list_result_graphs), stages) 
	merged.stages = len(stages)
	merged.steps = sum(map(len, stages))
	merged.subgraphs = [subgraph for subgraphs in map(lambda graph: graph.subgraphs, list_result_graphs) for subgraph in subgraphs] # one-level flatten
	return merged


def assemble_list_result_graphs(list_result_graphs, list_goal_parts=None):
	assembled = []
	for graph in list_result_graphs:
		assembled.append(_assemble_result_graph(graph))
	for seq in assembled:
		assert len(seq) == 1	# seq expected to be (sequence,) 1-tuple
	assembled = [seq[0] for seq in assembled]
	if list_goal_parts is not None:
		set_goal_parts = set(map(tuple, list_goal_parts))
		unassembled = set_goal_parts.difference(assembled)
		assert len(unassembled) == 0 # some goal parts remain unassembled 
	return assembled

def _assemble_result_graph(graph):
	seq = hashable(graph.part)
	list_intermediate_subparts = list(seq)
	for steps in graph._stages:
#	for i, steps in enumerate(graph._stages):
#		print 'stage', i+1,
		list_intermediate_subparts = _assemble_result_graph_stage(list_intermediate_subparts, steps)
	return list_intermediate_subparts
	
def _assemble_result_graph_stage(list_intermediate_subparts, steps):
	''' adaptation of pairing_optimised.assemble, 
	WARNING: expects all stage step pairs to universally applicable in any 
		order, so only works for an individual graph and stage
	'''
#	print 'steps:', steps
#	print 'before:', list_intermediate_subparts
	for pair in steps:
		l, r = pair
		joined, empty = l + r, ()
		for i in positions(pair, list_intermediate_subparts):
			list_intermediate_subparts[i:i+2] = joined, empty
	list_intermediate_subparts = [tup for tup in list_intermediate_subparts if len(tup) != 0]
#	print 'after:', list_intermediate_subparts
#	print
	return list_intermediate_subparts

def positions(small, big):
	''' Return start indices of sequence 'small' in longer sequence 'big'
	Adapted from http://stackoverflow.com/a/3847585/97790
	'''
	positions = []
	for i in xrange(len(big)-len(small)+1):
		for j in xrange(len(small)):
			if big[i+j] != small[j]:
				break
		else:
			positions.append(i)
	return positions

