''' pairing is an efficient implementation of algorithm Anew from [2]. 

@author: Jonathan Blakes <jonathan.blakes@nottingham.ac.uk>
@author: Ofir Raz <ofir.raz@weizmann.ac.il>
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


from collections import Counter, defaultdict
from random import Random
from operator import itemgetter
from time import time


def pairing(sequences, runs, 
		random=Random(), 
		sort=True, 							# disable heuristic if False
		sort_most_to_least_frequent=True,	# worse than random if False
		use_counting_sort=False, 			# marginally slower if True but better asymptotically faster
		shuffle_first=False, 			# for single runs (with seed=None)
		conserve_memory=False):				# slightly slower if True
	
	_start = time()
	
	# assert components of targets are hashable and singular (1-tuples)
	if len(sequences) > 0:
		for seq in sequences:
			for part in seq:
				assert isinstance(part, tuple)
				assert len(part) == 1

	# retain reference to initial sequences in order to repopulate between runs
	_sequences = sequences
	
	# find all overlapping pairs and trios
	_pairs, _trios = _pairs_and_trios(sequences, conserve_memory)

	# associate pairs with trios and trios with pairs for exclusion later
	_pairToTrios, _trioToPairs = _setmultimaps(_trios)
	
	# to return
	all_stages = {}
	time_runs = []
	
	for run in xrange(runs):

		start = time()

		# deep copy initial sequences
		sequences = [seq[:] for seq in _sequences]
		
		# reinstate initial variables 
		pairs, trios = _pairs, _trios
		pairToTrios, trioToPairs = _pairToTrios, _trioToPairs

		stages = []		# multi-stage plan to populate in run
		
		# simulate concatenating common pairs at expense of rarer ones
		# solution is the steps simulated at each stage
		while len(pairs) > 0: # done if no pairs in any sequence as all whole
			
			# plan assembly stage
			
			# score pairs
			if conserve_memory:
				scored = [(pair, sum(pairs[pair].values())) for pair in pairs.keys()]
			else:
				scored = [(pair, sum(map(len, pairs[pair].values()))) for pair in pairs.keys()]
	
			# shuffle list to randomize order of equal scoring pairs in stable sort  
			if run > 1 or shuffle_first:
				random.shuffle(scored)	# Fisher-Yates O(n)
			
			# stably sort pairs by most to least frequent
			if sort:
				if use_counting_sort:
					# Counting sort  O(n+r) average and worst; r=range(maxval)
					_inplace_counting_sort_scored(scored, 
						maxval=max(map(itemgetter(1), scored)), 
						reverse=sort_most_to_least_frequent)
				else:
					# Timsort O(n) best, O(n log n) average, O(n log n) worst
					scored.sort(key=itemgetter(1), 
						reverse=sort_most_to_least_frequent)
			
			# add pairs as concatenation steps in sorted order 
			# exclude subsequent pairs that overlap with steps
			steps = [] # in this stage
			excluded = set()
			for pair, _ in scored: # _ is unused score
				
#				# a pair can be excluded by assigning a score <= 0 above 
#				if score <= 0:
#					excluded.add(pair)	# sequences will not be fully assembled
				
				# exclude other pairs that are in a trio with this pair
				if not pair in excluded:
					overlapping = set()
					for trio in pairToTrios[pair]:
						overlapping.update(trioToPairs[trio])
					excluded.update(overlapping)
					
					# pair is next step in stage
					steps.append(pair)
			
			# done if no steps with score > 0
			if len(steps) == 0:
				break
			
			# stage plan complete 
			stages.append(steps)			
			
			# inplace update of sequences simulating assembly stage
			assemble(sequences, steps, pairs, conserve_memory)
			
			# find new overlapping pairs and trios
			pairs, trios = _pairs_and_trios(sequences, conserve_memory)
			
			# compute new multimaps
			pairToTrios, trioToPairs = _setmultimaps(trios)
			
			# loop to plan next stage
		
#		return stages
		stop = time()
		time_runs.append(stop - start)
		
		efficacy = len(stages), sum(map(len, stages))
		num_stages, sum_steps = efficacy
		
		all_stages.setdefault(num_stages, dict()).setdefault(sum_steps, []).append(tuple(frozenset(steps) for steps in stages))
	
	_stop = time()
	time_total = _stop - _start
	
	return all_stages, time_runs, time_total


def _inplace_counting_sort_scored(scored, maxval, reverse=False):
	m = maxval + 1
	count = [0] * m					# init with zeros
	bins = {}
	for pair_score in scored:
		score = pair_score[1]
		bins.setdefault(score, []).append(pair_score)
		count[score] += 1
	indices = range(m)
	if reverse:
		indices = indices[::-1]
	i = 0
	for a in indices:
		if count[a] > 0:
			j = i + count[a]
			scored[i:j] = bins[a]
			i = j
	

def _pairs_and_trios(sequences, conserve_memory=False):
	# if conserve_memory:
	#	only count pairs per sequence
	# else:
	# 	indices used later to retrieve only those sequence containing pair
	# trios used later exclude sets of pairs which overlap pair in step 
	if conserve_memory:
		# only count pairs per sequence
		pairs = defaultdict(Counter)
	else:
		# list left indices of pairs in sequence for counting *and* editing
		pairs = defaultdict(dict)
	trios = set()
	for index, seq in enumerate(sequences):
		penultimate = len(seq) - 2			# left index of rightmost pair
		for i in xrange(0, penultimate + 1):
			pair = tuple(seq[i:i+2])		# two consecutive parts
			if conserve_memory:
				pairs[pair][index] += 1	# increment pair count
			else:
				pairs[pair].setdefault(index, []).append(i)	# store left index of pair
			if i < penultimate:
				trio = tuple(seq[i:i+3])	# three consecutive parts
				trios.add(trio)
	return pairs, trios


def _setmultimaps(trios):
	pairToTrios = defaultdict(set)
	trioToPairs = defaultdict(set)
	for trio in trios:
		lp = trio[:2]	# recreate left and right overlapping pairs of trio
		rp = trio[1:]
		pairToTrios[lp].add(trio)	# use pairs to lookup trio later
		pairToTrios[rp].add(trio)
		trioToPairs[trio].add(lp)	# use trio to lookup pairs later
		trioToPairs[trio].add(rp)
	return pairToTrios, trioToPairs


#@profile
def assemble(sequences, steps, pairs=None, conserve_memory=False):
	''' Concatenate each pair in steps, updating sequences in-place '''
		
	if pairs is None: # when called from assemble_stages (i.e. not pairing)
		pairs, _ = _pairs_and_trios(sequences, conserve_memory)

	# execute planned steps in-place on occurrences in containing sequences
	# move elements of right tuple into left tuple
	# leave right empty so as to not create new pairs by removing too soon 
	changed = set()		# avoid cleaning unchanged later
	empty = () 			# the empty tuple
	for pair in steps:
		l, r = pair
		joined = l + r	# assembly intermediate once
		
		if conserve_memory:
			counter = pairs[pair]	# keys are indices of containing sequences 
			for index in counter.keys():
				seq = sequences[index]								# crucial! 
				mut = list(seq)				# mutable copy
				n = len(mut) - 1	# number of pairs in seq
				h = -2
				i = mut.index(l)	# fast-forward to first pair
				while i < n:
					if mut[i+1] == r and not i == h + 1:
						mut[i:i+2] = joined, empty		# replace pair
						h = i
					try:
						i = mut.index(l, i + 1)	# fast-forward to next pair
					except ValueError:
						i = n
				sequences[index] = mut	# replace sequence
				changed.add(index)		# remember changed
		else:
			for index, indices in pairs[pair].iteritems():
				seq = sequences[index]							# crucial!
				mut = list(seq)				# mutable copy
				h = -2
				for i in indices:
					# avoid consecutive indices: 1,1,1 => (1,1),(),1 => (1,1),(1,1)   
					if not i == h + 1:
						mut[i:i+2] = joined, empty		# replace pair
						h = i
				sequences[index] = mut	# replace sequence
				changed.add(index)		# remember changed
	
	# remove empty tuples (former right pairs)
	for index in changed:
		seq = sequences[index]
		sequences[index] = [tup for tup in seq if len(tup) > 0]

	# in-place so implicitly returns None


def assemble_stages(sequences, stages, conserve_memory=False):
	''' pairing specific '''
	assembled = sequences[:]
	for steps in stages:
		assemble(assembled, steps, None, conserve_memory)
	for seq in assembled:
		try:
			assert len(seq) == 1	# seq expected to be (sequence,) 1-tuple
		except AssertionError:
			unassembled = [seq for seq in assembled if len(seq) > 1]
			for seq in unassembled:
				print 'unassembled', seq
	return map(tuple, [tup[0] for tup in assembled])

