#!/usr/bin/env python2.7
''' Comparison of algorithms Anew [2] and A1 [1]; also acting as usage example.

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
from operator import itemgetter
from random import Random
from collections import defaultdict

from Anew import pairing, assemble_stages
from A1 import (min_stages, create_asm_graph_mgp_optimised_for_multiple_runs,
	all_list_result_graphs_to_all_stages, assemble_list_result_graphs)
from datasets import (example,
	small_artificial, phagemid, iGEM_2008,
	unott_azurin, ukb, uh1, uh2, ethz,
	remove_duplicates_retaining_order) 

import numpy as np


def compare(libraries, runs_Anew=10, runs_A1=10, seed=1, 
		write_dot=False, ext='pdf', view=False, viewer='evince'):
	''' Runs Anew and A1 on each library, checks correctness of all assembly 
	plans and prints efficacy and efficiency results for each set of runs. 
	Returns {library: (algorithm, plan)} dict, e.g.: {
		'uh1': ('A1', A1_list_result_graphs_best),
		'uh2': ('Anew', stages_best_Anew)}
	'''
	
	# library header
	print 'library k n_max stages_min' 
	# results header
	print 'algorithm runs', '\t', 
	print 'best mode worst', 
	print '(time_total time_once time_runs_sum time_runs_mean time_runs_std)'
	print '================================================================='
	
	best_plans = {} # to return
	
	for func in libraries:
		
		library = str(func).split()[1]
		targets = remove_duplicates_retaining_order(func())
		k = len(targets)
		n_max = max(map(len, targets))
		stages_min = max(map(min_stages, targets))#max(map(lambda target: int(math.ceil(math.log(len(target), 2))), targets))

		print library, k, n_max, stages_min


		# Anew

		print 'Anew', runs_Anew, '\t',

		all_stages_Anew, time_runs_Anew, time_total_Anew = \
			pairing(targets, 
				runs=runs_Anew, 
				random=Random(seed),
#				sort=False, 						# disable heuristic if False
#				sort_most_to_least_frequent=False,	# worse than random if False
#				use_counting_sort=True, 			# marginally slower if True but better asymptotically faster
#				shuffle_first=True,			# for single runs (with seed=None) 
#				conserve_memory=True,				# slightly slower if True
				shuffle_first=True)
				
		# check all assembly plans 
		all_stages_Anew_list = []
		for _, d in all_stages_Anew.iteritems():
			for list_stages in d.values():
				all_stages_Anew_list.extend(list_stages)
		# remove 1-tuple wrappers needed by Anew and make target hashable
		_targets = [tuple([part[0] for part in target]) for target in targets]
		for stages in all_stages_Anew_list:
			# assemble_stages raises AssertionError if assembly fails
			assembled = assemble_stages(targets, stages, conserve_memory=False)
			# double check
			assert set(assembled) == set(_targets)

		best_Anew, mode_Anew, worst_Anew = compute_efficacies(all_stages_Anew)
		print best_Anew, mode_Anew, worst_Anew, 
		print compute_times(time_runs_Anew, time_total_Anew) 


		# A1
		print 'A1  ', runs_A1, '\t',

		# remove 1-tuple wrappers needed by Anew 
		list_goal_parts = [[part[0] for part in target] for target in targets]

		all_list_result_graphs, time_runs_A1, time_total_A1 = \
			create_asm_graph_mgp_optimised_for_multiple_runs(
				list_goal_parts, 
				hash_part_library={}, 
				extra_slack=0, 
				runs=runs_A1, 
				random=Random(seed),
				shuffle_first=True)
	
		all_stages_A1, A1_stages_to_list_result_graphs = \
			all_list_result_graphs_to_all_stages(all_list_result_graphs)

		# check all assembly plans 
		for list_result_graphs in A1_stages_to_list_result_graphs.values():
			# assemble_stages raises AssertionError if assembly fails
			assembled = assemble_list_result_graphs(list_result_graphs, list_goal_parts)
			# double check	
			assert set(assembled) == set(map(tuple, list_goal_parts)) 

		best_A1, mode_A1, worst_A1 = compute_efficacies(all_stages_A1)
		print best_A1, mode_A1, worst_A1, 
		print compute_times(time_runs_A1, time_total_A1) 


		# choose best plan

		# get representative best stages 
		stages_best_Anew = all_stages_Anew[best_Anew[0]][best_Anew[1]][0]
		# and corresponding A1_list_result_graphs
		stages_best_A1 = all_stages_A1[best_A1[0]][best_A1[1]][0]
		A1_list_result_graphs_best = A1_stages_to_list_result_graphs[stages_best_A1]
		
		if best_A1 < best_Anew:
			best_plans[library] = ('A1', A1_list_result_graphs_best)
		else:
			best_plans[library] = ('Anew', stages_best_A1) # easier to assemble
		
		
		# vizualise best assembly plans of each using stages 
		# (requires GraphViz) 
		if write_dot:
			import subprocess
			for algorithm in 'A1', 'Anew':
				if algorithm == 'A1':
					runs = runs_A1
					best = best_A1
					dot = stages_dot(stages_best_A1)
				elif algorithm == 'Anew':
					runs = runs_Anew
					best = best_Anew
					dot = stages_dot(stages_best_Anew)
				best_stages, best_steps = best
				dotf = '%s-%s-%s-%s_%s.dot' % (library, algorithm, runs, best_stages, best_steps)
				imgf = '%s.%s' % (dotf, ext)
				cmd_graphviz = 'dot -T%s -o %s %s' % (ext, imgf, dotf)
				graphviz = subprocess.Popen(cmd_graphviz, shell=True)
				with open(dotf, 'w') as f:
					f.writelines(dot)
				if view and graphviz.wait() == 0:
					subprocess.Popen('%s %s' % (viewer, imgf), shell=True)
	
		print '-----------------------------------------------------------------'
	
	print
	
	return best_plans



# node numbering
def posints():
	i = 0
	while True:
		i += 1
		yield i


def stages_dot(stages, ids=False, tooltips=False, id2part=None, part2id=None, id2step=None, step2id=None):
	if tooltips and id2part is None:
		pass #TODO vtooltips
	
	v_ids = posints()
	vs = defaultdict(v_ids.next)
	vtooltips = {} 
	es = '''
	
	'''
	ranks = defaultdict(set)
	targets = set()
	for stage, steps in enumerate(stages):
		steps = tuple(steps)
		for step in steps:
			l, r = step
			i = vs[l]
			j = vs[r]
			try:
				targets.remove(l)
			except:
				pass
			try:
				targets.remove(r)
			except:
				pass
			p = vs[step]
			t = tuple(l + r)
			targets.add(t)
			n = vs[t]
			es += '%s->%s [penwidth=2.0];\n' % (i, p)
			es += '%s->%s;\n' % (j, p)
			es += '%s->%s [weight=2.0];\n' % (p, n)
			ranks[stage].add(t)

			vtooltips[i] = '(%s)' % ', '.join(['\\"%s\\"' % _ for _ in l])
#			print i, vtooltips[i], l
			vtooltips[j] = '(%s)' % ', '.join(['\\"%s\\"' % _ for _ in r])
#			print j, vtooltips[j], r
			vtooltips[p] = '%s + %s' % (vtooltips[i], vtooltips[j])
#			print p, vtooltips[p], step
			vtooltips[n] = '(%s)' % ', '.join(['\\"%s\\"' % _ for _ in t])
#			print n, vtooltips[n], t

			
	inputs = []
	vs_ = []
	ts = []
	ss = {}
	for v, i in sorted(vs.items(), key=itemgetter(1)):
		if isinstance(v[0], tuple):
			if ids:
				s = '%s [shape=circle, width=0.1, style=filled, fillcolor="#000000"];' % i
			else:
				s = '%s [shape=point, fillcolor="#000000"];' % i
			if tooltips:
				s = s.replace('];',', tooltip="%s"];' % vtooltips[i])
			vs_.append(s)
		else:
			s = '%s [shape=square, width=0.4];' % i
			if len(v) == 1:
				s = s.replace('];',', style=filled, fillcolor="#0000FF"];')
				if tooltips:
					s = s.replace('];',', tooltip="%s"];' % vtooltips[i])
				inputs.append(s)
			else:
				if v in targets:
					s = s.replace('];',', style=filled, fillcolor="green3"];')
					if tooltips:
						s = s.replace('];',', tooltip="%s"];' % vtooltips[i])
					ts.append(s)
				else:
					if tooltips:
						s = s.replace('];',', tooltip="%s"];' % vtooltips[i])
					vs_.append(s)
		ss[v] = s
		
	ranked = sorted(ranks.items())
	auxs = []
	for _, vs in ranked[:-1]:
		auxs_ = []
		for v in vs:
			if v not in targets:
				auxs_.append(ss[v].replace('];',', style=filled, fillcolor="red3"];'))
		auxs.append(auxs_)
		
	header = 'digraph {'

	if ids:
		header += 'node [label="", color=none, fillcolor="#000000EE"];'
	else:
		header += 'node [label="", fillcolor="#000000EE"];'  
		
	return header + '''

// dot-specific attributes
//size="7.75,10.25";	
//orientation="landscape";
//ratio="compress";
ranksep="0.8 equally";
nodesep=0.15;
edge [arrowhead=none, color="#00000022"];

{rank=min;
%s
}

%s

%s

{rank=max;
%s
}

%s
}
''' % (
	'\n'.join(inputs), '\n'.join(vs_), '\n\n'.join(
		['{rank=same;\n%s\n}' % r for r in [
			'\n'.join(bs) for bs in auxs]]), '\n'.join(ts), es)



def compute_efficacies(all_stages):
	''' determine best, mode and worst efficacy '''
	
	_efficacies_counter = Counter()
	for num_stages, d in all_stages.items():
		for sum_steps, stages in d.items():
			_efficacies_counter.update([(num_stages, sum_steps)] * len(stages))
	_efficacies_ascending = sorted(_efficacies_counter.elements())
	
	# best and worst are now easy
	best = _efficacies_ascending[0]
	worst = _efficacies_ascending[-1]
	
	# mode is a bit more complication 
	_by_mode_descending = sorted(_efficacies_counter.most_common(), 
		key=itemgetter(1), reverse=True)
	# [((4, 213), 2), ((4, 209), 2), ((4, 208), 1)]
	_mode_count = _by_mode_descending[0][1] 
	# 2
	_has_mode_count = filter(lambda e: e[1] == _mode_count, _by_mode_descending)
	# [((4, 213), 2), ((4, 209), 2)]
	_best_of_mode_count = sorted(_has_mode_count, key=itemgetter(0))
	# [((4, 209), 2)]
	mode = _best_of_mode_count[0][0]
	# (4, 209)
	
	return best, mode, worst


def compute_times(time_runs, time_total):
	time_runs_sum = np.sum(time_runs)
	time_once = time_total - time_runs_sum 
	time_runs_mean = np.mean(time_runs)
	time_runs_std = np.std(time_runs)
	# to microsecond accuracy
	time_total = round(time_total, 6)
	time_once = round(time_once, 6)
	time_runs_sum = round(time_runs_sum, 6)
	time_runs_mean = round(time_runs_mean, 6)
	time_runs_std = round(time_runs_std, 6)
	return time_total, time_once, time_runs_sum, time_runs_mean, time_runs_std


def fast_comparison():
	''' 200 runs of Anew and 2 runs of A1 on phagemid dataset used in [1] & [2] 
	'''
	compare([phagemid], runs_Anew=200, runs_A1=2, seed=1)


def full_comparison():
	''' 1 and then 1000 runs each of Anew and A1 on all datasets used in [2] '''
	libraries = [ # by A1 running time ascending
		example,
		small_artificial,
		phagemid,
		ukb,
		unott_azurin,
		uh1,
		iGEM_2008,
		ethz,
		uh2,
	]
	compare(libraries, 1, 1)		# to estimate running time 
	compare(libraries, 1000, 1000) 	# to validate results in [2]


if __name__ == '__main__':
#	fast_comparison()
	full_comparison()
#	compare([ukb], 1, 1, write_dot=True, view=True)
