#!/usr/bin/env python2.7
''' iGEM submissions by year 2007-2012, of [2] Figures 10 & 11 and Table 4. 

iGEM_submissions reads files of BioBrick parts scraped from partsregistry.org, 
removing duplicates to obtain a set of 4242 unique composite parts.
Those were grouped by #parts (n) and randomly sampled to yield libraries of 
fixed size (k) and n (where n = n_max for all parts).

@author: Jonathan Blakes <jonathan.blakes@nottingham.ac.uk>
@copyright (c) 2013-present Jonathan Blakes
@date August 31th 2013
@license: LGPLv3, http://www.gnu.org/licenses/lgpl-3.0.txt

References:

[2]	Blakes, J., Raz, O., Feige, U., Bacardit, J. Ben-Yehezkel, T., Shapiro, E. 
	and Krasnogor, N. A heuristic for maximizing DNA reuse in synthetic 
	DNA library assembly. ACS Synthetic Biology (submitted 2013).

'''


from collections import defaultdict, OrderedDict 
from itertools import groupby
from random import sample


def groupedby(iterable, key):
	groupedby = OrderedDict()
	for k, group in groupby(sorted(iterable, key=key), key=key):
		groupedby[k] = list(group)
	return groupedby


years = (2007, 2008, 2009, 2010, 2011, 2012) # 2013 only had 2 so far


def iGEM_submissions(years=years):
	''' Returns sequences of sub-parts mapped to composite-parts, by year,
	individually: 
		[{((sub),(parts)): set(part1, part2), # 2007 
		  ((other),(sub),(parts)): set(part3), ...},
	     {...}, # 2008
	     ...]   # and so on; *in input order*
	and cumulatively:
		[2007, 2007+2008, 2007+2008+2009, ...]
		             
	'''
	individual = []
	cumulative = []
	years = map(str, years)
	for year in years:
		dnald = '%s.dnald' % year
		with open(dnald) as f:
			string = f.read()
			multimap = defaultdict(set)
			sequences = [map(lambda s: s.split(' '), line.strip().split(':=')) 
			             for line in string.strip().split('\n')]
			for name, seq in sequences:
				seq = tuple([(i,) for i in seq])
				assert len(name) == 1
				multimap[seq].add(name[0])
			individual.append(dict((seq[:], set(list(names))) for seq, names in multimap.items()))
			cumulative.append(multimap)
			if len(cumulative) > 1:
				this = cumulative[-1]
				prev = cumulative[-2]
				for seq, names in prev.items():
					this[seq].update(names)
	return individual, cumulative


def check():
	individual, cumulative = iGEM_submissions()
	
	print 'years sum_parts unique_parts'
	print '============================'
	
	def len2nd(duo):
		return len(duo[1])
	
	print '----------------- individual'
#	individual_unique = 0
	for i, year in enumerate(years):
		multimap = individual[i]
		sum_names = sum(map(len2nd, multimap.items()))
		print year, sum_names, len(multimap)
#		individual_unique += len(multimap)
#	assert individual_unique == 4353
	
	print '----------------- cumulative'
	for i, year in enumerate(years):
		multimap = cumulative[i]
		sum_names = sum(map(len2nd, multimap.items()))
		if not year == 2007: # is same as individual[0]
			print '2007-%s' % year, sum_names, len(multimap)
	
	all_parts = cumulative[-1].keys()
	cumulative_unique = len(all_parts) 
	assert cumulative_unique == 4242

	grouped = groupedby(all_parts, len)
	
	print
	print 'n k_max'
	print '======='
	for n, parts in grouped.items():
		print n, len(parts)

	for k in range(60,224,60):
		for i in range(10):
			_ = sample(grouped[10], k)


if __name__ == '__main__':
	check()
	
