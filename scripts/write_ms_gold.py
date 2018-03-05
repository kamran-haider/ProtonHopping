"""
A quick fix to generate ms_gold file, will 
"""
with open("../data/gramicidin/prep_structures/ms_gold", "w") as f:
	for i in xrange(1, 529):
    	f.write("HOHW%04d\n" % i)
