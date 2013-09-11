from nipype.utils.filemanip import load_json, save_json
from copy import deepcopy
import sys

def get_strand_idx(stap):
    """
    stap is a list of staples from vstrands[idx]["stap"]
    """
    i = [idx for idx,s in enumerate(stap) if not s==[-1,-1,-1,-1] ]
    return i

def calculate_shift():
    """
    We want to copy staples from the idx of the beginning scaffold strand up to the last staple idx, n times by 21 or until the last idx of the added staple pattern is > the last scaffold idx.
    
    """
    all_shifts = []
    for i,vs in enumerate(vstrands):
        stap = get_strand_idx(vs["stap"])
        scaf = get_strand_idx(vs["scaf"])
        first_stap, last_stap = stap[0], stap[-1]
        first_scaf, last_scaf = scaf[0], scaf[-1]
        #find the closest multiple of 21 past staple length 
        staple_length = last_stap - first_stap
        if staple_length % 21:
            shift = staple_length + (21 - staple_length % 21)
        else:
            shift = staple_length
        all_shifts.append(shift)

    return max(all_shifts)

def add_except_n1(l,val):
    """
    Add val to the 1'st and 3rd number in the list except if there is a -1
    """
    L = deepcopy(l)

    if L[1] != -1:
        L[1] += val
    if L[3] != -1:
        L[3] += val

    return L

def do_shift(shift,N=1):
    """
    apply shift to dictionary
    """

    # remember to check we aren't going to shift past the scaffold strand
    for i,vs in enumerate(vstrands):
        stapidx = get_strand_idx(vs["stap"])
        scafidx = get_strand_idx(vs["scaf"])
        for n in range(1,N+1):
            for stidx in stapidx:
                if stidx+shift*n > scafidx[-1]:
                    raise Exception("Shift value is too high!")
                vs["stap"][stidx+shift*n] = add_except_n1(vs["stap"][stidx],shift*n)

        vs["stap_colors"] = vs["stap_colors"]*(N+1) 


help = "USAGE: python <filename> <number of repeats> <out filename>"
if len(sys.argv) != 4:
	raise Exception(help)
else:
	print sys.argv
	f = load_json(sys.argv[1])
	vstrands = f["vstrands"]
	num_v = len(vstrands)
	s = calculate_shift()
	print "shifting by ", s
	do_shift(s,int(sys.argv[2]))
	save_json(sys.argv[3],f)


