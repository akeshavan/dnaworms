from nipype.utils.filemanip import load_json, save_json
from copy import deepcopy

f = load_json("test2.json")
vstrands = f["vstrands"]
num_v = len(vstrands)

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

s = calculate_shift()
print "shift value = ", s

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


def do_shift(shift,n=1):
    """
    apply shift to dictionary
    """

    # remember to check we aren't going to shift past the scaffold strand
    for i,vs in enumerate(vstrands):
        stapidx = get_strand_idx(vs["stap"])
	scafidx = get_strand_idx(vs["scaf"])
	# TODO: Loop this over n to repeat
	for stidx in stapidx:
	    if stidx+shift > scafidx[-1]:
                raise Exception("Shift value is too high!")
	    vs["stap"][stidx+shift] = add_except_n1(vs["stap"][stidx],shift)

        vs["stap_colors"] = vs["stap_colors"]*2 ## This will change according to n!
    
print "result"
do_shift(s)

save_json("test3.json",f)
