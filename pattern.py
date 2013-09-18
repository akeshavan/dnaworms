from copy import deepcopy
import sys
 
def load_json(filename):
    """Load data from a json file

    Parameters
    ----------
    filename : str
        Filename to load data from.

    Returns
    -------
    data : dict

    """

    fp = file(filename, 'r')
    data = json.load(fp)
    fp.close()
    return data

def save_json(filename, data):
    """Save data to a json file

    Parameters
    ----------
    filename : str
        Filename to save data in.
    data : dict
        Dictionary to save in json file.

    """

    fp = file(filename, 'w')
    json.dump(data, fp, sort_keys=True, indent=4)
    fp.close()
  
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
        #scaf = get_strand_idx(vs["scaf"])
        if len(stap):
            first_stap, last_stap = stap[0], stap[-1]
            staple_length = last_stap - first_stap
            if staple_length % 21:
                shift = staple_length + (21 - staple_length % 21)
            else:
                shift = staple_length + 21
            all_shifts.append(shift)
        else:
            continue
        #if len(scaf):
        #    first_scaf, last_scaf = scaf[0], scaf[-1]
        #find the closest multiple of 21 past staple length 
       
    return max(all_shifts)

def get_max_shifts(shift):
    import numpy as np
    Ns = []
    for i,vs in enumerate(vstrands):
        stap = get_strand_idx(vs["stap"])
        scaf = get_strand_idx(vs["scaf"])
        N = np.floor(np.int(scaf[-1]-stap[-1])/shift)
        Ns.append(N)
    return int(min(Ns))

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
        for stidx in stapidx:
            for n in range(1,N+1):
                if stidx+shift*n > scafidx[-1]:
                    print "Shift value is too high! Stopping at n = ", n
                    break
                else:
                    vs["stap"][stidx+shift*n] = add_except_n1(vs["stap"][stidx],shift*n)
                    
        L = len(vs["stap_colors"])
        for j in range(1,n+1):
            vs["stap_colors"] += [[s[0]+shift*j,s[1]] for s in vs["stap_colors"][:L]]   
      

help = "USAGE: python <filename> <number of repeats> <out filename>"
if len(sys.argv) != 4:
    raise Exception(help)
else:
    print sys.argv
    f = load_json(sys.argv[1])
    vstrands = f["vstrands"]
    num_v = len(vstrands)
    s = calculate_shift()
    N = get_max_shifts(s)
    n = int(sys.argv[2])
    if n>N:
        print "number of copies too high, adjusting to N = ", N
        n = N
    print "shifting by ", s
    do_shift(s,n)
    save_json(sys.argv[3],f)


