from copy import deepcopy
import sys
import json 
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

def calculate_shift(vstrands,key="stap"):
    """
    We want to copy staples from the idx of the beginning scaffold strand up to the last staple idx, n times by 21 or until the last idx of the added staple pattern is > the last scaffold idx.
    
    """
    all_shifts = []
    for i,vs in enumerate(vstrands):
        stap = get_strand_idx(vs[key])
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
        if len(scaf) and len(stap):
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

def connector(vstrands,key="scaf"):
    """
    Connect staple strands if there is a break anywhere
    """
    for vs in vstrands:
        idx = [i for i,s in enumerate(vs[key][:-2]) if (s[2]  == -1 and s[3]==-1 and s[1] != -1 and s[0] != -1)\
                and (vs[key][i+1][0]==-1 and  vs[key][i+1][1]==-1 and vs[key][i+1][2]!=-1 and  vs[key][i+1][3]!=-1)]    
        idx2 = [i for i,s in enumerate(vs[key][:-2]) if (s[2]  != -1 and s[3]!=-1 and s[1] == -1 and s[0] == -1)\
                and (vs[key][i+1][0]!=-1 and  vs[key][i+1][1]!=-1 and vs[key][i+1][2]==-1 and  vs[key][i+1][3]==-1)]    
        for id in idx:
            a,b,_,_ = vs[key][id]
            _,_,c,d = vs[key][id+1]
            print vs[key][id], a,b,c,d
            vs[key][id] = [a,b,c,d-1]
            vs[key][id+1] = [a,b+1,c,d]

        for id in idx2:
            a,b,_,_ = vs[key][id+1]
            _,_,c,d = vs[key][id]
            print vs[key][id], a,b,c,d
            vs[key][id] = [a,b-1,c,d]
            vs[key][id+1] = [a,b,c,d+1]


def do_shift(vstrands,shift,N=1,key="stap"):
    """
    apply shift to dictionary
    """

    # remember to check we aren't going to shift past the scaffold strand
    for i,vs in enumerate(vstrands):
        stapidx = get_strand_idx(vs[key])
        #scafidx = get_strand_idx(vs["scaf"])
        if key == "scaf":
            vs["scaf"] += [[-1,-1,-1,-1]]*shift*N # add shift amount of spaces
            vs["stap"] += [[-1,-1,-1,-1]]*shift*N
            vs["loop"] += [0]*shift*N
            vs["skip"] += [0]*shift*N
        for stidx in stapidx:
            for n in range(1,N+1):
                #if stidx+shift*n > scafidx[-1]:
                #    print "Shift value is too high! Stopping at n = ", n
                #    break
                #else:
                vs[key][stidx+shift*n] = add_except_n1(vs[key][stidx],shift*n)
    
        L = len(vs["stap_colors"])
        if key=="stap":
            for j in range(1,n+1):
                vs["%s_colors"%key] += [[s[0]+shift*j,s[1]] for s in vs["%s_colors"%key][:L]]   
      

if __name__ == "__main__":
    help = "USAGE: python <filename> <number of repeats> <out filename> <scaf or stap> <Connect True(1) or False(0)>"
    if len(sys.argv) != 6:
        raise Exception(help)
    else:
        f = load_json(sys.argv[1])
        vstrands = f["vstrands"]
        if sys.argv[-1] == "stap":
            s = calculate_shift(vstrands)
            N = get_max_shifts(s)
            n = int(sys.argv[2])
            if n>N:
                print "number of copies too high, adjusting to N = ", N
                n = N
            print "shifting by ", s
            do_shift(vstrands,s,n)
            if bool(sys.argv[-1]):
                connector(vstrands,"stap")
            save_json(sys.argv[3],f)
        else: 
            s = calculate_shift(vstrands,"scaf")
            do_shift(vstrands,s,int(sys.argv[2]),"scaf")
            connector(vstrands)
            save_json(sys.argv[3],f)

