from pattern import *
import numpy as np
import sys

"""
We need to define a set of staples and then fill the scaffold strand
"""


"""Complimentary bases dict"""
B = {"A":"T","T":"A","C":"G","G":"C"}

def find_num(a,col):
    # doesn't go in order, need to navigate to a["vstrands"][num]
    return [idx for idx,x in enumerate(a["vstrands"]) if x["num"]==col][0]

def stap_length(a,row,col,idx=0):
    """
    Calculate staple length given a vstrand index(vstart) and index
    """

    row_new, col_new = a["vstrands"][row]["stap"][col][2], a["vstrands"][row]["stap"][col][3]
    if row_new == -1 and col_new == -1:
        return idx+1
    else:
        idx +=1
        return stap_length(a,find_num(a,row_new),col_new,idx)

"""
Let's fill in the complimentary to the staples on the scaffold, 
unless the scaffold isn't there at that spot
"""

def fill_scaf(a,stap_row,stap_col,stap_seq,idx=0):
    scaf = a["vstrands"][stap_row]["scaf"][stap_col]
    row_new, col_new = a["vstrands"][stap_row]["stap"][stap_col][2], a["vstrands"][stap_row]["stap"][stap_col][3]
    if not scaf == [-1,-1,-1,-1]:
        a["vstrands"][stap_row]["scaf_seq"][stap_col] = B[stap_seq[idx]]
    idx +=1
    if row_new == -1 and col_new == -1:
        return None
    return fill_scaf(a,find_num(a,row_new), col_new, stap_seq, idx)

def find_beg(a):

    return [[i,j] for i, strand in enumerate(a["vstrands"]) for j, row in enumerate(strand["scaf"]) if row[:2] == [-1,-1] and not row[2:]==[-1,-1]]

def walk_scaf(a,col,row,seq=[]):
    seq.append(a["vstrands"][col]["scaf_seq"][row])
    col_new, row_new = a["vstrands"][col]["scaf"][row][2:]
    if col_new ==-1 and row_new == -1:
        return seq
    return walk_scaf(a,find_num(a,col_new),row_new,seq)

def scaf_seq(a):
    beg_col, beg_row = find_beg(a)[0]
    seq  = walk_scaf(a,beg_col,beg_row)
    print ''.join(seq)
    return seq

if __name__ == "__main__":
    helpmsg = """
USAGE: python sequence.py <input json file> <output json file>
"""

    if not len(sys.argv) == 3:
        raise Exception(helpmsg)
    
    a = load_json(sys.argv[1])
    #Initialize new keys
    for i, strand in enumerate(a["vstrands"]):
        strand["stap_seq"]=[]
        strand["scaf_seq"] = ['?']*len(strand["scaf"])

    # Generate staples and fill scaffold
    for i, strand in enumerate(a["vstrands"]):
        for j, col in enumerate(strand["stap_colors"]):
            L = stap_length(a,i,col[0])
            strand["stap_seq"] +=[[col[0], ''.join([B.keys()[x] for x in np.random.random_integers(0,3,L)])]]
            fill_scaf(a,i,col[0],strand["stap_seq"][-1][-1])

    #print total scaffold
    scaf_seq(a)
    save_json(sys.argv[2],a)














