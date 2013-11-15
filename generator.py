import numpy as np
from pattern import load_json, save_json
import networkx as nx
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt


class NoodleBase(object):
    def __init__(self, length=84):
        self._length = length+1
        self.helices = {}
        self.graph = None
        self.Xs = {}
        self.Ks = {}
        self.stap = {}
        self.scaf_Xs = {}
        self.scaf = {}

        self._init_helices()
        self._init_crossovers()
        self._init_kinks()
        self._init_scaffold()

    def _init_helices(self):
        for i in range(6):
            self.helices[i] = -np.ones(self._length)
        self._helix_graph()
        for u, v, d in self.graph.edges(data=True):
            num = (self._length-d["start_idx"])/21
            idx = np.asarray([d["start_idx"]+21*n for n in range(num+1)])
            self.helices[u][idx] = v
            self.helices[v][idx] = u
            self.helices[u][idx[np.nonzero(idx)]-1] = v
            self.helices[v][idx[np.nonzero(idx)]-1] = u

    def _helix_graph(self):
        G=nx.Graph()
        G.add_nodes_from(range(6))
        G.node[0]["direction"]=3
        G.node[2]["direction"]=3
        G.node[4]["direction"]=3
        G.node[1]["direction"]=5
        G.node[3]["direction"]=5
        G.node[5]["direction"]=5
        G.add_edge(0,1,label="0 + 21n",start_idx=0)
        G.add_edge(1,2,label="14 + 21n",start_idx=14)
        G.add_edge(2,3,label="7 + 21n",start_idx=7)
        G.add_edge(3,4,label="0 + 21n",start_idx=0)
        G.add_edge(4,5,label="14 +21n",start_idx=14)
        G.add_edge(0,5,label="7 + 21n",start_idx=7)
        self.graph = G

    def draw_helix_graph(self):
        nx.draw_circular(self.graph, node_color="lightyellow", node_size=1000)
        pos = nx.circular_layout(self.graph)
        edge_labels = dict([((u,v,),d['label'])
                            for u,v,d in self.graph.edges(data=True)])
        nx.draw_networkx_edge_labels(self.graph, pos, edge_labels=edge_labels)
        plt.title("Crossover Locations for Staples")

    def _init_crossovers(self):
        for i in range(6):
            self.Xs[i] = np.zeros(self._length).astype(bool)

    def validate_crossovers(self):
        for helix, data in self.helices.iteritems():
            if -1 in data[self.Xs[helix]]:
                raise Exception("Invalid!")
            #TODO: add in part for reciprocity
        return True

    def _init_kinks(self):
        for i in range(6):
            self.Ks[i] = np.zeros(self._length).astype(bool)

    def validate_kinks(self):
        for helix, data in self.Xs.iteritems():
            if -1 in data[self.Ks[helix]]:
                raise Exception("Invalid!")
            kidx = np.nonzero(self.Ks[helix])[0]
            kdist = np.diff(kidx)
            if (kdist<6).any():
                raise(Exception("Kinks are too close together!"))
        return True

    def _init_cadnano(self,helix, direction, length):
        foo = []
        for i in range(self._length):
            if direction == 5:
                foo.append([helix,i-1,helix,i+1])
            elif direction == 3:
                foo.append([helix,i+1,helix,i-1])
            if i==0:
                if direction==5: foo[i] = [-1,-1,helix,i+1]
                elif direction==3: foo[i] = [helix,i+1,-1,-1]
            elif i==length-2:
                if direction==5: foo[i] = [helix,i-1,-1,-1]
                elif direction==3: foo[i] = [-1,-1,helix,i-1]
        return foo

    def _draw_Xovers_and_kinks(self,cad,direction,X,K,data,helix):
        Xidx = np.nonzero(X)[0]
        Kidx = np.nonzero(K)[0]

        for idx in Xidx:
            if direction==5 and not idx%7: cad[idx] = [int(data[idx]),idx,helix,idx+1];
            elif direction==5 and idx%7: cad[idx] = [helix,idx-1,int(data[idx]),idx];

            elif direction==3 and idx%7: cad[idx] = [int(data[idx]),idx,helix,idx-1];
            elif direction==3 and not idx%7: cad[idx] = [helix,idx+1,int(data[idx]),idx];

        for idx in Kidx:
            if direction==5: cad[idx] = [helix,idx-1,-1,-1]
            elif direction==3: cad[idx] = [-1,-1,helix,idx-1]

        return cad

    def _stap_cadnano(self):

        for helix, d in self.graph.nodes(data=True):
            cad = self._init_cadnano(helix,d['direction'],self._length)
            X = self.Xs[helix][:-1]
            K = self.Ks[helix][:-1]
            data = self.helices[helix][:-1]
            self.stap[helix] = self._draw_Xovers_and_kinks(cad,d['direction'],X,K,data,helix)

    def plot_staples(self):
        plt.figure(figsize=(12,4))
        plt.hold(True)
        for i, data in self.helices.iteritems():
            plt.plot([0,self._length],[i,i],color="b")
            idx = np.nonzero(data !=-1)[0]
            y = data[idx]
            plt.plot(idx,y,'o',color="grey",alpha=0.5)

            for ii,d in enumerate(y):
                plt.annotate('%d'%d,[idx[ii],i])

            Xdata = self.Xs[i]
            Xidx = np.nonzero(Xdata)[0] # <-- list of indices where crossovers happen!
            for Xi in Xidx:
                plt.plot([Xi,Xi],[i,data[Xi]])

            Kdata = self.Ks[i]
            Kidx = np.nonzero(Kdata)[0]
            Ky = np.ones(len(Kidx))*i
            plt.plot(Kidx,Ky,'x',color="black")

        plt.axis([-2,self._length,-1,6])
        plt.yticks(range(6))
        plt.xticks(range(0,self._length,7))

    def plot_scaffold(self):
        plt.figure(figsize=(12,4))
        plt.hold(True)
        for i,data in self.scaf_helices.iteritems():
            plt.plot([0, self._length],[i, i], color="b")
            idx = np.nonzero(data != -1)[0]
            y = data[idx]
            plt.plot(idx, y, 'o', color="grey", alpha=0.5)

            for ii,d in enumerate(y):
                plt.annotate('%d'%d,[idx[ii],i])

            Xdata = self.scaf_Xs[i]
            Xidx = np.nonzero(Xdata)[0] # <-- list of indices where crossovers happen!
            for Xi in Xidx:
                plt.plot([Xi,Xi],[i,data[Xi]])

        plt.axis([-2,self._length,-1,6])
        plt.yticks(range(6))
        plt.xticks(range(0,self._length,7))

    def randomize_staple(self):
        for helix, data in self.helices.iteritems():
            idx = np.nonzero(data !=-1)[0]
            self.Xs[helix][idx] = np.random.random_integers(0, 1, len(idx)).astype(bool)

        for helix, data in self.helices.iteritems():
            neighbors = np.unique(data)[1:].astype(int)
            for n in neighbors:
                idx = np.nonzero(data.astype(int)==n)[0]
                self.Xs[n][idx] += self.Xs[helix][idx]

        for helix, data in self.helices.iteritems():
            Xdata = self.Xs[helix]
            self.Ks[helix][np.multiply(Xdata == False,data != -1)] = True

    def _init_scaffold(self):
        self.scaf_helices = {}
        for helix, data in self.helices.iteritems():
            self.scaf_helices[helix] = -np.ones(self._length+5)
            neighbors = np.unique(data)[1:]
            for N in neighbors:
                idx = np.nonzero(data==N)[0]
                self.scaf_helices[helix][idx+5] = N
                self.scaf_helices[helix][idx-5] = N
            self.scaf_helices[helix] = self.scaf_helices[helix][:-5]


        for helix, data in self.scaf_helices.iteritems():
            self.scaf_Xs[helix] = np.zeros(len(data))

    def randomize_scaffold(self):
        for helix, data in self.scaf_helices.iteritems():
            idx = np.nonzero(data !=-1)[0]
            self.scaf_Xs[helix][idx] = np.random.binomial(1, 0.05, len(idx)).astype(bool)

            for i,value in enumerate(self.scaf_Xs[helix]):
                if value:
                    if self.scaf_helices[helix][i-1] != -1:
                        self.scaf_Xs[helix][i-1] = True
                    elif self.scaf_helices[helix][i+1] != -1:
                        self.scaf_Xs[helix][i+1] = True

        for helix, data in self.scaf_helices.iteritems():
            neighbors = np.unique(data)[1:].astype(int)
            for n in neighbors:
                idx = np.nonzero(data.astype(int)==n)[0]
                self.scaf_Xs[n][idx] += self.scaf_Xs[helix][idx]
            self.scaf_Xs[helix] = self.scaf_Xs[helix].astype(bool)

    def _init_cadnano_scaf(self,helix, direction, length):
        foo = []
        for i in range(length):
            if direction == 3:
                foo.append([helix,i-1,helix,i+1])
            elif direction == 5:
                foo.append([helix,i+1,helix,i-1])
            if i==0:
                if direction==3: foo[i] = [-1,-1,helix,i+1]
                elif direction==5: foo[i] = [helix,i+1,-1,-1]
            elif i==length-2:
                if direction==3: foo[i] = [helix,i-1,-1,-1]
                elif direction==5: foo[i] = [-1,-1,helix,i-1]
        return foo

    def _draw_Xovers_and_kinks_SCAF(self,cad,direction,X,data,helix):
        Xidx = np.nonzero(X)[0]

        for idx in Xidx:
            if direction==5:
                if (idx+5)%7==6 or (idx-5)%7==6:
                    cad[idx] = [int(data[idx]),idx,helix,idx-1]
                elif (idx+5)%7==0 or (idx-5)%7==0:
                    cad[idx] = [helix,idx+1,int(data[idx]),idx]

            elif direction==3:
                if (idx+5)%7==0 or (idx-5)%7==0:
                    cad[idx] = [int(data[idx]),idx,helix,idx+1]
                elif (idx+5)%7 ==6 or (idx-5)%7 == 6:
                    cad[idx] = [helix,idx-1,int(data[idx]),idx]

        return cad

    def _scaf_cadnano(self):

        for helix, da in self.graph.nodes(data=True):
            cad = self._init_cadnano_scaf(helix,da['direction'],self._length)
            X = self.scaf_Xs[helix][:-1]
            data = self.scaf_helices[helix][:-1]
            self.scaf[helix] = self._draw_Xovers_and_kinks_SCAF(cad,da['direction'],X,data,helix)

    def to_cadnano(self,outfile="noodle.json"):
        self._stap_cadnano()
        self._scaf_cadnano()
        json = load_json('template.json')
        for i,Helix in enumerate(json['vstrands']):
            #for key,item in Helix.iteritems():
            Helix['skip'] = [[0]]*(self._length-1)
            Helix['loop'] = [[0]]*(self._length-1)
            Helix['scaf'] = self.scaf[i][:-1]
            Helix['stap'] = self.stap[i][:-1]
        save_json(outfile,json)

    def randomize(self):
        self.randomize_scaffold()
        self.randomize_staple()
