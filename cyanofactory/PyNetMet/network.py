#    Copyright (C) 2012 Daniel Gamermann <daniel.gamermann@ucv.es>
#
#    This file is part of PyNetMet
#
#    PyNetMet is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    PyNetMet is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with PyNetMet.  If not, see <http://www.gnu.org/licenses/>.
#
#    
#    Please, cite us in your reasearch!
#


import Image
import ImageDraw

def stats(lista):
    """
    Calculates the average and standard deviation of the numbers in lista.
    """
    N = float(len(lista))
    xmed = sum(lista)/N
    s2 = [(xmed-ele)**2 for ele in lista]
    s = (sum(s2)/(N-1.))**.5
    return [xmed, s]

class Network:

    def __init__(self, M, names=[]):
        """ This defines a directed network network object, or undirected if M
        is symmetrical.
        M_ij=1 means the ith node has a directed link to the jth node."""
        if not names:
            nn = len(M)
            names = [str(ele) for ele in xrange(nn)]
        self.nodesnames = names
        nnodes = len(names)
        self.nnodes = nnodes
        links = []
        linksin = [[] for ii in xrange(nnodes)]
        linksout = [[] for ii in xrange(nnodes)]
        neigbs = [[] for ii in xrange(nnodes)]
        directed = [M[i][j] == M[j][i] for i in xrange(nnodes) for j in xrange(i+1,nnodes)]
        directed = reduce(lambda x, y: x and y, directed)
        self.directed = not directed
        for ii in xrange(nnodes):
            for jj in xrange(nnodes):
                if M[ii][jj]:
                    linksin[jj].append(ii)
                    linksout[ii].append(jj)
                    if ii not in neigbs[jj]:
                        neigbs[jj].append(ii)
                    if jj not in neigbs[ii]:
                        neigbs[ii].append(jj)
                    links.append((ii, jj))
        self.linksin = linksin
        self.linksout = linksout
        self.neigbs = neigbs
        self.links = links
        nlinks = len(links)
        self.nlinks = nlinks
        CCs = [[list(set(a) & set(b)) for a in self.neigbs]
                for b in self.neigbs]
        uCCs = [[list(set(a) | set(b)) for a in self.neigbs]
                 for b in self.neigbs]
        noneigbs = [ii for ii in xrange(nnodes) if len(self.neigbs[ii]) == 0]
        if len(noneigbs) != 0:
            self.no_neigbs = True
            print "Warning: nodes "+", ".join([names[ele] for ele in \
                                      noneigbs])+" have no connections at all."
            print "nCCs could not be calculated (divison by zero)."
            nCCs = []
        else:
            self.no_neigbs = False
            nCCs = [[float(len(CCs[ii][jj])+(ii in self.neigbs[jj]))/
                   (min(len(self.neigbs[ii]), len(self.neigbs[jj])))
                   for ii in xrange(nnodes)] for jj in xrange(nnodes)]
        self.uCCs = uCCs
        self.CCs = CCs
        self.nCCs = nCCs
        weight = []
        sd_wei = []
        for ele in nCCs:
            [ww, sd] = stats(ele)
            weight.append(ww)
            sd_wei.append(sd)
        self.weight = weight
        self.sd_wei = sd_wei
        Cis = [0. for i in range(self.nnodes)]
        kis = [len(ele) for ele in self.neigbs]
        for ii, neigs in enumerate(self.neigbs):
            ki = kis[ii]
            Ei = 0
            for jj, ele in enumerate(neigs):
                for kk, ele2 in enumerate(neigs[jj:]):
                    if ele in self.neigbs[ele2] or ele2 in self.neigbs[ele]:
                        Ei += 1
            if ki == 1 or ki == 0:
                Ci = 0.
            else:
                Ci = 2.*Ei/((ki-1)*ki)
            Cis[ii] = Ci
        self.Cis = Cis
        self.kis = kis

    def __repr__(self):
        if self.directed:
            dire = "Directed"
            nlinks = self.nlinks
        else:
            dire = "Undirected"
            nlinks = self.nlinks / 2
        repre = "%s Network with:\n %f Nodes \n %f Links" % (dire, self.nnodes,
                                                          nlinks)
        return repre

    def plot_nCCs(self, size=[1000,1000], grid=50, output="clust_CC.jpeg"):
        '''Orders and plot matrix with nCCs '''
        matr = [ele[:] for ele in self.nCCs]
        names = self.nodesnames[:]
        nnodes = self.nnodes
        orde = []
        nl = nnodes
        nums = range(nl)
        start = self.weight.index(min(self.weight))
        orde.append(nums.pop(start))
        old = start
        while len(nums) > 0:
            #hh = [((matr[start][ii]+matr[old][ii])/2.) for ii in xrange(nl)]
            hh=matr[start][:]
            tosearch = [sum([((hh[ii]-matr[ele][ii])/(hh[ii]+matr[ele][ii]))**2/max(0.00001,self.Cis[ii])
                        for ii in nums if matr[ele][ii]>0.]) for ele in nums]
            maax = min(tosearch)
            old = start
            start = tosearch.index(maax)
            start = nums.pop(start)
            orde.append(start)
        CCord = [[matr[n1][n2] for n2 in orde] for n1 in orde]
        names_ord = [names[n1] for n1 in orde]
        plot = Image.new("RGBA", size, (255, 255, 255, 255))
        draw = ImageDraw.Draw(plot)
        nodx = float(size[0])/nnodes
        nody = float(size[1])/nnodes
        for i in xrange(nnodes):
            for j in xrange(nnodes):
                stax = i * nodx
                stay = j * nody
                finx = (i+1) * nodx
                finy = (j+1) * nody
                if CCord[i][j] >= 0:
                    coll = int(256*(1.-CCord[i][j]))
                    lumi=(coll, coll, coll, 255)
                else:
                    lumi = (255, 25, 50, 255)
                draw.rectangle((stax, stay, finx, finy), outline=lumi,
                                fill=lumi)
        # The commands below draw a grid
        for i in xrange(0,nnodes, grid):
            for j in xrange(0, nnodes, grid):
                stax = i * nodx
                stay = j * nody
                finx = (i+1) * nodx * grid
                finy = (j+1) * nody * grid
                draw.line((stax, stay, finx, stay), fill=0)
                draw.line((stax, stay, stax, finy), fill=0)
                draw.line((finx, stay, finx, finy), fill=0)
                draw.line((stax, finy, finx, finy), fill=0)
        del draw
        plot.save(output, "JPEG")
        self.CCord = CCord
        self.names_ord = names_ord

    def plot_matr(self, matr, orde, size=[1000,1000], grid=50, 
                output="clust_CC.jpeg", returning = False):
        '''Plot matrix matr. '''
        CCord = [[matr[n1][n2] for n2 in orde] for n1 in orde]
        nnodes = len(CCord)
        plot = Image.new("RGBA", size, (255, 255, 255, 255))
        draw = ImageDraw.Draw(plot)
        nodx = float(size[0])/nnodes
        nody = float(size[1])/nnodes
        for i in xrange(nnodes):
            for j in xrange(nnodes):
                stax = i * nodx
                stay = j * nody
                finx = (i+1) * nodx
                finy = (j+1) * nody
                if CCord[i][j] >= 0:
                    coll = int(256*(1.-CCord[i][j]))
                    lumi = (coll, coll, coll, 255)
                else:
                    lumi = (255, 25, 50, 255)
                draw.rectangle((stax, stay, finx, finy), outline = lumi,
                                fill = lumi)
        # The commands below draw a grid
        for i in xrange(0, nnodes, grid):
            for j in xrange(0, nnodes, grid):
                stax = i * nodx
                stay = j * nody
                finx = (i+1) * nodx * grid
                finy = (j+1) * nody * grid
                draw.line((stax, stay, finx, stay), fill=0)
                draw.line((stax, stay, stax, finy), fill=0)
                draw.line((finx, stay, finx, finy), fill=0)
                draw.line((stax, finy, finx, finy), fill=0)
        del draw
        plot.save(output, "JPEG")
        if returning:
            return CCord

    def kruskal(self, matr, minimo=True):
        """
        Kruskal algorith for matr.
        """
        ncomps = len(matr)
        already = []
        ramas = []
        dists = []
        matrix = [ele[:] for ele in matr]
        dic_posis = [-(ele+1) for ele in xrange(ncomps)]
        comps = [[] for ele in xrange(ncomps)]
        nrama = 0
        mins = []
        for ii in xrange(ncomps):
            for jj in xrange(ii+1, ncomps):
                if matrix[ii][jj] >= 0:
                    mins.append((ii, jj, matrix[ii][jj]))
        mins.sort(key=lambda x:x[2], reverse=not minimo)
        for edge in mins:
            i = edge[0]
            j = edge[1]
            dist = edge[2]
            if dic_posis[i] != dic_posis[j]:
                dists.append(dist)
                if i not in already and j not in already:
                    ramas.append((i, j))
                    already.append(i)
                    already.append(j)
                    comps[i].append(j)
                    comps[j].append(i)
                    dic_posis[i] = nrama
                    dic_posis[j] = nrama
                    nrama += 1
                elif i not in already:
                    ramas.append((i, ramas[dic_posis[j]]))
                    comps[i].append(j)
                    already.append(i)
                    for ele in comps[j]:
                        comps[i].append(ele)
                        comps[ele].append(i)
                        dic_posis[ele] = nrama
                    comps[j].append(i)
                    dic_posis[i] = nrama
                    dic_posis[j] = nrama
                    nrama += 1
                elif j not in already:
                    ramas.append((j, ramas[dic_posis[i]]))
                    comps[j].append(i)
                    already.append(j)
                    for ele in comps[i]:
                        comps[j].append(ele)
                        comps[ele].append(j)
                        dic_posis[ele] = nrama
                    comps[i].append(j)
                    dic_posis[i] = nrama
                    dic_posis[j] = nrama
                    nrama += 1
                else:
                    ramas.append((ramas[dic_posis[i]], ramas[dic_posis[j]]))
                    for ele in comps[j]:
                        dic_posis[ele] = nrama
                        if i not in comps[ele]:
                            comps[ele].append(i)
                        if ele not in comps[i]:
                            comps[i].append(ele)
                        for ele2 in comps[i]:
                            dic_posis[ele2] = nrama
                            if j not in comps[ele2]:
                                comps[ele2].append(j)
                            if ele2 not in comps[j]:
                                comps[j].append(ele2)
                            if ele2 not in comps[ele]:
                                comps[ele].append(ele2)
                            if ele not in comps[ele2]:
                                comps[ele2].append(ele)
                    if j not in comps[i]:
                        comps[i].append(j)
                    if i not in comps[j]:
                        comps[j].append(i)
                    dic_posis[i] = nrama
                    dic_posis[j] = nrama
                    nrama += 1
        orde = [int(ele) for ele in str(ramas[nrama-1]).replace(" ", "").\
                                  replace("(", "").replace(")", "").split(",")]
        self.tree = ramas[nrama-1]
        self.dists = dists
        self.branches = ramas
        self.krusk_ord = orde

    def distance(self, enz_start, distances, net, to_visit, curdist):
        ''' Calculates distances of all other nodes from enz_start.
        '''
        while len(to_visit) > 0:
            # updates distances
            curdist += 1
            for ele in net[enz_start]:
                if distances[ele] > curdist:
                    distances[ele] = curdist
            # where to go from here
            bla = [distances[ii] for ii in to_visit]
            mindist = min(bla)
            tolook = to_visit[bla.index(mindist)]
            topop = to_visit.index(tolook)
            new_enz = to_visit.pop(topop)
            curdist = distances[new_enz]
            if curdist == "X": return None
            #print new_enz
            self.distance(new_enz, distances, net, to_visit, curdist)

    def calc_dist(self, net, ii):
        '''distance of one node to all others'''
        to_visit = range(self.nnodes)
        distances = ["X" for ele in range(self.nnodes)]
        distances[ii] = 0
        currdist = 0
        to_visit.pop(ii)
        self.distance(ii, distances, net, to_visit, currdist)
        return distances

    def calc_all_dists(self, net):
        """
        Calculates all distances between all nodes.
        """
        dist = []
        for ii in range(self.nnodes):
            dist.append(self.calc_dist(net, ii))
        return dist

    def distance_wp(self, enz_start, distances, net, to_visit, curdist, paths):
        ''' Calculates distances of all other nodes from enz_start.
        wp means with path to nodes '''
        while len(to_visit):
            # updates distances
            curdist += 1
            for ele in net[enz_start]:
                if distances[ele] > curdist:
                    distances[ele] = curdist
                    paths[ele] = [ele2 for ele2 in paths[enz_start]]
            # where to go from here
            bla = [distances[ii] for ii in to_visit]
            mindist = min(bla)
            tolook = to_visit[bla.index(mindist)]
            topop = to_visit.index(tolook)
            new_enz = to_visit.pop(topop)
            curdist = distances[new_enz]
            if curdist == "X": return None
            paths[new_enz].append(new_enz)
            self.distance_wp(new_enz, distances, net, to_visit, curdist, paths)

    def calc_dist_wp(self, net, ii):
        '''distance of one node to all others'''
        to_visit = range(self.nnodes)
        distances = ["X" for ele in range(self.nnodes)]
        distances[ii] = 0
        currdist = 0
        to_visit.pop(ii)
        paths = [[] for ele in range(self.nnodes)]
        paths[ii].append(ii)
        self.distance_wp(ii, distances, net, to_visit, currdist, paths)
        return [paths, distances]

    def calc_all_dists_wp(self, net):
        '''distance of all nodes to all others'''
        dist = []
        paths = []
        for ii in range(self.nnodes):
            [pat, dis] = self.calc_dist_wp(net,ii)
            dist.append(dis)
            paths.append(pat)
        return [paths, dist]

    def components(self):
        """
        Finds disconected components of the network.
        """
        [net, dis] = self.calc_all_dists_wp(self.neigbs)
        self.distances = dis
        self.paths = net
        comps = [reduce(lambda x,y:set(x) | set(y), ele) for ele in net]
        disc_comps = []
        for ele in comps:
            if ele not in disc_comps:
                disc_comps.append(ele)
        disc_comps = [list(ele) for ele in disc_comps]
        self.disc_comps = disc_comps

