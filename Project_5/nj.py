import argparse
import sys
import phy_parser
import numpy as np

parser = argparse.ArgumentParser(description = "Implementation of the Neighbour-Joining algorithm; accepts phylip-format distance matrices and outputs newick-format trees.")

parser.add_argument("file", help = "File should be a phylip-format distance-matrix")

args = parser.parse_args()


class NeighborJoining:
    def __init__(self, disMatrix, n):
        adj = self.runNeighborJoining(disMatrix, n)
        self.printGraph(adj)
        self.saveResult(adj)

    def saveResult(self, adj):
        f = open('result.txt', 'w')
        for i, nodes in enumerate(adj):
            for d, w in nodes:
                f.write(str(i)+'->'+str(d)+':'+'%0.3f' % w+'\n')

    def printDistMatrix(self, distMatrix):
        for d in distMatrix:
            print(' '.join([str(i) for i in d]))

    def printGraph(self, adj):
        for i, nodes in enumerate(adj):
            for d, w in nodes:
                print(str(i)+'->'+str(d)+':'+'%0.3f' % w)

    def runNeighborJoining(self, disMatrix, n):
        D = np.array(disMatrix, dtype = float)
        clusters = [i for i in range(n)]
        adj = [[] for i in range(n)]
        if len(D) <= 1:
            return adj
        while True:
            if 2 == n:
                adj[len(adj)-1].append((len(adj)-2, D[0][1]))
                adj[len(adj)-2].append((len(adj)-1, D[0][1]))
                break
            totalDist = np.sum(D, axis = 0)
            D1 = (n-2) * D
            D1 = D1 - totalDist
            D1 = D1 - totalDist.reshape((n, 1))
            np.fill_diagonal(D1, 0.)

            index = np.argmin(D1)
            i = index // n
            j = index % n
            delta = (totalDist[i] - totalDist[j])/(n-2)
            li = (D[i, j]+delta)/2
            lj = (D[i, j]-delta)/2
            d_new = (D[i, :]+D[j, :]-D[i, j])/2
            D = np.insert(D, n, d_new, axis = 0)
            d_new = np.insert(d_new, n, 0., axis = 0)
            D = np.insert(D, n, d_new, axis = 1)
            D = np.delete(D, [i, j], 0)
            D = np.delete(D, [i, j], 1)

            m = len(adj)
            adj.append([])
            adj[m].append((clusters[i], li))
            adj[clusters[i]].append((m, li))
            adj[m].append((clusters[j], lj))
            adj[clusters[j]].append((m, lj))
            if i < j:
                del clusters[j]
                del clusters[i]
            else:
                del clusters[i]
                del clusters[j]
            clusters.append(m)

            n -= 1

        return adj

if __name__ == "__main__":
    disMatrix = phy_parser.parser(args.file)
    n = len(disMatrix)

    NeighborJoining(disMatrix, n)