import argparse
import phy_parser
import numpy as np
import time

parser = argparse.ArgumentParser(description = "Implementation of the Neighbour-Joining algorithm; accepts phylip-format distance matrices and outputs newick-format trees.")

parser.add_argument("file", help = "File should be a phylip-format distance-matrix")

parser.add_argument("-r", "--runtime",
                    help = "For runtime analysis purposes",
                    action = "store_true")

args = parser.parse_args()


class NeighborJoining:
    def __init__(self, disMatrix, n):
        self.taxa_names = [taxon.name for taxon in disMatrix]
        adj = self.runNeighborJoining(disMatrix, n)
        self.saveResult(adj)

    def saveResult(self, adj):
        f = open('result.newick', 'w')
        f.write(self.buildTreeString(adj, len(adj)-1))
        f.write("\n")
        f.close()

    def buildTreeString(self, adj, root):
        visited = set()
        return self.buildTreeStringHelper(adj, visited, root)

    def buildTreeStringHelper(self, adj, visited, node):
        if node in visited:
            return ""
        visited.add(node)
        children = []
        for neighbor, distance in adj[node]:
            if neighbor not in visited:
                child_str = self.buildTreeStringHelper(adj, visited, neighbor)
                if child_str:
                    children.append(child_str + ":" + "{:.2f}".format(distance))
        if len(children) == 0:
            return self.taxa_names[node]
        else:
            return "(" + ",".join(children) + ")"

    def runNeighborJoining(self, disMatrix, n):
        D = np.array(disMatrix, dtype = float)
        clusters = [i for i in range(n)]
        adj = [[] for _ in range(n)]
        if len(D) <= 1:
            return adj
        while len(clusters) > 1:
            totalDist = np.sum(D, axis = 0)
            D1 = (n-2) * D
            D1 = D1 - totalDist - totalDist.reshape((n, 1))
            np.fill_diagonal(D1, 0.)

            i, j = np.unravel_index(np.argmin(D1, axis = None), D1.shape) # Get row (i) and column (j) of the minimum value

            if n-2 == 0:
                delta = 0
            else:
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

    if args.runtime:
        st = time.time()

    NeighborJoining(disMatrix, n)

    if args.runtime:
        et = time.time()
        print(et - st)