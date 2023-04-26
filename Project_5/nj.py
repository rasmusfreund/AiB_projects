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
    """
    Neighbor Joining class for creating a tree from distance matrix
    """

    def __init__(self, disMatrix, n):
        """
        Initializes the NeighborJoining object

        Args:
        - disMatrix: A list of distance between taxa
        - n: The number of taxa
        """
        self.taxa_names = [taxon.name for taxon in disMatrix]
        adj = self.runNeighborJoining(disMatrix, n)
        self.saveResult(adj)

    def saveResult(self, adj):
        """
        Save the resulting tree to a file in Newick format

        Args:
        - adj: Adjacency list of the resulting tree
        """
        f = open('result.newick', 'w')
        f.write(self.buildTreeString(adj, len(adj)-1))
        f.write("\n")
        f.close()

    def buildTreeString(self, adj, root):
        """
        Build a Newick-format string for the resulting tree

        Args:
        - adj: Adjacency list of the resulting tree
        - root: The root node of the resulting tree
        """
        visited = set()
        return self.buildTreeStringHelper(adj, visited, root)

    def buildTreeStringHelper(self, adj, visited, node):
        """
        Recursive helper function for building Newick-format string for the resulting tree

        Args:
        - adj: Adjacency list of the resulting tree
        - visited: A set of visited nodes
        - node: The current node to process
        """
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
        """
        Runs the neighbor-joining algorithm on the input distance matrix.

        Args:
        - disMatrix (list of lists of floats): The distance matrix to run the algorithm on.
        - n (int): The number of taxa in the distance matrix.

        Returns:
        - adj (list of lists of tuples): The adjacency list representation of the resulting tree.
        """
        # Distance matrix as a numpy array
        D = np.array(disMatrix, dtype = float)

        # List of clusters
        clusters = [i for i in range(len(D))]

        # Adjacency list of the resulting tree
        adj = [[] for _ in range(len(D))]

        # If there is only one taxon in the distance matrix, return an empty adjacency list
        if len(D) == 1:
            return adj

        while len(clusters) > 1:
            # Compute the total distance for each taxon
            totalDist = np.sum(D, axis = 0)

            # Compute the matrix D' based on the current distance matrix D
            D1 = (n-2) * D
            D1 = D1 - totalDist - totalDist.reshape((n, 1))
            np.fill_diagonal(D1, 0.)

            # Find the indices of the minimum value in D'
            i, j = np.unravel_index(np.argmin(D1, axis = None), D1.shape)

            # Compute the branch lengths of the two new branches
            if n-2 == 0:
                delta = 0
            else:
                delta = (totalDist[i] - totalDist[j])/(n-2)
            li = (D[i, j]+delta)/2
            lj = (D[i, j]-delta)/2

            # Compute the distances from the new node to all other nodes
            d_new = (D[i, :]+D[j, :]-D[i, j])/2

            # Update the distance matrix D
            D = np.insert(D, n, d_new, axis = 0)
            d_new = np.insert(d_new, n, 0., axis = 0)
            D = np.insert(D, n, d_new, axis = 1)
            D = np.delete(D, [i, j], 0)
            D = np.delete(D, [i, j], 1)

            # Add the new node to the adjacency list
            m = len(adj)
            adj.append([])
            adj[m].append((clusters[i], li))
            adj[clusters[i]].append((m, li))
            adj[m].append((clusters[j], lj))
            adj[clusters[j]].append((m, lj))

            # Remove the two merged clusters from the list of clusters
            if i < j:
                del clusters[j]
                del clusters[i]
            else:
                del clusters[i]
                del clusters[j]

            # Add the new cluster to the list of clusters
            clusters.append(m)

            # Decrement the number of taxa
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