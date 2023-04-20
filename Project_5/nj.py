import argparse
import phy_parser
import numpy as np

parser = argparse.ArgumentParser(description = "Implementation of the Neighbour-Joining algorithm; accepts phylip-format distance matrices and outputs newick-format trees.")

parser.add_argument("file", help = "File should be a phylip-format distance-matrix")

args = parser.parse_args()


def neighbour_joining(dist_matrix):
    # Convert distance matrix to numpy array
    D = np.array(dist_matrix)
    # Create a dictionary to hold the labels of the nodes in the tree
    node_labels = {i: str(i) for i in range(len(D))}
    # Create a dictionary to hold the distances between nodes in the tree
    node_distances = {}

    # Keep iterating until we've reduced the matrix to a single element
    while len(D) > 1:
        # Calculate the Q matrix
        n = len(D)
        Q = np.zeros((n, n))
        for i in range(n):
            for j in range(i+1, n):
                total_distance_i = np.sum(D[i])
                total_distance_j = np.sum(D[j])
                Q[i,j] = (n-2)*D[i,j] - total_distance_i - total_distance_j
        # Find the minimum value in the Q matrix
        min_index = np.argmin(Q)
        i, j = np.unravel_index(min_index, Q.shape)
        delta = (np.sum(D[i]) - np.sum(D[j])) / (n - 2)
        # Calculate the new distance between the nodes
        new_distance = (D[i,j] + delta, D[i,j] - delta)
        # Create a new node to represent the cluster
        new_label = "(" + node_labels[i] + ":" + str(new_distance[0]) + "," + node_labels[j] + ":" + str(new_distance[1]) + ")"
        node_labels[n] = new_label
        node_distances[new_label] = new_distance
        # Update the distance matrix
        new_row = np.zeros(n-1)
        new_D = np.zeros((n-1, n-1))
        for k in range(n):
            if k == i or k == j:
                continue
            new_row[k if k < j else k-1] = (D[i,k] + D[j,k] - D[i,j]) / 2
            for l in range(n):
                if l == i or l == j:
                    continue
                new_D[k if k < j else k-1, l if l < j else l-1] = D[k,l]
        new_D[-1] = new_row
        new_D[:,-1] = new_row
        D = new_D

    # Return the final tree in Newick format
    return node_labels[0] + ";"


def main():

    dist_matrix = phy_parser.parser(args.file)
    print(dist_matrix)


if __name__ == "__main__":
    main()