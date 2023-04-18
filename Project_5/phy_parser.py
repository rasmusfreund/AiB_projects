from typing import TextIO

class DistanceMatrix(object):
    def __init__(self):
        self.matrix = []
        
    def add(self, taxa):
        self.matrix.append(taxa)

    def __getitem__(self, item):
        return self.matrix[item]

    def __len__(self):
        return len(self.matrix)

    def __repr__(self):
        print(self.matrix)

    def __str__(self):
        taxa = ["\0"]
        for i in range(len(self.matrix)):
            taxa.append(self.matrix[i].name)
            for j in self.matrix[i]:
                taxa.append(str(j))
            while i < len(self.matrix) - 1:
                taxa.append("\n")
                break
        return " ".join(taxa)


class NamedList(list):
    def __init__(self, name, *args):
        super().__init__(*args)
        self.name = name


def parser(file: TextIO) -> list:
    matrix = DistanceMatrix()
    with open(file, 'r') as f:
        nr_taxa = int(f.readline())
        for _ in range(nr_taxa):
            taxa = f.readline().rstrip("\n").split(maxsplit = 1)
            matrix.add(NamedList(taxa[0], map(float, taxa[1].split())))

    return matrix
