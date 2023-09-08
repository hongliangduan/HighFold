import numpy as np


def get_opt_path(matrix):
    """
    Floyd algorithm to find the shortest path
    :param matrix:
    :return:
    """
    path = np.zeros_like(matrix)
    for i in range(matrix.shape[0]):
        path[i] = [j for j in range(matrix.shape[0])]

    for m in range(matrix.shape[0]):
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[0]):
                if matrix[i][m] + matrix[m][j] < matrix[i][j]:
                    matrix[i][j] = matrix[i][m] + matrix[m][j]
                    path[i][j] = m

    return matrix


def calc_offset_matrix(seq_cycpep, c1, c2):
    """
    :param seq_cycpep: string of amino acid residues
    :param c1: list of Cys residues
    :param c2: list of the corresponding Cys residues to c1, the same length with c1
    :return:
    """
    if len(c1) != len(c2):
        return []

    # init adjacency matrix
    n_aa = len(seq_cycpep)
    matrix = np.zeros((n_aa, n_aa)) + n_aa
    for i in range(n_aa):
        matrix[i][i] = 0

    # linear peptide connection
    for i in range(n_aa - 1):
        matrix[i][i + 1] = 1
        matrix[i + 1][i] = 1

    # nc connection
    matrix[0][n_aa - 1] = 1
    matrix[n_aa - 1][0] = 1

    # ss connection
    for i in range(len(c1)):
        matrix[c1[i]][c2[i]] = 1
        matrix[c2[i]][c1[i]] = 1

    # get the shortest path
    matrix = get_opt_path(matrix)

    return matrix

