import multiprocessing as mp
import time
import pandas as pd
from math import sqrt, inf


def next_permutation(a):
    """Generate the lexicographically next permutation inplace.
    https://en.wikipedia.org/wiki/Permutation#Generation_in_lexicographic_order
    Return false if there is no next permutation.
    """
    # Find the largest index i such that a[i] < a[i + 1]. If no such
    # index exists, the permutation is the last permutation
    for i in reversed(range(len(a) - 1)):
        if a[i] < a[i + 1]:
            break  # found
    else:  # no break: not found
        return False  # no next permutation

    # Find the largest index j greater than i such that a[i] < a[j]
    j = next(j for j in reversed(range(i + 1, len(a))) if a[i] < a[j])

    # Swap the value of a[i] with that of a[j]
    a[i], a[j] = a[j], a[i]

    # Reverse sequence from a[i + 1] up to and including the final element a[n]
    a[i + 1:] = reversed(a[i + 1:])
    return True


def get_path_length(path, matrix, start):
    ret_val = 0.0
    prev = start
    for p in path:
        ret_val += matrix[prev][p]
        prev = p
    ret_val += matrix[0][path[-1]]
    return ret_val


def task(args):
    pid = args[0]
    matrix = args[1]
    start_length = matrix[0][pid + 1]
    result = inf
    path = []
    result_path = []
    for i in range(1, len(matrix)):
        if i != (pid + 1):
            path.append(i)

    while True:
        r = get_path_length(path, matrix, pid + 1)
        if r < result:
            result = r
            result_path = path.copy()
        if not next_permutation(path):
            break
    return (result + start_length, [0, pid + 1] + result_path)


if __name__ == '__main__':
    test = pd.read_csv("tests/test_14.csv", header=None, delimiter=";").values
    N = test.shape[0]

    matrix = []
    for i in range(N):
        row = []
        for j in range(N):
            row.append(sqrt(pow(test[j, 0] - test[i, 0], 2) + pow(test[j, 1] - test[i, 1], 2)))
        matrix.append(row)
    #print(matrix)

    with mp.Pool(N - 1) as p:
        res = p.map(task, [(i, matrix) for i in range(N - 1)])
        rr = inf
        pp = []
        for r in res:
            if r[0] < rr:
                pp = r[1]
                rr = r[0]
        print(rr)
        print(pp)