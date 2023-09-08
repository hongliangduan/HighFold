from itertools import combinations
import math


def rankby(x):
    r = []
    for i in range(len(x)):
        r.append(x[i][0])
        r.append(x[i][1])
    return tuple(r)


def ss_generation(pos_cys):
    n_cys = len(pos_cys)
    if n_cys <= 1:
        return []

    pair_cys = []
    for i in range(n_cys):
        for j in range(i + 1, n_cys, 1):
            pair_cys.append((i, j))

    n_ss = math.floor(n_cys / 2)
    ss_combination_candidate = list(combinations(pair_cys, n_ss))

    ss_combination = []
    for item in ss_combination_candidate:
        list_cys = []
        for i in range(n_ss):
            list_cys.append(item[i][0])
            list_cys.append(item[i][1])
        if len(list_cys) == len(list(set(list_cys))):
            ss_combination.append(item)

    ss_combination.sort(key=lambda x: rankby(x))
    ss_combination = list(set(ss_combination))
    ss_combination.sort(key=lambda x: rankby(x))

    ss_pos = []
    for i in range(len(ss_combination)):
        ss_pos_i = list(ss_combination[i])
        for j in range(len(ss_pos_i)):
            ss_pos_i[j] = (pos_cys[ss_pos_i[j][0]], pos_cys[ss_pos_i[j][1]])
        ss_pos.append(ss_pos_i)

    return ss_pos


def comb(iterable, r):
    pool = tuple(iterable)
    n = len(pool)
    if r > n:
        return
    indices = list(range(r))
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i + 1, r):
            indices[j] = indices[j - 1] + 1
        yield tuple(pool[i] for i in indices)


def n_choose_k(res_list, n, n_list, k, k_index):
    if k > len(n_list) or k < 0:
        return

    if k == 0 or len(n_list) == k:
        r = [i for i in k_index]
        for i in range(k):
            r[n_list[i]] = 1
        res_list.append(r)
        return

    k_index[n-len(n_list)] = 1
    n_choose_k(res_list, n, n_list[1:], k - 1, k_index)
    k_index[n-len(n_list)] = 0
    n_choose_k(res_list, n, n_list[1:], k, k_index)


# print('--------------n_choose_k----------------')
n_constant = 4
total = 0
for k in range(2, n_constant+1, 2):
    res = []
    n_choose_k(res, n_constant, list(range(n_constant)), k, [0]*n_constant)
    for item in res:
        index_cys = [i+1 for i in range(len(item)) if item[i] == 1]
        ss = ss_generation(index_cys)
        total += len(ss)
