import os
import numpy as np
from pdbfixer.pdbfixermodel.pdbfixer.pdbfixer import PDBFixer
from openmm.app import PDBFile
import pandas as pd
import time
from DockQ import calc_DockQ


FLAG_TO_REMOVE = ['ANISOU', 'HETATM', 'TER']
ALPHABETA = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
DICT_F2S = {'G': 'GLY', 'A': 'ALA', 'V': 'VAL', 'L': 'LEU', 'I': 'ILE', 'M': 'MET', 'F': 'PHE', 'W': 'TRP', 'P': 'PRO',
            'S': 'SER', 'T': 'THR', 'C': 'CYS', 'Y': 'TYR', 'N': 'ASN', 'Q': 'GLN',
            'D': 'ASP', 'E': 'GLU', 'K': 'LYS', 'R': 'ARG', 'H': 'HIS'}
DICT_S2F = {'GLY': 'G', 'ALA': 'A', 'VAL': 'V', 'LEU': 'L', 'ILE': 'I', 'MET': 'M', 'PHE': 'F', 'TRP': 'W', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'CYS': 'C', 'TYR': 'Y', 'ASN': 'N', 'GLN': 'Q',
            'ASP': 'D', 'GLU': 'E', 'LYS': 'K', 'ARG': 'R', 'HIS': 'H'}
MAX_VALUE = 99999


def add_missing_residues(pdb_in, pdb_out):
    fixer = PDBFixer(pdb_in)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    pdb_tmp = pdb_out[:-4]+'_tmp.pdb'
    PDBFile.writeFile(fixer.topology, fixer.positions, open(pdb_tmp, "w"), keepIds=True)

    content_out = ''
    with open(pdb_in, 'r') as fr:
        for line in fr:
            if line[0:6] == 'SEQRES':
                content_out += line
    with open(pdb_tmp, 'r') as fr:
        for line in fr:
            content_out += line
    with open(pdb_out, 'w') as fw:
        fw.write(content_out)

    return 0


def line_clean(lines):
    """
    Remove the lines that start with items in FLAG_TO_REMOVE.
    ATOM   1251  CA ACYS A 710      -4.082  81.855  -0.752  0.50 10.97           C
    ANISOU 1251  CA ACYS A 710     1728   1009   1430   -263   -176     25       C
    :param lines:
    :return:
    """
    lines_new = []
    for line in lines:
        if line[0:6].strip() not in FLAG_TO_REMOVE:
            lines_new.append(line)
    return lines_new


def altloc_unique(lines):
    """
    ATOM   1251  CA ACYS A 710      -4.082  81.855  -0.752  0.50 10.97           C
    ATOM   1252  CA BCYS A 710      -4.304  81.737  -0.649  0.50 11.02           C
    :param lines:
    :return:
    """
    lines_new = []
    line_1 = ' '*80
    for line in lines:
        if line[0:6].strip() != 'ATOM':
            lines_new.append(line)
        if line[0:6].strip() == 'ATOM' and not (line_1[12:16] == line[12:16] and line_1[17:26] == line[17:26]):
            lines_new.append(line[:16]+' '+line[17:])
        line_1 = line
    return lines_new


def atom_and_amino_renumbering(lines):
    """
    This function is to update two types of line which starts with SEQRES or ATOM.
    For the case of SEQRES, renumber the chain id in the order of 'A,B,C...'.
    For the case of ATOM, renumber the atom index and amino acid index in the order of '1,2,3...'.
    :param lines:
    :return:
    """
    lines_new = []
    cnt_chain_seqres = 0
    chain_id_seqres = ''
    chain_id_pos_seqres = 11
    cnt_chain_atom = 0
    chain_id_atom = ''
    chain_id_pos_atom = 21
    atom_index = 0
    res_index = 0
    for line in lines:
        # case 'SEQRES', update chain_id
        if line[0:6].strip() == 'SEQRES':
            if cnt_chain_seqres == 0 or (cnt_chain_seqres > 0 and chain_id_seqres != line[chain_id_pos_seqres]):
                chain_id_seqres = line[chain_id_pos_seqres]
                cnt_chain_seqres += 1
            line = line[:chain_id_pos_seqres] + ALPHABETA[cnt_chain_seqres-1] + line[chain_id_pos_seqres+1:]
            #
            lines_new.append(line)

        # case 'ATOM' and 'TER', update chain_id, atom index and amino acid index
        if line[0:6].strip() == 'ATOM' and line[76:78].strip() != 'H' and not line[12:16].strip().endswith('XT'):
            # chain_id
            if cnt_chain_atom == 0 or (cnt_chain_atom > 0 and chain_id_atom != line[chain_id_pos_atom]):
                chain_id_atom = line[chain_id_pos_atom]
                cnt_chain_atom += 1
            line = line[:chain_id_pos_atom] + ALPHABETA[cnt_chain_atom-1] + line[chain_id_pos_atom+1:]
            # atom index
            atom_index += 1
            line = line[:6] + ('%5d' % atom_index) + line[11:]
            # residue index
            if line[12:16].strip() == 'N':
                res_index += 1
            line = line[0:22] + ('%4d' % res_index) + line[26:]
            #
            lines_new.append(line)

    return lines_new


def is_atom_ca(line):
    return True if line[0:4] == 'ATOM' and line[12:16].strip() == 'CA' else False


def polymer2monomer(lines):
    """
    find the ligand and receptor, and merge the receptor into one sequence.
    assume that the ligand is just one sequence.
    :param lines:
    :return:
    """
    # statistic the length of chains
    chain_len = {}
    for line in lines:
        if is_atom_ca(line):
            chain_id = line[21]
            chain_len[chain_id] = chain_len.get(chain_id, 0) + 1
#    print(chain_len)

    # find the ligand
    ligand_len_last = MAX_VALUE
    ligand_id_last = ''
    ligand_len_first = MAX_VALUE
    ligand_id_first = ''
    for k, v in chain_len.items():
        if v <= ligand_len_last:
            ligand_len_last = v
            ligand_id_last = k
        if v < ligand_len_first:
            ligand_len_first = v
            ligand_id_first = k

    chain_id_list = list(chain_len.keys())
    chain_index_last = -1
    chain_index_first = -1
    for i in range(len(chain_id_list)):
        if ligand_id_last == chain_id_list[i]:
            chain_index_last = i
        if ligand_id_first == chain_id_list[i]:
            chain_index_first = i
    chain_index = chain_index_last
    if chain_index_first in (0, len(chain_id_list)-1):
        chain_index = chain_index_first

    # renumber the chain_id
    pdb_new = []
    for line in lines:
        if line[0:4].strip() not in ('ATOM', 'TER'):
            pdb_new.append(line)
        elif line[0:4].strip() == 'ATOM':
            chain_id_new = 'B' if line[21] == chain_id_list[chain_index] else 'A'
            pdb_new.append(line[:21] + chain_id_new + line[22:])

    return pdb_new


def pdb_standardization(pdb_in, pdb_out):
    #
    lines_raw = []
    with open(pdb_in, 'r') as fr:
        for line in fr:
            lines_raw.append(line)

    lines_tmp = line_clean(lines_raw)
    lines_tmp = altloc_unique(lines_tmp)
    lines_tmp = atom_and_amino_renumbering(lines_tmp)
    lines_tmp = polymer2monomer(lines_tmp)

    pdb_content_std = ''.join(lines_tmp)
    with open(pdb_out, 'w') as fw:
        fw.write(pdb_content_std)

    return 0


def umeyama_alignment(x: np.ndarray, y: np.ndarray, with_scale: bool = False):
    """
    Computes the least squares solution parameters of an Sim(m) matrix
    that minimizes the distance between a set of registered points.
    Umeyama, Shinji: Least-squares estimation of transformation parameters
                     between two point patterns. IEEE PAMI, 1991
    :param x: mxn matrix of points, m = dimension, n = nr. of data points
    :param y: mxn matrix of points, m = dimension, n = nr. of data points
    :param with_scale: set to True to align also the scale (default: 1.0 scale)
    :return: r, t, c - rotation matrix, translation vector and scale factor
    """
    if x.shape != y.shape:
        print("data matrices must have the same shape")
        return 0, 0, 0

    # m = dimension, n = nr. of data points
    m, n = x.shape

    # means, eq. 34 and 35
    mean_x = x.mean(axis=1)
    mean_y = y.mean(axis=1)

    # variance, eq. 36
    # "transpose" for column subtraction
    sigma_x = 1.0 / n * (np.linalg.norm(x - mean_x[:, np.newaxis]) ** 2)

    # covariance matrix, eq. 38
    outer_sum = np.zeros((m, m))
    for i in range(n):
        outer_sum += np.outer((y[:, i] - mean_y), (x[:, i] - mean_x))
    cov_xy = np.multiply(1.0 / n, outer_sum)

    # SVD (text betw. eq. 38 and 39)
    u, d, v = np.linalg.svd(cov_xy)
    if np.count_nonzero(d > np.finfo(d.dtype).eps) < m - 1:
        print("Degenerate covariance rank, Umeyama alignment is not possible")
        return 0, 0, 0

    # S matrix, eq. 43
    s = np.eye(m)
    if np.linalg.det(u) * np.linalg.det(v) < 0.0:
        # Ensure a RHS coordinate system (Kabsch algorithm).
        s[m - 1, m - 1] = -1

    # rotation, eq. 40
    r = u.dot(s).dot(v)

    # scale & translation, eq. 42 and 41
    c = 1 / sigma_x * np.trace(np.diag(d).dot(s)) if with_scale else 1.0
    t = mean_y - np.multiply(c, r.dot(mean_x))

    return r, t, c


def pdb_align_model(x, y):
    """
    y = c * R * x + t
    :param x: info of pep_1, dict
    :param y: info of pep_2, sorted according to pep_1, dict
    :param mode: 'ca' for C_alpha only, 'ncc' for N-C_alpha-C, 'all' for all atoms, 'all_but_h' for all atoms but hydrogen
    :return:
    """
    x_ca_index = x.get('ca_index')
    y_ca_index = y.get('ca_index')
    x_coord = []
    y_coord = []
    for i in range(len(x_ca_index)):
        x_coord.append(x.get('coord')[x_ca_index[i]])
        y_coord.append(y.get('coord')[y_ca_index[i]])

    r, t, c = umeyama_alignment(np.array(x_coord).T, np.array(y_coord).T)

    return r, t, c


def pdb_parse(file_pdb):
    """
    extract coordinate info from pdb file
    :param file_pdb:
    :return:
    """
    # id = ''
    seqres = []
    uid = []  # 1_CA
    amino_id = []  # 1
    amino_id_by_chain = {}
    atom_name = []  # CA
    coord = []  #
    element = []  # C
    ca_index = []
    ncc_index = []
    non_hydrogen_index = []
    with open(file_pdb, 'r') as f:
        for line in f:
            if line[0:6].strip() == 'ATOM':  # if line[0:6].strip() in ('ATOM', 'HETATM'):
                uid.append(line[22:26].strip() + '_' + line[12:16].strip())
                amino_id.append(int(line[22:26].strip()))
                aas = amino_id_by_chain.get(line[21],[])
                if int(line[22:26].strip()) not in aas:
                    aas.append(int(line[22:26].strip()))
                amino_id_by_chain[line[21]] = aas
                atom_name.append(line[12:16].strip())
                coord.append([float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())])
                element.append(line[76:78].strip())
                if line[12:16].strip() == 'CA':
                    ca_index.append(len(atom_name) - 1)
                    seqres.append(line[17:20].strip())
                if line[12:16].strip() in ('N', 'CA', 'C'):
                    ncc_index.append(len(atom_name) - 1)
                if line[76:78].strip() != 'H' and not line[12:16].strip().endswith('XT'):
                    non_hydrogen_index.append(len(atom_name) - 1)
    fasta = ''.join([DICT_S2F.get(x, 'X') for x in seqres])
    info_pep = {'fasta': fasta, 'seqres': seqres, 'uid': uid, 'amino_id': amino_id, 'atom_name': atom_name,
                'coord': coord, 'element': element, 'ca_index': ca_index, 'ncc_index': ncc_index,
                'non_hydrogen_index': non_hydrogen_index, 'amino_id_by_chain': amino_id_by_chain}
    return info_pep


def align_pdb_model2native(file_native, file_model, file_model_align):
    index_type = 'ca_index'
    non_hydrogen_index = 'non_hydrogen_index'
    info_native = pdb_parse(file_native)
    info_model = pdb_parse(file_model)

    if len(info_native.get(index_type)) != len(info_model.get(index_type)):
        return info_native.get(non_hydrogen_index), np.array([])

    pam_rot, pam_trm, pam_sca = pdb_align_model(info_model, info_native)
    coord_model = np.array(info_model.get('coord'))[info_model.get(non_hydrogen_index)].T
    # coord_model = np.array(info_model.get(non_hydrogen_index)).T
    coord_model_aligned = pam_rot.dot(coord_model) + np.tile(pam_trm, [len(info_model.get(non_hydrogen_index)), 1]).T

    update_aligned_coord(file_model, list(coord_model_aligned.T), file_model_align)

    return info_native.get(non_hydrogen_index), coord_model_aligned.T


def update_aligned_coord(pdb_file, coord, pdb_file_align):
    """
    update aligned coord, chain index and residue index
    :param pdb_file:
    :param coord:
    :param pdb_file_align:
    :return:
    """
    pdb_new = ''
    index = 0
    with open(pdb_file, 'r') as fr:
        for line in fr:
            if line[0:6].strip() == 'ATOM' and line[76:78].strip() != 'H' and not line[12:16].strip().endswith('XT'):
               line = line[0:30] + ('%8.3f' % coord[index][0]) + \
                      ('%8.3f' % coord[index][1]) + ('%8.3f' % coord[index][2]) + line[54:]
               index += 1
            pdb_new += line

    with open(pdb_file_align, 'w') as fw:
        fw.write(pdb_new)
    return 0


def calculate_rmsd_ca_ncc_all(f1, f2):
    info1 = pdb_parse(f1)
    # print(info1)
    info2 = pdb_parse(f2)
    # print(info2)

    fasta = [info1.get('fasta'), info2.get('fasta')]
    rmsd = [-1, -1, -1]
    rmsd_chain = []
    if fasta[0] != fasta[1]:
        print('there are the different amino acids in the two input peptides')
        return rmsd, rmsd_chain

    index_chain = dict(info1.get('amino_id_by_chain'))
    n_chain = len(info1.get('amino_id_by_chain'))
    for i in range(n_chain):
        rmsd_chain.append([])
    index_type = 'ca_index'
    if len(info1.get(index_type)) != len(info2.get(index_type)):
        return rmsd, rmsd_chain
    coord1 = np.array(info1.get('coord'))[info1.get(index_type)].T
    coord2 = np.array(info2.get('coord'))[info2.get(index_type)].T
    diff = coord1 - coord2
    rmsd[0] = np.sqrt(np.sum(diff * diff) / coord1.shape[1])
    i = 0
    for k in index_chain.keys():
        index1 = np.array(info1.get(index_type))[[x-1 for x in index_chain.get(k)]]
        coord1 = np.array(info1.get('coord'))[index1].T
        index2 = np.array(info2.get(index_type))[[x-1 for x in index_chain.get(k)]]
        coord2 = np.array(info2.get('coord'))[index2].T
        diff = coord1 - coord2
        rmsd_chain[i].append(np.sqrt(np.sum(diff * diff) / coord1.shape[1]))
        i += 1

    index_type = 'ncc_index'
    if len(info1.get(index_type)) != len(info2.get(index_type)):
        return rmsd, rmsd_chain
    coord1 = np.array(info1.get('coord'))[info1.get(index_type)].T
    coord2 = np.array(info2.get('coord'))[info2.get(index_type)].T
    diff = coord1 - coord2
    rmsd[1] = np.sqrt(np.sum(diff * diff) / coord1.shape[1])
    i = 0
    for k in index_chain.keys():
        index1 = np.array(info1.get(index_type))[[x-1 for x in index_chain.get(k)]]
        coord1 = np.array(info1.get('coord'))[index1].T
        index2 = np.array(info2.get(index_type))[[x-1 for x in index_chain.get(k)]]
        coord2 = np.array(info2.get('coord'))[index2].T
        diff = coord1 - coord2
        rmsd_chain[i].append(np.sqrt(np.sum(diff * diff) / coord1.shape[1]))
        i += 1

    index_type = 'non_hydrogen_index'
    if len(info1.get(index_type)) != len(info2.get(index_type)):
        return rmsd, rmsd_chain
    coord1 = np.array(info1.get('coord'))[info1.get(index_type)].T
    coord2 = np.array(info2.get('coord'))[info2.get(index_type)].T
    diff = coord1 - coord2
    rmsd[2] = np.sqrt(np.sum(diff * diff) / coord1.shape[1])
    i = 0
    for k in index_chain.keys():
        index1 = np.array(info1.get(index_type))[[x-1 for x in index_chain.get(k)]]
        coord1 = np.array(info1.get('coord'))[index1].T
        index2 = np.array(info2.get(index_type))[[x-1 for x in index_chain.get(k)]]
        coord2 = np.array(info2.get('coord'))[index2].T
        diff = coord1 - coord2
        rmsd_chain[i].append(np.sqrt(np.sum(diff * diff) / coord1.shape[1]))
        i += 1

    return rmsd, rmsd_chain


def get_all_files(folder):
    names = []
    names_only = []
    for root, dirs, files in os.walk(folder):
        if root != folder:
            break
        for file in files:
            path = os.path.join(root, file)
            names.append(path)
            names_only.append(file)
    return names, names_only


def get_chain_info_from_header(pdb_with_header):
    chain_info = {'chain_len': [], 'chain_id': [], 'ligand_index': -1, 'fasta': []}
    fasta = ''
    with open(pdb_with_header, 'r') as fr:
        for line in fr:
            if line[0:6].strip() == 'SEQRES':
                if line[11] not in chain_info['chain_id']:
                    chain_info['chain_id'].append(line[11])
                    chain_info['chain_len'].append(int(line[12:18].strip()))
                fasta += ''.join([DICT_S2F.get(x, 'X') for x in line[19:].strip().split(' ')])
    # fasta
    chain_fasta = []
    index = 0
    for i in range(len(chain_info['chain_len'])):
        chain_fasta.append(fasta[index:index+chain_info['chain_len'][i]])
        index += chain_info['chain_len'][i]
    chain_info['fasta'] = chain_fasta
    # ligand_index
    pep_len = chain_info['chain_len'][0]
    chain_info['ligand_index'] = 0
    for i in range(len(chain_info['chain_len'])):
        if pep_len > chain_info['chain_len'][i]:
            pep_len = chain_info['chain_len'][i]
            chain_info['ligand_index'] = i

    return chain_info


def run(pdb_native_raw, pdb_native_fix, pdb_native_std, pdb_model_raw, pdb_model_std, pdb_model_std_align):
    # add the missing coord of residues to pdb file
    add_missing_residues(pdb_native_raw, pdb_native_fix)
    # pdb file preprocess, incuding renumbering atom index, renumbering residue index and removing the noisy line
    pdb_standardization(pdb_native_fix, pdb_native_std)
    pdb_standardization(pdb_model_raw, pdb_model_std)
    # coord align between the predicted pdb and the native pdb
    align_pdb_model2native(pdb_native_std, pdb_model_std, pdb_model_std_align)
    # extract chain infor from header
    chain_info = {}
    try:
        chain_info = get_chain_info_from_header(pdb_native_std)
    except Exception:
        print('error format in header')

    rmsd = [-1,-1,-1]
    rmsd_chain = []
    info = {'DockQ':-1, 'irms':-1, 'Lrms':-1, 'fnat':-1}
    if os.path.exists(pdb_model_std_align):
        rmsd, rmsd_chain = calculate_rmsd_ca_ncc_all(pdb_native_std, pdb_model_std_align)
        try:
            info = calc_DockQ(pdb_model_std_align, pdb_native_std)
        except Exception as e:
            print('in dockq:', e)

    # result show
    result = {}
    result['pid'] = pdb_native_raw.split('/')[-1][:4]
    result['pdb_native'] = pdb_native_raw
    result['pdb_model'] = pdb_model_raw
    for k, v in chain_info.items():
        result[k] = v
    for k,v in info.items():
        result[k.lower()] = v
    result['rmsd_ca_complex'] = rmsd[0]
    result['rmsd_ncc_complex'] = rmsd[1]
    result['rmsd_all_complex'] = rmsd[2]
    result['rmsd_ca_ligand'] = -1
    result['rmsd_ncc_ligand'] = -1
    result['rmsd_all_ligand'] = -1

    return result


# run batch
def run_batch(path_common, folder_native, folder_model, model_suffix):
    current_time = time.strftime('%Y%m%d%H%M%S', time.localtime())
    path_common_tmp = path_common + 'tmp_' + current_time + '/'
    if not os.path.exists(path_common_tmp):
        os.makedirs(path_common_tmp)

    files, names = get_all_files(path_common + folder_native + '/')

    columns = ['pid', 'pdb_native', 'pdb_model', 'chain_len', 'chain_id', 'ligand_index', 'fasta',
               'dockq', 'irms', 'lrms', 'fnat', 'rmsd_ca_complex', 'rmsd_ncc_complex', 'rmsd_all_complex',
               'rmsd_ca_ligand', 'rmsd_ncc_ligand', 'rmsd_all_ligand', 'rmsd_chain']
    data = [[] for _ in range(len(columns))]
    for i in range(len(files)):
        pdb_native_raw = path_common + folder_native + '/' + names[i]
        pdb_native_fix = path_common_tmp + names[i][:-4] + '_fix.pdb'
        pdb_native_std = path_common_tmp + names[i][:-4] + '_std.pdb'
        pdb_model_raw = path_common + folder_model + '/' + names[i][:-4] + model_suffix + '.pdb'
        pdb_model_std = path_common_tmp + names[i][:-4] + '_model_std.pdb'
        pdb_model_std_align = path_common_tmp + names[i][:-4] + '_model_align.pdb'
        res = run(pdb_native_raw, pdb_native_fix, pdb_native_std, pdb_model_raw, pdb_model_std, pdb_model_std_align)

        for k in range(len(columns)):
            data[k].append(res.get(columns[k]))

    df = pd.DataFrame(np.array(data).T, columns=columns)
    df.to_csv(path_common + 'result_'+current_time+'.csv')


path_common = './examples/'
folder_native = 'native'
folder_model = 'af'
model_suffix = '_af'
run_batch(path_common, folder_native, folder_model, model_suffix)
