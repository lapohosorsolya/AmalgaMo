import os
import numpy  as np
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd
import logomaker
import matplotlib.pyplot as plt


def read_individual_motif_files(input_dir, fmtcode):
    '''
    Read and collect all motifs from a directory.
    
    Parameters
    ----------
    input_dir : string
        path to the directory
    fmtcode : string
        motif format code; one of ["hocomoco", "jaspar", "meme"]

    Returns
    -------
    motif_dict : dict with numpy arrays
        collection of motif PFMs
    '''
    files = sorted(os.listdir(input_dir))
    motif_dict = {}
    for file in files:
        name, motif = _read_motif(os.path.join(input_dir, file), fmtcode)
        motif_dict[name] = { 'pfm': motif }
    return motif_dict


def _read_motif(filename, fmtcode):
    '''
    Read a motif file and convert it to a numpy array. 
    Before calling this function, check if the format code is valid.

    Parameters
    ----------
    filename : string
        path to the motif file
    fmtcode : string
        motif format code; one of ["hocomoco", "jaspar", "meme"]

    Returns
    -------
    name : string
        motif id
    motif : numpy array
        numpy array with rows corresponding to nucleotides [A, C, G, T]
    '''
    with open(filename, 'r') as f:
        lines = f.readlines()
        lines = [ line.rstrip() for line in lines ]
    name, arr = eval('_{}_lines_to_numpy'.format(fmtcode))(lines)
    return name, arr


def _parse_multi_motif_file(filename, fmtcode):
    '''
    Split a file containing multiple motifs.

    Parameters
    ----------
    filename : string
        path to the motif file
    fmtcode : string
        motif format code; one of ["hocomoco", "jaspar", "meme"]
    
    Returns
    -------
    motif_lines : list
        lines split by empty newlines
    '''
    with open(filename, 'r') as f:
        lines = f.readlines()
        lines = [ line.rstrip() for line in lines ]
    # get rid of header (for MEME files)
    if fmtcode == 'meme':
        for i in range(len(lines)):
            tokens = lines[i].split()
            if len(tokens) > 0:
                if tokens[0] == 'MOTIF':
                    start_idx = i
                    break
        lines = lines[start_idx:]
    # split lines by whitespace
    motif_lines = []
    current_lines = []
    for line in lines:
        if not line.strip():
            motif_lines.append(current_lines)
            current_lines = []
        else:
            current_lines.append(line)
    return motif_lines


def read_single_motif_file(filename, fmtcode):
    '''
    Read all motifs from a single file.
    
    Parameters
    ----------
    filename : string
        path to the input file
    fmtcode : string
        motif format code; one of ["hocomoco", "jaspar", "meme"]

    Returns
    -------
    motif_dict : dict with numpy arrays
        collection of motif PFMs
    '''
    list_of_lines = _parse_multi_motif_file(filename, fmtcode)
    motif_dict = {}
    for lines in list_of_lines:
        name, motif = eval('_{}_lines_to_numpy'.format(fmtcode))(lines)
        motif_dict[name] = { 'pfm': motif }
    return motif_dict


def _hocomoco_lines_to_numpy(lines):
    '''
    Convert lines read from a HOCOMOCO motif file to numpy format.
    '''
    name = lines[0].split('\t')[0].split('>')[-1]
    rows = []
    l = len(lines) - 1
    for i in range(1, l+1):
        rows.append(lines[i].split('\t'))
    arr = np.zeros((4, l), dtype = float)
    for i in range(4):
        for j in range(l):
            arr[i][j] = float(rows[j][i])
    return name, arr


def _jaspar_lines_to_numpy(lines):
    '''
    Convert lines read from a JASPAR motif file to numpy format.
    '''
    id = lines[0].split('\t')[0].split('>')[-1]
    tf = lines[0].split('\t')[1]
    name = id + '_' + tf
    rows = [ lines[i].split('[')[-1].split(']')[0].split() for i in range(1, 5) ]
    l = len(rows[0])
    arr = np.zeros((4, l))
    for i in range(4):
        for j in range(l):
            arr[i][j] = float(rows[i][j])
    return name, arr


def _meme_lines_to_numpy(lines):
    '''
    Convert lines read from a MEME motif file to numpy format.
    '''
    for i in range(len(lines)):
        tokens = lines[i].split()
        if len(tokens) > 0:
            if tokens[0] == 'MOTIF':
                name = tokens[1]
                start_idx = i+2
                break
    good_lines = lines[start_idx:]
    cols = []
    for line in good_lines:
        tokens = line.split()
        if len(tokens) == 4:
            cols.append(tokens)
        else:
            break
    l = len(cols)
    arr = np.zeros((4, l))
    for j in range(l):
        for i in range(4):
            arr[i][j] = float(cols[j][i])
    return name, arr


def write_motif(filename, name, motif, fmtcode):
    '''
    Write a motif to file, converting it to the specified format. 
    Before calling this function, check if the format code is valid.

    Parameters
    ----------
    filename : string
        path to the motif file
    motif : numpy array
        numpy array with rows corresponding to nucleotides [A, C, G, T]
    fmtcode : string
        motif format code; one of ["hocomoco", "jaspar", "meme"]
    '''
    lines = eval('_numpy_to_{}_lines'.format(fmtcode))(name, motif)
    with open(filename, 'w') as f:
        f.writelines(lines)


def write_multiple_motifs(filename, motif_dict, fmtcode):
    '''
    Write a motif to file, converting it to the specified format. 
    Before calling this function, check if the format code is valid.

    Parameters
    ----------
    filename : string
        path to the motif file
    motif_dict : dict with numpy arrays
        collection of PPMs to write to file
    fmtcode : string
        motif format code; one of ["hocomoco", "jaspar", "meme"]
    '''
    if fmtcode == 'meme':
        lines = ['MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n']
    else:
        lines = []
    with open(filename, 'w') as f:
        f.writelines(lines)
    for name, motif in motif_dict.items():
        lines = eval('_numpy_to_{}_lines'.format(fmtcode))(name, motif, no_header = True)
        lines.append('\n')
        with open(filename, 'a') as f:
            f.writelines(lines)


def _numpy_to_hocomoco_lines(name, motif, **kwargs):
    '''
    Convert a numpy array to HOCOMOCO motif file lines.
    '''
    lines = ['>{}\n'.format(name)]
    for j in range(motif.shape[1]):
        lines.append('\t'.join([ '{0:.6f}'.format(i) for i in motif[:, j] ]) + '\n')
    return lines


def _numpy_to_jaspar_lines(name, motif, **kwargs):
    '''
    Convert a numpy array to JASPAR motif file lines.
    '''
    lines = ['>{}\t{}\n'.format(name, name)]
    for i, base in enumerate(['A', 'C', 'G', 'T']):
        lines.append(base + '  [\t' + '\t'.join([ '{0:.6f}'.format(j) for j in motif[i, :].tolist() ]) + '  ]\n')
    return lines


def _numpy_to_meme_lines(name, motif, **kwargs):
    '''
    Convert a numpy array to MEME motif file lines.
    '''
    lines = ['MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n']
    if len(kwargs) > 0:
        if kwargs['no_header'] == True:
            lines = []
    lines.append('MOTIF {}\nletter-probability matrix: alength= 4 w= {} nsites= 100\n'.format(name, motif.shape[1]))
    for i in range(motif.shape[1]):
        lines.append('  '.join([ '{0:.6f}'.format(j) for j in motif[:, i].tolist() ]) + '\n')
    return lines


def get_motif_reverse_complement(motif):
    '''
    Get the reverse complement of a motif matrix.

    Parameters
    ----------
    motif : numpy array
        motif weights; first dimension should represent DNA bases, second dimension should represent positions

    Returns
    -------
    rc : numpy array
        reversed and complemented motif matrix
    '''
    complement = { 0: 3, 1: 2, 2: 1, 3: 0 }
    rc = np.zeros(motif.shape)
    l = motif.shape[1]
    for i in range(motif.shape[0]):
        for j in range(motif.shape[1]):
            rc[complement[i]][l - j - 1] = motif[i][j]
    return rc


def make_info_content_logo(motif):
    '''
    Convert a motif matrix (PFM) into a logo using information content.

    Parameters
    ----------
    motif : numpy array
        motif weights; first dimension should represent DNA bases, second dimension should represent positions

    Returns
    -------
    ic_logo : numpy array
        the information content logo matrix
    '''
    seq_len = motif.shape[1]
    H_l = np.zeros(seq_len)
    R_l = np.zeros(seq_len)
    f_bl = motif / np.sum(motif, axis = 0)
    for l in range(seq_len):
        H_l[l] = - sum([ f_bl[b][l]*np.log2(f_bl[b][l] + 0.01) for b in range(4) ])
    R_l = 2 - (H_l + 0.01)
    ic_logo = np.zeros((4, seq_len))
    for l in range(seq_len):
        for b in range(4):
            ic_logo[b, l] = f_bl[b][l] * R_l[l]
    return ic_logo


def plot_ic_logo(ppm, ax):
    '''
    Plot the information content logo of a PPM.

    Parameters
    ----------
    ppm : numpy array
        position-probability matrix; first dimension should represent DNA bases, second dimension should represent positions
    '''
    color_scheme = {'A': '#487b43', 'G': 'orange', 'C': 'royalblue', 'T': '#da2358'}
    motif = make_info_content_logo(ppm)
    df = pd.DataFrame(motif.T)
    df.columns = ['A', 'C', 'G', 'T']
    logomaker.Logo(df, ax = ax, color_scheme = color_scheme)


def write_logo(filepath, motif):
    '''
    Write a motif logo to file.
    '''
    fig, ax = plt.subplots(figsize = (10, 1))
    plot_ic_logo(motif, ax)
    ax.set_ylim(0, 2)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.savefig(filepath, bbox_inches = 'tight')
    plt.close(fig)


def get_siwcs_score(ppm1, ppm2, min_overlap = 0.9):
    '''
    Calculate shared information-weighted cosine similarity scores between two PPMs. 
    Returns the best possible score from the best possible alignment.
    When calculating a score for an overhanging alignment, the longer motif (or ppm1) is padded with uniform probabilities.

    Parameters
    ----------
    ppm1, ppm2 : numpy arrays
        motif position probability matrices; first dimension should represent DNA bases, second dimension should represent positions
    min_overlap : float
        minimum proportion of the shorter motif that must overlap with the longer motif; default 7/8 bases must overlap

    Returns
    -------
    best_score : numpy array
        the best length-normalized score obtained by comparing the motif positions 
    best_pos : numpy array
        the alignment position resulting in the best score (relative to the longer motif, or ppm1 if they have the same length); this needs adjustment based on max_overhang
    max_overhang : int
        the maximum overhang allowed for alignment; should be used to adjust position
    '''
    l1 = ppm1.shape[1]
    l2 = ppm2.shape[1]
    # calculate the max overhang allowed
    max_overhang = int(min(l1, l2) * (1 - min_overlap))
    # pad the longer motif on each end, by max overhang
    if l1 >= l2:
        len_1 = l1 + 2*max_overhang
        len_2 = l2
        padded = np.zeros((4, len_1))
        padded[:] = 0.25
        padded[:, max_overhang:l1+max_overhang] = ppm1
        ppm_1 = padded
        ppm_2 = ppm2
    else:
        len_1 = l1
        len_2 = l2 + 2*max_overhang
        padded = np.zeros((4, len_2))
        padded[:] = 0.25
        padded[:, max_overhang:l2+max_overhang] = ppm2
        ppm_1 = ppm1
        ppm_2 = padded
    similarities = cosine_similarity(ppm_1.T, ppm_2.T) # the longer motif is padded here
    # walk motif 2 along motif 1
    if len_1 >= len_2:
        uniform_sims = cosine_similarity(ppm_1.T, np.array([0.25, 0.25, 0.25, 0.25]).reshape(1, -1)).flatten()
        n = len_1 - len_2 + 1
        scores = np.zeros(n)
        for i in range(n):
            score_by_pos = uniform_sims.copy()
            score_by_pos[i:i+len_2] = np.array([ similarities[i+j][j] for j in range(len_2) ])
            padded_ppm_2 = np.zeros((4, len_1))
            padded_ppm_2[:] = 0.25
            padded_ppm_2[:, i:i+len_2] = ppm_2.copy()
            si = calculate_shared_information(ppm_1, padded_ppm_2)
            si_weights = si / np.sum(si)
            weighted_score_by_pos = score_by_pos * si_weights
            scores[i] = np.sum(weighted_score_by_pos)
    # walk motif 1 along motif 2
    else:
        uniform_sims = cosine_similarity(ppm_2.T, np.array([0.25, 0.25, 0.25, 0.25]).reshape(1, -1)).flatten()
        n = len_2 - len_1 + 1
        scores = np.zeros(n)
        for i in range(n):
            score_by_pos = uniform_sims.copy()
            score_by_pos[i:i+len_1] = np.array([ similarities[j][i+j] for j in range(len_1) ])
            padded_ppm_1 = np.zeros((4, len_2))
            padded_ppm_1[:] = 0.25
            padded_ppm_1[:, i:i+len_1] = ppm_1.copy()
            si = calculate_shared_information(padded_ppm_1, ppm_2)
            si_weights = si / np.sum(si)
            weighted_score_by_pos = score_by_pos * si_weights
            scores[i] = np.sum(weighted_score_by_pos)
    # return the best score and the corresponding alignment position
    return np.max(scores), np.argmax(scores), max_overhang


def calculate_shared_information(ppm1, ppm2):
    '''
    Calculate shared information between two PPMs.

    Parameters
    ----------
    ppm1, ppm2 : numpy arrays
        motif position probability matrices; first dimension should represent DNA bases, second dimension should represent positions

    Returns
    -------
    si : numpy array
        the shared information at each position
    '''
    p_xy = ppm1*ppm2
    si = -2*np.sum(p_xy * np.log2(p_xy, out = np.zeros_like(p_xy, dtype = np.float64), where = (p_xy!=0)), axis = 0)
    return 2 - si


def find_high_info_bounds(positional_ic, high_info_cutoff = 1):
    '''
    Get the left and right bounds of the high information window in a motif.
    If no positions reach the high information cutoff, the bounds will both be returned as 0.

    Parameters
    ----------
    positional_ic : numpy array
        information content at each position of the motif
    high_info_cutoff : int
        the minimum bits to qualify as "high-information" (default 1 bit)

    Returns
    -------
    left_bound : int
        the first high-information position
    right_bound : int
        the last high-information position
    '''
    if np.max(positional_ic) >= 1:
        left_bound = 0
        for pos in range(positional_ic.shape[0]):
            info = positional_ic[pos]
            if info >= high_info_cutoff:
                left_bound = pos
                break
        right_bound = positional_ic.shape[0]
        for pos in range(positional_ic.shape[0] - 1, -1, -1):
            info = positional_ic[pos]
            if info >= high_info_cutoff:
                right_bound = pos
                break
    else:
        left_bound = 0
        right_bound = 0
    return left_bound, right_bound


def select_representatives(initial_ordered_names, initial_similarities, merge_map, motif_dict):
    '''
    Select a representative motif for each merged set.
    
    Parameters
    ----------
    initial_ordered_names : numpy array
        list of motif names corresponding to the order of indices in the initial_similarities array
    initial_similarities: numpy array
        pairwise similarity scores for all input motifs
    merge_map : dict
        mapping of merged motifs to original motif names

    Returns
    -------
    remap : dict
        a re-mapping of merged motif names to selected representatives from the original motifs
    '''
    merged_names = [ i for i in merge_map.keys() if 'merge_' in i ]
    # make a new mapping with representatives
    remap = {}
    for merge in merged_names:
        # get original motifs in the cluster
        cluster_motifs = np.sort(np.array(merge_map[merge]))
        if len(cluster_motifs) > 2: # find medioid for clusters with more than 2 motifs
            # get their indices
            indices = np.searchsorted(initial_ordered_names, cluster_motifs)
            # subset the similarities and get the medioid
            clust_sim = initial_similarities[indices][:, indices]
            medioid_idx = np.argmin(np.sum(clust_sim, axis = 0))
            rep_name = cluster_motifs[medioid_idx]
        else: # select the motif with the greatest core information content
            info_contents = []
            for motif in cluster_motifs:
                ic = motif_dict[motif]['logo']
                pos_ic = np.sum(ic, axis = 0)
                left_bound, right_bound = find_high_info_bounds(pos_ic)
                core_ic = np.sum(pos_ic[left_bound:right_bound+1])
                info_contents.append(core_ic)
            rep_name = cluster_motifs[np.argmax(info_contents)]
        remap[merge] = rep_name
    return remap