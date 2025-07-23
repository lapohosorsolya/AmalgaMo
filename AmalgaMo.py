import os, sys, getopt
import numpy as np
from motifs import *
from tqdm import tqdm
import json

'''
USAGE NOTES

Required Inputs
---------------
    -i (input directory)
        the directory containing the input motif matrices (individual files)
    -o (output directory)
        the directory where the output should be written; will make this if it doesn't exist
    -f (motif format code)
        the file format of the input motifs; one of ("hocomoco", "jaspar", "meme")

Optional Inputs
---------------
    -s (float between 0.5 and 1)
        the similarity cutoff for merging; default is 0.9
    -a (float between 0.5 and 1)
        the minimum alignment overlap; default is 0.9
    -m (int)
        the maximum allowed length difference between two merged motifs; default is 2
    -t (float between 0 and 1)
        the total information ratio requirement for merging; default is 0.85
    -r (int)
        the difference in number of positions that can be tolerated when comparing high-information windows; default is 0

'''


def main(argv):
    try:
        opts, _ = getopt.getopt(argv, 'm:t:r:s:i:a:o:f:')
    except getopt.GetoptError:
        print('\n::: Error: cannot parse command line inputs')
        sys.exit(2) 
    for opt, arg in opts:
        if opt == '-m':
            global max_allowed_len_diff
            max_allowed_len_diff = int(arg)
        elif opt == '-t':
            global tir_cutoff
            tir_cutoff = float(arg)
        elif opt == '-r':
            global window_relax
            window_relax = int(arg)
        elif opt == '-i':
            global input_dir
            input_dir = arg
        elif opt == '-s':
            global sim_cutoff
            sim_cutoff = float(arg)
        elif opt == '-a':
            global min_overlap
            min_overlap = float(arg)
        elif opt == '-o':
            global output_dir
            output_dir = arg
        elif opt == '-f':
            global form
            form = arg


if __name__ == "__main__":

    ext = { 'hocomoco': 'pfm', 'jaspar': 'jaspar', 'meme': 'meme' }

    # default parameters
    sim_cutoff = 0.9
    min_overlap = 0.9
    tir_cutoff = 0.85
    max_allowed_len_diff = 2
    window_relax = 0
    form = None
    main(sys.argv[1:])

    # check if motif format code is valid
    print('\nMotif file format:')
    if form not in ['hocomoco', 'jaspar', 'meme']:
        print('Please indicate the format of the input motif matrices. Valid formats include:\n\t- hocomoco\n\t- jaspar\n\t- meme')
        sys.exit(2)
    else:
        print('\t{}'.format(form))

    # check if parameters are valid
    print('\nParameters:')
    if sim_cutoff >= 1 or sim_cutoff < 0.5:
        print('Please select a similarity score cutoff between 0.5 and 1.')
        sys.exit(2)
    else:
        print('\ts = {}'.format(sim_cutoff))
    if min_overlap >= 1 or min_overlap < 0.5:
        print('Please select a minimum alignment overlap between 0.5 and 1.')
        sys.exit(2)
    else:
        print('\ta = {}'.format(min_overlap))
    if tir_cutoff >= 1 or tir_cutoff < 0:
        print('Please select a total information ratio cutoff between 0 and 1.')
        sys.exit(2)
    else:
        print('\tt = {}'.format(tir_cutoff))
    if max_allowed_len_diff < 0:
        print('Please select a maximum length difference greater than or equal to 0.')
        sys.exit(2)
    else:
        print('\tm = {}'.format(max_allowed_len_diff))
    if window_relax < 0:
        print('Please select a maximum core length difference greater than or equal to 0.')
        sys.exit(2)
    else:
        print('\tr = {}'.format(window_relax))

    # check if input and output directories exist
    print('\nPaths:')
    if os.path.isdir(input_dir) == False:
        print('Input directory does not exist: {}'.format(input_dir))
        sys.exit(2)
    else:
        print('\tinput directory = {}'.format(input_dir))
    if os.path.isdir(output_dir) == False:
        os.mkdir(output_dir)
    print('\toutput directory = {}'.format(output_dir))

    # get motifs
    print('\nReading motifs and calculating information content...')
    files = sorted(os.listdir(input_dir))
    n_motifs = len(files)
    motif_dict = {}
    names = []
    for file in files:
        name, motif = read_motif(os.path.join(input_dir, file), form)
        motif_dict[name] = { 'pfm': motif }
        names.append(name)
    for motif in motif_dict.keys():
        arr = motif_dict[motif]['pfm']
        motif_dict[motif]['length'] = arr.shape[1]
        ppm = arr / np.sum(arr, axis = 0)
        motif_dict[motif]['ppm'] = ppm
        ic = make_info_content_logo(ppm)
        motif_dict[motif]['logo'] = ic
        left_bound, right_bound = find_high_info_bounds(np.sum(ic, axis = 0))
        motif_dict[motif]['high_info_region_length'] = right_bound - left_bound
        motif_dict[motif]['total_info'] = np.sum(ic)

    # calculate matrix of total information ratios
    tirs = np.zeros((n_motifs, n_motifs))
    for i in tqdm(range(n_motifs), total = n_motifs):
        total_info_1 = motif_dict[names[i]]['total_info']
        for j in range(i):
            total_info_2 = motif_dict[names[j]]['total_info']
            r = min(total_info_1, total_info_2) / max(total_info_1, total_info_2)
            tirs[i, j] = r
            tirs[j, i] = r

    # best pairwise similarities; mirrored across diagonal
    pairwise_similarities = np.zeros((n_motifs, n_motifs))

    # best alignment orientations for merging; 1 if the forward versions align, -1 if one motif needs to be reverse complemented; mirrored across diagonal
    alignment_orientation = np.zeros((n_motifs, n_motifs))

    # best alignment positions for merging; alignment position is always relative to the longer motif, according to the saved alignment orientation
    # filled in such that the alignment of "antisense" motif pairs uses the forward version of the motif indexed by the row
    alignment_position = np.zeros((n_motifs, n_motifs), dtype = int)
    alignment_position[:] = -1

    # overhangs for alignment positions; 
    pairwise_overhangs = np.zeros((n_motifs, n_motifs), dtype = int)

    # need motif lengths for alignment position matrix
    motif_lengths = np.zeros(n_motifs, dtype = int)

    # calculate all pairwise similarities
    print('\nCalculating initial pairwise similarities...')
    for i in tqdm(range(n_motifs), total = n_motifs):

        # this is the "reference" motif - only need to use positive version of this one
        pos_motif_1 = np.copy(motif_dict[names[i]]['ppm'])
        l1 = pos_motif_1.shape[1]
        motif_lengths[i] = l1
        hir_1 = motif_dict[names[i]]['high_info_region_length']

        for j in range(i):

            l2 = motif_lengths[j]
            hir_2 = motif_dict[names[j]]['high_info_region_length']

            # only do all this if the length difference and total information ratio qualifies
            if abs(hir_1 - hir_2) <= window_relax and abs(l1 - l2) <= max_allowed_len_diff and tirs[i, j] > tir_cutoff:

                # read another motif
                pos_motif_2 = np.copy(motif_dict[names[j]]['ppm'])
                neg_motif_2 = get_motif_reverse_complement(pos_motif_2)

                # find out which orientation aligns best with the reference motif
                pos_neg_sims = np.zeros(2)
                pos_neg_alignment = np.zeros(2)
                pos_neg_overhang = np.zeros(2)
                for d, motif in enumerate([pos_motif_2, neg_motif_2]):
                    score, position, overhang = get_siwcs_score(pos_motif_1, motif, min_overlap = min_overlap)
                    pos_neg_sims[d] = score
                    pos_neg_alignment[d] = position
                    pos_neg_overhang[d] = overhang

                # fill in values
                if pos_neg_sims[0] >= pos_neg_sims[1]:
                    pairwise_similarities[i, j] = pairwise_similarities[j, i] = pos_neg_sims[0]
                    alignment_orientation[i, j] = alignment_orientation[j, i] = 1
                    # alignment position for a pair of forward motifs can be the same at [i,j] and [j,i], always relative to the longer motif
                    x = pos_neg_alignment[0]
                    o = pos_neg_overhang[0]
                    alignment_position[i, j] = x
                    pairwise_overhangs[i, j] = pairwise_overhangs[j, i] = o
                    if l1 == l2:
                        alignment_position[j, i] = 2 * o - x
                    else:
                        alignment_position[j, i] = x
                else:
                    pairwise_similarities[i, j] = pairwise_similarities[j, i] = pos_neg_sims[1]
                    alignment_orientation[i, j] = alignment_orientation[j, i] = -1
                    # motif j is reverse complemented to the alignment positions are different at [i,j] and [j,i], always relative to the longer motif
                    x = pos_neg_alignment[1]
                    o = pos_neg_overhang[1]
                    alignment_position[i, j] = x
                    pairwise_overhangs[i, j] = pairwise_overhangs[j, i] = o
                    if l1 == l2:
                        alignment_position[j, i] = x
                    else:
                        alignment_position[j, i] = abs(l1 - l2) + 2 * o - x

    # save initial data
    np.save(os.path.join(output_dir, 'AmalgaMo_initial_similarities.npy'), pairwise_similarities)
    np.save(os.path.join(output_dir, 'AmalgaMo_initial_names.npy'), np.array(names))

    # iterative merging
    print('\nIteratively merging motifs...')

    # when motifs i and j are merged:
        # add their names to the end of an ordered list and associate the list with a new key in the motif dictionary
        # make a new pfm by summing the merged PFMs (using their alignment orientation and position)
        # save the merged pfm to the dictionary
        # delete data on motifs i and j from all arrays
        # then, add a new row and column for the merged motif with newly calculated similarity scores, etc.

    # list of motif names to be updated during iterations; indexing will always match arrays
    ordered_motif_names = names.copy()

    # dictionary of PFMs
    for motif in motif_dict:
        motif_dict[motif]['merge'] = []

    # get the maximum similarity score
    max_score = pairwise_similarities.max()
    max_score_by_iter = [max_score]

    count = 1

    # iterative merging
    while max_score > sim_cutoff:

        # find the most similar pair of motifs
        i1, i2 = np.unravel_index(pairwise_similarities.argmax(), pairwise_similarities.shape)
        idx1, idx2 = max(i1, i2), min(i1, i2)
        
        # get the names of these motifs
        name1, name2 = ordered_motif_names[idx1], ordered_motif_names[idx2]
        l1, l2 = motif_dict[name1]['length'], motif_dict[name2]['length']

        # find out how many PPMs have already been merged to obtain these two, and get their names
        prev_merges_1 = len(motif_dict[name1]['merge'])
        prev_merges_2 = len(motif_dict[name2]['merge'])
        if prev_merges_1 == 0:
            merge_names_1 = [name1]
        else:
            merge_names_1 = motif_dict[name1]['merge']
        if prev_merges_2 == 0:
            merge_names_2 = [name2]
        else:
            merge_names_2 = motif_dict[name2]['merge']

        # get the PPMs
        ppm1 = np.copy(motif_dict[name1]['ppm'])
        align_orient = alignment_orientation[idx1, idx2]
        if align_orient == 1:
            ppm2 = np.copy(motif_dict[name2]['ppm'])
        else:
            ppm2 = get_motif_reverse_complement(np.copy(motif_dict[name2]['ppm']))

        # pad motifs with overhang allowance
        o = pairwise_overhangs[idx1, idx2]
        x = alignment_position[idx1, idx2]
        if l1 >= l2:
            # pad the longer motif
            padded1 = np.zeros((4, l1 + 2 * o))
            padded1[:] = 0.25
            padded1[:, o:o+l1] = ppm1
            # pad the shorter motif
            padded2 = np.zeros_like(padded1)
            padded2[:] = 0.25
            padded2[:, x:x+l2] = ppm2
        else:
            # pad the longer motif
            padded2 = np.zeros((4, l2 + 2 * o))
            padded2[:] = 0.25
            padded2[:, o:o+l2] = ppm2
            # pad the shorter motif
            padded1 = np.zeros_like(padded2)
            padded1[:] = 0.25
            padded1[:, x:x+l1] = ppm1

        # make a merged PPM, using the number of previous merges to weigh the average
        merged_ppm = ((prev_merges_1 + 1)*padded1 + (prev_merges_2 + 1)*padded2) / (prev_merges_1 + prev_merges_2 + 2)

        # trim the merged PPM
        n_unique = np.array([ np.unique(merged_ppm[:, col]).shape[0] for col in range(merged_ppm.shape[1]) ])
        uniform_idx = np.flip(np.where(n_unique == 1)[0])
        for idx in uniform_idx:
            if idx >= o + max(l1, l2) or idx < o:
                merged_ppm = np.delete(merged_ppm, idx, 1)
            
        # add the merged PPM to the dictionary
        merged_original_motif_names = merge_names_1 + merge_names_2
        merged_length = merged_ppm.shape[1]
        new_name = '-'.join(merged_original_motif_names)
        if new_name in motif_dict.keys():
            new_name = 'merge_{}'.format(i)
        ic = make_info_content_logo(merged_ppm)
        merged_tic = np.sum(ic)
        left_bound, right_bound = find_high_info_bounds(np.sum(ic, axis = 0))
        merged_hir = right_bound - left_bound
        motif_dict[new_name] = { 'merge': merged_original_motif_names, 'length': merged_length, 'ppm': np.copy(merged_ppm), 'logo': ic, 'total_info': merged_tic, 'high_info_region_length': merged_hir }
        ordered_motif_names.append(new_name)

        # delete data from arrays, starting with the greater index
        sorted_idx = sorted([idx1, idx2], reverse = True)
        for idx in sorted_idx:
            tirs = np.delete(tirs, idx, 0)
            tirs = np.delete(tirs, idx, 1)
            pairwise_similarities = np.delete(pairwise_similarities, idx, 0)
            pairwise_similarities = np.delete(pairwise_similarities, idx, 1)
            alignment_orientation = np.delete(alignment_orientation, idx, 0)
            alignment_orientation = np.delete(alignment_orientation, idx, 1)
            alignment_position = np.delete(alignment_position, idx, 0)
            alignment_position = np.delete(alignment_position, idx, 1)
            pairwise_overhangs = np.delete(pairwise_overhangs, idx, 0)
            pairwise_overhangs = np.delete(pairwise_overhangs, idx, 1)
            del ordered_motif_names[idx]

        # calculate new data for merged motif
        new_tirs_row = np.zeros(len(ordered_motif_names) - 1)
        new_tirs_col = np.zeros(len(ordered_motif_names))
        new_sim_row = np.zeros(len(ordered_motif_names) - 1)
        new_sim_col = np.zeros(len(ordered_motif_names))
        new_orient_row = np.zeros(len(ordered_motif_names) - 1, dtype = int)
        new_orient_col = np.zeros(len(ordered_motif_names), dtype = int)
        new_pos_row = np.zeros(len(ordered_motif_names) - 1, dtype = int)
        new_pos_col = np.zeros(len(ordered_motif_names), dtype = int)
        new_pos_row[:] = -1
        new_pos_col[:] = -1
        new_overhang_row = np.zeros(len(ordered_motif_names) - 1, dtype = int)
        new_overhang_col = np.zeros(len(ordered_motif_names), dtype = int)

        for j, motif_2 in enumerate(ordered_motif_names[:-1]):

            # get info content anc calculate total information ratio
            total_info_2 = motif_dict[motif_2]['total_info']
            r = min(merged_tic, total_info_2) / max(merged_tic, total_info_2)
            new_tirs_row[j] = r
            new_tirs_col[j] = r
            motif_2_len = motif_dict[motif_2]['length']
            hir_2 = motif_dict[motif_2]['high_info_region_length']

            # check if the length difference and total information ratio qualify for alignment and calculate similarity
            if abs(merged_hir - hir_2) <= window_relax and abs(merged_length - motif_2_len) <= max_allowed_len_diff and r > tir_cutoff:

                # read motif ppm
                pos_motif_2 = np.copy(motif_dict[motif_2]['ppm'])
                neg_motif_2 = get_motif_reverse_complement(np.copy(pos_motif_2))

                # find out which orientation aligns best with the reference motif
                pos_neg_sims = np.zeros(2)
                pos_neg_alignment = np.zeros(2)
                pos_neg_overhang = np.zeros(2)
                for d, motif in enumerate([pos_motif_2, neg_motif_2]):
                    score, position, overhang = get_siwcs_score(merged_ppm, motif, min_overlap = min_overlap)
                    pos_neg_sims[d] = score
                    pos_neg_alignment[d] = position
                    pos_neg_overhang[d] = overhang

                # fill in values
                if pos_neg_sims[0] >= pos_neg_sims[1]:
                    new_sim_row[j] = new_sim_col[j] = pos_neg_sims[0]
                    new_orient_row[j] = new_orient_col[j] = 1
                    # alignment position for a pair of forward motifs can be the same at [i,j] and [j,i], always relative to the longer motif
                    x = pos_neg_alignment[0]
                    o = pos_neg_overhang[0]
                    new_pos_row[j] = x
                    new_overhang_row[j] = new_overhang_col[j] = o
                    if l1 == l2:
                        new_pos_col[j] = 2 * o - x
                    else:
                        new_pos_col[j] = x
                else:
                    new_sim_row[j] = new_sim_col[j] = pos_neg_sims[1]
                    new_orient_row[j] = new_orient_col[j] = -1
                    # motif j is reverse complemented to the alignment positions are different at [i,j] and [j,i], always relative to the longer motif
                    x = pos_neg_alignment[1]
                    o = pos_neg_overhang[1]
                    new_pos_row[j] = x
                    new_overhang_row[j] = new_overhang_col[j] = o
                    if l1 == l2:
                        new_pos_col[j] = x
                    else:
                        new_pos_col[j] = abs(l1 - l2) + 2 * o - x

        tirs = np.concatenate([tirs, new_tirs_row.reshape((-1, len(ordered_motif_names)-1))], axis = 0)
        tirs = np.concatenate([tirs, new_tirs_col.reshape(len(ordered_motif_names), 1)], axis = 1)

        pairwise_similarities = np.concatenate([pairwise_similarities, new_sim_row.reshape((-1, len(ordered_motif_names)-1))], axis = 0)
        pairwise_similarities = np.concatenate([pairwise_similarities, new_sim_col.reshape(len(ordered_motif_names), 1)], axis = 1)

        alignment_orientation = np.concatenate([alignment_orientation, new_orient_row.reshape((-1, len(ordered_motif_names)-1))], axis = 0)
        alignment_orientation = np.concatenate([alignment_orientation, new_orient_col.reshape(len(ordered_motif_names), 1)], axis = 1)

        alignment_position = np.concatenate([alignment_position, new_pos_row.reshape((-1, len(ordered_motif_names)-1))], axis = 0)
        alignment_position = np.concatenate([alignment_position, new_pos_col.reshape(len(ordered_motif_names), 1)], axis = 1)

        pairwise_overhangs = np.concatenate([pairwise_overhangs, new_overhang_row.reshape((-1, len(ordered_motif_names)-1))], axis = 0)
        pairwise_overhangs = np.concatenate([pairwise_overhangs, new_overhang_col.reshape(len(ordered_motif_names), 1)], axis = 1)

        max_score = pairwise_similarities.max()
        max_score_by_iter.append(max_score)

        print(count, end = '\r')
        count += 1

    # assess result
    all_merges = []
    max_merge = 0
    for motif in ordered_motif_names:
        mergelist = motif_dict[motif]['merge']
        n = len(mergelist)
        if n > 0:
            all_merges.append(motif)
            if n > max_merge:
                max_merge = n
    print('\n\nGenerated {} merged motifs, with the largest merged set consisting of {} original motifs'.format(len(all_merges), max_merge))
    print('Original motifs left unmerged: {}'.format(len(ordered_motif_names) - len(all_merges)))
    print('Resulting non-redundant motifs: {}'.format(len(ordered_motif_names)))

    # save results
    print('\nSaving results...')
    output_dict = {}
    output_ppms = {}
    merge_counter = 1
    for name in ordered_motif_names:
        entry = motif_dict[name]
        names = entry['merge']
        if len(names) > 0:
            new_name = 'merge_{}'.format(merge_counter)
            merge = names
            merge_counter += 1
        else:
            new_name = name
            merge = [name]
        output_dict[new_name] = merge
        output_ppms[new_name] = entry['ppm']

    with open(os.path.join(output_dir, 'AmalgaMo_results.json'), 'w') as f:
        json.dump(output_dict, f, indent = 4)

    np.savez(os.path.join(output_dir, 'AmalgaMo_PPMs.npz'), **output_ppms)

    np.save(os.path.join(output_dir, 'AmalgaMo_final_similarities.npy'), pairwise_similarities)
    np.save(os.path.join(output_dir, 'AmalgaMo_final_names.npy'), np.array(ordered_motif_names))

    param_dict = { 'tir_cutoff': tir_cutoff, 'max_length_diff': max_allowed_len_diff, 'max_core_length_diff': window_relax, 'min_overlap': min_overlap, 'sim_cutoff': sim_cutoff, 'input_dir': input_dir, 'output_dir': output_dir, 'motif_format': form }
    with open(os.path.join(output_dir, 'AmalgaMo_params.json'), 'w') as f:
        json.dump(param_dict, f, indent = 4)

    print('\nSaving merged motif logos...')
    merged_motif_names = [ i for i in output_ppms.keys() if 'merge_' in i ]
    ppm_dir = os.path.join(output_dir, 'merged_ppms')
    os.mkdir(ppm_dir)
    logo_dir = os.path.join(output_dir, 'merged_logos')
    os.mkdir(logo_dir)
    for name in tqdm(merged_motif_names, total = len(merged_motif_names)):
        write_motif(os.path.join(ppm_dir, '{}.{}'.format(name, ext[form])), name, output_ppms[name], form)
        write_logo(os.path.join(logo_dir, '{}.pdf'.format(name)), output_ppms[name])

    print('\nSaved results to {}\n'.format(output_dir))

