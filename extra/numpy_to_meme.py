import os, sys, getopt
import numpy as np

'''
USAGE NOTES

Required Inputs
---------------
    -i (input file)
        the numpy archive (.npz) containing the input motif matrices or the input directory containing individual (.npy) files, or simply a numpy file (.npy)
    -o (output file)
        the full path and new name of the output file

Optional Inputs
---------------
    -s (subset file)
        path to a file containing the names of motifs to include in the output; if provided, it will be used to filter the output, if not provided, all motifs will be included

'''


def main(argv):
    try:
        opts, _ = getopt.getopt(argv, 'i:o:s:')
    except getopt.GetoptError:
        print('\n::: Error: cannot parse command line inputs')
        sys.exit(2) 
    for opt, arg in opts:
        if opt == '-i':
            global input_file
            input_file = arg
        elif opt == '-o':
            global output_file
            output_file = arg
        elif opt == '-s':
            global subset_file
            subset_file = arg


if __name__ == "__main__":

    # by default, all motifs are in the output
    subset_file = None
    main(sys.argv[1:])

    # check if input and output directories exist
    if os.path.isfile(input_file) == False and os.path.isdir(input_file) == False:
        sys.exit(2)

    # get motifs
    if os.path.isfile(input_file):
        motifs = np.load(input_file)
        if isinstance(motifs, np.ndarray):
            names = ['.'.join(input_file.split('/')[-1].split('.')[:-1])]
            motifs = { names[0]: motifs }
        else:
            names = sorted(list(motifs.files))
    else:
        files = sorted(os.listdir(input_file))
        names = [ i.split('.npy')[0] for i in files ]
        motifs = {}
        for i, name in enumerate(names):
            motifs[name] = np.load(os.path.join(input_file, files[i]))

    # subset motifs if needed
    if subset_file is not None:
        with open(subset_file, 'r') as f:
            lines = f.readlines()
            subset_names = [ line.rstrip() for line in lines ]
        names = subset_names
        print('Subsetting motifs to {}'.format(', '.join(names)))
        subset_motifs = { name: motifs[name] for name in names }

    # generate header
    text = 'MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n' 

    # convert motifs
    for name in names:
        arr = motifs[name]
        arr = arr / np.sum(arr, axis = 0)
        unit_text = 'MOTIF {}\nletter-probability matrix: alength= 4 w= {} nsites= 100\n'.format(name, arr.shape[1])
        for i in range(arr.shape[1]):
            probas = arr[:, i].tolist()
            unit_text += '  '.join([ '{0:.6f}'.format(j) for j in probas ]) + '\n'
        unit_text += '\n'
        text += unit_text

    # write file
    with open(output_file, 'w') as f:
        f.write(text)

    print('Done')
