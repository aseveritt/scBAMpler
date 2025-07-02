import argparse, sys
from build_dict import main as main_BuildInput
from perform_sampling import main as main_Sampling

def main(argv=None):

    #Clean up arguments in case they're passed with tabs or newlines on accident.
    #section via ChatGPT
    raw_args = sys.argv[1:] if argv is None else argv
    cleaned_args = []
    for arg in raw_args:
        if '\t' in arg or '\n' in arg or ' ' in arg:
            cleaned_args.extend(shlex.split(arg))
        else:
            cleaned_args.append(arg)

    parser = argparse.ArgumentParser(prog="scBAMpler")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: create-dictionary
    parser_buildInput = subparsers.add_parser("create-dictionary", help="Create Input Dictionary for scBAMpler:sampling")
    parser_buildInput.add_argument('-b', '--bam_file', 
                                   help='Path to the coordinate sorted, input BAM file.', required=True, type=str)
    parser_buildInput.add_argument('-p', '--peak_file', 
                                   help='Path to peak file in BED6 format', required=True, type=str)
    parser_buildInput.add_argument('-o', '--output_file', 
                                   help='Name and location of final cell type dictionary stored as a pickle file', required=True, type=str)
    parser_buildInput.add_argument('-i', '--intersect_file', 
                                   help='optional, if the cell barcode:peak read mapping (*.bed.gz) already exists, you can provide it here to skip that initial step.', required=False, default=None, type=str)
    parser_buildInput.add_argument('--delete_intersect', 
                                   help='By default, cell barcode:peak read mapping file is saved. Setting this flag will delete it.', required=False, action="store_true")
    
    # Subcommand: sampler
    parser_sampler = subparsers.add_parser("sampler", 
                                           help="Downsample input BAM as instructed")
    parser_sampler.add_argument('-i', '--input_file', 
                                help='name of input pickle file', required=True, type=str)
    parser_sampler.add_argument('-o', '--output_file', 
                                help='name of output file', required=True, type=str)
    parser_sampler.add_argument('-e', '--edit', 
                                help='type of downsampling to be performed', required=True, type=str, 
                                choices=["cells", "reads", "frip", "peakreads"])
    parser_sampler.add_argument('-v', '--value', 
                                help='specified value (e.g. ncells, nreads, frip)', required=True, type=float)
    parser_sampler.add_argument('--seed', 
                                help='type of downsampling to be performed', required=False, default=42, type=int)
    parser_sampler.add_argument('--nproc', 
                                help='type of downsampling to be performed', required=True, type=int)
    parser_sampler.add_argument('-b', '--bam_file', 
                                help='name of original bam file', required=True, type=str)

    args = parser.parse_args(cleaned_args)
        
    if args.command == "create-dictionary":
        main_BuildInput(args)
        
    elif args.command == "sampler":
        main_Sampling(args)

if __name__ == "__main__":
    main()
