import argparse, sys
from scBAMpler import build_dict, perform_sampling, downsampling_functions, generateBAM

def main(argv=None):

    #Clean up arguments in case they're passed with tabs or newlines on accident.
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
    parser_buildInput.add_argument('-v', '--verbose',
                                   help='Print update messages.', required=False, action="store_true")

    
    # Subcommand: sampler
    parser_sampler = subparsers.add_parser("sampler", 
                                           help="Downsample input BAM as instructed")
    parser_sampler.add_argument('-i', '--input_pickle', 
                                help='Path to the pickle file', required=True, type=str)
    parser_sampler.add_argument('-b', '--input_bam', 
                                help='Path to the coordinate-sorted input BAM file. Reads will be directly extracted from this file.', 
                                required=True, type=str)
    parser_sampler.add_argument('-o', '--output_prefix', 
                                help='Prefix for all output files.', required=True, type=str)
    parser_sampler.add_argument('--downsample_by', dest="edit",
                                help='Type of downsampling to perform.', required=True, type=str, 
                                choices=["cells", "reads", "frip", "peakreads"])
    parser_sampler.add_argument('--downsample_to', dest="value",
                                help='Target value for the downsampling operation. ', required=True, type=float)
    parser_sampler.add_argument('--seed', 
                                help='Random seed for reproducibility.', required=False, default=42, type=int)
    parser_sampler.add_argument('--nproc', 
                                help='Number of processors to use.', required=True, type=int)
    parser_sampler.add_argument('--output_fragment', 
                                help='If set, will also output a `fragment.tsv.bgz` file in addition to the BAM file. ', 
                                required=False, action="store_true")
    parser_sampler.add_argument('-v', '--verbose',
                                help='Print update messages.', required=False, action="store_true")

    

    # Subcommand: generateBAM
    parser_justBAM = subparsers.add_parser("generateBAM", 
                                           help="Given a list of input reads, downsample BAM file accordingly")
    parser_justBAM.add_argument('--input_bam', 
                            help='Path to the coordinate-sorted input BAM file. Reads will be directly extracted from this file.', 
                            required=True, type=str)
    parser_justBAM.add_argument('--output_bam', 
                            help='Desired name of output BAM file', 
                            required=True, type=str)
    parser_justBAM.add_argument('--selected_reads', 
                            help='Plain-file text with new line separated reads to keep.', 
                            required=True, type=str)
    parser_justBAM.add_argument('--output_fragment', 
                            help='If set, will also output a `fragment.tsv.bgz` file in addition to the BAM file. ', 
                            required=False, action="store_true")
    parser_justBAM.add_argument('-v', '--verbose',
                                help='Print update messages.', required=False, action="store_true")

    
    args = parser.parse_args(cleaned_args)
        
    if args.command == "create-dictionary":
        build_dict.main(args)
        
    elif args.command == "sampler":
        perform_sampling.main(args)
        
    elif args.command == "generateBAM":
        main_genBam(args)
        

if __name__ == "__main__":
    main()
