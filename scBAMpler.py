import argparse
from build_dict import main as main_BuildInput
from perform_sampling import main as main_Sampling

def main():
    parser = argparse.ArgumentParser(prog="scBAMpler")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: create-dictionary
    parser_buildInput = subparsers.add_parser("create-dictionary", help="Create Input Dictionary for scBAMpler:sampling")
    parser_buildInput.add_argument('-b', '--bam_file', 
                                   help='path to bam file', required=True, type=str)
    parser_buildInput.add_argument('-p', '--peak_file', 
                                   help='path to bed file of called peaks', required=True, type=str)
    parser_buildInput.add_argument('-o', '--output_file', 
                                   help='name of output pickle file', required=True, type=str)
    parser_buildInput.add_argument('-i', '--intersect_file', 
                                   help='optional, provide mapping of peak-to-reads', required=False, default=None, type=str)
    parser_buildInput.add_argument('--delete_intersect', 
                                   help='delete output of bedtools', required=False, action="store_true")
    
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

    args = parser.parse_args()

    if args.command == "create-dictionary":
        main_BuildInput(args)
    elif args.command == "sampler":
        main_Sampling(args)

if __name__ == "__main__":
    main()
