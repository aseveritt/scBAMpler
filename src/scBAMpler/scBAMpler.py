import argparse, sys, shlex
from scBAMpler import build_dict, perform_sampling, generateBAM, CellSim_make, CellSim_mix, CellSim_extract

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
                                choices=["cells", "reads", "frip"]) #"peakreads" is not implemented in perform_sampling; re-add here when it is
    parser_sampler.add_argument('--downsample_to', dest="value",
                                help='Target value for the downsampling operation. ', required=True, type=float)
    parser_sampler.add_argument('--seed', 
                                help='Random seed for reproducibility.', required=False, default=42, type=int)
    parser_sampler.add_argument('--nproc', 
                                help='Number of processors to use.', required=True, type=int, default=1)
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
    parser_justBAM.add_argument('--nproc', 
                                help='Number of processors to use.', required=True, type=int, default=1)

    parser_justBAM.add_argument('-v', '--verbose',
                                help='Print update messages.', required=False, action="store_true")



    # Subcommand: make-pseudobulks
    parser_pseudo = subparsers.add_parser("make-pseudobulks", 
                                           help="Generate pseudo-bulk ATAC-seq profiles by constrained k-means for future population mixing.")
    parser_pseudo.add_argument("--input", 
                               required=True, metavar="FILE",
                               help="Path to input HDF5 file (peakmat_input.h5)")
    parser_pseudo.add_argument("--output", 
                               required=True, metavar="FILE", 
                               help="Path for output pickle file (e.g. medoids_s5000.pickle)")
    parser_pseudo.add_argument("--dimred", 
                               default="umap", choices=["umap", "tsne"], 
                               help="Embedding to use for clustering")
    parser_pseudo.add_argument("--label-col", 
                               default="CellLine", metavar="COL", 
                               help="Name of the grouping column in the embedding (column 3 of the H5 embedding table)")
    parser_pseudo.add_argument("--cluster-size", 
                               type=int, default=500, metavar="N",
                               help="Target number of cells per pseudo-bulk cluster")
    parser_pseudo.add_argument("--nproc",
                               type=int, default=8, metavar="N",
                               help="Number of parallel processes for clustering")
    parser_pseudo.add_argument("--seed",
                               type=int, default=42, metavar="N",
                               help="Random seed for reproducibility")

    # Subcommand: mix-pseudobulks
    parser_mix = subparsers.add_parser("mix-pseudobulks", 
                                       help="Generate and score random pseudo-bulk combinations")
    parser_mix.add_argument("--input", 
                        required=True, metavar="FILE",
                        help="Pickle file from make-pseudobulks")
    parser_mix.add_argument("--output", 
                        required=True, metavar="FILE",
                        help="Output CSV file path")
    parser_mix.add_argument("--groups", 
                        nargs="+", default=["all"], metavar="LABEL",
                        help="Labels to restrict cluster pool to (dominant label per cluster), "
                        "or 'all' for unbiased sampling. E.g.: --groups K562 HEPG2")
    parser_mix.add_argument("--n-combos", 
                        type=int, default=2000, metavar="N",
                        help="Number of random combinations to generate per ft-size")
    parser_mix.add_argument("--cluster-size", 
                        type=int, default=500, metavar="N",
                        help="Cells per cluster (from make-pseudobulks), used to compute r=ft_size/cluster_size")
    parser_mix.add_argument("--ft-sizes", 
                        type=int, nargs="+", metavar="N", default=[500, 1000, 2000, 5000, 10000, 15000, 20000],
                        help="Target footprint sizes in cells. One round of sampling per size.")
    parser_mix.add_argument("--label-col", 
                        default="CellLine", metavar="COL",
                        help="Name of the grouping column in embedding_df")
    parser_mix.add_argument("--nproc", 
                        type=int, default=8, metavar="N",
                        help="Number of parallel scoring processes")
    parser_mix.add_argument("--seed", 
                        type=int, default=42,
                        help="Random seed for reproducibility")
    

    # Subcommand: extract-populations
    parser_extract = subparsers.add_parser("extract-populations",
                                       help="Write bash scripts that extract each selected population into a BAM")
    parser_extract.add_argument("--barcode-dir",
                        required=True, metavar="DIR",
                        help="Directory of barcode files named combo_<ID>.<label>.barcodes.csv "
                             "(as written by notebooks/inspect-select-combos.ipynb)")
    parser_extract.add_argument("--output",
                        required=True, metavar="DIR",
                        help="Output directory for scripts and BAMs")
    parser_extract.add_argument("--bam",
                        required=True, action="append", metavar="LABEL=PATH",
                        help="BAM for one label, repeatable. One is required for every label "
                             "present in the barcode files. "
                             "e.g. --bam HEPG2=test_data/HEPG2_subset.bam")
    parser_extract.add_argument("--combo-ids",
                        type=int, nargs="+", default=None, metavar="ID",
                        help="Only write scripts for these population IDs. Default: all found.")
    parser_extract.add_argument("--nproc",
                        type=int, default=8, metavar="N",
                        help="Processors passed to sinto and samtools")
    parser_extract.add_argument("--output_fragment",
                        action="store_true",
                        help="Also produce a .frags.tsv.bgz per population")
    
    
    args = parser.parse_args(cleaned_args)
        
    if args.command == "create-dictionary":
        build_dict.main(args)
        
    elif args.command == "sampler":
        perform_sampling.main(args)
        
    elif args.command == "generateBAM":
        generateBAM.main(args)
        
    elif args.command == "make-pseudobulks":
        CellSim_make.main(args)
        
    elif args.command == "mix-pseudobulks":
        CellSim_mix.main(args)
        
    elif args.command == "extract-populations":
        CellSim_extract.main(args)
        
        
        

if __name__ == "__main__":
    main()
