#perform_sampling.py
import os, pysam, time, pickle, shutil
import downsampling_functions as dsfs
from datetime import timedelta
from downsampling_functions import Cells

def main(args):    

    ########################
    #USER CHECKS
    def validateFile(arg):
        if not os.path.isfile(arg): print(f'ERROR: The file "{arg}" does not exist!')
        else: return 
    def validateTools(tool_name):
        if shutil.which(tool_name) is None:
            print(f"ERROR: Required tool '{tool_name}' not found in PATH.")
            sys.exit(1)

    validateFile(args.input_pickle)
    validateTools("samtools"); validateTools("bedtools")
    
    out_dir, _ = os.path.split(args.output_prefix)
    if out_dir == "": out_dir = "."
    if not os.path.exists(out_dir): os.mkdir(out_dir)
    
    ########################
    start_time = time.time()
    with open(args.input_pickle, "rb") as f:
        cb_dict, cb_encoder, qname_encoder = pickle.load(f)
    if args.verbose: print("--- %s h:m:s to read in pickle ---" % str(timedelta(seconds=(time.time() - start_time))))

    ########################
    rp_pre = dsfs.TotalReadPairs(cb_dict)
    cells_pre = len(cb_dict.keys())   

    ########################
    if args.edit == "cells":
        cb_dict = dsfs.DownsampleCells(cb_dict, N_cells=args.value, seed=args.seed, verbose=args.verbose)
    elif args.edit == "reads":
        dsfs.DownsampleReads(cb_dict, N_desired_reads=args.value, seed=args.seed, verbose=args.verbose)
    elif args.edit == "frip":
        dsfs.DownsampleFRIP(cb_dict, frip=args.value, seed=args.seed, verbose=args.verbose)
    else:
        print("ERROR: Invalid downsampling type"); sys.exit(1)

    read_file = args.output_prefix+".txt"
    dsfs.OutputDict(cb_dict, qname_encoder, read_file, verbose=args.verbose)

    ########################
    rp_post = dsfs.TotalReadPairs(cb_dict)
    cells_post = len(cb_dict.keys())   
    
    dsfs.WriteLog(args.output_prefix+".summary.txt", 
                  {'input_file':args.input_pickle, 'sampling_type':args.edit, 'value': args.value, 'seed':args.seed}, 
                  {'N_cells_removed': cells_pre-cells_post, 'N_reads_removed':rp_pre-rp_post}, 
                  dsfs.Summary(cb_dict, output_as = "dict"), verbose=args.verbose)

    output_bam = args.output_prefix+".bam"
    returncode = dsfs.GenerateOutputBam(args.input_bam, read_file, args.nproc, output_bam, verbose=args.verbose)
    
    if (returncode == 0 and args.output_fragment == True):
        validateTools("sinto")
        dsfs.GenerateOuputFragment(output_bam, args.output_prefix+".frags.tsv.bgz", args.nproc, verbose=args.verbose)

    
if __name__ == '__main__':
    main()

