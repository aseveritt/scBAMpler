#perform_sampling.py
import argparse, os, pysam, time, pickle
import downsampling_functions as dsfs
from datetime import timedelta
from downsampling_functions import Cells

def main():    
    def validatefile(arg):
        if not os.path.isfile(arg): parser.error('The file "{}" does not exist!'.format(arg))
        else: return 
    validatefile(args.input_file)
    
    start_time = time.time()
    with open(args.input_file, "rb") as f:
        cb_dict, cb_encoder, qname_encoder = pickle.load(f)
    print("--- %s h:m:s to read in pickle ---" % str(timedelta(seconds=(time.time() - start_time))))
    
    out_dir, filename = os.path.split(args.output_file)
    if out_dir == "": out_dir = "."
    rp_pre = dsfs.TotalReadPairs(cb_dict)
    cells_pre = len(cb_dict.keys())   

    if args.edit == "cells":
        start_time = time.time()
        cb_dict = dsfs.DownsampleCells(cb_dict, N_cells=args.value, seed=args.seed)
        print("--- %s h:m:s to downsample cells ---" % str(timedelta(seconds=(time.time() - start_time))))
        
        start_time = time.time()
        dsfs.OutputDict(cb_dict, cb_encoder, args.output_file, "cells", os.path.splitext(filename)[0])
        print("--- %s h:m:s to write reads to file ---" % str(timedelta(seconds=(time.time() - start_time))))
    
    elif args.edit == "reads":
        start_time = time.time()
        dsfs.DownsampleReads(cb_dict, N_desired_reads=args.value, seed=args.seed)
        print("--- %s h:m:s to downsample reads ---" % str(timedelta(seconds=(time.time() - start_time))))

        start_time = time.time()
        dsfs.OutputDict(cb_dict, qname_encoder, args.output_file, "reads")
        print("--- %s h:m:s to write reads to file ---" % str(timedelta(seconds=(time.time() - start_time))))
    
    elif args.edit == "frip":
        start_time = time.time()
        dsfs.DownsampleFRIP(cb_dict, frip=args.value, seed=args.seed)
        print("--- %s h:m:s to downsample reads ---" % str(timedelta(seconds=(time.time() - start_time))))
    
        start_time = time.time()
        dsfs.OutputDict(cb_dict, qname_encoder, args.output_file, "reads")
        print("--- %s h:m:s to write reads to file ---" % str(timedelta(seconds=(time.time() - start_time))))
    
    elif args.edit == "peakreads":
        start_time = time.time()
        dsfs.DownsamplePeakReads(cb_dict,  N_desired_reads=args.value, seed=args.seed)
        print("--- %s h:m:s to downsample peak reads ---" % str(timedelta(seconds=(time.time() - start_time))))
    
        start_time = time.time()
        dsfs.OutputDict(cb_dict, qname_encoder, args.output_file, "reads")
        print("--- %s h:m:s to write reads to file ---" % str(timedelta(seconds=(time.time() - start_time))))
        
    else:
        print("ERROR: invalid argument type, dont even think I can reach this part of the code tbh. "); sys.exit(1)
    
    
    rp_post = dsfs.TotalReadPairs(cb_dict)
    cells_post = len(cb_dict.keys())   
    
    start_time = time.time()
    dsfs.WriteLog(args.output_file, 
                  {'input_file':args.input_file, 'sampling_type':args.edit, 'value': args.value, 'seed':args.seed}, 
                  {'N_cells_removed': cells_pre-cells_post, 'N_reads_removed':rp_pre-rp_post}, 
                  dsfs.Summary(cb_dict, output_as = "dict"))
    print("--- %s h:m:s to write log file ---" % str(timedelta(seconds=(time.time() - start_time))))

    start_time = time.time()
    output_bam = os.path.splitext(filename)[0]+".bam"
    returncode = dsfs.GenerateOutputBam(args.edit, args.bam_file, args.output_file, args.nproc, out_dir, output_file=output_bam)
    print("--- %s h:m:s to write bam file ---" % str(timedelta(seconds=(time.time() - start_time))))
    
    if (returncode == 0):
        start_time = time.time()
        dsfs.GenerateOuputFragment("%s/%s" % (out_dir, output_bam), args.nproc)
        print("--- %s h:m:s to write fragment file ---" % str(timedelta(seconds=(time.time() - start_time))))
    else: print("bam didn't generate correctly, not generating fragments")

    
if __name__ == '__main__':
    main()

