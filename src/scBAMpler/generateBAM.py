#generateBAM.py
import os, shutil
import downsampling_functions as dsfs

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

    validateTools("samtools"); validateTools("bedtools")
    validateFile(args.input_bam)
    validateFile(args.selected_reads)
    
    returncode = dsfs.GenerateOutputBam(args.input_bam, args.selected_reads, args.nproc, args.output_bam, verbose=args.verbose)
    
    if (returncode == 0 and args.output_fragment == True):
        validateTools("sinto")
        sample_name = os.path.splitext(os.path.basename(args.output_bam))[0]
        dsfs.GenerateOuputFragment(args.output_bam, sample_name+".frags.tsv.bgz", args.nproc, verbose=args.verbose)
    
    return
    
if __name__ == '__main__':
    main()
    