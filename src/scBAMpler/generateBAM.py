#generateBAM.py
import os, shutil, sys
from scBAMpler import downsampling_functions as dsfs

def main(args):    

    ########################
    #USER CHECKS
    def validateFile(arg):
        if not os.path.isfile(arg):
            print(f'ERROR: The file "{arg}" does not exist!')
            sys.exit(1)
        return
    def validateTools(tool_name):
        if shutil.which(tool_name) is None:
            print(f"ERROR: Required tool '{tool_name}' not found in PATH.")
            sys.exit(1)

    validateTools("samtools"); validateTools("bedtools")
    validateFile(args.input_bam)
    validateFile(args.selected_reads)

    out_dir = os.path.dirname(args.output_bam)
    if out_dir: os.makedirs(out_dir, exist_ok=True)

    returncode = dsfs.GenerateOutputBam(args.input_bam, args.selected_reads, args.nproc, args.output_bam, verbose=args.verbose)
    if returncode != 0:
        print(f"ERROR: failed to write '{args.output_bam}'.")
        sys.exit(1)

    if args.output_fragment:
        validateTools("sinto")
        #build the fragment path from output_bam so it lands beside the BAM rather than
        #in the current working directory, matching what `sampler` does.
        output_fragment = os.path.splitext(args.output_bam)[0] + ".frags.tsv.bgz"
        frag_status = dsfs.GenerateOuputFragment(args.output_bam, output_fragment, args.nproc, verbose=args.verbose)
        if frag_status != 0:
            print(f"ERROR: failed to write '{output_fragment}'.")
            sys.exit(1)

    return
    
if __name__ == '__main__':
    main()
    