import os, pysam, pickle, sys, shutil, subprocess
import downsampling_functions as dsfs

def main(args):     

    ########################
    #USER CHECKS
    def validateFile(arg):
        if not os.path.isfile(arg): print(f'ERROR: The file "{arg}" does not exist!')
        else: return 

    def validateBAM(bam_path):
        try:
            subprocess.check_output(['samtools', 'quickcheck', '-v', bam_path], stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print(f"ERROR: BAM file '{bam_path}' failed samtools quickcheck:\n{e.output.decode().strip()}")
            sys.exit(1)
    
    def validateTools(tool_name):
        if shutil.which(tool_name) is None:
            print(f"ERROR: Required tool '{tool_name}' not found in PATH.")
            sys.exit(1)
    
    validateBAM(args.bam_file)
    validateFile(args.peak_file)
    if args.intersect_file is not None: validateFile(args.intersect_file)
    validateTools("samtools"); validateTools("bedtools")
    
    #Check if bam index exists. 
    if not os.path.exists("%s.bai" % args.bam_file): 
        pysam.index(args.bam_file)
        print("--- Created BAM index file")

    
    ########################
    #Generate peak file if not provided
    if args.intersect_file is None:
        intersect_file = "%s.peaks.bed.gz" % os.path.splitext(args.output_file)[0]
        if os.path.isfile(intersect_file): parser.error('The file "{}" exists! Cannot overwrite'.format(intersect_file))
        dsfs.IntersectPeaks(args.bam_file, args.peak_file, intersect_file)
    else:
        intersect_file = args.intersect_file

    ########################
    #make dictionary of Cells Objects
    cb_dict, cb_encoder, qname_encoder = dsfs.BuildCellDict(args.bam_file)

    ########################
    #Add peak information to dictionary
    dsfs.AddPeakInfo(cb_dict, intersect_file, cb_encoder, qname_encoder, args.delete_intersect)

    ########################
    #Output dictionary as pickle which can be read in for downsampling. 
    with open(args.output_file, "wb") as f: # "wb" because we want to write in binary mode
        pickle.dump([cb_dict, cb_encoder, qname_encoder], f) 
    
    mylog = {'bam':args.bam_file, 'peak':args.peak_file}
    mylog.update(dsfs.Summary(cb_dict, output_as = "dict"))
    with open(os.path.splitext(args.output_file)[0]+'.summary.txt', 'w') as f:
        for i in mylog: f.write(i + "\t" + str(mylog[i]) + "\n")


if __name__ == '__main__':
    main()

