import os, pysam, pickle #time
import downsampling_functions as dsfs
#from datetime import timedelta

def main(args):     
    #Check if bam index exists. 
    if not os.path.exists("%s.bai" % args.bam_file): 
        print("--- Indexing bam file")
        pysam.index(args.bam_file) 

    #Generate peak file if not provided
    if args.intersect_file is None:
        intersect_file = "%s.peaks.bed.gz" % os.path.splitext(args.output_file)[0]
        if os.path.isfile(intersect_file): parser.error('The file "{}" exists! Cannot overwrite'.format(intersect_file))
       
        #start_time = time.time()
        dsfs.IntersectPeaks(args.bam_file, args.peak_file, intersect_file)
        #print("--- %s h:m:s to generate intersect file ---" % (str(timedelta(seconds=(time.time() - start_time)))))
    else:
        intersect_file = args.intersect_file
    
    #make dictionary of Cells Objects
    #start_time = time.time()
    cb_dict, cb_encoder, qname_encoder = dsfs.BuildCellDict(args.bam_file)
    #print("--- %s h:m:s to build cell dictionary ---" % str(timedelta(seconds=(time.time() - start_time))))

    #Add peak information to dictionary
    #start_time = time.time()
    dsfs.AddPeakInfo(cb_dict, intersect_file, cb_encoder, qname_encoder, args.delete_intersect)
    #print("--- %s h:m:s to add peak information ---" % (str(timedelta(seconds=(time.time() - start_time)))))

    #Output dictionary as pickle which can be read in for downsampling. 
    #start_time = time.time()
    with open(args.output_file, "wb") as f: # "wb" because we want to write in binary mode
        pickle.dump([cb_dict, cb_encoder, qname_encoder], f) 
    #print("--- %s h:m:s to output pickle ---" % str(timedelta(seconds=(time.time() - start_time))))
    
    mylog = {'bam':args.bam_file, 'peak':args.peak_file}
    mylog.update(dsfs.Summary(cb_dict, output_as = "dict"))
    with open(os.path.splitext(args.output_file)[0]+'.summary.txt', 'w') as f:
        for i in mylog: f.write(i + "\t" + str(mylog[i]) + "\n")

if __name__ == '__main__':
    main()

