#downsampling_functions.py

import pysam, os, subprocess, sys, itertools, functools, shlex
import pandas as pd
import numpy as np
from collections import Counter
from datetime import datetime

def internal_timer(func):
    @functools.wraps(func)
    def wrapper_decorator(*args, **kwargs):
        verbose = kwargs.get("verbose", True)
        start_time = datetime.now()
        result = func(*args, **kwargs)  # No extra reference
        end_time = datetime.now()
        if verbose:
            print(f"--- {end_time - start_time} h:m:s to run '{func.__name__}' ---")
        return result
    return wrapper_decorator
    

#FUNCTIONS FOR CREATING THE DICTIONARY + OBJECT
#################################################

class Cells(object):
    #initiate a class to store all of the relevant information. 
    #I read it was more memory efficient if you specify slots so it doesn't create a dictionary (?)
    __slots__ = ('cb', 'readslist', 'readcount', 'peaklist', "peakcount", "nonpeakcount", "n_edits")
    def __init__(self, cb):
        self.cb = cb
        self.readslist = set()
        self.peaklist = []
        self.n_edits = 0
        self.peakcount = 0
        self.nonpeakcount = 0

@internal_timer
def BuildCellDict(bam_file, verbose=True):

    #read in bam file using pysam
    cb_dict = {}
    cb_encoder = {}
    qname_encoder = {}
    c = 0
    q = 0
    no_cb = 0
    curr_bam = pysam.AlignmentFile(bam_file, 'rb')

    #for read in tqdm(curr_bam, desc="Progress adding cell barcodes", total=curr_bam.mapped):
    for read in curr_bam:
        if not read.is_read1: continue    #we only want read pairs, so select only the first read.

        #not every alignment is guaranteed to carry a barcode, so tally and skip these
        #rather than raising a KeyError partway through the file.
        if not read.has_tag("CB"):
            no_cb += 1
            continue

        cb = read.get_tag("CB")
        qname = read.query_name #in new files

        #check if we've seen this cell barcode before. regardless, get its integer encoding
        if cb not in cb_encoder:
            cb_int = c; c+=1
            cb_encoder[cb] = cb_int #add new entry to encoder
            cb_dict[cb_int] = Cells(cb_int) # add new entry to dictionary. 
        else:
            cb_int = cb_encoder[cb]

        #check if we've seen this read before. regardless, get its integer encoding
        if qname not in qname_encoder:
            qname_int = q; q+=1
            qname_encoder[qname] = qname_int
        else:
            qname_int = qname_encoder[qname]

        cb_dict[cb_int].readslist.add(qname_int) #append readname to list of reads
        #purposely using a set here so the reads aren't redundant
        

    if no_cb > 0 and verbose:
        print(f"--- Skipped {no_cb} read1 alignments with no CB tag")

    if not cb_dict:
        print(f"ERROR: No reads with a CB tag were found in '{bam_file}'. "
              "scBAMpler expects cell barcodes in the CB:Z: tag.")
        sys.exit(1)

    #when all reads are done, quickly go back and add a tally of how many read-pairs per cell there are.
    for item in cb_dict.keys():
        cb_dict[item].readslist = list(cb_dict[item].readslist) #new
        a = len(cb_dict[item].readslist)
        cb_dict[item].readcount = a
        cb_dict[item].peaklist = np.zeros(a) #initialize vector of 0s for peak info. 
        
    return cb_dict, cb_encoder, qname_encoder

#################################################

@internal_timer
def IntersectPeaks(bam_file, peak_file, intersect_file, timeout = 21600, verbose=True):
    #timeout in 6hrs. 
    
    awk_statement = '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){sub(/^CB:Z:/, "", $i); print $i, $1 }}}'
    cmd = ("set -o pipefail; "
       "bedtools intersect "
       "-abam {bam} "
       "-b <(bedtools sort -i {bed} -faidx <(samtools view -H {bam} | grep '^@SQ' | sed 's/.*SN://' | cut -f1)) "
       "-sorted -f 0.75 -ubam | "
       "samtools view -h - | awk {awk} | sort | uniq | gzip > {out}").format(
    bam=shlex.quote(bam_file), bed=shlex.quote(peak_file),
    awk=shlex.quote(awk_statement), out=shlex.quote(intersect_file))

    #exit on failure rather than returning: downstream steps would otherwise read an
    #empty intersect file and fail with a confusing pandas error instead of this one.
    #remove the partial output first, otherwise the next run refuses to overwrite it.
    try:
        subprocess.check_output(cmd, shell=True, executable='/bin/bash', stderr=subprocess.STDOUT, timeout=timeout)
    except subprocess.CalledProcessError as e:
        stderr_output = e.output  # This contains the stderr output
        print("ERROR: Bedtools command failed. stderr:", stderr_output.decode())
        _remove_partial(intersect_file)
        sys.exit(1)
    except BaseException as e: #BaseException so Ctrl-C also cleans up
        print("ERROR: An error occurred:", str(e))
        _remove_partial(intersect_file)
        sys.exit(1)

    return


def _remove_partial(path):
    if os.path.exists(path):
        os.remove(path)
        print(f"--- Removed partial output '{path}'")

    
@internal_timer
def AddPeakInfo(cb_dict, intersect_file, cb_encoder, qname_encoder, delete, verbose=True):
    rip_df  = pd.read_csv(intersect_file, compression='gzip', header=None, sep=' ', names=["cb", "qname"])

    unmatched = []

    def custom_function(group):
        l = []
        for i in group['qname']:
            try:
                l.append(qname_encoder[i])
            except KeyError:
                l.append(None)
                #collect rather than print per-read: on a full-size BAM this can be
                #millions of lines of output.
                unmatched.append(i)
        return(l)

    #select the column explicitly after groupby: the grouping column is excluded from
    #apply() in pandas 3, and being explicit here silences that FutureWarning without
    #changing behaviour, since custom_function only ever reads 'qname'.
    grouped_results = rip_df.groupby('cb')[['qname']].apply(custom_function)
    grouped_results_filt = grouped_results.dropna()
    q_dict = grouped_results_filt.to_dict()

    if unmatched and verbose:
        print(f"--- {len(unmatched)} peak reads had no matching read1 in the BAM "
              f"(possible unmapped R2), e.g. {unmatched[0]}")

    for curr_cb in q_dict.keys(): #for all cb we need to update
        curr_cb_int = cb_encoder[curr_cb]
        idxs = np.where(np.isin(cb_dict[curr_cb_int].readslist, q_dict[curr_cb]))[0]
        cb_dict[curr_cb_int].peaklist[idxs] = 1 

    for cb_int in cb_dict.keys():
        m = cb_dict[cb_int].peaklist
        ones_count = np.count_nonzero(m==1)
        cb_dict[cb_int].peakcount = ones_count
        cb_dict[cb_int].nonpeakcount = len(m) - ones_count
    
    if (delete):
        try:
            os.remove(intersect_file)
        except OSError as e:
            print("ERROR: An error occurred while removing intersect file:", str(e))
        
    return

#################################################

def CalculateFRIP(cb_dict):
    peak = np.sum([cb_dict[item].peakcount for item in cb_dict.keys()])
    nonpeak = np.sum([cb_dict[item].nonpeakcount for item in cb_dict.keys()])
    
    if peak == 0: frip="NA"
    else: frip = round(peak/(peak+nonpeak),3)
    
    return frip, peak, nonpeak


def TotalReadPairs(cb_dict):
    tot = np.sum([cb_dict[item].readcount for item in cb_dict.keys()])
    return tot


def Summary(cb_dict):
    #Summarize cb_dict as a dictionary, written to the .summary.txt logs.
    edit_list = np.array([int(cb_dict[i].n_edits) for i in cb_dict.keys()])
    curr_frip, peakPairs, nonpeakPairs = CalculateFRIP(cb_dict)

    output_dict = {}
    output_dict["Ncells"] = len(cb_dict.keys())
    output_dict["Nreadpairs"] = TotalReadPairs(cb_dict)
    output_dict["Npeakpairs"] = peakPairs
    output_dict["Nnonpeakpairs"] = nonpeakPairs
    output_dict["FRIP"] = curr_frip
    output_dict["Nedits"] = np.sum(edit_list)
    output_dict["Ncells_with_edits"] = np.sum(edit_list > 0)
    return output_dict

#################################################
## GENERIC DOWNSAMPLING FUNCTIONS

def ChooseCells(cb_dict, N, sample_case, seed):
    #sample_case: "random", "peaks", "nonpeaks" -- what we're drawing from. 
    #N: int -- how many reads are we REMOVING

    cells = list(cb_dict.keys()) 
    
    #What is the maximum number of draws per cell we can do?
    match sample_case:
        case "random":
            replace_limits = [cb_dict[i].readcount for i in cb_dict] #number of max reads
        case "peaks":
            replace_limits = [cb_dict[i].peakcount for i in cb_dict] 
        case "nonpeaks":
            replace_limits = [cb_dict[i].nonpeakcount for i in cb_dict]        
    
    #Weight our sampling based on how many draws we could do. 
    weights = replace_limits/np.sum(replace_limits)
    
    #repeat cells and weights the exact number of times we could draw it so we don't have to sample w/replacement
    #e.g CellA: 4 reads, CellB: 1 reads -> [CellA, CellA, CellA, CellA, CellB]
    repeated_cells = np.repeat(cells, replace_limits) 
    tmp = np.repeat(weights, replace_limits)
    repeated_weights = tmp/np.sum(tmp)
    
    #set seed and sample - return a redundant list of cell barcode names.
    np.random.seed(seed)
    chosen_cells = np.random.choice(repeated_cells, size=int(N), replace=False, p=repeated_weights)

    #Count how many times the CB was selected. {"CB":5, "CB2":10}
    chosen_cells_count = dict(Counter(chosen_cells)) #count how may times cell barcode appears
    
    #for all cells, save how many edits they're getting in the Cells object. 
    for i in chosen_cells_count: cb_dict[i].n_edits = chosen_cells_count[i] #update class
    
    return


def RemoveReads(cell_object, seed, sample_case):
    if cell_object.n_edits == 0:
        return #quick bail to save time
    
    idx_range = np.array(range(cell_object.readcount))
    
    match sample_case:
        case "random":
            idx_select_from = range(cell_object.readcount) #can choose any of the indexes
            np.random.seed(seed);
            idx_rem = np.random.choice(idx_select_from, size=cell_object.n_edits, replace=False)
            
        case "peaks":
            idx_select_from = np.where(cell_object.peaklist == 1)[0] #can only choose indices in readlist where peaklist =1
            np.random.seed(seed);
            idx_rem = np.random.choice(idx_select_from, size=cell_object.n_edits, replace=False)
    
        case "nonpeaks":
            idx_select_from = np.where(cell_object.peaklist == 0)[0] #can only choose indices in readlist where peaklist =0
            np.random.seed(seed);
            idx_rem = np.random.choice(idx_select_from, size=cell_object.n_edits, replace=False)
    
    #update cell_object in place so the dictionary is accurate. 
    idx_keep = idx_range[~np.isin(idx_range,idx_rem)] #np.in1d is deprecated in numpy 2.x
    cell_object.readslist = np.array(cell_object.readslist)[idx_keep]
    cell_object.peaklist = cell_object.peaklist[idx_keep]
    cell_object.readcount = len(cell_object.readslist)
    cell_object.peakcount = np.count_nonzero(cell_object.peaklist == 1)
    cell_object.nonpeakcount = np.count_nonzero(cell_object.peaklist == 0)        
    return



def CleanCells(cb_dict):
    #sometimes we can remove all reads from a cell barcode. 
    #this function just removes those from the dict so we can see how many are completely lost. 
    
    rem_list = []
    for item in cb_dict.keys():
        if cb_dict[item].readcount == 0:
            rem_list.append(item)
    #cant del immediately in loop bc dict size changes
    for i in rem_list:
        del cb_dict[i]
    return


@internal_timer
def OutputDict(cb_dict, encoder, out_readfile, verbose):
    
    tmp = [v.readslist for v in cb_dict.values()] #numeric list of reads needed. 
    results_int = list(itertools.chain.from_iterable(tmp)) 
    
    inv_encoder = {v: k for k, v in encoder.items()} #invert the cb_encoder dictionary.
    results_str = [inv_encoder[i] for i in results_int]
    results_str = sorted(results_str) #easier later, but can remove depending on how long it takes

    with open(out_readfile, mode='wt', encoding='utf-8') as f:
        for i in results_str: f.write(i+"\n")
    
    return

#################################################
## FUNCTIONS SPECIFIC TO DOWNSAMPLING CELLS

@internal_timer
def DownsampleCells(cb_dict, N_cells, seed, verbose):
    
    #set seed, randomly choose N cell barcodes without replacement. 
    #either return as a (1) dictionary -- for consistency with other functions
    #or (2) list of cell names, obv faster. 
    
    if N_cells < 0 or N_cells > len(cb_dict):
        print("ERROR: Requested %d cells but must be between 0 and %d (current total)" % (N_cells, len(cb_dict)))
        sys.exit(1)

    np.random.seed(seed)
    chosen_cells = np.random.choice(list(cb_dict.keys()), size=int(N_cells), replace=False)
    cb_dict_sub = {cb:cb_dict[cb] for cb in chosen_cells}
    return cb_dict_sub
    
    
#################################################
## FUNCTIONS SPECIFIC TO DOWNSAMPLING READS

@internal_timer
def DownsampleReads(cb_dict, N_desired_reads, seed, verbose):

    total_reads = TotalReadPairs(cb_dict)
    if N_desired_reads < 0 or N_desired_reads > total_reads:
        print("ERROR: Requested %d read pairs but must be between 0 and %d (current total)" % (N_desired_reads, total_reads))
        sys.exit(1)
    Nreads_to_remove = total_reads-N_desired_reads
    
    #choose which cells are going to get downsampled and by how much. stored in {'cb1':{n_edits = 3}}
    ChooseCells(cb_dict, N=Nreads_to_remove, seed=seed, sample_case ="random")
    
    #Remove the read from that cb's Cells() object (in place so nothing to return)
    #iterating through all keys here, but that function will quickly bail if no edits need to be made. 
    for item in cb_dict.keys():
        RemoveReads(cb_dict[item], seed=seed, sample_case="random")

    CleanCells(cb_dict) 
    return


#################################################
## FUNCTIONS SPECIFIC TO DOWNSAMPLING FRIP

@internal_timer
def DownsampleFRIP(cb_dict, frip, seed, verbose):
    curr_frip, peakPairs, nonpeakPairs = CalculateFRIP(cb_dict)
    desired_frip = frip

    #CalculateFRIP returns "NA" when there are no peak reads, which cannot be compared below
    if peakPairs == 0:
        print("ERROR: No reads overlap peaks in this dictionary, so FRIP cannot be downsampled. "
              "Check that the peak file passed to create-dictionary matches the BAM."); sys.exit(1)
    #both formulas below divide by frip or 1-frip
    if not 0 < desired_frip < 1:
        print(f"ERROR: FRIP must be between 0 and 1 (exclusive), got {desired_frip}"); sys.exit(1)

    if curr_frip > desired_frip: #need to remove peak reads
        Npeak_to_remove = peakPairs - round((desired_frip*nonpeakPairs)/(1-desired_frip)) #p/n+p = frip, solved for p

        #user check that we're not decimating either peak/nonpeak counts too much (here, setting as 1000 read pairs)
        if Npeak_to_remove > peakPairs-1000:
            print("ERROR: Minimum amount of peak reads is 1000. Cannot satisfy this FRIP value"); sys.exit(1)
        
        #choose cells to remove peak reads from 
        ChooseCells(cb_dict, N=Npeak_to_remove, seed=seed, sample_case ="peaks")
        for item in cb_dict.keys():
            RemoveReads(cb_dict[item], seed=seed, sample_case="peaks") #remove them
    
    elif curr_frip < desired_frip: #need to remove nonpeak reads  
        Nnonpeak_to_remove = nonpeakPairs - round((peakPairs*(1-desired_frip))/desired_frip) #p/n+p = frip, solved for n
        
        if Nnonpeak_to_remove > nonpeakPairs-1000:
            print("ERROR: Minimum amount of nonpeak reads is 1000. Cannot satisfy this FRIP value"); sys.exit(1)
        
        ChooseCells(cb_dict, N=Nnonpeak_to_remove, seed=seed, sample_case ="nonpeaks")
        for item in cb_dict.keys():
            RemoveReads(cb_dict[item], seed=seed, sample_case="nonpeaks")
    
    else: print("ERROR current and desired frip are identical. Please fix"); sys.exit(1)
    
    CleanCells(cb_dict)
    return


def _submit_cmd(cmd, err = "ERROR"):
    #returns 0 on success, non-zero otherwise, so callers can stop rather than
    #running later steps against output an earlier step never produced.
    try:
        subprocess.check_output(f"set -o pipefail; {cmd}", shell=True, executable='/bin/bash', stderr=subprocess.STDOUT)
        return 0
    except subprocess.CalledProcessError as e:
        stderr_output = e.output  # This contains the stderr output
        print(err, stderr_output.decode())
        return e.returncode if e.returncode != 0 else 1
    except Exception as e:
        print(err, str(e))
        return 1



@internal_timer   
def GenerateOutputBam(input_bam, read_file, nproc, output_file, verbose):
        
    q = shlex.quote
    cmd = 'samtools view -N %s -o %s %s -@ %s' % (q(read_file), q(output_file), q(input_bam), str(nproc))
    status = _submit_cmd(cmd, "ERROR: in generate output bam")
    if status != 0: return status

    cmd2 = 'samtools index %s' % q(output_file)
    return _submit_cmd(cmd2, "ERROR: in indexing output bam")


@internal_timer
def GenerateOutputFragment(input_bam, output_fragment, nproc, verbose):
    tmp_output = output_fragment + "_tmp"
    q = shlex.quote

    #strip the whole suffix, not just the last extension: splitext only removes '.bgz',
    #which left '.frags.tsv' embedded in every cell name in the fragment file.
    #note we cannot split on '.' here -- prefixes legitimately contain them (e.g. f0.2_s33).
    sample_name = os.path.basename(output_fragment)
    for suffix in (".frags.tsv.bgz", ".frags.tsv", ".tsv.bgz", ".bgz"):
        if sample_name.endswith(suffix):
            sample_name = sample_name[:-len(suffix)]
            break

    cmd1 = "sinto fragments --collapse_within -p %s -b %s -f %s > /dev/null" % (nproc, q(input_bam), q(tmp_output))
    status = _submit_cmd(cmd1, "ERROR: in sinto fragment creation (step 1)")
    if status != 0:
        #bail out here: without the tmp file, steps 2 and 4 can only fail too, and
        #reporting three errors for one root cause makes the real problem harder to see.
        print("ERROR: skipping remaining fragment steps. No fragment file was produced.")
        return status

    #pound and dash do not work btw for archr.
    #pass the sample name to awk as a variable rather than splicing it into the program text
    awk_part = '{print $1, $2, $3, name":"$4, $5}'
    cmd2 = (f"bedtools sort -i {q(tmp_output)} | awk -v name={q(sample_name)} {q(awk_part)} "
            fr"| tr ' ' '\t' | bgzip -c > {q(output_fragment)}")
    status = _submit_cmd(cmd2, "ERROR: bedtools bgzipped (step 2)")
    if status != 0: return status

    #cmd3 = f"tabix {outfile}"
    #_submit_cmd(cmd3, "ERROR: in indexing bgzipped (step 3)")

    try:
        os.remove(tmp_output)
    except OSError as e:
        print("ERROR: in removing file (step 4)", str(e))
        return 1
    return 0


@internal_timer
def WriteLog(output_file, a, b, c, verbose):
    logfile = open(output_file, 'w')
    
    logfile.write("## params ##\n")
    for k,v in a.items(): logfile.write(f"{k}\t{v}\n")
    
    logfile.write("\n## edit info ##\n")
    for k,v in b.items(): logfile.write(f"{k}\t{v}\n")
        
    logfile.write("\n## resulting object ##\n")
    for k,v in c.items(): logfile.write(f"{k}\t{v}\n")
    
    logfile.close()
    return
