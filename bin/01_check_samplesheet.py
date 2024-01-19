#!/usr/bin/env python3

# modification of the nf-core/rnaseq pipe.
# At IEO user and sequencing facility are providing mostly rid, sid identifiers and the run_folder
# with these parameters would be sufficient to create "channel" per sample merging fastqs across multiple re-seq runs.

# sample test: file_in="/hpcnfs/data/GN2/fgualdrini/Master_batch_scripts/INUTILOMA/Library_testing/Sample_sheet.csv"

import os
import sys
import errno
import argparse
from pathlib import Path
import glob
import re

def parse_args(args=None):
    Description = "Reformat rnaseq samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = f"ERROR: Please check samplesheet -> {error}"
    if context != "" and context_str != "":
        error_str = f"ERROR: Please check samplesheet -> {error}\n{context.strip()}: '{context_str.strip()}'"
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure: use csv comma separation:
    SID,RID,PATH,PAIR,strandedness,name,replicate
    SID,RID,/../../,(forward,reverse,unstranded),SAMPLE_name,replicate in the form of rep1/rep2/rep3/...
    this design follow the logic of sample annotation at the IEO where each sample is assigned a unique SID id
    associated to each run of the specified sample is a RID
    """
    sample_mapping_dict = {}
    with open(file_in, "r", encoding='utf-8-sig') as fin:
        ## Check header
        MIN_COLS = 13
        HEADER = ["rid","sid","sample","replicate","path","lanes","grna_start","grnalength","r1_scaffold_trim","r1_scaffold_pos","umi_read","umi_pos","umi_length","r2_scaffold_trim","r2_scaffold_pos"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print(
                f"ERROR: Please check samplesheet header -> {','.join(header)} != {','.join(HEADER)}"
            )
            sys.exit(1)
        ## Check sample entries
        for line in fin:
            if line.strip():
                lspl = [x.strip().strip('"') for x in line.strip().split(",")]
                ## Check valid number of columns per row
                if len(lspl) < len(HEADER):
                    print_error(
                        f"Invalid number of columns (minimum = {len(HEADER)})!",
                        "Line",
                        line,
                    )
                num_cols = len([x for x in lspl if x])
                if num_cols < MIN_COLS:
                    print_error(
                        f"Invalid number of populated columns (minimum = {MIN_COLS})!",
                        "Line",
                        line,
                    )
                ## Check sample name entries
                rid,sid,sample,replicate,path,lanes,grna_start,grnalength,r1_scaffold_trim,r1_scaffold_pos,umi_read,umi_pos,umi_length,r2_scaffold_trim,r2_scaffold_pos = lspl[: len(HEADER)]
                sample=sample + "_" + replicate  # create a unique identifie
                if sample.find(" ") != -1:
                    print(
                        f"WARNING: Spaces have been replaced by underscores for sample: {sample}"
                    )
                    sample = sample.replace(" ", "_")
                if sample.find("-") != -1:
                    print(
                        f"WARNING: - have been replaced by underscores for sample: {sample}"
                    )
                    sample = sample.replace("-", "_")    
                if not sample:
                    print_error("Sample entry has not been specified!", "Line", line)
                if not r1_scaffold_pos in ["3p","5p"]:
                    print_error("r1_scaffold_pos must be one of 3p or 5p")
                if not r2_scaffold_pos in ["3p","5p"]:
                    print_error("r2_scaffold_pos must be one of 3p or 5p")
                # According to IEO file path policy samples are folders :
                # "..._rid_sid/" within samples are with extension: ".fastq.gz" and again rid_sid_..._R1...fastq.gz and rid_sid_..._R2...fastq.gz
                # the default is paired end sequencing at IEO
                # Given the way the crisper screen library is designed R1 and R2 must be present! Check:
                p=Path(path)
                # here we consider the "lanes" in case issues associated to the NOVASeq lanes are present.
                # lanes can have values in any, or the L001, L002, ... 
                # Illumina outputs include: "_L002_R1_001.fastq.gz"
                ## Create sample mapping dictionary = {sample: [[ single_end, fastq_1, fastq_2, input ]]}
                lanes=lanes.strip().split(":")
                if any(x!='any' for (x) in lanes) and any(x.startswith('L00') for (x) in lanes) and any(x[-1].isdigit() for (x) in lanes):
                    print("lane lenth",len(lanes))
                    # in here we need to split by common partial match:
                    # Important the Lanes must be in the folder specified!
                    
                    subfolders = list(p.glob("*"+rid+"*"+sid+"*"))
                    search = '*fastq.gz'
                    fastqs = list(subfolders[0].glob('./' + search))
                    
                    if len(fastqs) == 0:
                        print_error(
                                f"The path provided does not contain the files as RID+SID",
                                "path:",
                                path
                            )


                    for ll in lanes:
                        search = "*"+rid+"*"+sid+'*'+ll+'*R1*fastq.gz'
                        fastqs_1 = list(subfolders[0].glob(search))
                        search = "*"+rid+"*"+sid+'*'+ll+'*R2*fastq.gz'
                        fastqs_2 = list(subfolders[0].glob(search))
                        
                        if len(fastqs_1) ==0 and len(fastqs_2) ==0:
                            print(
                                f"ERROR: There must be both R1 and R2. R1 containing the gRNA and R2 the UMI"
                            )
                            sys.exit(1)
                        
                        sample_info = [str(fastqs_1[0]), str(fastqs_2[0]), grna_start,grnalength,r1_scaffold_trim,r1_scaffold_pos,umi_read,umi_pos,umi_length,r2_scaffold_trim,r2_scaffold_pos]
                        # for each we append:
                        if sample not in sample_mapping_dict:
                            sample_mapping_dict[sample] = [sample_info]
                        else:
                            if sample_info in sample_mapping_dict[sample]: 
                                print_error("Samplesheet contains duplicate rows!", "Line", line)
                            else:
                                print("append")
                                sample_mapping_dict[sample].append(sample_info)
                elif all(x=='any' for (x) in lanes):
                    print("procede with any lane")
                    # we need to get to the unique elements of "Lanes"
                    subfolders = list(p.glob("*"+rid+"*"+sid+"*"))
                    search = '*R1*fastq.gz'
                    fastqs = list(subfolders[0].glob('./' + search))
                    fastqs = [ re.sub(r'^.*?_L', 'L', str(x)).split('_')[0] for x in fastqs]
                    
                    for ll in fastqs:
                        search = "*"+rid+"*"+sid+'*'+ll+'*R1*fastq.gz'
                        fastqs_1 = list(subfolders[0].glob(search))
                        search = "*"+rid+"*"+sid+'*'+ll+'*R2*fastq.gz'
                        fastqs_2 = list(subfolders[0].glob(search))
                        
                        if len(fastqs_1) ==0 and len(fastqs_2) ==0:
                            print(
                                f"ERROR: There must be both R1 and R2. R1 containing the gRNA and R2 the UMI"
                            )
                            sys.exit(1)
                        sample_info = [str(fastqs_1[0]), str(fastqs_2[0]), grna_start,grnalength,r1_scaffold_trim,r1_scaffold_pos,umi_read,umi_pos,umi_length,r2_scaffold_trim,r2_scaffold_pos ]
                        # for each we append:
                        if sample not in sample_mapping_dict:
                            sample_mapping_dict[sample] = [sample_info]
                        else:
                            if sample_info in sample_mapping_dict[sample]: 
                                print_error("Samplesheet contains duplicate rows!", "Line", line)
                            else:
                                print('append')
                                sample_mapping_dict[sample].append(sample_info)
                else:
                    print(lanes)
                    print_error(
                        f"The lane format must by either any or in the from L00* with multiple lanes separate by : check sample sheet!"
                    )
    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(
                ",".join(["sample", "fastq_1", "fastq_2", "grna_start","grnalength","r1_scaffold_trim","r1_scaffold_pos","umi_read","umi_pos","umi_length","r2_scaffold_trim","r2_scaffold_pos" ])
                + "\n"
            )
            for sample in sorted(sample_mapping_dict.keys()):
                ## Check that multiple runs of the same sample have identical attributes:
                if not all(x[2:len(x)] == sample_mapping_dict[sample][0][2:len(sample_mapping_dict[sample][0])] for x in sample_mapping_dict[sample]):
                    print_error(
                        f"Multiple runs of a sample must have same attributes!!",
                        "Sample",
                        sample,
                    )
                for idx, val in enumerate(sample_mapping_dict[sample]):
                    fout.write(",".join([f"{sample}_T{idx+1}"] + val) + "\n")
    else:
        print_error(f"No entries to process!", "Samplesheet: {file_in}")


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())