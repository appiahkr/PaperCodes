import sys, os, subprocess  
from Bio import SeqIO      
import gzip                
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec  
import numpy as np                  
import collections                  
from collections import defaultdict


# ------------------------
# Main function
# ------------------------
def main():
    usage = 'usage:' + sys.argv[0] + '<sexData><keyDir><fastqDir>'
    
    # Check if the correct number of arguments are provided
    if len(sys.argv) != 4:
        print(usage)
        sys.exit(1)

    # Assign command-line arguments to variables
    sexData = sys.argv[1]   # CSV mapping accession IDs to sex
    keyDir = sys.argv[2]    # Directory containing barcode key files
    fastqDir = sys.argv[3]  # Directory or FASTQ file to process

    # Read sex metadata
    accToSex = getSexData(sexData)

    # Map barcodes to samples based on sex
    sexBarcodeToAcc, barcodeLengthList = getGenderFastq(keyDir, accToSex)

    # Extract reads from FASTQ files using the barcodes
    extractFastqUsingBarcodes(fastqDir, sexBarcodeToAcc, barcodeLengthList)

# ------------------------
# Extract FASTQ reads by barcode
# ------------------------
def extractFastqUsingBarcodes(fastqFile, sexBarcodeToAcc, barcodeLengthList):
    """
    Extract reads from a FASTQ file matching barcodes for each sample.
    
    Args:
        fastqFile (str): path to gzipped FASTQ file
        sexBarcodeToAcc (dict): mapping of sex -> (barcode, readName, lane) -> accession
        barcodeLengthList (list): list of barcode lengths (unused in this function)
    """
    
    # Parse read name and lane from the FASTQ file name
    readName1, lane1 = fastqFile.split('_')[2], fastqFile.split('_')[4]
    
    # Open gzipped FASTQ file for reading
    with gzip.open(fastqFile, "rt") as fastq_file:
        # Iterate over all reads in the FASTQ file
        for record in SeqIO.parse(fastq_file, "fastq"):
            seq = record.seq  # Sequence of the current read
            
            # Iterate over each sex category (male/female)
            for sex, accDict in sexBarcodeToAcc.items():
                # Iterate over all barcode/readName/lane combinations for this sex
                for (barcode, readName2, lane2) in accDict:
                    barcodeLen = len(barcode)
                    
                    # Check if this read matches the barcode, readName, and lane
                    if readName2 == readName1 and lane2 == lane1 and seq[:barcodeLen] == barcode:
                        acc = accDict[(barcode, readName2, lane2)]
                        
                        # Generate output file name
                        fileName = f"{sex}_{acc}_{readName1}_s_{lane1}.fastq"
                        
                        # Write matching read to the output FASTQ file
                        # NOTE: Currently opens file in append mode without closing properly
                        out = open(fileName, 'a')
                        SeqIO.write(record, out, 'fastq')
                        # Improvement: Use 'with open(fileName, 'a') as out:' to ensure file is closed

# ------------------------
# Read sex metadata
# ------------------------
def getSexData(sexData):
    """
    Read a CSV file mapping accession IDs to sex.
    
    Args:
        sexData (str): path to CSV file (acc,sex)
        
    Returns:
        dict: acc -> 'male' or 'female' assuming female is represented as 0
        and female is represented as 1
    """
    accToSex = {}
    for line in open(sexData):
        fields = line.strip().split(',')
        if len(fields) > 1:
            acc, sex = fields
            # Map '0' to female, any other value to male
            if sex == '0':
                accToSex[acc] = 'female'
            else:
                accToSex[acc] = 'male'
    return accToSex

# ------------------------
# Read barcode key files and associate barcodes with sex/accession
# ------------------------
def getGenderFastq(keysDir, accToSex):
    """
    Read key files and create a mapping from sex -> (barcode, readName, lane) -> accession.
    
    Args:
        keysDir (str): directory containing key files
        accToSex (dict): mapping accession -> sex
    
    Returns:
        tuple: 
            sexToAcc (dict): sex -> (barcode, readName, lane) -> accession
            sortedList (list): sorted list of barcode lengths
    """
    sexToAcc = defaultdict(dict)
    files = os.listdir(keysDir)
    barcodeLengthList = set()  # store all barcode lengths
    
    # Iterate over all files in the directory
    for filename in files:
        file_path = os.path.join(keysDir, filename)
        if 'key' in filename:  # only process files with 'key' in the name
            for line in open(file_path):
                fields = line.strip().split(',')
                readName, lane, barcode, acc = fields[0], fields[1], fields[2], fields[3]
                barcodeLengthList.add(len(barcode))
                acc = acc.lower().replace('_','')
                
                # Only include accessions that are in sex metadata
                if acc in accToSex:
                    sex = accToSex[acc]
                    sexToAcc[sex][(barcode, readName, lane)] = acc
    
    # Return sorted list of barcode lengths for reference
    sortedList = sorted(barcodeLengthList)
    return sexToAcc, sortedList

# ------------------------
# Run the main function
# ------------------------
main()
B
