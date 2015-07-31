### * Description

# Entry points for the command line scripts

### * Wishlist

# pygenes parse genbank-records/* --genes gene.table --records records.table
# pygenes hash genes.table --peptideField --hash md5 -o genes.table
# pygenes mergePeptides genes.table --maxDissimilarity 0.05 --mapping mergedPep.mapping

### * Set up

### ** Import

import sys
import argparse
import hashlib
import collections
from Bio import SeqIO
import pygenes as pygenes

### * Parser

def makeParser() :
    """Prepare the parser

    Returns:
        ArgumentParser: An argument parser

    """
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help = "One of the available actions")
    # Parse GenBank records
    sp_parse = subparsers.add_parser("parse",
                                     help = "Parse GenBank records to a table")
    sp_parse.add_argument("gbRecords", metavar = "GB_RECORD", nargs = "+",
                          type = str,
                          help = "GenBank records")
    sp_parse.add_argument("-g", "--genes", metavar = "FILE", type = str,
                          help = "Output file for gene table")
    sp_parse.add_argument("-r", "--records", metavar = "FILE", type = str,
                          help = "Output file for record table")
    sp_parse.set_defaults(action = "parse")
    # Calculate hash digests
    sp_hash = subparsers.add_parser("hash",
                                    help = "Produce hash digest (e.g. for peptides)")
    sp_hash.add_argument("input", metavar = "INPUT_TABLE", type = str,
                         help = "Input table")
    sp_hash.add_argument("-o", "--output", metavar = "FILE", type = str,
                         help = "Output file for gene table")
    sp_hash.add_argument("--hash", metavar = "HASH_ALGORITHM",
                         choices = ["md5", "sha1", "sha224", "sha256", "sha384",
                                    "sha512"],
                         default = "md5",
                         help = "Hash algorithm to use for unique sequence signature "
                         "(default: md5)")
    sp_hash.set_defaults(action = "hash")
    # Merge peptides (from gene table)
    sp_mergePep = subparsers.add_parser("mergePep",
                                     help = "Merge similar peptides")
    sp_mergePep.add_argument("input", metavar = "INPUT_TABLE", type = str,
                          help = "Input table")
    sp_mergePep.add_argument("-o", "--output", metavar = "FILE", type = str,
                          help = "Output file for gene table")
    sp_mergePep.add_argument("-d", "--dissim", metavar = "FLOAT", type = float,
                          help = "Maximum dissimilarity for merging "
                             "(between 0 and 1)")
    sp_mergePep.add_argument("-f", "--fasta", metavar = "FILE", type = str,
                      help = "Output fasta file for merged peptides")
    sp_mergePep.set_defaults(action = "mergePep")
    # Extract columns
    sp_extract = subparsers.add_parser("extract",
                                       help = "Extract some columns from a "
                                       "table. Output is sent to stdout.")
    sp_extract.add_argument("inputTable", metavar = "TABLE", type = str,
                            help = "Input table file")
    sp_extract.add_argument("columns", metavar = "COLNAME", type = str,
                            nargs = "+",
                            help = "Name of the column(s) to extract")
    sp_extract.add_argument("-t", "--type", choices = ["gene", "record"],
                            help = "Table type (default: gene)",
                            default = "gene")
    sp_extract.set_defaults(action = "extract")
    # Return
    return parser

### * Mains

### ** Main entry point (dispatch)

def main(args = None, stdout = None, stderr = None) :
    """Main entry point

    Args:
        args (namespace): Namespace with script arguments, parse the command 
          line arguments if None
        stdout (file): Writable stdout stream (if None, use `sys.stdout`)
        stderr (file): Writable stderr stream (if None, use `sys.stderr`)

    """
    if args is None :
        parser = makeParser()
        args = parser.parse_args()
    if stdout is None :
        stdout = sys.stdout
    if stderr is None :
        stderr = sys.stderr
    dispatch = dict()
    dispatch["parse"] = main_parse
    dispatch["hash"] = main_hash
    dispatch["mergePep"] = main_mergePep
    dispatch["extract"] = main_extract
    dispatch[args.action](args, stdout, stderr)
    
### ** Main parse

def main_parse(args, stdout, stderr) :
    if args.genes is None and args.records is None :
        msg = "You should use at least one of -g or -r. Use -h for help.\n"
        stdout.write(msg)
        sys.exit(0)
    if args.genes is not None :
        geneTable = pygenes.GeneTable()
        geneTable.parseGenBankRecord(args.gbRecords)
        stderr.write("Calculating geneId hash\n")
        geneTable.makeGeneId()
        stderr.write("Writing gene table\n")
        geneTable.writeTable(args.genes)
    if args.records is not None :
        stderr.write("Getting record items\n")
        recordTable = pygenes.RecordTable()
        for r in args.gbRecords :
            recordTable.addGenBankRecord(r)
        stderr.write("Writing record table\n")
        recordTable.writeTable(args.records)
    sys.exit(0)

### ** Main hash

def main_hash(args, stdout, stderr) :
    geneTable = pygenes.GeneTable()
    geneTable.loadTable(args.input)
    hashFunctions = { "md5" : hashlib.md5,
                      "sha1" : hashlib.sha1,
                      "sha224" : hashlib.sha224,
                      "sha256" : hashlib.sha256,
                      "sha384" : hashlib.sha384,
                      "sha512" : hashlib.sha512}
    args.hash = hashFunctions[args.hash]
    geneTable.hashPeptides(args.hash)
    geneTable.writeTable(args.output)
    sys.exit(0)
    
### ** Main mergePep

def main_mergePep(args, stdout, stderr) :
    geneTable = pygenes.GeneTable()
    geneTable.loadTable(args.input)
    geneTable.extractSimplifiedPeptides(args.dissim)
    geneTable.writeTable(args.output)
    geneTable.writeSimplifiedPeptides(args.fasta)

### ** Main extract

def main_extract(args, stdout, stderr) :
    if args.type == "gene" :
        table = pygenes.GeneTable()
    else :
        assert args.type == "record"
        table = pygenes.RecordTable()
    table.loadTable(args.inputTable)
    stdout.write("#" + "\t".join(args.columns) + "\n")
    columns = table.col(*args.columns)
    for x in columns :
        stdout.write("\t".join(x) + "\n")
