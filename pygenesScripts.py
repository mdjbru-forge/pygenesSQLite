### * Description

# Entry points for the command line scripts

DESCRIPTION = ( "Parse gene data present in GenBank records to produce a gene "
"table containing summary information about each gene, such as the record of "
"origin, the coding sequence, the protein sequence, hash values for the "
"sequences, location, unique gene identifier, ... This utility can also merge "
"similar sequences before running blastp, in order to reduce the run time when "
"many very similar sequences are present. The gene table can be used by tools "
"for gene clustering for example.")

### * Wishlist

# pygenes parse genbank-records/* --genes gene.table --records records.table
# pygenes hash genes.table --peptideField --hash md5 -o genes.table
# pygenes mergePeptides genes.table --maxDissimilarity 0.05 --mapping mergedPep.mapping

### * Set up

### ** Import

import sys
import argparse
import shutil
import hashlib
import os
import collections
import sqlite3 as sql
from Bio import SeqIO
import pygenes as pygenes

### * Parser

def makeParser() :
    """Prepare the parser

    Returns:
        ArgumentParser: An argument parser

    """
    parser = argparse.ArgumentParser(description = DESCRIPTION)
    subparsers = parser.add_subparsers()
### ** Parse GenBank records
    sp_parseGB = subparsers.add_parser("parseGB",
                                     description = "Parse GenBank records into "
                                     "a gene table.",
                                     help = "Parse GenBank records to a table")
    sp_parseGB.add_argument("gbRecords", metavar = "GB_RECORD", nargs = "+",
                          type = str,
                          help = "GenBank records")
    sp_parseGB.add_argument("-g", "--genes", metavar = "FILE", type = str,
                          help = "Output file for gene table")
    sp_parseGB.add_argument("-r", "--records", metavar = "FILE", type = str,
                          help = "Output file for record table")
    sp_parseGB.set_defaults(action = "parseGB")
### ** Parse EMBL records
    sp_parseEMBL = subparsers.add_parser("parseEMBL",
                                     description = "Parse EMBL records into "
                                     "a gene table.",
                                     help = "Parse EMBL records to a table")
    sp_parseEMBL.add_argument("emblRecords", metavar = "EMBL_RECORD", nargs = "+",
                          type = str,
                          help = "EMBL records")
    sp_parseEMBL.add_argument("-g", "--genes", metavar = "FILE", type = str,
                          help = "Output file for gene table")
    sp_parseEMBL.add_argument("-r", "--records", metavar = "FILE", type = str,
                          help = "Output file for record table")
    sp_parseEMBL.set_defaults(action = "parseEMBL")
### ** Build SQL Genomes, Cds and Records tables from EMBL files
    sp_SQL_EMBL = subparsers.add_parser("parseEMBLtoSQL",
                                           description = "Parse EMBL files into an SQLite "
                                           "database",
                                           help = "Parse EMBL files into a new SQLite "
                                           "database")
    sp_SQL_EMBL.add_argument("-o", "--outDb", type = str,
                                help = "Output database (tables will be "
                                "deleted in the database if it already exists)")
    sp_SQL_EMBL.add_argument("emblFiles", metavar = "EMBL_FILE", nargs = "+",
                                type = str,
                                help = "EMBL file(s) (not compressed)")
    sp_SQL_EMBL.set_defaults(action = "SQL_EMBL")
### ** # Calculate hash digests
    # sp_hash = subparsers.add_parser("hash",
    #                                 help = "Produce hash digest (e.g. for peptides)")
    # sp_hash.add_argument("input", metavar = "INPUT_TABLE", type = str,
    #                      help = "Input table")
    # sp_hash.add_argument("-o", "--output", metavar = "FILE", type = str,
    #                      help = "Output file for gene table")
    # sp_hash.add_argument("--hash", metavar = "HASH_ALGORITHM",
    #                      choices = ["md5", "sha1", "sha224", "sha256", "sha384",
    #                                 "sha512"],
    #                      default = "md5",
    #                      help = "Hash algorithm to use for unique sequence signature "
    #                      "(default: md5)")
    # sp_hash.set_defaults(action = "hash")
### ** Merge peptides (from gene table)
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
### ** Extract columns
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
### ** Count unique sequences in a fasta file
    sp_countUniqFasta = subparsers.add_parser("count",
                                              help = "Count unique sequences in "
                                              "a fasta file (using hash)")
    sp_countUniqFasta.add_argument("inputFiles", metavar = "FASTA_FILE",
                                   type = str, nargs = "+",
                                   help = "Input fasta file(s)")
    sp_countUniqFasta.set_defaults(action = "count")
### ** Return parser
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
    dispatch["parseGB"] = main_parseGB
    dispatch["parseEMBL"] = main_parseEMBL
    dispatch["SQL_EMBL"] = main_SQL_EMBL
    dispatch["hash"] = main_hash
    dispatch["mergePep"] = main_mergePep
    dispatch["extract"] = main_extract
    dispatch["count"] = main_count
    dispatch[args.action](args, stdout, stderr)
    
### ** Main parseGB

def main_parseGB(args, stdout, stderr) :
    if args.genes is None and args.records is None :
        msg = "You should use at least one of -g or -r. Use -h for help.\n"
        stdout.write(msg)
        sys.exit(0)
    if args.genes is not None :
        geneTable = pygenes.GeneTable()
        geneTable.parseGenBankRecord(args.gbRecords)
        stderr.write("Calculating geneId hash\n")
        geneTable.writeTable(args.genes)
    if args.records is not None :
        stderr.write("Getting record items\n")
        recordTable = pygenes.RecordTable()
        for r in args.gbRecords :
            recordTable.addGenBankRecord(r)
        stderr.write("Writing record table\n")
        recordTable.writeTable(args.records)
    sys.exit(0)

### ** Main parseEMBL

def main_parseEMBL(args, stdout, stderr) :
    # if args.genes is None and args.records is None :
    #     msg = "You should use at least one of -g or -r. Use -h for help.\n"
    #     stdout.write(msg)
    #     sys.exit(0)
    if args.genes is not None :
        geneTable = pygenes.GeneTable()
        headers = geneTable.getHeaders()
        with open(args.genes, "w") as fo :
            fo.write("#" + "\t".join(headers) + "\n")
            geneTable.streamEMBLRecord(args.emblRecords, outFile = fo,
                                       headers = headers)
    if args.records is not None :
        msg = "Parsing records into record table"
        stderr.write(msg + "\n")
        recordTable = pygenes.RecordTable()
        headers = recordTable.getHeaders()
        with open(args.records, "w") as fo :
            fo.write("#" + "\t".join(headers) + "\n")
            for r in args.emblRecords :
                recordTable.streamEMBLRecord(r, fo, headers)
    sys.exit(0)

### ** Main SQL_EMBL

def main_SQL_EMBL(args, stdout, stderr) :
    # Create the table
    dbConnection = sql.connect(args.outDb)
    cursor = dbConnection.cursor()
    cursor.execute("DROP TABLE IF EXISTS Cds")
    cursor.execute("CREATE TABLE Cds ("
                   "record_id TEXT, "
                   "pepSeq TEXT, "
                   "nucSeq TEXT, "
                   "pepLen INTEGER, "
                   "location TEXT, "
                   "translationTable INTEGER, "
                   "geneName TEXT, "
                   "productName TEXT, "
                   "productAccNum TEXT)")
    cursor.execute("DROP TABLE IF EXISTS Records")
    cursor.execute("CREATE TABLE Records ("
                   "id TEXT UNIQUE, "
                   "name TEXT, "
                   "description TEXT, "
                   "seq TEXT, "
                   "seqLen INTEGER, "
                   "genome_filename TEXT)")
    cursor.execute("DROP TABLE IF EXISTS Genomes")
    cursor.execute("CREATE TABLE Genomes ("
                   "filename TEXT UNIQUE, "
                   "biosample TEXT UNIQUE, "
                   "organism TEXT, "
                   "nRecords INTEGER, "
                   "refs TEXT)")
    # Go through the EMBL files
    total = str(len(args.emblFiles))
    for (i, f) in enumerate(args.emblFiles) :
        stderr.write("Processing file " + str(i) + "/" + total + " - " +
                     os.path.basename(f) + "\n")
        # Genomes table
        d = pygenes.EMBLFileInfo(f)
        if d is not None  :
            cursor.execute("INSERT INTO Genomes (filename, biosample, "
                           "organism, nRecords, refs) "
                           "VALUES (\"{filename}\", \"{biosample}\", "
                           "\"{organism}\", {nRecords}, \"{references}\" "
                           ")".format(filename = d["filename"],
                                      biosample = d["biosample"],
                                      organism = d["organism"],
                                      nRecords = str(d["nRecords"]),
                                      references = d["references"]))
            dbConnection.commit()
        # Cds and Records tables
        pygenes.parseEMBLtoSQL(f, cursor)
        dbConnection.commit()
    # Close the connection
    dbConnection.close()
    
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
    msg = "Building the hash index"
    stderr.write(msg + "\n")
    hashIndex = pygenes.buildGeneTableFileIndex(args.input)
    msg = "Building the length index"
    stderr.write(msg + "\n")
    lengthIndex = pygenes.buildLengthFileIndex(hashIndex)
    hashToMerged = dict()
    with open(args.fasta, "w") as fo :
        for l in lengthIndex.keys() :
            msg = "Merging sequences of length " + str(l)
            stderr.write(msg + "\n")
            sequences = pygenes.gatherSequences(args.input, lengthIndex, l)
            mergedSequences = pygenes.mergeSequences(sequences, args.dissim)
            for (k, v) in mergedSequences.items() :
                originalHash = pygenes.md5hash(k)
                mergedHash = pygenes.md5hash(v)
                assert not hashToMerged.get(originalHash, False)
                hashToMerged[originalHash] = mergedHash
            newMerged = set(mergedSequences.values())
            for seq in newMerged :
                fo.write(">" + pygenes.md5hash(seq) + "\n")
                fo.write(seq + "\n")
    with open(args.input, "r") as fi :
        with open("tmp." + args.output, "w") as fo :
            headers = fi.readline()
            headerElements = headers.lstrip("#").strip().split("\t")
            fo.write(headers)
            for line in fi :
                content = dict(zip(headerElements, line.strip().split("\t")))
                content["mergedPeptideHash"] = hashToMerged[content["peptideHash"]]
                fo.write("\t".join([content[x] for x in headerElements]) + "\n")
    shutil.move("tmp." + args.output, args.output)
                
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

### ** Main count

def main_count(args, stdout, stderr) :
    def h(string) :
        h = hashlib.md5()
        h.update(string)
        return h.hexdigest()
    cont = True
    for f in args.inputFiles :
        if cont :
            try :
                a = SeqIO.parse(f, "fasta")
                uniqHash = set()
                for s in a :
                    uniqHash.add(h(str(s.seq)))
                stdout.write(f + "\t" + str(len(uniqHash)) + "\n")
            except IOError :
                cont = False

