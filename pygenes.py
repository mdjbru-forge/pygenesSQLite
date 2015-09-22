### * Description

# Complementary module for genbank.py

### * Setup

### ** Import

import sys
import collections
import itertools
import warnings
import os
import gzip
import hashlib
from Bio import SeqIO

### * Functions

### ** seqDistance(seq1, seq2)

def seqDistance(seq1, seq2) :
    """Calculate the dissimilarity between two sequences.
    Compare pairs of characters, and count number of mismatches and pairs 
    with at least one X (also considered mismatch).
    Note: matching gaps would be considered matches here! No gaps should be 
    present in the sequences.

    Args:
        seq1 (str): First sequence
        seq2 (str): Second sequence (same length as the first sequence)

    Returns:
        float: Dissimilarity (0-1) between the two sequences
    """
    mismatch = 0
    assert len(seq1) == len(seq2)
    for (x, y) in zip(seq1, seq2) :
        if x != y or "x" in (x + y) or "X" in (x + y) :
            mismatch += 1
    return mismatch * 1. / len(seq1)

### ** seqDistances(seqList, distances = None)

def seqDistances(seqList, distances = None, minDistance = None) :
    """Calculate the distances between all pairs of sequences in seqList

    Args:
        seqList (list of str): List of sequences
        distances (dict): A previous output from :func:`seqDistances`, from
          which precalculated distances can be taken
        minDistance (float): Previous minimum distance valuex

    Returns:
        tuple (dict, float): Dictionary mapping frozen sets {s1, s2} to the 
          distance between s1 and s2, where {s1, s2} are all the possible sets 
          of distinct sequences; minimum distance observed, 1 if no 
          pairs of sequences.

    """
    distances_old = distances
    if distances_old is None :
        distances_old = dict()
    if minDistance is None :
        minDistance = 1
    distances = dict()
    for i in seqList :
        for j in seqList :
            if i != j :
                if not distances.get(frozenset([i, j]), False) :
                    if distances_old.get(frozenset([i, j]), False) :
                        distances[frozenset([i, j])] = distances_old[frozenset([i, j])]
                    else :
                        d = seqDistance(i, j)
                        distances[frozenset([i, j])] = d
                        if d < minDistance :
                            minDistance = d
    return (distances, minDistance)

### ** SeqDistanceCompo(compo1, compo2)

def seqDistanceCompo(compo1, compo2) :
    """Calculate the distance between compositions of protein sequences.  The
    distance is calculated as the sum of the absolute differences in count for
    each character, divided by 2 * (seq length). This is not a correct to
    measure distance between sequences, but it should give a proxy of a lower
    boundary for the actual distance.
    
    Args:
        compo1 (Counter): Composition of sequence 1
        compo2 (Counter): Composition of sequence 2

    Return:
        float
    
    """
    assert sum(compo1.values()) == sum(compo2.values())
    characters = set(compo1.keys() + compo2.keys())
    d = 0
    for c in characters :
        d += abs(compo1.get(c, 0) - compo2.get(c, 0))
    return d * 1. / (2 * sum(compo1.values()))

### ** compoNodeDistances(nodes, distances = None)

def compoNodeDistances(nodes, distances = None, minDistance = None) :

    """Calculate the distances between all pairs of sequences in seqList
    This is a fast implementation, looking not at the sequence itself but
    at the composition of each sequence. It is extremely rough and can merge
    sequences which are very different but have similar compositions.

    Args:
        nodes (list of Node): List of Node namedtuples (defined in 
          mergeSequencesCompo)
        distances (dict): A previous output from :func:`compoNodeDistances`, 
          from which precalculated distances can be taken
        minDistance (float): Previous minimum distance valuex

    Returns:
        tuple (dict, float): Dictionary mapping frozen sets {s1, s2} to the 
          distance between s1 and s2, where {s1, s2} are all the possible sets 
          of distinct node ids; minimum distance observed, 1 if no 
          pairs of sequences.

    """
    distances_old = distances
    if distances_old is None :
        distances_old = dict()
    if minDistance is None :
        minDistance = 1
    distances = dict()
    k = 0
    for i in nodes :
        for j in nodes :
            k += 1
            if i.id != j.id :
                if not distances.get(frozenset([i.id, j.id]), False) :
                    if distances_old.get(frozenset([i.id, j.id]), False) :
                        distances[frozenset([i.id, j.id])] = distances_old[frozenset([i.id, j.id])]
                    else :
                        characters = set(i.counter.keys() + j.counter.keys())
                        d = sum([abs(i.counter.get(x, 0) - j.counter.get(x, 0)) for x in characters])
                        d = d * 1. / (2 * sum(i.counter.values()))
                        assert d <= 1
                        distances[frozenset([i[0], j[0]])] = d
                        if d < minDistance :
                            minDistance = d
            if k % 100000 == 0 :
                sys.stderr.write("Distance no " + str(k) + "\n")
    return (distances, minDistance)

### ** nodeAverageCompo(nodeList)

def nodeAverageCompo(nodeList) :
    """Determine the merged Node from a list of Node

    Args:
        modeList (list of Node): List of Node namedtuples (defined in 
          mergeSequencesCompo)

    Return:
        Node: Merged Node with id = "NA", counter = average of input Node 
          counters and seqs = merged list of seqs from input Node

    """
    Node = collections.namedtuple("Node", ["id", "counter", "seqs"])
    averageCounter = dict()
    characters = list()
    seqs = list()
    for n in nodeList :
        characters += list(n.counter.keys())
        seqs += n.seqs
    characters = set(characters)
    seqs = list(set(seqs))
    for c in characters :
        averageCount = sum([x.counter.get(c, 0) for x in nodeList]) * 1. / len(nodeList)
        averageCounter[c] = averageCount
    return Node(id = "NA", counter = averageCounter, seqs = seqs)        
    
### ** seqConsensus(seq1, seq2)

def seqConsensus(seq1, seq2) :
    """Determine the consensus between two sequences, replacing mismatches by X

    Args:
        seq1 (str): First sequence
        seq2 (str): Second sequence (same length as the first sequence)
    """
    o = ""
    for (x, y) in zip(seq1, seq2) :
        if x == y :
            o += x
        else :
            o += "X"
    return o

### ** seqConsensusList(seqList)

def seqConsensusList(seqList) :
    """Determine the consensus between sequences, replacing mismatches by X

    Args:
        seqList (list of str): Sequences

    Returns:
        str: The consensus sequence

    """
    o = ""
    assert len(seqList) > 0
    length = len(seqList[0])
    assert all([len(x) == length for x in seqList])
    for i in range(length) :
        chars = set([x[i] for x in seqList])
        if len(chars) == 1 :
            o += seqList[0][i]
        else :
            o += "X"
    return o

### ** groupSets(setList)

def groupSets(setList) :
    """Merge sets into larger groups when elements are shared between them.

    Args:
        setList (list of set): List of sets

    Returns:
        list of set: List of sets with sets with shared elements merged
    """
    bag = set([frozenset(x) for x in setList])
    groups = set()
    while len(bag) > 0 :
        previousSize = 0
        element = set(bag.pop())
        while previousSize != len(element) :
            previousSize = len(element)
            overlapping = [x for x in bag if len(element & x) > 0]
            [bag.remove(x) for x in overlapping]
            [element.update(x) for x in overlapping]
        group = frozenset(element)
        groups.add(group)
    return list(groups)

### ** mergeSequences(sequences, maxDistance, stderr = None)

def mergeSequences(sequences, maxDistance, stderr = None) :
    """We can imagine a dendrogram based on the distance between sequences.

    Terminal nodes are the original sequences, and intermediates nodes are
    consensus sequences of their children nodes (parent to the root, children
    towards the leaves) where mismatches between children are replaced by "X".

    One intermediate node can have more than two children.

    We call "simplified sequences" a set of nodes so that every leaf has
    exactly one of these nodes among its parents. At every instant, there is a
    set of current "simplified nodes" and each leaf is assigned to one of them.

    The starting state is with "simplified sequences" being the leaves.

    The target state is with "simplified sequences" being only one or with
    distance between them greater than the maxDistance value. We then return
    the mapping between the leaves and the simplified nodes.

    We proceed iteratively:
    1. Check that there is more than one "simplified node"
    2. Calculate all the distances between the current "simplified nodes"
    3. Determine the minimum distance and check that it is not greater than the
       threshold
    4. Select all the pairs with this minimum distance
    5. Make groups of connected pairs among those
    6. For each independent group of connected pairs, merge them:
         a. Determine the simplified sequence from those simplified sequences.
            This is a newly determined node in the dendrogram.
         b. Map the leaves which were mapped to the simplified sequences of the
            group to this new node
         c. Remove the simplified sequences of the group from the current set 
            of simplified sequences and add the new simplified sequence

    Args:
        sequences (iterable of str): Original sequences
        maxDistance (float): Maximum dissimilarity between sequences to be
          merged
        stderr (file stream): Optional, stderr for verbosity

    Returns:
        dict: Dictionary mapping (original sequence, simplified sequence)
    
    """
    closeStderr = False
    if stderr is None :
        stderr = open(os.devnull, "w")
        closeStderr = True
    stderr.write("Prepare nodes. ")
    leaves = list(sequences)
    stderr.write(".")
    simpleNodes = list(sequences)
    stderr.write(".")
    mapping = dict(zip(leaves, leaves))
    stderr.write(".\n")
    stderr.write("Calculating distances. ")
    (distances, minDistance) = seqDistances(simpleNodes)
    stop = (len(simpleNodes) < 2) or (minDistance > maxDistance)
    while not stop :
        stderr.write("Merging " + str(len(simpleNodes)) + " sequences. ")
        bestPairs = [x for x in distances.keys() if distances[x] == minDistance]
        mergingGroups = groupSets(bestPairs)
        for group in mergingGroups :
            simpleNode = seqConsensusList(list(group))
            for leaf in leaves :
                if mapping[leaf] in group :
                    mapping[leaf] = simpleNode
            [simpleNodes.remove(x) for x in group]
            simpleNodes.append(simpleNode)
        stderr.write("Calculating distances. ")
        (distances, minDistance) = seqDistances(simpleNodes, distances)
        stop = (len(simpleNodes) < 2) or (minDistance > maxDistance)
    stderr.write("\n")    
    if closeStderr :
        stderr.close()
    return mapping

### ** mergeSequencesCompo(sequences, maxDistance, stderr = None)

def mergeSequencesCompo(sequences, maxDistance, stderr = None) :
    """NOTE: THIS FUNCTION USES THE COMPOSITION OF EACH SEQUENCE TO CALCULATE
    DISTANCES, NOT THE ACTUAL SEQUENCE ITSELF. USE ONLY IF YOU KNOW WHAT YOU
    ARE DOING!

    YOU SHOULD PROBABLY USE THE mergeSequences FUNCTION, NOT THIS ONE.

    THIS FUNCTION CALCULATES DISTANCES BASED ONLY ON SEQUENCE COMPOSITION, NOT
    ON THE ACTUAL ARRANGEMENT OF CHARACTERS IN THE STRING! THIS MEANS THAT TWO
    SEQUENCES WHICH ARE VERY DIFFERENT BUT HAVE THE SAME COMPOSITION WILL BE
    MERGED TOGETHER.

    We can imagine a dendrogram based on the distance between composition of
    sequences.

    Terminal nodes are the original sequences, and intermediates nodes are just
    described by the average composition of their children nodes (parent to the
    root, children towards the leaves).

    One intermediate node can have more than two children.

    One node is described by its composition (a collections.Counter object) and
    by a list of sequences it contains. This seems the downstream structure of
    a given node is not known, except for the leaves (sequences) it contains.

    We call "simplified nodes" a set of nodes so that every leaf has
    exactly one of these nodes among its parents. At every instant, there is a
    set of current "simplified nodes" and each leaf is assigned to one of them.

    The starting state is with "simplified nodes" being the leaves.

    The target state is with "simplified nodes" being only one or with
    distance between them greater than the maxDistance value. We then return
    the mapping between the leaves and the simplified nodes.

    We proceed iteratively:
    1. Check that there is more than one "simplified node"
    2. Calculate all the distances between the current "simplified nodes"
    3. Determine the minimum distance and check that it is not greater than the
       threshold
    4. Select all the pairs with this minimum distance
    5. Make groups of connected pairs among those
    6. For each independent group of connected pairs, merge them:
         a. Determine the simplified node composition from those simplified nodes.
            This is a newly determined node in the dendrogram.
         b. Map the leaves (sequences) which were mapped to the simplified nodes
            of the group to this new node
         c. Remove the simplified nodes of the group from the current set 
            of simplified nodes and add the new simplified node

    Args:
        sequences (iterable of str): Original sequences
        maxDistance (float): Maximum value for the estimate of minimum 
          dissimilarity between sequences to be merged
        stderr (file stream): Optional, stderr for verbosity

    Returns:
        dict: Dictionary mapping (original sequence, simplified node id)

    """
    warnings.warn("Using mergeSequencesCompo, which does NOT give correct "
                  "estimates for the distances between sequences. YOU "
                  "PROBABLY SHOULDN'T DO THAT! Use mergeSequences "
                  "instead.")
    closeStderr = False
    if stderr is None :
        stderr = open(os.devnull, "w")
        closeStderr = True
    stderr.write("Prepare nodes")
    Node = collections.namedtuple("Node", ["id", "counter", "seqs"])
    nodeId = 0
    nodes = dict()
    for seq in sequences :
        newNode = Node(nodeId, collections.Counter(seq), [seq])
        nodes[nodeId] = newNode
        nodeId += 1
    stderr.write(".\n")
    stderr.write("Calculating distances (" + str(len(nodes)) + "seqs)\n")
    (distances, minDistance) = compoNodeDistances(nodes.values())
    stop = (len(nodes) < 2) or (minDistance > maxDistance)
    while not stop :
        stderr.write("Merging\n")
        stderr.write("Min distance: " + str(minDistance) + "\n")
        bestPairs = [x for x in distances.keys() if distances[x] == minDistance]
        print(bestPairs)
        mergingGroups = groupSets(bestPairs)
        for group in mergingGroups : # group is a set of node ids
            print(group)
            newNode = nodeAverageCompo([nodes[x] for x in group]) # Return a Node namedtuple
            newNode2 = Node(id = nodeId, counter = newNode.counter, seqs = newNode.seqs)
            [nodes.pop(x) for x in group]
            nodes[nodeId] = newNode2
            nodeId += 1
        stderr.write("Calculating distances (" + str(len(nodes)) + "seqs)\n")
        (distances, minDistance) = compoNodeDistances(nodes.values(), distances)
        stop = (len(nodes) < 2) or (minDistance > maxDistance)
    if closeStderr :
        stderr.close()
    mapping = dict()
    for node in nodes.values() :
        for seq in node.seqs :
            mapping[seq] = node.id
    return mapping

### ** md5hash(string)

def md5hash(string) :
    """Calculate the md5 hash for a string

    Args:
        string (str): Input string

    Returns:
        str: Md5 hash for the input string

    """
    h = hashlib.md5()
    h.update(string)
    return h.hexdigest()

### ** invertDict(inputDict)

def invertDict(inputDict) :
    """Build the reciprocal dictionary from an input dictionary, i.e. from a
    mapping (k, v) where all k are unique a mapping (v, [k]) where [k] is the 
    list of keys mapping to the same value v in the original dictionary

    Args:
        inputDict (dict): Input dictionary (k, v)

    Returns:
        dict: Dictionary (v, [k])

    """
    o = dict()
    for (k, v) in inputDict.items() :
        o[v] = o.get(v, [])
        o[v].append(k)
    return o

### ** mergeSequencesOld(sequences, maxDistance)

def mergeSequencesOld(sequences, maxDistance) :

    """Merge biological sequences of same length based on their similarity

    Args:
        sequences (iterable of str): List or set of strings
        maxDistance (float): Maximum distance allowed to merge sequences

    Returns:
        dictionary: Mapping between the original sequences and the merged 
          sequences
    """
    # https://en.wikipedia.org/wiki/Hierarchical_clustering
    # www.ijetae.com/files/Volume2Issue5/IJETAE_0512_48.pdf (Rambabu 2012, IJETAE,
    # "Clustering Orthologous Protein Sequences thru Python Based Program")
    #
    # Initialization
    sequences = list(sequences)
    clusters = set(sequences)
    assert len(clusters) > 2
    ancestors = dict()
    descendants = dict()
    for c in clusters :
        ancestors[c] = c
    for c in clusters :
        descendants[c] = c
    distances = dict()
    for i in clusters :
        for j in clusters :
            if i != j and not distances.get(frozenset([i, j]), False) :
                distances[frozenset([i, j])] = seqDistance(i, j)
    # Make clusters
    done = False
    while (not done) :
        sortedDistances = sorted(distances.items(), key = lambda x: x[1])
        if sortedDistances[0][1] > maxDistance :
            done = True
        else :
            # Merge
            toMerge = [x[0] for x in sortedDistances if x[1] == sortedDistances[0][1]]
            for pair in toMerge :
                distances.pop(pair)
                pair = list(pair)
                try  :
                    clusters.remove(descendants[pair[0]])
                except :
                    pass
                try :
                    clusters.remove(descendants[pair[1]])
                except :
                    pass
                newCluster = seqConsensus(descendants[pair[0]], descendants[pair[1]])
                clusters.add(newCluster)
                for ancestor in ancestors[pair[0]] :
                    descendants[ancestor] = newCluster
                for ancestor in ancestors[pair[1]] :
                    descendants[ancestor] = newCluster
                descendants[newCluster] = newCluster
                ancestors[newCluster] = ancestors.get(newCluster, [])
                ancestors[newCluster] += list(pair)
            # Refresh distances
            for i in clusters :
                for j in clusters :
                    if i != j and not distances.get(frozenset([i, j]), False) :
                        distances[frozenset([i, j])] = seqDistance(i, j)
            if len(clusters) < 2 :
                done = True
    # Return
    return dict(zip(sequences, [descendants[x] for x in sequences]))
    
### ** extractCodingSeqFast(CDS, seqRecord)

def extractCodingSeqFast(CDS, seqRecord) :
    """Helper function to get the CDS sequence faster than with the extract
    method of a CDS object.
    Note: This will NOT cope with complex locations. It will just extract the 
    nucleotide sequence from the start to the end positions of the CDS location
    and possibly reverse-complement it if needed. Use extractCodingSeqReliable 
    for a safer extraction.
    
    Return "None" if CDS has no "translation" qualifier.

    Args:
        CDS (Bio.SeqFeature.SeqFeature): CDS object
        seqRecord (Bio.SeqRecord.SeqRecord): Original record to extract the 
          nucleotide from
    """
    if not CDS.qualifiers.get("translation", False) :
        return "None"
    if len(CDS.location.parts) > 1 :
        warnings.warn("CDS with complex location, using "
                      "extractCodingSeqReliable()\n" +
                      str(CDS.location))
        seq = extractCodingSeqReliableSeq(CDS, seqRecord)
    else :
        seq = seqRecord.seq[CDS.location.start:CDS.location.end]
        if CDS.location.strand == -1 :
            seq = seq.reverse_complement()
    seq += "NN" # to enable translation of final codon if unambiguous
    # Test for translation
    checked = False
    i = 0
    while (not checked and i < 10) :
        translated = str(seq[i: ].translate(table = 11))
        i += 1
        expected = CDS.qualifiers["translation"][0]
        exp_len = len(expected) * 3
        if compareExpPredProt(translated, expected) :
            checked = True
            seq = seq[(i-1):]
    if not checked :
        sys.stderr.write("Expected: " + expected + "\n")
        sys.stderr.write("Predictd: " + translated.strip("*") + "\n")
    seq = seq[0:exp_len]
    if not str(seq.translate(table= 11)[1:]).replace("*", "U") == expected[1:] :
        print(seq.translate(table= 11))
        print(expected)
        sys.exit(1)
    return str(seq)

### ** compareExpPredProt(exp, pred)

def compareExpPredProt(exp, pred) :
    """Compare two protein sequences

    Args:
        exp (str): Expected sequence (e.g. from qualifiers["translation"][0])
        pred (str): Predicted sequence (e.g. from Seq.translate())

    Returns:
        bool: True if sequence identical after removing trailing stop symbols,
          removing the first aa and checking for selenocysteins in the sequence

    """
    exp = exp.rstrip("*X")
    pred = pred.rstrip("*X")
    if exp[1:] == pred[1:] :
        return True
    elif len(exp) != len(pred) :
        return False
    else :
        mismatches = [(i, j) for (i, j) in zip(exp, pred) if i != j]
        ok = True
        for m in mismatches :
            if m != ("U", "*") and m != ("*", "U") :
                ok = False
        if ok :
            return True
        else :
            return False

### ** extractCodingSeqReliable(CDS, seqRecord)

def extractCodingSeqReliable(CDS, seqRecord) :
    """Helper function to get the CDS sequence for a CDS object.
    Note: This will cope with complex locations, but is slow. See 
    extractCodingSeqFast for an alternative for simple locations.
    
    Args:
        CDS (Bio.SeqFeature.SeqFeature): CDS object
        seqRecord (Bio.SeqRecord.SeqRecord): Original record to extract the 
          nucleotide from
    """
    return str(CDS.extract(seqRecord).seq)

### ** extractCodingSeqReliableSeq(CDS, seqRecord)

def extractCodingSeqReliableSeq(CDS, seqRecord) :
    """Helper function to get the CDS sequence for a CDS object.
    Note: This will cope with complex locations, but is slow. See 
    extractCodingSeqFast for an alternative for simple locations.
    
    This function returns a Seq object, not a string.

    Args:
        CDS (Bio.SeqFeature.SeqFeature): CDS object
        seqRecord (Bio.SeqRecord.SeqRecord): Original record to extract the 
          nucleotide from
    """
    return CDS.extract(seqRecord).seq

### ** revComp(DNAstring)

def revComp(DNAstring) :
    """Determine the reverse complement of a DNA sequence

    Args:
        DNAstring (Str): string with A, T, G, C, X, -

    Returns:
        str: Reverse complement

    """
    rev = dict(zip(["A", "T", "G", "C", "X", "-"],
                   ["T", "A", "C", "G", "X", "-"]))
    return("".join([rev[x] for x in reversed(DNAstring)]))

### ** locStr2int(locationString, ignoreFuzzy = False)

def locStr2int(locationString, ignoreFuzzy = False) :
    """Convert a location string from a gene entry into a sequence of
    integers (positions in the original record, zero-based like lists)

    Args:
        locationString (str): Location string (e.g. "[2093106:2093455](-)")
        ignoreFuzzy (bool): If True, do not take into account "<" and ">" in 
          fuzzy location (see Biopython FeatureLocation)

    Returns:
        list of (int, strand): List of positions so that 
          [record.seq[x[0]] for x in list] produces the corresponding coding 
          sequence and strand ("+" or "-")

    """
    l = locationString
    if l.startswith("join{") and l.endswith("}") :
        # Split the string into elemental locations
        l = l[5:-1].split(", ")
        o = []
        for sublist in l :
            o += locStr2int(sublist, ignoreFuzzy)
        return o
    if ignoreFuzzy :
        start = int(locationString.split(":")[0].lstrip("[").lstrip("<>"))
        end = int(locationString.split(":")[1].split("]")[0].lstrip("<>"))
    else :
        start = int(locationString.split(":")[0].lstrip("["))
        end = int(locationString.split(":")[1].split("]")[0])
    sense = locationString[-2]
    indices = list(range(start, end))
    if sense == "-" :
        indices.reverse()
    return [(x, sense) for x in indices]

### ** buildGeneTableFileIndex(filename)

def buildGeneTableFileIndex(filename) :
    """Make one pass through a gene table file to build an index mapping 
    peptide hash, last file position of corresponding record and length

    Args:
        filename (str): Path to the input file

    Returns:
        dict: Dictionary (peptideHash, (lastFilePos, peptideLength))

    """
    o = dict()
    with open(filename, "r") as fi :
        headers = fi.readline().lstrip("#").strip().split("\t")
        hashI = headers.index("peptideHash")
        lengthI = headers.index("peptideLength")
        (pos, line) = (fi.tell(), fi.readline().strip())
        while line :
            elements = line.split("\t")
            if (len(elements) > 0) :
                o[elements[hashI]] = (pos, elements[lengthI])
            (pos, line) = (fi.tell(), fi.readline().strip())
    return o

### ** buildGeneTableFileIndexField(filename, field)

def buildGeneTableFileIndexField(filename, field) :
    """Make one pass through a gene table file to build an index mapping 
    one of the field with the file positions of corresponding gene entries

    Args:
        filename (str): Path to the input file
        field (str): Header of the column to use as a field (e.g. "mergedPeptideHash")

    Returns:
        dict: Dictionary (peptideHash, [filePositions])

    """
    o = dict()
    with open(filename, "r") as fi :
        headers = fi.readline().lstrip("#").strip().split("\t")
        fieldI = headers.index(field)
        (pos, line) = (fi.tell(), fi.readline().strip())
        while line :
            elements = line.split("\t")
            if (len(elements) > 0) :
                o[elements[fieldI]] = o.get(elements[fieldI], [])
                o[elements[fieldI]].append(pos)
            (pos, line) = (fi.tell(), fi.readline().strip())
    return o

### ** buildGeneTableMergedHashDict(filename)

def buildGeneTableMergedHashDict(filename) :
    """Build a dictionary mapping merged peptide hash and protein sequences

    Args:
        filename (str): Path to the input file

    Returns:
        dict: Dictionary (peptideHash, [(geneId, protSeq)])

    """
    o = dict()
    with open(filename, "r") as fi :
        headers = fi.readline().lstrip("#").strip().split("\t")
        mergedHashI = headers.index("mergedPeptideHash")
        geneIdI = headers.index("geneId")
        protSeqI = headers.index("peptideSeq")
        (pos, line) = (fi.tell(), fi.readline().strip())
        while line :
            elements = line.split("\t")
            if (len(elements) > 0) :
                o[elements[mergedHashI]] = o.get(elements[mergedHashI], [])
                o[elements[mergedHashI]].append((elements[geneIdI],
                                            elements[protSeqI]))
            (pos, line) = (fi.tell(), fi.readline().strip())
    return o

### ** buildLengthFileIndex(hashFileIndex)

def buildLengthFileIndex(hashFileIndex) :
    """Build a dictionary mapping length to a list of peptide hash and file
    position from the output of buildGeneTableFileIndex

    Args:
        hashFileIndex (dict): Mapping (peptideHash, (lastFilePos, 
          peptideLength)), output from buildGeneTableFileIndex

    Returns:
        dict: Dictionary (length, [(hash, pos)])
    
    """
    o = dict()
    for (k, v) in hashFileIndex.items() :
        length = int(v[1])
        o[length] = o.get(length, [])
        o[length].append((k, v[0]))
    return o

### ** gatherSequences(filename, lengthIndex, targetLength)

def gatherSequences(filename, lengthIndex, targetLength) :
    """Collect the peptide sequences of a given length from a gene table
    file.

    Args:
        filename (str): Path to the gene table file
        lengthIndex (dict): Dictionary (length, [(hash, pos)]), output from
          buildLengthFileIndex
        targetLength (int): Length for which to extract sequences

    Returns:
        list of str: List of sequences of target length

    """
    o = []
    with open(filename, "r") as fi :
        headers = fi.readline().lstrip("#").strip().split("\t")
        sequenceI = headers.index("peptideSeq")
        targets = lengthIndex.get(targetLength, [])
        for target in targets :
            fi.seek(target[1])
            line = fi.readline().strip().split("\t")
            o.append(line[sequenceI])
    return o

### ** loadLightGeneTable(inputFile, fields)

def loadLightGeneTable(inputFile, fields) :
    """Load a light version of a gene table, as a dict indexed by geneId

    Args:
        inputFile (str): Gene table file name
        fields (list of str): List of fields to extract from the table

    Returns:
        dict: Mapping between geneId and a list of values corresponding to the
          requested fields
    """
    o = dict()
    with open(inputFile, "r") as fi :
        headers = fi.next().strip("#").strip().split("\t")
        headers_i = [headers.index(f) for f in fields]
        geneId_i = headers.index("geneId")
        for l in fi :
            l = l.strip().split("\t")
            o[l[geneId_i]] = [l[x] for x in headers_i]
    return o

### ** refToStr

def refToStr(refList) :
    """Convert an EMBL reference list to a string

    Args:
        refList (list): record.annotations["references"] data from an EMBL 
          record

    Returns:
        str: A string containing all the reference information
    """
    return "<REFSEP>\n".join([str(x) for x in refList])

### ** dbxrefsToBiosample

def dbxrefsToBiosample(dbxrefs, defaultBiosample = None) :
    """Extract the biosample information from a list of database cross references
    of an EMBL record

    Args:
        dbxrefs (list): record.dbxrefs from an EMBL record
        defaultBiosample (str): Return value if not exactly one BioSample value 
          is found. If None, throw an exception.

    Returns:
        str: The biosample identity
    """
    if defaultBiosample is None :
        assert sum([x.startswith("BioSample:") for x in dbxrefs]) == 1
    if sum([x.startswith("BioSample:") for x in dbxrefs]) == 1 :
        return [x.split(":")[1] for x in dbxrefs if x.startswith("BioSample:")][0]
    else :
        print("Not exactly one BioSample identity for " + str(dbxrefs))
        return defaultBiosample

### ** EMBLFileInfo

def EMBLFileInfo(emblFile, defaultBiosample = "NA") :

    """Extract basic information from an EMBL file. Returns None if there
    was a problem with the file parsing (e.g. empty file).

    Args:
        emblFile (str): Path to an EMBL file (not compressed)
        defaultBiosample (str): Default string when no BioSample information is 
          present. If set to None, throw an error if BioSample is missing.

    Returns:
        dict: A dictionary with filename, organism, biosample, references,
          nRecords
    """
    entry = SeqIO.parse(emblFile, "embl")
    try :
        record = entry.next()
    except :
        print("Error with " + emblFile)
        return None
    nRecords = 1
    filename = os.path.basename(emblFile)
    organism = record.annotations["organism"]
    references = refToStr(record.annotations["references"])
    biosample = dbxrefsToBiosample(record.dbxrefs, defaultBiosample)
    for r in entry :
        r_organism = record.annotations["organism"]
        r_references = refToStr(record.annotations["references"])
        r_biosample = dbxrefsToBiosample(record.dbxrefs, defaultBiosample)
        assert all([organism == r_organism,
                    references == r_references,
                    biosample == r_biosample])
        nRecords += 1
    o = dict()
    o["organism"] = organism
    o["references"] = references
    o["biosample"] = biosample
    o["nRecords"] = nRecords
    o["filename"] = filename
    return o    

### ** openEMBLfile(emblFilename)

def openEMBLfile(emblFilename) :
    """Open an EMBL file (can be compressed)

    Args:
        emblFilename (str): Path to an EMBL file (can be gzipped)

    Returns:
        list: A list of the SeqRecord objects from the EMBL file
    """
    if not emblFilename.endswith(".gz") :
        EMBLRecord = SeqIO.parse(emblFilename, "embl")
        return list(EMBLRecord)
    else :
        with gzip.open(emblFilename, "r") as EMBLRecordGz :
            EMBLRecord = SeqIO.parse(EMBLRecordGz, "embl")
            return list(EMBLRecord)

### ** parseEMBLtoSQL(emblFile, SQLiteCursor)

def parseEMBLtoSQL(emblFile, SQLiteCursor) :
    """Load the content of an EMBL file into an SQLite database

    Args:
        emblFile (str): Path to the EMBL file (compressed or not)
        SQLiteCursor (sqlite3.Cursor): Cursor to a connected SQLite database

    """
    EMBLRecord = openEMBLfile(emblFile)
    filename = os.path.basename(emblFile)
    if filename.endswith(".gz") :
        filename = filename[-3:]
    rows = list()
    # Collect the data
    for record in EMBLRecord :
        allCDS = [x for x in record.features if x.type == "CDS"]
        for CDS in allCDS :
            recordId = record.id
            assert len(CDS.qualifiers.get("translation", ["null"])) == 1
            pepSeq = CDS.qualifiers.get("translation", ["null"])[0]
            if pepSeq == "null" :
                pepLen = "null"
            else :
                pepLen = len(pepSeq)
                          
            location = str(CDS.location)
            row = (recordId,
                   pepSeq,
                   extractCodingSeqFast(CDS, record),
                   pepLen,
                   location,
                   CDS.qualifiers.get("transl_table", ["null"])[0],
                   CDS.qualifiers.get("gene", ["null"])[0],
                   CDS.qualifiers.get("product", ["null"])[0],
                   CDS.qualifiers.get("protein_id", ["null"])[0])
            rows.append(row)
    # Store into the database
    SQLiteCursor.executemany("INSERT INTO Cds ("
                             "record, pepSeq, nucSeq, pepLen, location, "
                             "translationTable, geneName, productName, "
                             "productAccNum) VALUES ("
                             "?, ?, ?, ?, ?, ?, ?, ?, ?" 
                             ")", rows)

### * Named tuples

# How to set default values for a named tuple:
# http://stackoverflow.com/questions/11351032/named-tuple-and-optional-keyword-arguments

Gene = collections.namedtuple("Gene", ["recordId", "peptideSeq", "codingSeq",
                                       "peptideLength", "location",
                                       "translationTable", "gene", "product",
                                       "proteinId", "function", "essentiality",
                                       "peptideHash", "mergedPeptideHash",
                                       "geneId"])
Gene.__new__.__defaults__ = ("None", ) * 14

Record = collections.namedtuple("Record", ["recordId", "strainName", "database",
                                           "sequence", "organism", "description",
                                           "references", "length"])
Record.__new__.__defaults__ = ("None", ) * 8

AlnPos = collections.namedtuple("AlnPos", ["geneId", "pepAlnPos", "aa",
                                           "nucAlnPos", "base", "recordPos",
                                           "strand"])
AlnPos.__new__.__defaults__ = ("None", ) * 7

### * Classes

### ** ObjectTable()

class ObjectTable(object) :

    """Parent class for more specific table classes"""

### *** __init__(self)

    def __init__(self) :
        self.items = []
        self.itemType = None

### *** __len__(self)

    def __len__(self) :
        return len(self.items)

### *** __getitem__(self, key)

    def __getitem__(self, key) :
        return self.items[key]

### *** __repr__(self)

    def __repr__(self) :
        return ("<ObjectTable (" + str(len(self)) + " items) for " +
                self.itemType.__doc__ + ">")

### *** _oneCol(self, colName)

    def _oneCol(self, colName) :
        """Return an iterator on the "column" with name colName

        Args:
            colName (str): One of the named attributes of the item type

        Returns:
            iterator
        """
        for x in self :
            yield x.__getattribute__(colName)

### *** col(self, *colNames)

    def col(self, *colNames) :
        """Return an iterator on the "columns" with names colNames

        Args:
            colNames (str): One of the named attributes of the item type

        Returns:
            iterator
        """
        # http://stackoverflow.com/questions/243865/how-do-i-merge-two-python-iterators
        return itertools.izip(*[self._oneCol(x) for x in colNames])

### *** loadTable(self, path)

    def loadTable(self, path) :
        """Load gene information from a tabular file. The first line
        contains the headers.

        Args:
            path (str): Path to the file
        """
        with open(path, "r") as fi :
            headers = fi.readline().strip("\n").strip("#").split("\t")
            for line in fi :
                line = line.strip("\n").split("\t")
                data = dict(zip(headers, line))
                self.items.append(self.itemType(**data))
            
### *** writeTable(self, path)

    def writeTable(self, path) :
        """Write gene information to a tabular file

        Args:
            path (str): Path to the file
        """
        with open(path, "w") as fo :
            headers = list(self.itemType()._asdict().keys())
            fo.write("#" + "\t".join(headers) + "\n")
            for item in self.items :
                itemDict = item._asdict()
                fo.write("\t".join([str(itemDict[x]) for x in headers]) + "\n")

### *** getHeaders(self)

    def getHeaders(self) :
        """Return the headers of the table

        """
        headers = list(self.itemType()._asdict().keys())
        return headers

### ** AlnPosTable() ObjectTable

class AlnPosTable(ObjectTable) :
    """Store a table containing alignment information"""

### *** __init__(self)
    
    def __init__(self) :
        ObjectTable.__init__(self)
        self.itemType = AlnPos

### *** addAlnPos(self, alnPos)

    def addAlnPos(self, alnPos) :
        """Add the information about an alignment position

        Args:
            alnPos (AlnPos): Alignment position data

        """
        self.items.append(alnPos)

### *** addAlnRow(self, geneId, pepAln, geneEntry, ignoreFuzzy = False)

    def addAlnRow(self, geneId, pepAln, geneEntry, ignoreFuzzy = False) :
        """Parse the alignmnent information from a peptide alignment row and 
        related Gene and Record entries and add it to the AlnPos table

        Args:
            geneId (str): gene identifier (used to search Gene entry)
            pepAln (str): Peptide aligned sequences (with gaps if any)
            geneEntry (Gene named tuple): Gene entry
            ignoreFuzzy (bool): If True, ignore fuzzy locations in locStr2int
              calls

        """
        assert pepAln.replace("-", "") == geneEntry.peptideSeq
        ntSeq = geneEntry.codingSeq
        ntAln = ""
        k = 0
        for l in pepAln :
            if l == "-" :
                ntAln += "---"
            else :
                ntAln += ntSeq[k:k+3]
                k += 3
        recordPos = locStr2int(geneEntry.location, ignoreFuzzy = ignoreFuzzy)
        assert len(recordPos) == len(ntSeq)
        gaps = 0
        for i in range(len(ntAln)) :
            p = dict()
            p["geneId"] = geneId
            p["pepAlnPos"] = int(i / 3)
            p["aa"] = pepAln[int(i / 3)]
            p["nucAlnPos"] = i
            p["base"] = ntAln[i]
            if ntAln[i] != "-" :
                p["recordPos"] = recordPos[i - gaps][0]
                p["strand"] = recordPos[i - gaps][1]
            else :
                p["recordPos"] = "NA"
                p["strand"] = "NA"
                gaps += 1
            self.items.append(AlnPos(**p))

### *** writeCompactFile(self, outFile)

    def writeCompactFile(self, outFile) :
        """Write an alignment detailled table to a compact file

        Args:
            outFile (str): Output file name

        """
        assert len(self) % 3 == 0
        with open(outFile, "w") as fo :
            currentGeneId = None
            for i in xrange(len(self) / 3) :
                j = 3 * i
                (a,b,c) = (self[j], self[j+1], self[j+2])
                assert a.geneId == b.geneId == c.geneId
                assert a.pepAlnPos == b.pepAlnPos == c.pepAlnPos
                if (a.geneId != currentGeneId) :
                    currentGeneId = a.geneId
                    fo.write(">" + currentGeneId + "\n")
                line = [str(a.pepAlnPos), a.aa]
                if a.aa != "-" :
                    line += [",".join([a.base, b.base, c.base])]
                    line += [",".join([str(a.recordPos), str(b.recordPos),
                                       str(c.recordPos)])]
                    line += [",".join([a.strand, b.strand, c.strand])]
                fo.write("\t".join(line) + "\n")
            
### ** RecordTable() ObjectTable

class RecordTable(ObjectTable) :
    """Store a table containing record information"""
    
### *** __init__(self)

    def __init__(self) :
        ObjectTable.__init__(self)
        self.itemType = Record
        self._recordIdDict = dict()
        self.nRecords = 0
        
### *** refreshRecordIdDict(self)

    def refreshRecordIdDict(self) :
        self._recordIdDict = dict()
        for (i, r) in enumerate(self.items) :
            self._recordIdDict[r.recordId] = i

### *** recordId(self, Id)

    def recordId(self, Id) :
        """Return the Record corresponding to a given id

        """
        try :
            record = self[self._recordIdDict[Id]]
        except KeyError :
            self.refreshRecordIdDict()
            record = self[self._recordIdDict[Id]]
        if not record.recordId == Id :
            self._recordIdDict = dict()
            self.refreshRecordId()
        return self[self._recordIdDict[Id]]

### *** addGenBankRecord(self, gbRecord)

    def addGenBankRecord(self, gbRecord) :
        """Add a GenBank record data

        Args:
            gbRecord (Bio.SeqRecord.SeqRecord): GenBank Record object
        """
        if isinstance(gbRecord, str) :
            gbRecord = SeqIO.read(gbRecord, "genbank")
        d = dict()
        d["recordId"] = "GI:" + gbRecord.annotations["gi"]
        if gbRecord.description.endswith(", complete genome.") :
            d["strainName"] = gbRecord.description[:-18]
        else :
            d["strainName"] = "NA"
        d["database"] = "GenBank"
        d["sequence"] = str(gbRecord.seq)
        d["organism"] = gbRecord.annotations["organism"]
        d["description"] = gbRecord.description
        d["references"] = "<REFSEP>".join([str(x) for x in gbRecord.annotations["references"]]).replace("\n", "<FIELDSEP>")
        d["length"] = len(d["sequence"])
        self.items.append(Record(**d))

### *** streamEMBLRecord(self, EMBLRecord, outFile, headers)

    def streamEMBLRecord(self, EMBLRecord, outFile, headers) :
        """Add an EMBL record data to an output file

        Args:
            EMBLRecord (Bio.SeqRecord.SeqRecord): EMBL Record object
            outFile (file): File handle to output file
            headers (list of str): Headers of the ouput file

        """
        self.nRecords += 1
        msg = "(stream) Parsing record " + str(self.nRecords)
        sys.stderr.write(msg + "\n")
        if isinstance(EMBLRecord, str) :
            if not EMBLRecord.endswith(".gz") :
                EMBLRecordGzToClose = False
                EMBLRecord = SeqIO.parse(EMBLRecord, "embl")
            else :
                EMBLRecordGzToClose = True
                EMBLRecordGz = gzip.open(EMBLRecord, "r")
                EMBLRecord = SeqIO.parse(EMBLRecordGz, "embl")
        for record in EMBLRecord :
            d = dict()
            d["recordId"] = record.id
            if record.description.endswith(", complete genome.") :
                d["strainName"] = record.description[:-18]
            else :
                d["strainName"] = "NA"
            d["database"] = "Ensembl/EMBL"
            d["sequence"] = str(record.seq)
            d["organism"] = record.annotations["organism"]
            d["description"] = record.description
            d["references"] = "<REFSEP>".join([str(x) for x in record.annotations["references"]]).replace("\n", "<FIELDSEP>")
            d["length"] = len(d["sequence"])
            outFile.write("\t".join([str(d.get(x, "None")) for x in headers]) + "\n")
        if EMBLRecordGzToClose :
            EMBLRecordGz.close()
            
### ** GeneTable() ObjectTable

class GeneTable(ObjectTable) :
    """Store a table containing bacterial gene information"""

### *** __init__(self)

    def __init__(self) :
        ObjectTable.__init__(self)
        self.itemType = Gene
        self.stderr = sys.stderr
        self.nParsedRecords = 0
        self.simplifiedSeqs = None
        self._geneIdDict = dict()
        
### *** parseGenBankRecord(self, gbRecord, hashConstructor = hashlib.md5)

    def parseGenBankRecord(self, gbRecord, hashConstructor = hashlib.md5) :
        """Parse the content of a GenBank record

        Args:
            gbRecord (Bio.SeqRecord.SeqRecord): GenBank Record object. It can 
              also be a list of GenBank Record objects, a path to a GenBank
              record or a list of paths to GenBank records.
            hashConstructor (from hashlib module): If not None, used to 
              calculate hash for each peptide sequence
        """
        # TODO: Set simplifiedSeqs to None when new sequences are parsed?
        if isinstance(gbRecord, list) :
            for r in gbRecord :
                self.parseGenBankRecord(r, hashConstructor)
        else :
            self.nParsedRecords += 1
            msg = "Parsing GenBank record " + str(self.nParsedRecords)
            self.stderr.write(msg + "\n")
            if isinstance(gbRecord, str) :
                gbRecord = SeqIO.read(gbRecord, "genbank")
            allCDS = [x for x in gbRecord.features if x.type == "CDS"]
            for CDS in allCDS :
                peptideSeq = ";".join(CDS.qualifiers.get("translation", ["None"]))
                if hashConstructor is not None :
                    h = hashConstructor()
                    h.update(peptideSeq)
                    peptideHash = h.hexdigest()
                else :
                    peptideHash = "None"
                recordId = "GI:" + gbRecord.annotations["gi"]
                location = str(CDS.location)
                h = hashlib.md5()
                h.update(recordId + location)
                geneId = h.hexdigest()
                gene = self.itemType(recordId = recordId,
                  peptideSeq = peptideSeq,
                  peptideHash = peptideHash,
                  peptideLength = str(len(";".join(CDS.qualifiers.get("translation", ["None"])))),
                  codingSeq = extractCodingSeqFast(CDS, gbRecord),
                  location = location,
                  geneId = geneId,
                  translationTable = ";".join(CDS.qualifiers.get("transl_table", ["None"])),
                  gene = ";".join(CDS.qualifiers.get("gene", ["None"])),
                  product = ";".join(CDS.qualifiers.get("product", ["None"])),
                  proteinId = ";".join(CDS.qualifiers.get("protein_id", ["None"])))
                self.items.append(gene)

### *** parseEMBLRecord(self, EMBLRecord, hashConstructor = hashlib.md5)

    def parseEMBLRecord(self, EMBLRecord, hashConstructor = hashlib.md5) :
        """Parse the content of an EMBL record

        Args:
            EMBLRecord (Bio.SeqRecord.SeqRecord): EMBL Record object. It can 
              also be a list of EMBL Record objects, a path to a EMBL
              record or a list of paths to EMBL records.
            hashConstructor (from hashlib module): If not None, used to 
              calculate hash for each peptide sequence

        """
        # TODO: Set simplifiedSeqs to None when new sequences are parsed?
        if isinstance(EMBLRecord, list) :
            for r in EMBLRecord :
                self.parseEMBLRecord(r, hashConstructor)
        else :
            self.nParsedRecords += 1
            msg = "Parsing EMBL record " + str(self.nParsedRecords)
            self.stderr.write(msg + "\n")
            if isinstance(EMBLRecord, str) :
                if not EMBLRecord.endswith(".gz") :
                    EMBLRecordGzToClose = False
                    EMBLRecord = SeqIO.parse(EMBLRecord, "embl")
                else :
                    EMBLRecordGzToClose = True
                    EMBLRecordGz = gzip.open(EMBLRecord, "r")
                    EMBLRecord = SeqIO.parse(EMBLRecordGz, "embl")
            for record in EMBLRecord :
                allCDS = [x for x in record.features if x.type == "CDS"]
                for CDS in allCDS :
                    peptideSeq = ";".join(CDS.qualifiers.get("translation", ["None"]))
                    if hashConstructor is not None :
                        h = hashConstructor()
                        h.update(peptideSeq)
                        peptideHash = h.hexdigest()
                    else :
                        peptideHash = "None"
                    recordId = record.id
                    location = str(CDS.location)
                    h = hashlib.md5()
                    h.update(recordId + location)
                    geneId = h.hexdigest()
                    gene = self.itemType(recordId = recordId,
                      peptideSeq = peptideSeq,
                      peptideHash = peptideHash,
                      peptideLength = str(len(";".join(CDS.qualifiers.get("translation", ["None"])))),
                      codingSeq = extractCodingSeqFast(CDS, record),
                      location = location,
                      geneId = geneId,
                      translationTable = ";".join(CDS.qualifiers.get("transl_table", ["None"])),
                      gene = ";".join(CDS.qualifiers.get("gene", ["None"])),
                      product = ";".join(CDS.qualifiers.get("product", ["None"])),
                      proteinId = ";".join(CDS.qualifiers.get("protein_id", ["None"])))
                    self.items.append(gene)
            if EMBLRecordGzToClose :
                EMBLRecordGz.close()

### *** streamEMBLRecord(self, EMBLRecord, outFile, headers, hashConstructor = hashlib.md5)

    def streamEMBLRecord(self, EMBLRecord, outFile, headers, hashConstructor = hashlib.md5) :
        """Parse the content of an EMBL record and write it to a file instead 
        of storing it.

        Args:
            EMBLRecord (Bio.SeqRecord.SeqRecord): EMBL Record object. It can 
              also be a list of EMBL Record objects, a path to a EMBL
              record or a list of paths to EMBL records.
            hashConstructor (from hashlib module): If not None, used to 
              calculate hash for each peptide sequence
            outFile (file): File handle for output
            headers (list of str): Headers of the output file

        """
        # TODO: Set simplifiedSeqs to None when new sequences are parsed?
        if isinstance(EMBLRecord, list) :
            for r in EMBLRecord :
                self.streamEMBLRecord(r, outFile, headers, hashConstructor)
        else :
            self.nParsedRecords += 1
            msg = "(stream) Parsing EMBL record " + str(self.nParsedRecords)
            self.stderr.write(msg + "\n")
            if isinstance(EMBLRecord, str) :
                EMBLRecord = openEMBLfile(EMBLRecord)
            for record in EMBLRecord :
                allCDS = [x for x in record.features if x.type == "CDS"]
                for CDS in allCDS :
                    peptideSeq = ";".join(CDS.qualifiers.get("translation", ["None"]))
                    if hashConstructor is not None :
                        h = hashConstructor()
                        h.update(peptideSeq)
                        peptideHash = h.hexdigest()
                    else :
                        peptideHash = "None"
                    recordId = record.id
                    location = str(CDS.location)
                    h = hashlib.md5()
                    h.update(recordId + location)
                    geneId = h.hexdigest()
                    gene = self.itemType(recordId = recordId,
                      peptideSeq = peptideSeq,
                      peptideHash = peptideHash,
                      peptideLength = str(len(";".join(CDS.qualifiers.get("translation", ["None"])))),
                      codingSeq = extractCodingSeqFast(CDS, record),
                      location = location,
                      geneId = geneId,
                      translationTable = ";".join(CDS.qualifiers.get("transl_table", ["None"])),
                      gene = ";".join(CDS.qualifiers.get("gene", ["None"])),
                      product = ";".join(CDS.qualifiers.get("product", ["None"])),
                      proteinId = ";".join(CDS.qualifiers.get("protein_id", ["None"])))
                    geneDict = gene._asdict()
                    outFile.write("\t".join([str(geneDict[x]) for x in headers]) + "\n")
                    
### *** makeGeneId(self)

    def makeGeneId(self) :
        """Calculate the gene id for each row. Gene id is the md5 hash of the 
        record id concatenated with the location. At the same time, populate
        self._geneIdDict with a mapping between gene Id and row.

        """
        for (i, g) in enumerate(self.items) :
            h = hashlib.md5()
            h.update(g.recordId + g.location)
            hStr = h.hexdigest()
            geneData = g._asdict()
            geneData["geneId"] = hStr
            self.items[i] = self.itemType(**geneData)
            self._geneIdDict[hStr] = i

### *** refreshGeneIdDict(self)

    def refreshGeneIdDict(self) :
        self._geneIdDict = dict()
        for (i, g) in enumerate(self.items) :
            self._geneIdDict[g.geneId] = i
            
### *** geneId(self, hashId)

    def geneId(self, hashId) :
        """Return the Gene corresponding to a given id

        """
        try :
            gene = self[self._geneIdDict[hashId]]
        except KeyError :
            self.refreshGeneIdDict()
            gene = self[self._geneIdDict[hashId]]
        if not gene.geneId == hashId :
            self._geneIdDict = dict()
            self.refreshGeneId()
        return self[self._geneIdDict[hashId]]

### *** hashPeptides(self, hashConstructor)

    def hashPeptides(self, hashConstructor) :
        """Calculate hash value for each gene peptide sequence

        Args:
            hashConstructor (function): Hash algorithm to be used (from the 
              ``hashlib`` module)
        """
        for (i, g) in enumerate(self.items) :
            h = hashConstructor()
            h.update(g.peptideSeq)
            hStr = h.hexdigest()
            geneData = g._asdict()
            geneData["peptideHash"] = hStr
            self.items[i] = self.itemType(**geneData)
            
### *** extractUniquePeptides(self)

    def extractUniquePeptides(self) :
        """Extract the hash and sequences of unique peptides
        """
        uniquePep = []
        uniqueHash = set([])
        for g in self.items :
            assert g.peptideHash is not None
            if g.peptideHash not in uniqueHash :
                uniqueHash.add(g.peptideHash)
                uniquePep.append((g.peptideHash, g.peptideSeq))
        return uniquePep

### *** writeUniquePeptides(self, path)

    def writeUniquePeptides(self, path) :
        """Write the unique peptides to a fasta file

        Args:
            path (str): Path to the fasta file
        """
        uniquePep = self.extractUniquePeptides()
        with open(path, "w") as fo :
            for pep in uniquePep :
                fo.write(">" + pep[0] + "\n")
                fo.write(pep[1] + "\n")

### *** extractSimplifiedPeptides(self, maxDissimilarity) :

    def extractSimplifiedPeptides(self, maxDissimilarity) :
        """From the unique peptide sequences, produce simplified sequences which 
        result from merging similar sequences together. Only sequences of same 
        length can be merged, based on their dissimilarity.

        This function returns the mapping and sets the ``simplifiedSeqs`` 
        attribute.

        TODO replace md5 by argument for hash constructor

        Args:
            maxDissimilarity (float): Comprised between 0 and 1, maximum 
              dissimilarity for merging

        Returns:
            dict: Dictionary mapping (original peptide, simplified peptide)

        """
        originalSeqs = set([x.peptideSeq for x in self.items])
        originalSeqsByLen = dict()
        for x in originalSeqs :
            l = len(x)
            originalSeqsByLen[l] = originalSeqsByLen.get(l, [])
            originalSeqsByLen[l].append(x)
        mapping = dict()
        for seqs in originalSeqsByLen.values() :
            print("Processing length " + str(len(seqs[0])) +
                  "(" + str(len(seqs)) + " sequences)")
            mapping.update(mergeSequences(seqs, maxDissimilarity))
        self.simplifiedSeqs = mapping
        # Update gene data
        for (i, g) in enumerate(self.items) :
            h = hashlib.md5()
            h.update(mapping[g.peptideSeq])
            geneData = g._asdict()
            geneData["mergedPeptideHash"] = h.hexdigest()
            self.items[i] = self.itemType(**geneData)
        return mapping
            
### *** writeSimplifiedPeptides(self, path)

    def writeSimplifiedPeptides(self, path) :
        """Write the simplified peptides to a fasta file (i.e. consensus 
        sequences with X at polymorphic positions, obtained by merging 
        sequences of equal length if not too dissimilar). The simplified 
        sequences must have been produced by the
        :func:`extractSimplifiedPeptides` beforehand.

        This is a good input for a blastp run.

        Args:
            path (str): Path to the fasta file

        """
        if self.simplifiedSeqs is None :
            raise Exception("The method to simplify sequences must be called "
                            "first (extractSimplifiedPeptides)")
        with open(path, "w") as fo :
            for pep in set(self.simplifiedSeqs.values()) :
                h = hashlib.md5()
                h.update(pep)
                hStr = h.hexdigest()
                fo.write(">" + hStr + "\n")
                fo.write(pep + "\n")

