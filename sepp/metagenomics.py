#from itertools import izip
# from itertools import izip
import os
import re
from sepp.exhaustive_tipp import *
from sepp.alignment import MutableAlignment
from sepp.alignment import _write_fasta
from sepp.config import options
import tempfile
import sys, time, datetime


"""
Collection of functions for metagenomic pipeline for taxonomic classification
Created on June 3, 2014

@author: namphuon
"""

"""
Changes

+ 

Edited by Erin Molloy on January 21, 2017.
"""

global taxon_map, level_map, key_map

global levels
levels = ["species", "genus", "family", "order", "class", "phylum"]
coglist = ['COG0012', 'COG0016', 'COG0018', 'COG0048', 'COG0049',
             'COG0052', 'COG0080', 'COG0081', 'COG0087', 'COG0088',
             'COG0090', 'COG0091', 'COG0092', 'COG0093', 'COG0094',
             'COG0096', 'COG0097', 'COG0098', 'COG0099', 'COG0100',
             'COG0102', 'COG0103', 'COG0124', 'COG0172', 'COG0184',
             'COG0185', 'COG0186', 'COG0197', 'COG0200', 'COG0201',
             'COG0202', 'COG0215', 'COG0256', 'COG0495', 'COG0522',
             'COG0525', 'COG0533', 'COG0541', 'COG0552', 'COG0085']

#TODO Fix parameter passing
#TODO Make taxonomy loading a class
def load_taxonomy(taxonomy_file, lower=True):
    global taxon_map, level_map, key_map
    f = open(taxonomy_file, 'r')
  
    #First line is the keywords for the taxonomy, need to map the keyword to the positional 
    #index of each keyword
    results = f.readline().lower().replace('"','').strip().split(',')
    key_map = dict([(results[i],i) for i in range(0,len(results))])
    
    #Now fill up taxonomy, level maps keep track of what taxa exist at each level, taxon_map  
    #keep track of entire taxonomy
    taxon_map = {}
    level_map = {"species":{}, "genus":{}, "family":{}, "order":{}, "class":{}, "phylum":{}}

    for line in f:
        results = line.replace('"','').strip()
        if (lower):
            results.lower()
        results = results.split(',')
        #insert into taxon map
        taxon_map[results[0]] = results
    
        #insert into level map
        for level in levels:      
            if (results[key_map[level]] == ''):
                continue
            else:
                if (results[key_map[level]] not in level_map[level]):
                    level_map[level][results[key_map[level]]] = {}
                level_map[level][results[key_map[level]]][results[0]]=results[0]
    return (taxon_map, level_map, key_map)


def distribution(classification_files, output_dir):
    global taxon_map, level_map, key_map, levels, level_names
    distribution = {"species": {}, "genus": {}, "family": {}, "order": {}, "class": {}, "phylum": {}}
    total_frags = 0
    for class_input in classification_files:
        class_in = open(class_input, 'r')
        frag_info = {"species": {'unclassified': 1}, "genus": {'unclassified': 1}, "family": {'unclassified': 1},
                     "order": {'unclassified': 1}, "class": {'unclassified': 1}, "phylum": {'unclassified': 1}}
        (old_name, old_rank, old_probability, old_line) = ("", None, -1, "")
        for line in class_in:
            results = line.strip().split(',')
            if (len(results) > 5):
                results = [results[0], results[1], results[2], results[-2], results[-1]]
            (name, id, rank, probability) = (results[0], results[1], results[3], float(results[4]));
            if (rank not in distribution):
                continue
            if (old_name == ""):
                (old_name, old_rank, old_probability, old_line) = (name, rank, probability, line)
            if (name != old_name):
                total_frags += 1
                assert frag_info['phylum']['unclassified'] != 1
                for clade, cladeval in frag_info.items():
                    for clade_name, cnc in cladeval.items():
                        if (clade_name, cnc not in distribution[clade]):
                            distribution[clade][clade_name] = 0
                        distribution[clade][clade_name] += cnc
                frag_info = {"species": {'unclassified': 1}, "genus": {'unclassified': 1},
                             "family": {'unclassified': 1}, "order": {'unclassified': 1}, "class": {'unclassified': 1},
                             "phylum": {'unclassified': 1}}
                (old_name, old_rank, old_probability, old_line) = (name, rank, probability, line)
            if (id not in frag_info[rank]):
                frag_info[rank][id] = 0
            frag_info[rank][id] += probability
            frag_info[rank]['unclassified'] -= probability
        total_frags += 1
        assert frag_info['phylum']['unclassified'] != 1
        for clade, cladeval in frag_info.items():
            for clade_name, cnc in cladeval.items():
                if (clade_name not in distribution[clade]):
                    distribution[clade][clade_name] = 0
                distribution[clade][clade_name] += cnc

    level_names = {1: 'species', 2: 'genus', 3: 'family', 4: 'order', 5: 'class', 6: 'phylum'}
    for level in level_names:
        f = open(output_dir + "/abundance.distribution.%s.csv" % level_names[level], 'w');
        f.write('taxa\tabundance\n')
        lines = []
        for clade, value in distribution[level_names[level]].items():
            name = clade
            if (name != 'unclassified'):
                name = taxon_map[clade][key_map['tax_name']]
            lines.append('%s\t%0.4f\n' % (name, float(value) / total_frags))
        lines.sort()
        f.write(''.join(lines))
        f.close()
    return distribution


def remove_unclassified_level(classifications, level=6):
    global taxon_map, level_map, key_map, levels
    frags = list(classifications.keys())
    for frag in frags:
        if classifications[frag][level] == 'NA':
            del classifications[frag]


def write_classification(class_input, output):
    '''Writes a classification file
    '''
    class_out = open(output, 'w')
    class_out.write("fragment\tspecies\tgenus\tfamily\torder\tclass\tphylum\n");
    keys = list(class_input.keys())
    keys.sort()
    for frag in keys:
        class_out.write("%s\n" % "\t".join(class_input[frag]));
    class_out.close()
  
#Fix problem with NA being unclassified  
def write_abundance(classifications, output_dir, labels=True, remove_unclassified=True):
    global taxon_map, level_map, key_map, levels

    level_abundance = {1: {'total': 0}, 2: {'total': 0}, 3: {'total': 0}, 4: {'total': 0}, 5: {'total': 0},
                       6: {'total': 0}}

    level_names = {1: 'species', 2: 'genus', 3: 'family', 4: 'order', 5: 'class', 6: 'phylum'}
    for lineage in classifications.values():
        # insert into level map
        for level in range(1, 7):
            if (lineage[level] == 'NA'):
                if ('unclassified' not in level_abundance[level]):
                    level_abundance[level]['unclassified'] = 0
                level_abundance[level]['unclassified'] += 1
                level_abundance[level]['total'] += 1
                # continue
            else:
                if (lineage[level] not in level_abundance[level]):
                    level_abundance[level][lineage[level]] = 0
                level_abundance[level][lineage[level]] += 1
                level_abundance[level]['total'] += 1
    for level in level_names:
        f = open(output_dir + "/abundance.%s.csv" % level_names[level], 'w');
        f.write('taxa\tabundance\n')
        lines = []
        for clade in level_abundance[level]:
            if clade == 'total':
                continue
            name = clade
            if labels and name != 'unclassified':
                name = taxon_map[clade][key_map['tax_name']]
            lines.append('%s\t%0.4f\n' % (name, float(level_abundance[level][clade]) / level_abundance[level]['total']))
        lines.sort()
        f.write(''.join(lines))
        f.close()


def generate_classification(class_input, threshold):
    global taxon_map, level_map, key_map, levels
    class_in = open(class_input, 'r')
    level_map_hierarchy = {"species": 0, "genus": 1, "family": 2, "order": 3, "class": 4, "phylum": 5, "root": 6}
    # Need to keep track of last line so we can determine when we switch to new classification
    old_name = "";
    old_probability = 1;
    old_id = "";
    old_rank = "";

    # keep track of all fragment names
    names = {}
    classification = {}
    for line in class_in:
        results = line.strip().split(',')
        if (len(results) > 5):
            results = [results[0], results[1], results[2], results[-2], results[-1]]
        (name, id, rank, probability) = (results[0], results[1], results[3], float(results[4]));
        names[name] = name;
        if (name != old_name):
            # when we switch to new fragment, output last classification for old fragment
            if (old_name != ""):
                lineage = taxon_map[old_id];
                output_line = [old_name]
                for level in levels:
                    clade = lineage[key_map[level]];
                    if (clade == ""):
                        clade = "NA"
                    output_line.append(clade)
                classification[old_name] = output_line
            old_name = name;
            old_rank = "root";
            old_probability = 1;
            old_id = '1';

        # Switch to new rank if the new probability is higher than threshold
        # and our rank is more specific than our original rank
        if (rank in level_map_hierarchy and (level_map_hierarchy[old_rank] > level_map_hierarchy[rank]) and (
            probability > threshold)):
            old_rank = rank
            old_probability = probability
            old_id = id
        # Switch to new rank if the new rank matches old rank but has higher probability
        elif (rank in level_map_hierarchy and (level_map_hierarchy[old_rank] == level_map_hierarchy[rank]) and (
            probability > old_probability)):
            old_rank = rank
            old_probability = probability
            old_id = id

    if old_id in taxon_map:
        lineage = taxon_map[old_id];
        output_line = [old_name]
        for level in levels:
            clade = lineage[key_map[level]];
            if (clade == ""):
                clade = "NA"
            output_line.append(clade)
        classification[name] = output_line
    return classification

# NEED TO FIX ABOVE...


def read_mapping(input, header=False, delimiter='\t'):
    '''Read a mapping file

    Parameters
    ----------

    Returns
    -------
    '''
    d = {}
    with open(input) as f:
        for line in f:
            if (header == True):
                header = False
                continue
            results = line.strip().split(delimiter)
            d[results[0]] = results[1]  # use to be just result and taking whole list
    return d

# THIS IS FASTER...
# Read file
# split file with spaces or 
#>>> import re  # Will be splitting on: , <space> - ! ? :
#>>> filter(None, re.split("[, \-!?:]+", "Hey, you - what are you doing here!?"))
#['Hey', 'you', 'what', 'are', 'you', 'doing', 'here']
  # # %%timeit
#with open('words.txt', 'rb') as lines, open('counts.bin', 'rb') as counts:
#    words = lines.read().split('\n')
#    values = array('d')
#    values.fromfile(counts, 333333)
#    dict(zip(words, values))

def take_reverse_complement(sequence):
    """
    Takes the reverse complement of a sequence

    Input
    -----
    sequence : string

    Returns
    -------
    revcom : string

    NOTE: This function produces the same answer as the following:
    > from Bio.Alphabet import generic_dna
    > from Bio.Seq import Seq
    > x = "TTATATTATTTTTTTATAAAAGTAATACTCAATTTACTTTTCAAAATATTTGATA"
    > x_rc = str(Seq(x, generic_dna).reverse_complement())
    """
    cmap = {'A':'T', 'a':'t', 'T':'A', 't':'a',
            'C':'G', 'c':'g', 'G':'C', 'g':'c',
            '-':'-'}
    revcom = "".join([cmap.get(c, c) for c in sequence[::-1]])
    return revcom


def run_blast(blastn, blastdb, input_file, output_file, nthreads):
    """
    Runs BLAST.

    Parameters
    ----------
    blastn : string
             path/name of the blastn binary
    blastdb : string
               path/name of the BLAST database
    input_file : string
                 name/path of input fasta file
    output_file : string
                  path/name of output BLAST file
    nthreads : int
               number of threads (CPUs) to use in BLAST search

    Returns
    -------
    Nothing, Output is written by BLAST
    """

    # os.system("%s " % blastn + \
    #             "-db %s " % blastdb + \
    #             "-query %s " % input_file + \
    #             "-out %s " % output_file + \
    #             "-num_threads %d " % nthreads + \
    #             # "-max_target_seqs 1 -outfmt 6")
    #             "-max_target_seqs 1 -outfmt \" 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \"")

    #***  Nute removing max_target_seqs; 06/15/2018 ***#
    print("running blast temp: %s" % datetime.datetime.now())
    output_file_temp = output_file + '.temp'
    os.system("%s " % blastn + \
                "-db %s " % blastdb + \
                "-query %s " % input_file + \
                "-out %s " % output_file_temp + \
                "-num_threads %d " % nthreads + \
                "-outfmt \" 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \"")
    #pause for a bit until the file has been created, then continue
    ct=0
    while not os.path.isfile(output_file_temp) and ct <= 30:
        print("temp blast file doesn't exist yet, ct=%s" % ct)
        time.sleep(2)
        if ct == 30:
            raise FileNotFoundError
        else:
            ct += 1

    os.system('awk -F \'\\t\' \'!seen[$1]++\' %s > %s' % (output_file_temp,output_file))
    ct = 0
    while not os.path.isfile(output_file) and ct <= 30:
        time.sleep(2)
        if ct == 30:
            raise FileNotFoundError
        else:
            ct += 1

# Make another version where you treat blast data as
# dictionary and extract from streaming fasta file
# also time this stuff in TIPP...???

def blast_to_markers(args, genes, temp_dir):
    """
    Produce a fasta file for each marker on which TIPP will be run

    Uses BLAST to figure out
        1) which read goes in which marker fasta file
        2) which direction to write the read
        3) whether to trim a non-homologous section of the read (MN added 5/18)

    Parameters
    ----------
    args : 
    temp_dir : 

    Returns
    -------
    bins : python dictionary
           marker gene to number of reads in respective fasta file
    """
    input_file = args.fragment_file.name
    output_dir = args.outdir

    # Produce a BLAST file
    if args.blast_file is None:
        print("Blasting fragments against ALL markers!\n")
        blastn = args.blast.path # + "/blastn"
        # TO DO: Do *NOT* hard code this information!!
        blastdb = args.reference.path + \
                  "/blast/%s/alignment.fasta" % (args.genes)
                  # "/blast/%s/cogs_blast_db" % (args.genes)


        blast_file = output_dir + "/blast.out"
        nthreads = args.cpu
        run_blast(blastn, blastdb, input_file, blast_file, nthreads)
    else:
        blast_file = args.blast_file

    # NOTE: This is less readable but more efficient...
    # THIS SHOULD BE IT'S OWN FUNCTION....
    bins = dict(zip(genes, [0] * len(genes)))

    marker_map = read_mapping(args.reference.path + \
                              "/blast/%s/seq2marker.tab" % (args.genes))

    f = open(input_file, 'r')
    g = open(blast_file, 'r')

    gfiles = {}
    for gene in genes:
        gene_file = temp_dir + "/%s.frags.fas" % gene
        gfiles[gene] = open(gene_file, 'w')

    fforward = True
    gforward = True
    query_name_old = ""
    i = 0
    while(True):
        #print i

        if fforward:
            i = i + 1
            fasta_name = f.readline().replace('>', '').replace('\n', '')
        if gforward:
            blast_data = g.readline()

        # No more data in BLAST file
        if blast_data == "":
            print ("DONE!")
            break

        # [query_name, subject_name, percent_identity, alignment_length,
        #  mismatches, gap_opens, query_start, query_end, subject_start,
        #  subject_end, evalue, bit_score] = blast_data.split('\t')
        [query_name, subject_name, percent_identity, alignment_length,
        mismatches, gap_opens, query_start, query_end, subject_start,
        subject_end, evalue, bit_score, query_length, subject_length] = blast_data.split('\t')

        if fasta_name == query_name:
            query_name_old = query_name
            query_start = int(query_start)
            query_end = int(query_end)
            subject_start = int(subject_start)
            subject_end = int(subject_end)
            query_length = int(query_length)
            subject_length = int(subject_length)

            direction = (query_end - query_start) * \
                        (subject_end - subject_start)
            
            i = i + 1
            seq_data = f.readline().replace('\n', '')

            #trim sequence if sequence extends beyond reference domain:
            if query_start==1 and (subject_end==1 or subject_end==subject_length):
                seq_data = seq_data[:(query_end+1)]
            elif query_end==query_length and (subject_start==1 or subject_start==subject_length):
                seq_data = seq_data[(query_start-1):]

            # take the RC if needed
            if direction < 0:
                seq_data = take_reverse_complement(seq_data)
            
            # Find marker
            try:
                marker = marker_map[subject_name]
            except:
                print (blast_data)
                sys.exit(0)
            bins[marker] += 1

            gfile = gfiles[marker]
            gfile.write(">" + query_name + "\n")
            gfile.write(seq_data  + '\n')

            fforward = True
            gforward = True
        elif query_name_old == query_name:
            fforward = False
            gforward = True
        else:
            i = i + 1
            f.readline()
            gforward = False
            fforward = True

        # FIX THIS TO BE MORE GENERAL...
        if i > 20000000:
            print("PROBLEM - THEORETICAL END OF FASTA - all my files had 5M seqs...\n")
            f.close()
            g.close()
            sys.exit()
    
    f.close()
    g.close()
    for gfile in gfiles.values():
        gfile.close()

    return bins


def build_abundance_profile(args, genes):
    """


    Parameters
    ----------
    args :
    genes :

    Returns
    -------
    Nothing, files are written
    """
    input_file = args.fragment_file.name
    output_dir = args.outdir
    temp_dir = tempfile.mkdtemp(dir=args.tempdir)

    if True:
        # If using shotgun data,
        # use BLAST to bin sequences to specific markers
        bins = blast_to_markers(args, genes, temp_dir)

        if sum(bins.values()) == 0:
            print("Unable to bin any fragments!\n")
            return
        else:
            print("Finished binning the fragments!\n")
    else:
        # If using amplicon data (e.g., 16S),
        # you *still* need BLAST to determine direction of fragment (actually that's not true, amplicon data
        # usually comes pack oriented 3' --> 5')
        sys.exit("Not implemented yet...\n")
  
    # TO DO:
    # Make a class to do this?
    # Or use pandas...?
    #load up taxonomy for 30 marker genes
    # MN: I don't think we can do this just once. It's not necessarily the same taxonomy. I guess it probably was
    #      supposed to be but I screwed it up.
    if (args.genes == 'markers'):
        (taxon_map, level_map, key_map) = load_taxonomy(args.reference.path + '/refpkg/rpsB.refpkg/all_taxon.taxonomy')
    else:
        (taxon_map, level_map, key_map) = load_taxonomy(args.reference.path + '/refpkg/COG0012.refpkg/all_taxon.taxonomy')
    
    #all classifications stored here  
    classifications = {}
    classification_files = []
    
    #Now run TIPP on each fragment    
    gene_name = 'sate'
    
    if (args.genes == 'cogs'):
        gene_name = 'pasta'
        # gene_name = 'upp'   # FOR UPP RUN
    
    #for gene in bins.keys():
    # this is temporary, for our experiments (4/10)
    gene_set = []
    if args.genes == 'markers':
        gene_set = ["rplB", "rplD", "rplE"]
    elif args.genes == 'cogs':
        # gene_set = ["COG0088", "COG0090", "COG0094"]
        gene_set = coglist.copy()

    cmdfile=open(os.path.join(output_dir,'command_list.txt'),'w')
    for gene in gene_set:
        nfrags = bins[gene]
        if nfrags > 0:
            #Get size of each marker
            total_taxa = 0
            if gene_name=='upp':
                with open(args.reference.path + '/refpkg_upp/%s.refpkg/%s.size' % (gene, gene_name), 'r') as f:
                    total_taxa = int(f.readline().strip())
                    # FOR UPP RUN
            else:
                with open(args.reference.path + '/refpkg/%s.refpkg/%s.size' % (gene, gene_name), 'r') as f:
                    total_taxa = int(f.readline().strip())
        
            decomp_size = args.alignment_size
            if (decomp_size > total_taxa):
                decomp_size = int(total_taxa / 2)
        
            cpus = args.cpu
            if (nfrags < cpus):
                cpus = nfrags
            
            extra = ""
            if args.dist == True:
                extra = "-D"
    
            if args.max_chunk_size is not None:
                extra = extra + " -F %d" % args.max_chunk_size
    
            if args.cutoff != 0:
                extra = extra + " -C %f" % args.cutoff

            if gene_name=='upp':
                markerpath = args.reference.path + "/refpkg_upp/%s.refpkg/" % gene
            else:
                markerpath = args.reference.path + "/refpkg/%s.refpkg/" % gene

            cmd = "python3 /projects/tallis/nute/code/sepp-erin/sepp/run_tipp.py " + \
                  "-c %s " % args.config_file.name + \
                  "--cpu %s " % cpus + \
                  "-m %s " % args.molecule + \
                  "-f %s " % (temp_dir + "/%s.frags.fas" % gene) + \
                  "-t %s " % (markerpath + "%s.taxonomy" % gene_name) + \
                  "-a %s " % (markerpath + "%s.fasta" % gene_name) + \
                  "-r %s " % (markerpath + "%s.taxonomy.RAxML_info" % gene_name) + \
                  "-tx %s " % (markerpath + "all_taxon.taxonomy") + \
                  "-txm %s " % (markerpath + "species.mapping") + \
                  "-at %0.2f " % 0.0 + \
                  "-pt %0.2f " % 0.0 + \
                  "-A %d " % decomp_size + \
                  "-P %d " % args.placement_size + \
                  "-p %s " % (temp_dir + "/tipp_%s" % gene) + \
                  "-o %s " % ("tipp_%s" % gene) + \
                  "-d %s -tiny " % (output_dir + "/markers/") + \
                  "%s" % extra
            cmdfile.write(cmd + '\n\n*******\n\n')
            # "-at %0.2f " % args.alignment_threshold + \
            # "-pt %0.2f " % 0.0 + \
            # "-adt %s " % (markerpath + "%s.tree" % gene_name) + \
            print(cmd)
            os.system(cmd)

    # REMOVE ALL OF THIS AND TURN IT INTO A SECOND PYTHON SCRIPT
    # OTHERWISE YOU NEED TO RE-RUN TIPP EVERYTIME YOU JUST WANT TO USE A DIFFERENT
    # PLACEMENT THRESHOLD... 
    # MAYBE THIS SCRIPT IS MORE FOR RUNNING TIPP ON SHOTGUN DATA
    # AND THE OTHER SCRIPT CREATES TAXON IDS / ABUNDANCE PROFILES

    #if (not os.path.exists(output_directory+"/markers/tipp_%s_classification.txt" % gene)):
    #    continue

    #gene_classification = generate_classification(output_directory+"/markers/tipp_%s_classification.txt" % gene,options().placement_threshold)
    #classification_files.append(output_directory+"/markers/tipp_%s_classification.txt" % gene)
    
    #Now write individual classification and also pool classifications    
    #write_classification(gene_classification, output_directory + "/markers/tipp_%s.classification" % gene)    
    #classifications.update(gene_classification)    
    #remove_unclassified_level(classifications)
    #write_classification(classifications, output_directory+"/markers/all.classification")
    #write_abundance(classifications,output_directory)
  
    #if (options().dist == True):
    #    distribution(classification_files, output_directory)




def argument_parser():
    # Get the arguments from default SEPP
    parser = sepp.config.get_parser()

    tippGroup = parser.add_argument_group("TIPP OPTIONS", 
                                          "These arguments set settings specific to TIPP")                                 

    tippGroup.add_argument("-at", "--alignmentThreshold",
                           type=float, 
                           dest="alignment_threshold",
                           metavar="N", 
                           default = 0.95,
                           help="Enough alignment subsets are selected to reach a commulative probability of N. "
                                "This should be a number between 0 and 1 [default: 0.95].")                            

    tippGroup.add_argument("-pt", "--placementThreshold",
                           type=float, 
                           dest="placement_threshold",
                           metavar="N", 
                           default=0.95,
                           help="Enough placements are selected to reach a commulative probability of N. "
                                "This should be a number between 0 and 1 [default: 0.95].")    

    # Change this to be more flexible in two ways...
    # First denote either shotgun or amplicon data
    # if using shotgun then be allowed to look under just a specific gene... but that changes the database structure...
    # or you just need to delete other things from the dictionary before writting them...
    # allow people to do 16S from this file...
    # because otherwise you need to handle the directionality thing on your own
    # which sucks...

    #tippGroup.add_argument("-g",
    #                       "--gene",
    #                       type=str, 
    #                       dest="gene",
    #                       metavar="N", 
    #                       default=None,
    #                       help="Classify on only the specified gene.")    
                      
    tippGroup.add_argument("-G", "--genes",
                           type=str, 
                           dest="genes",
                           metavar="GENES", 
                           default="markers",
                           help="Use markers or cogs genes [default: markers].")

    tippGroup.add_argument("-b", "--blast_file",
                           type=str, 
                           dest="blast_file",
                           metavar="N", 
                           default=None,
                           help="Blast file with fragments already binned.")  
                      
    tippGroup.add_argument("-Di", "--dist",
                           dest="dist",
                           action='store_true',
                           default=False,
                           help="Treat fragments as distribution.")    
                      
    tippGroup.add_argument("-C", "--cutoff",
                           type=float, 
                           dest="cutoff",
                           metavar="N", 
                           default=0.0,
                           help="Placement probability requirement to count toward the distribution. "
                                "This should be a number between 0 and 1 [default: 0.0].")    
                      

def main():
    argument_parser()
    args = options()

    if args.alignment_size is None:
      args.alignment_size = 100

    if args.tempdir is None:
      print("A temporary directory must be specified by the user!\n")
      return

    # Could add extra options here to handle 16S data...
    # 

    if args.genes == "markers":
        genes = ['dnaG', 'frr', 'infC', 'nusA', 'pgk', 'pyrG1',
                 'pyrg', 'rplA', 'rplB', 'rplC', 'rplD', 'rplE',
                 'rplF', 'rplK', 'rplL', 'rplM', 'rplN', 'rplP',
                 'rplS', 'rplT', 'rpmA', 'rpsB', 'rpsC', 'rpsE',
                 'rpsI', 'rpsJ', 'rpsK', 'rpsM', 'rpsS', 'smpB']
    elif args.genes == "cogs":
        genes = ['COG0012', 'COG0016', 'COG0018', 'COG0048', 'COG0049',
                 'COG0052', 'COG0080', 'COG0081', 'COG0087', 'COG0088',
                 'COG0090', 'COG0091', 'COG0092', 'COG0093', 'COG0094',
                 'COG0096', 'COG0097', 'COG0098', 'COG0099', 'COG0100',
                 'COG0102', 'COG0103', 'COG0124', 'COG0172', 'COG0184',
                 'COG0185', 'COG0186', 'COG0197', 'COG0200', 'COG0201',
                 'COG0202', 'COG0215', 'COG0256', 'COG0495', 'COG0522',
                 'COG0525', 'COG0533', 'COG0541', 'COG0552', 'COG0085']
    else:
        print("Other gene options are not implemented yet...")
        return

    build_abundance_profile(args, genes)


if __name__ == '__main__':   
    main()
