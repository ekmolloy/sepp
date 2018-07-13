'''
Created on May 31, 2015

@author: namphuon
'''
import sys,random,argparse,os,shutil,datetime
from argparse import ArgumentParser, Namespace
from sepp import get_logger
from sepp.alignment import MutableAlignment, ExtendedAlignment,_write_fasta
from sepp.exhaustive import JoinAlignJobs, ExhaustiveAlgorithm
from sepp.jobs import PastaAlignJob
from sepp.filemgr import get_temp_file, check_or_make_dir_path
import dendropy
from sepp.tree import PhylogeneticTree
from sepp.config import options,valid_decomp_strategy, valid_dir_path
from sepp.jobs import HMMBuildJob
import sepp.config
from sepp.math_utils import lcm
from sepp.problem import SeppProblem, Problem
from sepp.scheduler import JobPool,Join
from multiprocessing import Pool, Manager
from sepp.alignment import ExtendedAlignment

_LOG = get_logger(__name__)
 
class EnsembleJoinSearchJobs(Join):
    '''
    After all search jobs have finished on tips, we need return the distribution
    of the bitscores for the search.  This join accomplishes this
    '''
    def __init__(self):
        Join.__init__(self)
        
    def setup_with_root_problem(self, root_problem):
        self.root_problem = root_problem            
        for p in root_problem.iter_leaves():
            self.add_job(p.jobs["hmmsearch"])

    def perform(self):
        '''
        A dummy join that waits for all the search results to complete     
        '''
        print("")
                
                
    def __str__(self):
        return "join search jobs for all tips of ", self.root_problem

class EnsembleExhaustiveAlgorithm(ExhaustiveAlgorithm):
    '''
    This implements the exhaustive algorithm where all alignments subsets
    are searched for every fragment. This is for UPP, meaning that no placement
    is performed, and that there is always only one placement subset (currently).
    '''
    def __init__(self):
        ExhaustiveAlgorithm.__init__(self)    
        self.symfrac = False
        self.elim = None
        self.filters = True 
               
    def check_options(self):
        options().info_file = "A_dummy_value"

        if options().tree_file is None or options().alignment_file is None:
            _LOG.error("Specify the backbone alignment and tree and query sequences")
            exit(-1)
        sequences = MutableAlignment()
        sequences.read_file_object(open(self.options.alignment_file.name))  
        return ExhaustiveAlgorithm.check_options(self)
            
    def check_and_set_sizes(self, total):
        assert (self.options().placement_size is None) or (
                self.options().placement_size >= total), \
                "currently eHMMs works with only one placement subset. Please leave placement subset size option blank."        
        ExhaustiveAlgorithm.check_and_set_sizes(self, total)
        self.options().placement_size = total
        
        
    def merge_results(self):
        ''' merges search results'''
        if "fragments.distribution.done" in self.root_problem.annotations:
            return
        sequence_scores = dict([(name, []) for name in self.root_problem.fragments.keys()])        
        for fragment_chunk_problem in self.root_problem.iter_leaves():
            align_problem = fragment_chunk_problem.get_parent()
            assert isinstance(align_problem, SeppProblem)
            '''For each subproblem start with an empty set of fragments, 
            and add to them as we encounter new best hits for that subproblem'''
            if align_problem.fragments is None: 
                align_problem.fragments = self.root_problem.fragments.get_soft_sub_alignment([])
            search_res = fragment_chunk_problem.get_job_result_by_name("hmmsearch")
            for key, val in search_res.items():
                 sequence_scores[key].append([val[1],val[0]])

                    
        # TODO: is the following efficient enough? Do we need to make lists
        # and then turn them to sets?
        notScored = []
        for key,v in sequence_scores.items():
            if len(v) == 0:
                notScored.append(key)
                    
        self.root_problem.annotations["fragments.distribution.done"] = 1

        ''' Make sure all fragments are in at least one subproblem. 
        TODO: what to do with those that are not?  For now, only output warning message'''
        #notScored = [k for k, v in max_evalues.iteritems() if v[1] is None]
        _LOG.warning("Fragments %s are not scored against any subset" %str(notScored))
        #assert len(notScored) == 0, "Fragments %s are not scored against any subset" %str(notScored)  
        self.results = sequence_scores      
        
    def connect_jobs(self):
        ''' a callback function called after hmmbuild jobs are finished'''
        def enq_job_searchfragment(result, search_job):
            search_job.hmmmodel = result 
            JobPool().enqueue_job(search_job)        
        assert isinstance(self.root_problem, SeppProblem)
        for placement_problem in self.root_problem.get_children():
            '''For each alignment subproblem, ...'''
            for alg_problem in placement_problem.children:
                assert isinstance(alg_problem, SeppProblem)                
                ''' create the build model job'''
                bj = alg_problem.jobs["hmmbuild"]            
                ''' create the search jobs'''
                for fc_problem in alg_problem.get_children():
                    sj = fc_problem.jobs["hmmsearch"]
                    ''' connect bulid and search jobs'''
                    bj.add_call_Back(lambda result, next_job = sj: enq_job_searchfragment(result, next_job))
        jsj = EnsembleJoinSearchJobs()
        jsj.setup_with_root_problem(self.root_problem)         
                
    def output_results(self):        
        search_results = self.results        
        _LOG.info("Generating csv of search results. ")
        outfilename = self.get_output_filename("scores.csv")
        not_matched = self.get_output_filename("unmatched.csv")
        f = open(outfilename, 'w')
        unmatched = open(not_matched, 'w')
        f.write("seq,bitscore,evalue\n")
        
        for key,value in search_results.items():
            if len(value) == 0:
                unmatched.write("%s " % key)
            else:
                for pair in value:
                    f.write("%s,%0.4f,%s\n" % (key,pair[0],"{:.3e}".format(pair[1])))
        f.close()
        unmatched.close()
                        
    def create_fragment_files(self):
        alg_subset_count = len(list(self.root_problem.iter_leaves()))
        frag_chunk_count = lcm(alg_subset_count,self.options.cpu)//alg_subset_count
        return self.read_and_divide_fragments(frag_chunk_count)

class EnsemblePackageBuildAlgorithm(EnsembleExhaustiveAlgorithm):
    next_serial_num = 1
    def __init__(self):
        EnsembleExhaustiveAlgorithm.__init__(self)

    def check_options(self):
        options().info_file = "A_dummy_value"

        if options().package_roster is None:
            _LOG.error("To build a package, specify a path to a roster of families.")
            exit(-1)
        self.read_package_build_roster()

        if options().package_location is None:
            _LOG.error("You must specify a folder for the package to be saved.")
            exit(-1)

        self.loc_package = os.path.abspath(options().package_location)
        self.loc_hmms = os.path.join(self.loc_package, 'hmm_files')
        self.loc_fastas = os.path.join(self.loc_package, 'hmm_input_fastas')
        check_or_make_dir_path(self.loc_hmms)
        check_or_make_dir_path(self.loc_fastas)

    def get_serial_number(self):
        newser = self.next_serial_num
        self.next_serial_num += 1
        return newser

    def check_and_set_sizes(self, total):
        if (options().alignment_size is None):
            options().alignment_size = int(total*.10)
        if (options().placement_size is not None):
            _LOG.info('Note that placement size is not currently used in this implementation of HIPPI.')
        if options().placement_size is not None and options().placement_size < options().alignment_size:
            _LOG.warning("alignment_size (%d) cannot be larger than placement_size (%d).   Setting alignment_size to be placement_size" %(options().alignment_size,options().placement_size))
            options().placement_size= options().alignment_size
        if options().placement_size is not None and options().placement_size > total:
            options().placement_size = total
            _LOG.warning("Input placement size larger than total number of sequences!  Setting placement size to %d" %(total))
        if options().alignment_size > total:
            options().alignment_size = total
            _LOG.warning("Input alignment size larger than total number of sequences!  Setting alignment size to %d" %(total))

        _LOG.info("Decomposition Sizes are set to alignment: %s, placement: %s" %(options().alignment_size, options().placement_size))



    def read_package_build_roster(self):
        self.gene_families={}
        for ln in options().package_roster:
            if len(ln)>1:
                a = ln.strip().split('\t')
                if os.path.isfile(a[1]) and os.path.isfile(a[2]):
                    self.gene_families[a[0]] = {
                        'name':a[0],
                        'tree_path':a[1],
                        'alignment_path':a[2],
                        'root_problem':None
                    }
                else:
                    _LOG.error("For gene family %s, either the alignment file or the tree file could not be found."
                               "It will be skipped for this package build.")

    def read_alignment_and_tree(self,tree_path,alignment_path):
        _LOG.info("Reading input alignment: %s" % (alignment_path))
        alignment = MutableAlignment()
        alignment_file = open(alignment_path,'r')
        alignment.read_file_object(alignment_file)

        _LOG.info("Reading input tree: %s" % (tree_path))
        tree_file = open(tree_path,'r')
        tree = PhylogeneticTree(dendropy.Tree.get_from_stream(tree_file,
                                                              schema="newick",
                                                              preserve_underscores=True))
        return (alignment, tree)

    def build_subproblems(self):
        self.root_problem = Problem(parent=None)
        self.alignment_problem_directory={} #keyed by labels
        self.alignment_problem_id_lookup={} #keyed by serial number, to match to label (redundant but helpful)


        # we have to iterate over the build_subproblems task once per family
        for fam in self.gene_families.keys():
            # This part has a lot of code duplication with the ExhaustiveAlgorithm.build_subproblems, but I think
            #   it has to happen.
            _LOG.info("Starting on gene family %s" % fam)
            (alignment,tree) = self.read_alignment_and_tree(self.gene_families[fam]['tree_path'],self.gene_families[fam]['alignment_path'])
            if options().distance != 1:
                self.compute_distances(alignment)

            assert isinstance(tree, PhylogeneticTree)
            assert isinstance(alignment, MutableAlignment)
            tree.get_tree().resolve_polytomies()
            tree.lable_edges()

            ''' Make sure size values are set, and are meaningful. '''
            self.check_and_set_sizes(alignment.get_num_taxa())

            '''Have to maintain a separate *root_problem* for each family'''
            self.gene_families[fam]['root_problem']=_create_pseudo_root_problem(fam, tree, alignment, self.root_problem)

            p_tree = PhylogeneticTree(dendropy.Tree(tree.den_tree))
            p_key = 0
            alignment_tree_map = PhylogeneticTree(dendropy.Tree(p_tree.den_tree)).decompose_tree(
                self.options.alignment_size,
                strategy=self.strategy,
                minSize=self.minsubsetsize,
                tree_map={},
                decomp_strategy=self.options.decomp_strategy,
                pdistance=options().distance,
                distances=self.distances,
                maxDiam=self.options.maxDiam)
            assert len(alignment_tree_map) > 0, ("Tree could not be decomposed"
                                                 " given the following settings; strategy:%s minsubsetsize:%s alignmet_size:%s"
                                                 % (self.strategy, self.minsubsetsize, self.options.alignment_size))

            for (a_key, a_tree) in alignment_tree_map.items():
                assert isinstance(a_tree, PhylogeneticTree)
                self.modify_tree(a_tree)
                assert isinstance(self.gene_families[fam]['root_problem'], SeppProblem), 'root problem for %s is not a SeppProblem' % fam
                alignment_problem  = SeppProblem(a_tree.leaf_node_names(),
                                                 parent=self.gene_families[fam]['root_problem'])
                alignment_problem.subtree = a_tree
                alignment_problem.label = "fam_%s_A_%s" %(str(fam),str(a_key))
                serial = self.get_serial_number()
                hmm_path = os.path.join(self.loc_hmms,'hmmbuild_%s_%s.hmm' % (alignment_problem.label, serial))
                input_fasta_path = os.path.join(self.loc_fastas, 'hmmbuild_input_%s_%s.fasta' % (alignment_problem.label, serial))
                n_seqs = alignment_problem.subalignment.get_num_taxa()
                self.alignment_problem_directory[alignment_problem.label]={'family':fam, 'label':alignment_problem.label,
                                                                           'serial_number':serial, 'hmm_path':hmm_path,
                                                                           'fasta_path':input_fasta_path, 'a_key':str(a_key),
                                                                           'num_seqs':n_seqs}
                self.alignment_problem_id_lookup[serial] = alignment_problem.label
                alignment_problem.write_subalignment_without_allgap_columns(input_fasta_path)


            _LOG.info("Breaking family %s into %d alignment subsets." % (fam,len(alignment_tree_map)))

        _LOG.debug("Subproblem structure: %s" % str(self.root_problem))
        return self.root_problem

    def document_package_build(self):
        package_index_path = os.path.join(self.loc_package,'alignment_problem_directory.tsv')
        packdir_f = open(package_index_path,'w')
        line_arg = '%(family)s\t%(label)s\t%(serial_number)s\t%(hmm_path)s\t%(fasta_path)s\t%(a_key)s\t%(num_seqs)s\n'
        # line_arg = '%(family)s\t%(label)s\t%(serial_number)s\t%(hmm_path)s\t%(fasta_path)s\t%(a_key)s\n'
        header = line_arg % {'family':'family', 'label':'label', 'serial_number':'serial_number',
                             # 'hmm_path': 'hmm_path', 'fasta_path': 'fasta_path', 'a_key': 'a_key'}
                             'hmm_path':'hmm_path', 'fasta_path':'fasta_path', 'a_key':'a_key', 'num_seqs':'num_seqs'}
        packdir_f.write(header)
        for k in self.alignment_problem_directory.keys():
            packdir_f.write(line_arg % self.alignment_problem_directory[k])
        packdir_f.close()

        package_metadata_path = os.path.join(self.loc_package,'package_metadata.txt')
        self.meta_f = open(package_metadata_path,'w')
        str_header = '''
        #*****************************************************
        #
        #  HIPPI: Build documentation for ensemble of HMMs on 
        #     multiple gene families.
        #
        #*****************************************************
        
        '''
        self.meta_f.write(str_header)
        self.meta_f.write('Date of Build: %s\n' % datetime.datetime.now())
        self.meta_f.write('Number of Gene Families: %s\n' % len(self.gene_families))
        self.meta_f.write('Total number of HMMs: %s\n' % len(self.alignment_problem_directory))
        self.meta_f.write('Path to final package: %s\n' % self.loc_package)
        self.meta_f.write('Path to the file containing the list of HMMs and properties: %s' % package_index_path)
        self.meta_f.write('\n')
        self.meta_f.write('Options and Settings at Build (dumped from argument parser):\n')
        for (k,v) in vars(options()).items():
            self.meta_f.write('\t%s: %s\n' % (k,v))

        self.meta_f.write('\n\n')
        self.meta_f.close()

        self.subprob_f = open(os.path.join(self.loc_package,'subproblem_structure.txt'))
        self.subprob_f.write('Subproblem Structure:\n')
        self.subprob_f.write(str(self.root_problem))
        self.subprob_f.close()

    def build_jobs(self):

        for fam_root in self.root_problem.get_children():
            fam_name = fam_root.label

            for alg_problem in fam_root.get_children():
                assert isinstance(alg_problem, SeppProblem)
                ''' create the build model job'''
                rec = self.alignment_problem_directory[alg_problem.label]
                bj = HMMBuildJob()
                bj.setup(infile=rec['fasta_path'],
                         outfile=rec['hmm_path'],
                         symfrac=self.symfrac,
                         informat='fasta',
                         molecule=self.molecule,
                         **vars(self.options.hmmbuild))
                alg_problem.add_job(bj.job_type, bj)


    def connect_jobs(self):
        '''
        We actually don't have muhc to connect here because we're just getting everything built.
        '''
        pass

    def merge_results(self):
        pass

    def output_results(self):
        self.document_package_build()

    def enqueue_firstlevel_job(self):
        for p in self.root_problem.children:
            for ap in p.children:
                JobPool().enqueue_job(ap.jobs["hmmbuild"])

def _create_pseudo_root_problem(family_name, tree, alignment, parent):
    ''' Create the root problem for this family'''
    root_problem = SeppProblem(tree.leaf_node_names(), parent=parent)
    root_problem.label = family_name
    root_problem.subalignment = alignment
    root_problem.subtree = tree
    return root_problem

def augment_parser():
    sepp.config.set_main_config_path(os.path.expanduser("~/.sepp/upp.config"))
    parser = sepp.config.get_parser()    
    parser.description = "This script runs the UPP algorithm on set of sequences.  A backbone alignment and tree can be given as input.  If none is provided, a backbone will be automatically generated."
    
    decompGroup = parser.groups['decompGroup']                                 
    decompGroup.__dict__['description'] = ' '.join(["These options",
        "determine the alignment decomposition size, backbone size, and how to decompose the backbone set."])
        
    
    decompGroup.add_argument("-A", "--alignmentSize", type = int, 
                      dest = "alignment_size", metavar = "N", 
                      default = 10,
                      help = "max alignment subset size of N "
                             "[default: 10]")    
    decompGroup.add_argument("-S", "--decomp_strategy", type = valid_decomp_strategy, 
                      dest = "decomp_strategy", metavar = "DECOMP",
                      default = "hierarchical", 
                      help = "decomposition strategy "
                             "[default: ensemble of HMMs (hierarchical)]")                              
                             
    inputGroup = parser.groups['inputGroup']                             
    inputGroup .add_argument("-s", "--sequence_file", type = argparse.FileType('r'),
                      dest = "sequence_file", metavar = "SEQ", 
                      default = None,
                      help = "Unaligned sequence file.  "
                             "If no backbone tree and alignment is given, the sequence file will be randomly split into a backbone set (size set to B) and query set (remaining sequences), [default: None]")                             
    inputGroup.add_argument("-c", "--config", 
                      dest = "config_file", metavar = "CONFIG",
                      type = argparse.FileType('r'), 
                      help = "A config file, including options used to run UPP. Options provided as command line arguments overwrite config file values for those options. "
                             "[default: %(default)s]")    
    inputGroup.add_argument("-t", "--tree", 
                      dest = "tree_file", metavar = "TREE",
                      type = argparse.FileType('r'), 
                      help = "Input tree file (newick format) "
                             "[default: %(default)s]")    
    inputGroup.add_argument("-a", "--alignment", 
                      dest = "alignment_file", metavar = "ALIGN",
                      type = argparse.FileType('r'), 
                      help = "Aligned fasta file "
                             "[default: %(default)s]")
    #Nute adding
    inputGroup.add_argument("-pa", "--package",
                            dest="package_location", metavar="PACKAGE", default=None,
                            type = valid_dir_path,
                            help = "Location of the HIPPI package to be used.")
    inputGroup.add_argument("-ro" ,"--package_roster", type=argparse.FileType('r'),
                            dest = "package_roster", metavar="ROSTER", default=None,
                            help = "this option is used if we wnat to have HIPPI precompute a package but not "
                                   "analyze it. This argument specifies a roster of families to become part of"
                                   "the package. The roster should be a tab-separated list (no header) of:"
                                   "<family_name>   <tree_file_path>    <alignment_file_path>"
                                   ""
                                   "This argument must be omitted unless you want to build a package of HMMs and"
                                   "save it to disk.")
                             
    uppGroup = parser.add_argument_group("UPP Options".upper(), 
                         "These options set settings specific to UPP")                                 
    
    seppGroup = parser.add_argument_group("SEPP Options".upper(), 
                         "These options set settings specific to SEPP and are not used for UPP.")                                 
    seppGroup.add_argument("-P", "--placementSize", type = int, 
                      dest = "placement_size", metavar = "N",
                      default = None, 
                      help = "max placement subset size of N "
                             "[default: 10%% of the total number of taxa]")                              
    seppGroup.add_argument("-r", "--raxml", 
                      dest = "info_file", metavar = "RAXML",
                      type = argparse.FileType('r'), 
                      help = "RAxML_info file including model parameters, generated by RAxML."
                             "[default: %(default)s]")    
    seppGroup.add_argument("-f", "--fragment",
                      dest = "fragment_file", metavar = "FRAG",
                      type = argparse.FileType('r'), 
                      help = "fragment file "
                             "[default: %(default)s]")          
                             
                                                   
def main():
    augment_parser()
    if options().package_roster is None:
        _LOG.info('No roster detected so proceeding with EnsembleExhaustiveAlgorithm.')
        EnsembleExhaustiveAlgorithm().run()
    else:
        _LOG.info('Gene family roster detected, so proceeding to construct a HIPPI package.')
        EnsemblePackageBuildAlgorithm().run()

if __name__ == '__main__':   
    main()
        
