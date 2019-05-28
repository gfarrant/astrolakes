#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import sys
import os
from optparse import OptionParser
from sys import stdout
import re
import kegg_toolbox

_author_ = 'Gregory K. FARRANT'
_project_ = 'AstroLakes - Marie Sklodowska-Curie Action No704956'
    
def read_blasts(input_dir):
    samples = set()
    blast_results = {}
    for file in os.listdir(input_dir):
        if not file[:6] == "BLASTN" and not "TOP95cov98ident.tab" in file:
            continue
            
        logging.info("    Reading file: %s"%file)
        sample_id = file.split('_')[1]
        samples.add(sample_id)
        
        file = input_dir+"/"+file
        for line in open(file):
            parsed = line.split()
            # geneID = "%s_%s"%(sample_id,parsed[1])
            geneID = parsed[1]
            
            if not sample_id in blast_results:
                blast_results[sample_id] = {}
            if not geneID in blast_results[sample_id]:
                blast_results[sample_id][geneID] = 0
            blast_results[sample_id][geneID] += 1
    
    ### Completion of the matrix ###
    # for sample in samples:
        # for geneID in blast_results[sample].keys():
            # if not geneID in blast_results[sample]:
                # blast_results[sample][geneID] = 0
                
    return (sorted(list(samples)),blast_results)

def read_diamonds(input_dir,identity_cutoff):
    gene_to_KEGGgeneID = {}
    
    for file in os.listdir(input_dir):
        if not file[:7] == "DIAMOND":
            continue
        logging.info("    Reading file: %s"%file)
        
        file = input_dir+"/"+file
        under_cutoff = 0
        over_cutoff = 0
        for line in open(file):
            parsed = line.split()
            # geneID = "%s_%s"%(sample_id,parsed[0])
            geneID = parsed[0]
            kegg_geneID = parsed[1]
            identity = float(parsed[2])
            
            if identity >= identity_cutoff:
                over_cutoff += 1
                if not geneID in gene_to_KEGGgeneID:
                    gene_to_KEGGgeneID[geneID] = kegg_geneID
                else:
                    logging.warn("Duplicated gene identifier: %s"%geneID)
            else:
                under_cutoff += 1
        logging.info("     -> Kept lines: %d / Skipped lines: %d"%(over_cutoff,under_cutoff))
    return gene_to_KEGGgeneID
    
def add_KEGG_orthogroup_and_taxonomy(samples,blast_results,taxM,gene_to_KEGGgeneID):

    ko_results = {}
    per_sample_dico = {}
    taxonomies = set()
    
    for sample in samples:
        success_gene = [0,0]
        success_ko = [0,0]
        for geneID in blast_results[sample].keys():
            ''' Get KEGGgeneID '''
            # geneID_short = '_'.join(geneID.split('_')[1:])
            if geneID in gene_to_KEGGgeneID:
                success_gene[0] += 1
                KEGGgeneID = gene_to_KEGGgeneID[geneID]
            else:
                success_gene[1] += 1
                KEGGgeneID = 'UNK_%s_%s'%(sample,geneID)
                
            ''' Get the corresponding KO '''
            if KEGGgeneID in taxM.gene_to_ko:
                KOs = taxM.gene_to_ko[KEGGgeneID]
                success_ko[0] += 1
            else:
                success_ko[1] += 1
                KOs = ([KEGGgeneID])
                
            ''' Get the corresponding strain taxonomy '''
            if not "UNK_" in KEGGgeneID:
                strainID = KEGGgeneID.split(":")[0]
                if not strainID in taxM.leaves:
                    logging.warn("Unknown strain identifier: %s"%strainID)
                    strain_taxonomy = "Unknown_ID_%s"%strainID
                else:
                    strain_taxonomy = ";".join(taxM.leaves[strainID].taxonomy)
            else:
                strain_taxonomy = "Unknown"
            taxonomies.add(strain_taxonomy)
            
            ''' Fill the KO dictionnary '''
            # for KO in KOs:
                # if not KO in ko_results:
                    # ko_results[KO] = {}
                # if not sample in ko_results[KO]:
                    # ko_results[KO][sample] = 0
                # ko_results[KO][sample] += blast_results[sample][geneID]
                
            for KO in KOs:
                if not sample in per_sample_dico:
                    per_sample_dico[sample] = {}
                if not KO in per_sample_dico[sample]:
                    per_sample_dico[sample][KO] = {}
                if not strain_taxonomy in per_sample_dico[sample][KO]:
                    per_sample_dico[sample][KO][strain_taxonomy] = 0
                per_sample_dico[sample][KO][strain_taxonomy] += blast_results[sample][geneID]
        
        logging.info("    [%s] GENE success %d - fail %d // KO success %d - fail %d"%(sample,success_gene[0],success_gene[1],success_ko[0],success_ko[1]))
    
    """ Matrix completion """
    for sample in samples:
        for ko in per_sample_dico[sample].keys():
            for taxonomy in taxonomies:
                if not taxonomy in per_sample_dico[sample][ko]:
                    per_sample_dico[sample][ko][taxonomy] = 0
    
    """ Sorted KOs list """
    tmp_KOs = {}
    for sample in samples:
        for ko in per_sample_dico[sample].keys():
            if not ko in tmp_KOs:
                tmp_KOs[ko] = 0
            for taxonomy in taxonomies:
                tmp_KOs[ko] += per_sample_dico[sample][ko][taxonomy]
    KOs = [elem[0] for elem in sorted(tmp_KOs.items(), key = lambda x:x[1], reverse = True)]
    
    return (ko_results,per_sample_dico,KOs)
    
def build_table(samples,blast_results,output_file):
    
    out = stdout
    if output_file is not None:
        out = open(output_file,'w')

    out.write("KO\t%s\n"%"\t".join(samples))
    for KO_md5 in blast_results.keys():
        out.write("%s\t%s\n"%(KO_md5,"\t".join([str(blast_results[KO_md5][sample]) for sample in samples])))
    if output_file is not None:
        out.close()
        
    
def convert_to_ontology_and_distribute_to_strains(samples,blast_results,conversion_table,md5_strain_table,tax_dico):
    per_sample_dico = {}
    taxonomies = set()
    warnings = set()
    for md5 in blast_results.keys():
        for sample in samples:
            if not sample in per_sample_dico:
                per_sample_dico[sample] = {}
            
            """ Get the KO """
            if md5 in conversion_table:
                ko = conversion_table[md5]
            else:
                ko = md5
            
            if not ko in per_sample_dico[sample]:
                per_sample_dico[sample][ko] = {}
                
            """ Get the strain & taxonomy """
            if md5[:4] == "UNK_":
                taxonomy = "Unknown strain"
            elif not md5 in md5_strain_table:
                logging.warn("Missing md5 - strain information (%s)"%md5)
                strain = md5
            elif not md5_strain_table[md5] in tax_dico:
                taxo_test = False
                for taxo_key in tax_dico.keys():
                    if " ".join(md5_strain_table[md5].split(" ")[:2]) in taxo_key:
                        if not md5_strain_table[md5] in warnings:
                            logging.warn("Approximative taxonomy found for strain: %s -> %s"%(md5_strain_table[md5],tax_dico[taxo_key]))
                        warnings.add(md5_strain_table[md5])
                        
                        taxonomy = tax_dico[taxo_key]
                        taxo_test = True
                        break
                if not taxo_test:
                    for taxo_key in tax_dico.keys():
                        if md5_strain_table[md5].split(" ")[0] in taxo_key:
                            if not md5_strain_table[md5] in warnings:
                                logging.warn("Very approximative taxonomy found for strain: %s -> %s"%(md5_strain_table[md5],tax_dico[taxo_key]))
                            warnings.add(md5_strain_table[md5])
                            
                            taxonomy = tax_dico[taxo_key]
                            taxo_test = True
                            break
                if not taxo_test:
                    logging.warn("No taxonomy found for strain: %s"%md5_strain_table[md5])
                    taxonomy = "Unknown taxonomy"
            else:
                taxonomy = tax_dico[md5_strain_table[md5]]
                
            """ Simplifying taxonomy by removing the species level """
            if ";" in taxonomy:
                taxonomy = ";".join(taxonomy.split(";")[:-1])
            
            if not taxonomy in per_sample_dico[sample][ko]:
                per_sample_dico[sample][ko][taxonomy] = 0
            per_sample_dico[sample][ko][taxonomy] += blast_results[md5][sample]
            taxonomies.add(taxonomy)
            
    """ Matrix completion """
    for sample in samples:
        for ko in per_sample_dico[sample].keys():
            for taxonomy in taxonomies:
                if not taxonomy in per_sample_dico[sample][ko]:
                    per_sample_dico[sample][ko][taxonomy] = 0
            
    return (per_sample_dico,sorted(list(taxonomies)))
    
def build_final_table(samples,final_results,hierarchy_table,md5_function_table,ko2descr_table,output_file):
    
    out = stdout
    if output_file is not None:
        out = open(output_file,'w')

    out.write("KO\tlvl1\tlvl2\tlvl3\tlvl4\t%s\n"%"\t".join(samples))
    for KO in final_results.keys():
        if KO in hierarchy_table:
            hierarchy = "\t".join(hierarchy_table[KO][0])
        elif KO in ko2descr_table:
            # logging.warn("Unknown KO: %s (Not found in hierarchy)"%KO)
            hierarchy = "unk\tunk\tunk\t%s"%ko2descr_table[KO]
        elif KO in md5_function_table: ## Case where KO = md5
            hierarchy = "unk\tunk\tunk\t%s"%md5_function_table[KO]
        else:
            # logging.warn("Unknown KO: %s (Not found nowhere)"%KO)
            hierarchy = "unk\tunk\tunk\tunk"
        counts = "\t".join([str(final_results[KO][sample]) for sample in samples])
        out.write("%s\t%s\t%s\n"%(KO,hierarchy,counts))
    if output_file is not None:
        out.close()
    
def build_final_table_with_taxonomy(samples,taxonomies,per_sample_dico,hierarchy_table,md5_function_table,ko2descr_table):
    
    tmp_KOs = {}
    for sample in samples:
        for ko in per_sample_dico[sample].keys():
            if not ko in tmp_KOs:
                tmp_KOs[ko] = 0
            for taxonomy in taxonomies:
                tmp_KOs[ko] += per_sample_dico[sample][ko][taxonomy]
    KOs = [elem[0] for elem in sorted(tmp_KOs.items(), key = lambda x:x[1], reverse = True)]
    
    tmp_taxonomies = {}
    for sample in samples:
        for ko in KOs: #per_sample_dico[sample].keys():
            for taxonomy in taxonomies:
                if not taxonomy in tmp_taxonomies:
                    tmp_taxonomies[taxonomy] = 0
                tmp_taxonomies[taxonomy] += per_sample_dico[sample][ko][taxonomy]
    taxonomies = [elem[0] for elem in sorted(tmp_taxonomies.items(), key = lambda x:x[1], reverse = True)]
    
    ko_dico = {}
    for sample in samples:
        with open("RESULT_%s.tab"%sample,'w') as out:
            out.write("Sample\tKO\tlvl1\tlvl2\tlvl3\tlvl4\t%s\n"%"\t".join(taxonomies))
            for KO in KOs: #per_sample_dico[sample].keys():
                if KO in hierarchy_table:
                    hierarchy = "\t".join(hierarchy_table[KO][0])
                elif KO in ko2descr_table:
                    # logging.warn("Unknown KO: %s (Not found in hierarchy)"%KO)
                    hierarchy = "unk\tunk\tunk\t%s"%ko2descr_table[KO]
                elif KO in md5_function_table: ## Case where KO = md5
                    hierarchy = "unk\tunk\tunk\t%s"%md5_function_table[KO]
                else:
                    # logging.warn("Unknown KO: %s (Not found nowhere)"%KO)
                    hierarchy = "unk\tunk\tunk\tunk"
                ko_dico[KO] = hierarchy
                counts = "\t".join([str(per_sample_dico[sample][KO][taxonomy]) for taxonomy in taxonomies])
                out.write("%s\t%s\t%s\t%s\n"%(sample,KO,hierarchy,counts))
                
    return (KOs,ko_dico)

def build_merged_taxonomy_small_dataset(per_sample_dico,samples,koM):
    logging.info("Merging the contigency table based on intelligently merged taxonomy")
        
    """Reversing the dictionnary"""
    # per_sample_dico[sample][ko][taxonomy] = Nb of reads
    for sample in samples:
        sub_results = {}
        for ko in per_sample_dico[sample].keys():
            for taxo in per_sample_dico[sample][ko].keys():
                if not taxo in sub_results:
                    sub_results[taxo] = {}
                if not ko in sub_results[taxo]:
                    sub_results[taxo][ko] = per_sample_dico[sample][ko][taxo]
                else:
                    logging.error("BUG ALERT HERE !!")
                    sys.exit(-1)
                    
        """Build taxonomic tree dictionnary"""
        taxonomic_tree_explorer = []
        for taxo in sub_results.keys():
            for y in range(1,8):
                if not ";".join(taxo.split(";")[:y]) in taxonomic_tree_explorer:
                    taxonomic_tree_explorer.append(";".join(taxo.split(";")[:y]))
        taxonomic_tree_explorer = sorted(taxonomic_tree_explorer)
            
        """Build the object structure for the leaves"""
        tree = LeavesManager()
        tree.import_leaves(sub_results)
        tree.set_threshold()
        for leaf in tree.leaves:
            leaf.test_and_set_significance(tree.threshold)
            
        """Merge the branches"""
        final_result = {}
        while len(tree.leaves) > 0:
            for taxonomic_level in taxonomic_tree_explorer:
                
                leaves_selection = tree.get_leaves_with_given_taxonomy(taxonomic_level)
                if len(leaves_selection) == 0:
                    continue
                
                
                significant_leaves = []
                for leaf in leaves_selection:
                    if leaf.is_significant:
                        significant_leaves.append(leaf)
                
                if len(significant_leaves) == 0:
                    """Sending the read count to the upper level"""
                    lvl = len(taxonomic_level.split(";"))
                    if lvl >= 1:
                        taxonomic_level = ";".join(taxonomic_level.split(";")[:lvl-1])+"_OTHERS"
                    else :
                        logging.warn("No significant branch from the root ! (Should not be possible)")
                        
                    while len(taxonomic_level.split(";")) <= 6:
                        taxonomic_level += ";*"
                    new_leaf = tree.merge_branches(taxonomic_level,leaves_selection)
                    if not new_leaf.taxonomy in final_result:
                        final_result[new_leaf.taxonomy] = {}
                        final_result[new_leaf.taxonomy] = new_leaf.abundances
                    else :
                        for k in final_result[new_leaf.taxonomy].keys():
                            final_result[new_leaf.taxonomy][k] += new_leaf.abundances[k]
                    tree.discard_leaves(leaves_selection)
                    
                elif len(significant_leaves) == 1:
                    """Merging minor branches into the major one"""
                    taxonomy = significant_leaves[0].taxonomy
                    new_leaf = tree.merge_branches(taxonomy,leaves_selection)
                    if not taxonomy in final_result:
                        final_result[new_leaf.taxonomy] = {}
                        final_result[new_leaf.taxonomy] = new_leaf.abundances
                    else :
                        logging.error("A leaf preexisted while building the final_result...")
                    tree.discard_leaves(leaves_selection)
                
                else:
                    continue
        
        reversed_final_results = {}
        taxo_set = set()
        for taxo in final_result.keys():
            taxo_set.add(taxo)
            for ko in final_result[taxo].keys():
                if not ko in reversed_final_results:
                    reversed_final_results[ko] = {}
                if not taxo in reversed_final_results[ko]:
                    reversed_final_results[ko][taxo] = final_result[taxo][ko]
        taxo_list = sorted(list(taxo_set))
        
        ko_dico = koM.kos

        tmp_KOs = {}
        for sample in samples:
            for ko in reversed_final_results.keys():
                if not ko in tmp_KOs:
                    tmp_KOs[ko] = 0
                for taxonomy in taxo_list:
                    tmp_KOs[ko] += reversed_final_results[ko][taxonomy]
        KOs = [elem[0] for elem in sorted(tmp_KOs.items(), key = lambda x:x[1], reverse = True)]
        
        with open("RESULTS_%s_with-taxonomy-merged.tab"%sample,'w') as out:
            out.write("KO\tlvl1\tlvl2\tlvl3\tlvl4\tlvl5\tlvl6\tlvl7\tlvl8\t%s\n"%"\t".join(taxo_list))
            for ko in KOs:
                counts = "\t".join([str(reversed_final_results[ko][taxo]) for taxo in taxo_list])
                if "UNK" in ko:
                    out.write("UNASSIGNED_GENE\t%s\tUnk\tUnk\tUnk\tUnk\tUnk\tUnk\tUnk\t%s\n"%(ko,counts))
                elif ":" in ko:
                    out.write("UNASSIGNED_KEGG_GENE\t%s\tUnk\tUnk\tUnk\tUnk\tUnk\tUnk\tUnk\t%s\n"%(ko,counts))
                else:
                    for category in ko_dico[ko].brite.keys():
                        hierarchies = ko_dico[ko].brite[category]
                        for hierarchy in hierarchies:
                            out.write("KO\t%s\t%s\t%s\t%s\n"%(ko,category,"\t".join(hierarchy),counts))
                
    return (reversed_final_results,taxo_list)

                # if "UNK" in ko or ":" in ko:
                    # hierarchy = "TUTU_UNK"
                # elif ko_dico[ko].brite is not None:
                    # hierarchy = list(ko_dico[ko].brite["KEGG Orthology (KO) [BR:ko00001]"])[0]
                # else :
                    # hierarchy = "TUTU"

class LeavesManager:
    
    def __init__(self):
        self.leaves = set()
        self.threshold = 0
    
    def import_leaves(self,taxo_dico):
        for (taxonomy,abundances) in taxo_dico.items():
            leaf = Leaf(taxonomy,abundances)
            leaf.set_total_abundance()
            self.leaves.add(leaf)
    
    def discard_leaves(self,leaves):
        for leaf in leaves:
            self.leaves.remove(leaf)
    
    def set_threshold(self):
        threshold = 0
        leaves_abundances = []
        for leaf in self.leaves:
            leaves_abundances.append(leaf.total_abundance)
            
        total = sum(leaves_abundances)
        sorted_abundances = sorted(leaves_abundances, reverse = True)
        
        sum_of_decreasing_abundances = 0
        for x in sorted_abundances:
            sum_of_decreasing_abundances += x
            if sum_of_decreasing_abundances > total * 0.95:
                self.threshold = x
                return
    
    def get_leaves_with_given_taxonomy(self,taxonomy):
        leaves_selection = []
        for leaf in self.leaves:
            if taxonomy.split(";") == leaf.taxonomy.split(";")[:len(taxonomy.split(";"))]:
                leaves_selection.append(leaf)
        return leaves_selection
    
    def merge_branches(self,taxonomy,leaves):
        new_abundances = {}
        for leaf in leaves:
            for sample in leaf.abundances.keys():
                if not sample in new_abundances:
                    new_abundances[sample] = 0
                new_abundances[sample] += leaf.abundances[sample]
                
        new_leaf = Leaf(taxonomy,new_abundances)
        return new_leaf
        
class Leaf:

    def __init__(self,taxonomy,abundances):
        self.taxonomy = taxonomy #text
        self.abundances = abundances #dictionnary
        self.total_abundance = 0
        self.is_significant = False
        
    def set_total_abundance(self):
        self.total_abundance = sum(self.abundances.values())
        
    def test_and_set_significance(self,threshold):
        if self.total_abundance >= threshold:
            self.is_significant = True
        
def build_EC_based_results(ko_dico,reversed_final_results,taxo_list,samples):
    logging.info("Building EC-based results with intelligently merged taxonomy")
    EC_dico = {}
    pattern = r"[0-9]+\.[0-9]+\.[0-9]+\.[0-9\-]+"
    for (ko,hierarchy) in ko_dico.items():
        hit = re.search(pattern,hierarchy)
        if hit:
            # print(type(hit),str(hit),hierarchy)
            
            EC = hit.group()
            if ko in EC_dico:
                EC_dico[ko] = EC
        else:
            EC_dico[ko] = "Unknown_EC"
                
    """ Merge results per EC """
    EC_based_results = {}
    for ko in reversed_final_results.keys():
        EC = EC_dico[ko]
        if not EC in EC_based_results:
            EC_based_results[EC] = {}
        for taxo in taxo_list:
            if not taxo in EC_based_results[EC]:
                EC_based_results[EC][taxo] = 0
            EC_based_results[EC][taxo] += reversed_final_results[ko][taxo]
            
    """ Build output file """
    sample = samples[0]
    with open("RESULTS_%s_EC_based_with-taxonomy-merged.tab"%sample,'w') as out:
            out.write("EC\t%s\n"%"\t".join(taxo_list))
            for EC in EC_based_results.keys():
                counts = "\t".join([str(EC_based_results[EC][taxo]) for taxo in taxo_list])
                out.write("%s\t%s\n"%(EC,counts))
        
if __name__ == '__main__':

    ##OPTION PARSER##
    parser = OptionParser(usage = "%prog -i input_file -o output_file (-v)")
    parser.add_option("-d","--dir",dest = "input_dir",help = "Name for directory containing input BLAST & DIAMOND files")
    parser.add_option("-i","--identity_cutoff",dest = "identity_cutoff",help = "Minimal identity to validate a functionnal assignment (default = 60%)",default=60.0)
    parser.add_option("-o","--output",dest = "output_file",help = "Name for output file")
    parser.add_option("-v","--verbose",dest="verbose",help="Enable verbose output.",action="store_true",default=False)
    (options,args) = parser.parse_args()

    if options.verbose is True :
        logging.basicConfig(level=logging.INFO)
    else :
        logging.basicConfig(level=logging.ERROR)
    
    if options.input_dir is not None:
        input_dir = options.input_dir
    else :
        logging.error("No directory provided (-d, --dir). Exiting.")
        sys.exit(-1)
    # coverage_cutoff = float(options.coverage_cutoff)
    identity_cutoff = float(options.identity_cutoff)
    
    ### READ DIAMOND & BLAST RESULTS ###
    logging.info("Starting reading BLAST files")
    (samples,blast_results) = read_blasts(input_dir)
    print(samples)
    # build_table(samples,blast_results,"raw_result.tab")
    
    logging.info("Starting reading DIAMOND files")
    gene_to_KEGGgeneID = read_diamonds(input_dir,identity_cutoff)
    
    ### ADD ONTOLOGY AND TAXONOMY ###
    logging.info("Adding KO information")
    # logging.info("Importing KEGG Taxonomy")
    taxM = kegg_toolbox.TaxonomyManager()

    (ko_results,per_sample_dico,KOs) = add_KEGG_orthogroup_and_taxonomy(samples,blast_results,taxM,gene_to_KEGGgeneID)
    # build_table(samples,ko_results,options.output_file)
    koM = kegg_toolbox.KoManager()
    (reversed_final_results,taxo_list) = build_merged_taxonomy_small_dataset(per_sample_dico,samples,koM)
