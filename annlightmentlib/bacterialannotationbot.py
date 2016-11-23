from collections import defaultdict
import csv
import logging
import pprint
import pywikibot
from pywikibot.data import api
import sys
import time
from annlightmentlib.gff3parser import Gff3Parser

class BacterialAnnotationBot():

    def __init__(self, site, annotation_file, interaction_file, strain_id,
                 databank):
        self.annotation_file = annotation_file
        self.interaction_file = interaction_file
        self.strain_id = strain_id
        self.site = site
        self.repo = self.site.data_repository()
        self.databank = databank
        self.locus_tag_dict = defaultdict(dict)
        self.id_locus_tag_dict = defaultdict(dict)
        self.number_of_uploaded_items = defaultdict(int)


    def _get_strain_name(self, strain_id):
        strain_item = pywikibot.ItemPage(self.repo, strain_id) 
        strain_item_dict = strain_item.get()
        strain_label_en = strain_item_dict['labels']['en']
        return(strain_label_en)

    def _make_log_entry_no_transcript(self, ncRNA_name, targets_locus_tag):
        timestr = time.strftime("%Y%m%d-%H_%M_%S")
        logging.basicConfig(filename='ANNlightment_' + timestr +
                            '.log', level=logging.DEBUG)
        logging.info("no matching transcript for ncRNA {} and targets locus "
                     "tag {}".format(ncRNA_name, targets_locus_tag))

    
    def _make_log_entry_interaction_exists(
            self, method, ncRNA_ID, transcript_ID):
        timestr = time.strftime("%Y%m%d-%H_%M_%S")
        logging.basicConfig(filename='ANNlightment_' + timestr +
                            '.log', level=logging.DEBUG)
        logging.info("interaction({}) between ncRNA {} and "
                     "transcript {} already exists".format(
                         method, ncRNA_ID, transcript_ID))
                     
    def _make_log_entry_create_interaction(
            self, method, ncRNA_ID, transcript_ID):
        timestr = time.strftime("%Y%m%d-%H_%M_%S")
        logging.basicConfig(filename='ANNlightment_' + timestr +
                            '.log', level=logging.DEBUG)
        logging.info("created interaction({}) between ncRNA {} and "
                     "transcript {}".format(method, ncRNA_ID, transcript_ID))
    
    def _make_log_entry(self, item_type, new_item_id):
        timestr = time.strftime("%Y%m%d-%H_%M_%S")
        logging.basicConfig(filename='ANNlightment_' + timestr +
                            '.log', level=logging.DEBUG)
        logging.info("created {} item with ID {}".format(item_type,
                                                         new_item_id))
    def _search_item_by_label(self, site, item_name):
        all_results = []
        last_continue = 0
        params = { 'action' :'wbsearchentities' , 'format' : 'json' ,
                   'language' : 'en', 'type' : 'item',
                   'search': item_name}
        while isinstance(last_continue, int):
            params["continue"] = last_continue
            request = api.Request(site=site,**params)
            result = request.submit()
            for item_entry in result['search']:
                all_results.append(item_entry)
            last_continue = result['search'
            '-continue'] if "search-continue" in result else None
        return(all_results)
   
    def _item_already_exists(self, name_entries, item_name, item_type,
                             strain_name):
        return_list = []        
        if name_entries != []:
            for list_item in name_entries:
                if (list_item['label'] == item_name and
                    list_item['description'] == "{} found in "
                    "{}".format(item_type, strain_name)):
                    return_list.append("yes")    
                elif (list_item['label'] == item_name
                      and list_item['description'] !="{} found in "
                      "{}".format(item_type, strain_name)):
                    return_list.append("yes, but different description")
                else:
                    return_list.append("no")
        elif name_entries == []:
            return_list.append("no")
        if "yes" in return_list:
            return("yes")
        elif ("yes" not in return_list and
              "no" not in return_list and
              "yes, but different description" in return_list):
            return("yes, but different description")
        elif ("yes" not in return_list and
              "no" in return_list):
            return("no")
        
    def all_features(self):        
        self.create_genes_products_and_claims()
        self.create_relating_claims()
        self.create_transcripts_and_claims()
        self.create_transcripts_for_parentless_genes()
        self.create_sRNA_interactions()
        self.print_overview()

    def print_overview(self):
        sum_of_items = (sum(self.number_of_uploaded_items.values()) -
                        self.number_of_uploaded_items['interactions'])     
        print("ANNlightment uploaded a total of {} items".format(sum_of_items))
        print(" {} genes, \n {} rRNAs, \n {} tRNAs, \n {} proteins, \n {} TSS,"
              "\n {} ncRNAs, \n {} transcripts derived from ANNOgesic \n "
              "{} transcripts derived from RefSeq and \n {} "
              "interactions".format(
                  self.number_of_uploaded_items['genes'],
                  self.number_of_uploaded_items['rRNAs'],
                  self.number_of_uploaded_items['tRNAs'],
                  self.number_of_uploaded_items['proteins'],
                  self.number_of_uploaded_items['TSS'],
                  self.number_of_uploaded_items['ncRNAs'],
                  self.number_of_uploaded_items['transcripts_ANNOgesic'],
                  self.number_of_uploaded_items['transcripts_RefSeq'],
                  self.number_of_uploaded_items['interactions']))
              
    def create_genes_products_and_claims(self):
        with open(self.annotation_file) as csv_file: 
            gff_parser = Gff3Parser()
            gff_iterator = gff_parser.entries(csv_file)
            for row in gff_iterator:
                if row.feature == "gene":
                    self._process_entry("gene", row)
                elif row.feature == "CDS":
                    self._process_entry("protein", row)
                elif row.feature == "rRNA":
                    self._process_entry("rRNA", row)
                elif row.feature == "tRNA":
                    self._process_entry("tRNA", row)
                elif row.feature == "TSS":
                    self._process_entry("transcription start site", row)
                elif row.feature == "ncRNA":
                    self._process_ncRNA_entry(row)
                        
    def _process_entry(self, entry_type, row):
        property_dict = self._get_property_dict()
        strain_name = self._get_strain_name(self.strain_id)
        genomic_start = row.start
        genomic_end = row.end
        if row.strand == "+":
            strand_direction = property_dict["Forward_Strand"]
        elif row.strand == "-":
            strand_direction = property_dict["Reverse_Strand"]        
        item_name = strain_name + " " + row.attributes['ID']
        entry_id = row.attributes['ID']
        if "Parent" in row.attributes.keys():
            parent = row.attributes["Parent"]
        else:
            parent = "NA"            
        item_type = "bacterial {}".format(entry_type)
        item_description = ("bacterial {} found in ".format(entry_type)
                            + strain_name)
        logging_expression = entry_type
        if entry_type in ["protein","rRNA","tRNA"]:
            item_alias = row.attributes['product']
            locus_tag = row.attributes['locus_tag']
        elif entry_type == "gene":
            item_alias = row.attributes['Name']
            locus_tag = row.attributes['locus_tag']
        elif entry_type == "transcription start site":
            item_alias = row.attributes['Name']
            locus_tag = "NA"
            logging_expression = "TSS"
        self._create_gene_or_product_item_if_non_existent(
            item_name, item_type, strain_name, item_description,
            logging_expression, strand_direction, genomic_start, genomic_end,
            locus_tag, parent, item_alias, entry_id)
        
        
    def _process_ncRNA_entry(self, row):
        property_dict = self._get_property_dict()
        strain_name = self._get_strain_name(self.strain_id)
        genomic_start = row.start
        genomic_end = row.end
        if row.strand == "+":
            strand_direction = property_dict["Forward_Strand"]
        elif row.strand == "-":
            strand_direction = property_dict["Reverse_Strand"]        
        item_name = strain_name + " " + row.attributes['ID']
        item_alias = row.attributes['Name']
        entry_id = row.attributes['ID']
        if "Parent" in row.attributes.keys():
            parent = row.attributes["Parent"]
        else:
            parent = "NA"                    
        locus_tag = "NA"
        item_type = "bacterial sRNA"
        item_description = ("bacterial sRNA found in " +
                            strain_name)
        logging_expression = "ncRNA"
        sRNA_ID = self._create_gene_or_product_item_if_non_existent(
                        item_name, item_type, strain_name, item_description,
                        logging_expression, strand_direction, genomic_start,
                        genomic_end, locus_tag, parent, item_alias, entry_id)
        g_item_type = "bacterial gene"
        g_item_description = ("bacterial gene found in " + strain_name)
        g_logging_expression = "gene"
        gene_ID = self._create_gene_or_product_item_if_non_existent(
                        item_name,
                        g_item_type, strain_name, g_item_description,
                        g_logging_expression, strand_direction, genomic_start,
                        genomic_end, locus_tag, parent, item_alias, entry_id)
        if sRNA_ID != None:
            self._add_claim_item(gene_ID, property_dict["encodes"],
                                 sRNA_ID)
            self._add_claim_item(sRNA_ID,
                                 property_dict["encoded by"],
                                 gene_ID)
            print("connected gene {} with product {}".format(
                gene_ID, sRNA_ID))
            
    def _create_gene_or_product_item_if_non_existent(
            self, item_name, item_type, strain_name, item_description,
            logging_expression, strand_direction, genomic_start, genomic_end,
            locus_tag, parent, item_alias, entry_id):
        name_entries = self._search_item_by_label(self.site, item_name)
        if "," in parent:
            parents = parent.split(',')
        elif parent == "NA":
            parents = []
            parents.append("NA")
        else:
            parents = []
            parents.append(parent)
        if self._item_already_exists(name_entries, item_name,
                                     item_type, strain_name) == "yes":
            print("item {} already exists".format(item_name))
        elif (self._item_already_exists(name_entries, item_name,
                                        item_type, strain_name) == "no" or
              self._item_already_exists(name_entries, item_name,
                                        item_type, strain_name) == ("yes, " 
                                        "but different description")):
            new_item_id = self._create_new_item(item_name, item_description)
            self._make_log_entry(logging_expression, new_item_id)
            if logging_expression in ["gene", "protein", "tRNA", "rRNA", "TSS"]:
                self._write_to_id_locus_tag_dict(
                    new_item_id, logging_expression, locus_tag, parents,
                    genomic_start, genomic_end, strand_direction, item_name,
                    entry_id)
            elif logging_expression == "ncRNA":
                self._write_sRNA_to_id_locus_tag_dict(
                    new_item_id, logging_expression, locus_tag, parents,
                    genomic_start, genomic_end, strand_direction, item_name,
                    item_alias, entry_id)    
            print("created {} item with ID {}".format(logging_expression,
                                                      new_item_id))
            if logging_expression == "gene":
                self._add_claims_to_new_gene_item(
                    new_item_id, strand_direction, genomic_start, genomic_end,
                    locus_tag, item_alias)
                self.number_of_uploaded_items['genes'] += 1
            elif logging_expression == "rRNA":
                self._add_claims_to_new_rRNA_item(new_item_id, item_alias)
                self.number_of_uploaded_items['rRNAs'] += 1
            elif logging_expression == "tRNA":
                self._add_claims_to_new_tRNA_item(new_item_id, item_alias)
                self.number_of_uploaded_items['tRNAs'] += 1
            elif logging_expression == "protein":
                self._add_claims_to_new_protein_item(
                    new_item_id, locus_tag, item_alias)
                self.number_of_uploaded_items['proteins'] += 1
            elif logging_expression == "TSS":
                self._add_claims_to_new_TSS_item(
                    new_item_id, strand_direction, genomic_start, genomic_end,
                    item_alias)
                self.number_of_uploaded_items['TSS'] += 1
            elif logging_expression == "ncRNA":
                self._add_claims_to_new_sRNA_item(new_item_id, item_alias)
                self.number_of_uploaded_items['ncRNAs'] += 1
            return(new_item_id)
                
    def _write_to_id_locus_tag_dict(
            self, item_id, feature_type, locus_tag, parents, genomic_start,
            genomic_end, strand_direction, item_name, entry_id):
        self.id_locus_tag_dict[item_id]['feature_type'] = feature_type
        self.id_locus_tag_dict[item_id]['locus_tag'] = locus_tag                
        self.id_locus_tag_dict[item_id]['parent'] = parents
        self.id_locus_tag_dict[item_id]['genomic_start'] = genomic_start
        self.id_locus_tag_dict[item_id]['genomic_end'] = genomic_end
        self.id_locus_tag_dict[item_id]['strand_direction'] = strand_direction
        self.id_locus_tag_dict[item_id]['item_name'] = item_name
        self.id_locus_tag_dict[item_id]['entry_id'] = entry_id

    def _write_sRNA_to_id_locus_tag_dict(
            self, item_id, feature_type, locus_tag, parents, genomic_start,
            genomic_end, strand_direction, item_name, item_alias, entry_id):
        self.id_locus_tag_dict[item_id]['feature_type'] = feature_type
        self.id_locus_tag_dict[item_id]['locus_tag'] = locus_tag                
        self.id_locus_tag_dict[item_id]['parent'] = parents
        self.id_locus_tag_dict[item_id]['genomic_start'] = genomic_start
        self.id_locus_tag_dict[item_id]['genomic_end'] = genomic_end
        self.id_locus_tag_dict[item_id]['strand_direction'] = strand_direction
        self.id_locus_tag_dict[item_id]['item_name'] = item_name
        self.id_locus_tag_dict[item_id]['item_alias'] = item_alias
        self.id_locus_tag_dict[item_id]['entry_id'] = entry_id
        
    def _add_claims_to_new_tRNA_item(self, new_item_id, item_alias):
        property_dict = self._get_property_dict()
        claim_item_dict = {property_dict["found_in_taxon"] : self.strain_id,
                           property_dict["instance_of"] :
                           property_dict["transfer RNA"]}
        for claim, target in claim_item_dict.items():
            self._add_claim_item(new_item_id, claim, target)
        self._add_alias_to_new_item(new_item_id, item_alias)
                
    def _add_claims_to_new_rRNA_item(self, new_item_id, item_alias):
        property_dict = self._get_property_dict()
        claim_item_dict = {property_dict["found_in_taxon"] : self.strain_id,
                           property_dict["instance_of"] :
                           property_dict["ribosomal RNA"]}
        for claim, target in claim_item_dict.items():
            self._add_claim_item(new_item_id, claim, target)    
        self._add_alias_to_new_item(new_item_id, item_alias)

    def _add_claims_to_new_gene_item(
            self, new_item_id, strand_direction, genomic_start, genomic_end,
            locus_tag, item_alias):
        property_dict = self._get_property_dict()
        claim_item_dict = {property_dict["found_in_taxon"] : self.strain_id,
                           property_dict["instance_of"] : property_dict["DNA"],
                           property_dict["strand_orientation"] :
                           strand_direction}
        claim_string_dict = {property_dict["genomic_start"] :
                             str(genomic_start),
                             property_dict["genomic_end"] :
                             str(genomic_end), property_dict["NCBI Locus tag"] :
                             locus_tag}
        additional_item_dict = {property_dict["instance_of"] :
                                property_dict["gene"]}                    
        for claim, target in claim_item_dict.items():
            self._add_claim_item(new_item_id, claim, target)
        for claim, target in additional_item_dict.items():
            self._add_claim_item(new_item_id, claim, target)                    
        for claim, target in claim_string_dict.items():
            self._add_claim_string(new_item_id, claim, target) 
        self._add_alias_to_new_item(new_item_id, item_alias)

    def _add_claims_to_new_protein_item(self, new_item_id, locus_tag,
                                        item_alias):
        property_dict = self._get_property_dict()
        claim_item_dict = {property_dict["found_in_taxon"] : self.strain_id,
                           property_dict["instance_of"] :
                           property_dict["protein"]}
        claim_string_dict = {property_dict["NCBI Locus tag"] : locus_tag}      
        for claim, target in claim_item_dict.items():
            self._add_claim_item(new_item_id, claim, target)              
        for claim, target in claim_string_dict.items():
            self._add_claim_string(new_item_id, claim, target)
        self._add_alias_to_new_item(new_item_id, item_alias)

    def _add_claims_to_new_TSS_item(self, new_item_id, strand_direction,
                                    genomic_start, genomic_end, item_alias):
        property_dict = self._get_property_dict()       
        claim_item_dict = {property_dict["found_in_taxon"] : self.strain_id,
                           property_dict["instance_of"] : property_dict["DNA"],
                           property_dict["strand_orientation"] :
                           strand_direction}
        claim_string_dict = {property_dict["genomic_start"] :
                             str(genomic_start),
                             property_dict["genomic_end"] :
                             str(genomic_end)}
        claim_additional_item_dict = {property_dict["instance_of"] :
                                      property_dict["TSS"]}                    
        for claim, target in claim_item_dict.items():
            self._add_claim_item(new_item_id, claim, target)
        for claim, target in claim_additional_item_dict.items():
            self._add_claim_item(new_item_id, claim, target)                    
        for claim, target in claim_string_dict.items():
            self._add_claim_string(new_item_id, claim, target)
        self._add_alias_to_new_item(new_item_id, item_alias)

    def _add_claims_to_new_sRNA_item(self, new_item_id, item_alias):
        property_dict = self._get_property_dict()
        claim_item_dict = {property_dict["found_in_taxon"] : self.strain_id,
                           property_dict["instance_of"] : property_dict["RNA"]}
        additional_item_dict = {property_dict["instance_of"] :
                                property_dict["small_RNA"]}                    
        for claim, target in claim_item_dict.items():
            self._add_claim_item(new_item_id, claim, target)
        for claim, target in additional_item_dict.items():
            self._add_claim_item(new_item_id, claim, target)
        self._add_alias_to_new_item(new_item_id, item_alias)
            
    def _add_alias_to_new_item(self, new_item_id, item_alias):
        alias = item_alias
        new_item = pywikibot.ItemPage(self.repo, new_item_id)
        newalias = {"en": [alias]}
        new_item.editAliases(newalias, summary="Setting new aliases.")
                
    def create_relating_claims(self):
        matching_IDs = self._get_matching_IDs_dict()
        property_dict = self._get_property_dict()
        for gene_ID, product_IDs in matching_IDs.items():
            for product_ID in product_IDs:
                self._add_claim_item(gene_ID, property_dict["encodes"],
                                     product_ID)
                self._add_claim_item(product_ID, property_dict["encoded by"],
                                 gene_ID)
                print("connected gene {} with product {}".format(gene_ID,
                                                                 product_ID)) 
            
    def _get_matching_IDs_dict(self):
        id_locus_tag_dict = self.id_locus_tag_dict
        matching_IDs = {}
        for item in id_locus_tag_dict.items():
            if (item[1]["feature_type"] == "gene" and
                item[1]["locus_tag"] != "NA"):
                cds_list = []
                for item2 in id_locus_tag_dict.items():
                    if (item2[1]["locus_tag"] != "NA" and
                        item2[1]["feature_type"] in ["protein", "rRNA", "tRNA"]): 
                        if (item2[1]["locus_tag"] == item[1]["locus_tag"] and
                            item2[0] != item[0]):
                            cds_list.append(item2[0])
                matching_IDs[item[0]] = cds_list
        return(matching_IDs)
    def create_transcripts_and_claims(self):
        strain_name = self._get_strain_name(self.strain_id)
        property_dict = self._get_property_dict()
        with open(self.annotation_file) as gff: 
            gff_parser = Gff3Parser()
            gff_iterator = gff_parser.entries(gff)
            for row in gff_iterator:
                genomic_start = row.start
                genomic_end = row.end
                if row.strand == "+":
                    strand_direction = property_dict["Forward_Strand"]
                elif row.strand == "-":
                    strand_direction = property_dict["Reverse_Strand"]
                if row.feature == "transcript":
                    locus_tag = "NA"
                    item_type = "bacterial transcript" 
                    item_description = ("bacterial transcript found in " +
                                        strain_name)
                    item_name = (strain_name + " " + row.attributes['ID'] +
                                 " " + str(row.start) + " " + str(row.end))
                    item_alias = row.attributes['Name']
                    logging_expression = "transcript"
                    ID = row.attributes['ID']
                    new_item_id = self._create_transcript_item_if_non_existent(
                        item_name, item_type, strain_name, item_description,
                        logging_expression, strand_direction, genomic_start,
                        genomic_end, ID)
                    determination_method = "ANNOgesic"
                    if new_item_id != None:
                        self._add_claims_to_new_transcript_item(new_item_id,
                            strand_direction, genomic_start, genomic_end,
                            determination_method, item_alias) 
                    
    def _create_transcript_item_if_non_existent(
            self, item_name, item_type, strain_name, item_description,
            logging_expression, strand_direction, genomic_start, genomic_end,
            ID):
        property_dict = self._get_property_dict()
        name_entries = self._search_item_by_label(self.site, item_name)
        new_item_id = None
        if self._item_already_exists(name_entries, item_name,
                                     item_type, strain_name) == "yes":
            print("item {} already exists".format(item_name))
        elif (self._item_already_exists(name_entries, item_name,
                                        item_type, strain_name) == "no" or
              self._item_already_exists(name_entries, item_name,
                                        item_type, strain_name) == ("yes, " 
                                        "but different description")):
            new_item_id = self._create_new_item(item_name, item_description)
            self._make_log_entry(logging_expression, new_item_id)
            print("created {} item with ID {}".format(logging_expression,
                                                      new_item_id))
            for child_ID in self._get_transcript_children(ID):
                self._add_claim_item(new_item_id, property_dict["has part"],
                                     child_ID)
                self._add_claim_item(child_ID, property_dict["part of"],
                                     new_item_id)
                print("created has part/part of connection "
                      "between {} and {}".format(new_item_id, child_ID))
        return(new_item_id)
                
    def _add_claims_to_new_transcript_item(self, new_item_id, strand_direction,
                                           genomic_start, genomic_end,
                                           determination_method, item_alias):
        property_dict = self._get_property_dict()
        if determination_method == "ANNOgesic":
            determination_method = property_dict["ANNOgesic"]
            self.number_of_uploaded_items['transcripts_ANNOgesic'] += 1
        elif determination_method == "RefSeq":
            determination_method = property_dict["RefSeq"]
            self.number_of_uploaded_items['transcripts_RefSeq'] += 1
        claim_item_dict = {property_dict["found_in_taxon"] : self.strain_id,
                           property_dict["instance_of"] : property_dict["RNA"],
                           property_dict["strand_orientation"] :
                           strand_direction,
                           property_dict["determination method"] :
                           determination_method}
        claim_string_dict = {property_dict["genomic_start"] :
                             str(genomic_start),
                             property_dict["genomic_end"] : str(genomic_end)}
        additional_claim_item_dict = {property_dict["instance_of"] :
                                      property_dict["transcript"]}
        for claim, target in claim_item_dict.items():
            self._add_claim_item(new_item_id, claim, target)            
        for claim, target in claim_string_dict.items():
            self._add_claim_string(new_item_id, claim, target)
        for claim, target in additional_claim_item_dict.items():
            self._add_claim_item(new_item_id, claim, target)
        self._add_alias_to_new_item(new_item_id, item_alias)
                
    def _get_transcript_children(self, ID):
        id_locus_tag_dict = self.id_locus_tag_dict
        IDs_of_children = []
        for item in id_locus_tag_dict.items():
            if (item[1]["feature_type"] == "gene" or
                item[1]["feature_type"] == "TSS"):
                if "parent" in item[1]:                   
                    for prnt in item[1]["parent"]:
                        if prnt == ID:
                            IDs_of_children.append(item[0])
        print("these are the children of the transcript {}".format(
            IDs_of_children))
        return(IDs_of_children)
    
    def create_transcripts_for_parentless_genes(self):
        property_dict = self._get_property_dict()
        strain_name = self._get_strain_name(self.strain_id)
        id_locus_tag_dict = self.id_locus_tag_dict
        item_description = "bacterial transcript found in " + strain_name
        determination_method = "RefSeq"
        item_type = "bacterial transcript"
        for item in id_locus_tag_dict.items():
            if item[1]["parent"][0] == "NA":
                item_name = (strain_name + " " "transcript" + " " +
                             str(item[1]["genomic_start"]) +
                             " " + str(item[1]["genomic_end"]))
                item_alias = "transcript by RefSeq"
                name_entries = self._search_item_by_label(self.site, item_name)
                if self._item_already_exists(name_entries, item_name,
                                             item_type, strain_name) == "yes":
                    print("item {} already exists".format(item_name))           
                elif (self._item_already_exists(name_entries, item_name,
                                                item_type, strain_name) == "no"
                      or self._item_already_exists(name_entries, item_name,
                                                   item_type, strain_name) ==
                      "yes, but different description"):
                    new_transcript_id = self._create_new_item(item_name,
                                                              item_description)
                    print("created transcript item with ID {}".format(
                        new_transcript_id))                    
                    self._make_log_entry("transcript", new_transcript_id)
                    self._add_claims_to_new_transcript_item(new_transcript_id,
                        item[1]["strand_direction"], item[1]["genomic_start"],
                        item[1]["genomic_end"], determination_method,
                        item_alias)
                    self._add_claim_item(new_transcript_id,
                                         property_dict["has part"],
                                         item[0])
                    self._add_claim_item(item[0], property_dict["part of"],
                                         new_transcript_id)
                    print("created has part/part of connection "
                          "between {} and {}".format(new_transcript_id,
                                                     item[0]))

    def create_sRNA_interactions(self):
        property_dict = self._get_property_dict()
        id_locus_tag_dict = self.id_locus_tag_dict
        strain_name = self._get_strain_name(self.strain_id)
        with open(self.interaction_file) as csvfile:
            sRNA_interactions_dict = csv.DictReader(csvfile, delimiter="\t")
            for row in sRNA_interactions_dict:
                specs = self._process_interactions_row(row)
                if (self._get_matching_transcripts(specs["sRNA_name"],
                                                   specs["targets_locus_tag"],
                                                   specs["sRNA_start_pos"],
                                                   specs["sRNA_end_pos"],
                                                   specs["entry_id"]) == (
                                                       None)):
                    print("no matching transcripts")
                    self._make_log_entry_no_transcript(
                        specs["sRNA_name"], specs["targets_locus_tag"])
                    continue
                sRNA_item_id, matching_transcript_IDs = (
                    self._get_matching_transcripts(specs["sRNA_name"],
                                                   specs["targets_locus_tag"],
                                                   specs["sRNA_start_pos"],
                                                   specs["sRNA_end_pos"],
                                                   specs["entry_id"]))
                specs = self._process_interactions_row(row)
                
                for transcript_ID in matching_transcript_IDs:
                    try:
                        self._add_claim_item_qualifier(
                            sRNA_item_id, specs["claim"], transcript_ID,
                            specs["quali_start"],
                            specs["start_pos_sRNA_plex"],
                            specs["quali_end"], specs["end_pos_sRNA_plex"],
                            specs["quali_method"], property_dict["RNAplex"])
                        self._add_claim_item_qualifier(
                            transcript_ID, specs["claim"], sRNA_item_id,
                            specs["quali_start"],
                            specs["start_pos_target_plex"],
                            specs["quali_end"],
                            specs["end_pos_target_plex"],
                            specs["quali_method"],
                            property_dict["RNAplex"])
                        print("created interaction (RNAplex) between "
                              "sRNA {} "
                              "and transcript {}".format(sRNA_item_id,
                                                         transcript_ID))
                        self._make_log_entry_create_interaction(
                            "plex", sRNA_item_id, transcript_ID)
                    except:
                        print("interaction (RNAplex) between sRNA {} "
                              "and transcript {} already exists".format(
                                  sRNA_item_id, transcript_ID))
                        self._make_log_entry_interaction_exists(
                            "plex", sRNA_item_id, transcript_ID)

                    try:
                        self._add_claim_item_qualifier(
                            sRNA_item_id, specs["claim"],
                            transcript_ID, specs["quali_start"],
                            specs["start_pos_sRNA_up"],
                            specs["quali_end"],
                            specs["end_pos_sRNA_up"],
                            specs["quali_method"],
                            property_dict["RNAup"])
                        self._add_claim_item_qualifier(
                            transcript_ID, specs["claim"],
                            sRNA_item_id, specs["quali_start"],
                            specs["start_pos_target_up"],
                            specs["quali_end"],
                            specs["end_pos_target_up"],
                            specs["quali_method"],
                            property_dict["RNAup"])
                        print("created interaction (RNAup) between "
                              "sRNA {} and transcript {}".format(
                                  sRNA_item_id, transcript_ID))
                        self._make_log_entry_create_interaction(
                            "up", sRNA_item_id, transcript_ID)
                        self.number_of_uploaded_items[
                            'interactions'] += 1
                    except:
                        print("interaction (RNAup) between sRNA {} "
                              "and transcript {} already exists".format(
                                  sRNA_item_id, transcript_ID))
                        self._make_log_entry_interaction_exists(
                            "up", sRNA_item_id, transcript_ID)
                        
    def _process_interactions_row(self, row):
        property_dict = self._get_property_dict()
        both_genome_positions_plex = row['sRNA_interacted_position_RNAplex']
        genome_pos_list_plex = both_genome_positions_plex.split("-")
        both_target_positions_plex = row['target_interacted_position_RNAplex']
        target_pos_list_plex = both_target_positions_plex.split("-")
        both_genome_positions_up = row['sRNA_interacted_position_RNAup']
        genome_pos_list_up = both_genome_positions_up.split("-")
        both_target_positions_up = row['target_interacted_position_RNAup']
        target_pos_list_up = both_target_positions_up.split("-")   
        sRNA_position_list = row['sRNA_position'].split('-')
        sRNA_start_position = sRNA_position_list[0]
        sRNA_end_position = sRNA_position_list[1]
        
        targets_locus_tag_and_short_name = row['target_locus_tag'].split("|")
        targets_locus_tag = targets_locus_tag_and_short_name[0]
        trgt_gene_entry_id = row["target_gene_ID"]
        specs = {
            "start_pos_sRNA_plex": genome_pos_list_plex[0],
            "end_pos_sRNA_plex": genome_pos_list_plex[1],
            "start_pos_target_plex": target_pos_list_plex[0],
            "end_pos_target_plex": target_pos_list_plex[1],
            "start_pos_sRNA_up": genome_pos_list_up[0],
            "end_pos_sRNA_up": genome_pos_list_up[1],
            "start_pos_target_up": target_pos_list_up[0],
            "end_pos_target_up": target_pos_list_up[1],
            "target_locus_position": row['target_position'],
            "target_strand": row['target_strand'],
            "claim": property_dict["physically interacts with"],
            "quali_start": property_dict["genomic_start"],
            "quali_end": property_dict["genomic_end"],
            "quali_method": property_dict["determination method"],
            "sRNA_name": row['sRNA'], "targets_locus_tag": targets_locus_tag,
            "sRNA_start_pos": sRNA_position_list[0],
            "sRNA_end_pos": sRNA_position_list[1],
            "entry_id": trgt_gene_entry_id}
        return(specs)
                                    
    def _get_matching_transcripts(self, sRNA_name, targets_locus_tag,
                                  sRNA_start_position, sRNA_end_position,
                                  entry_id):
        interacting_ncRNA_id = None
        real_transc = None
        id_locus_tag_dict = self.id_locus_tag_dict
        for item in id_locus_tag_dict.items():
            if (item[1]["locus_tag"] == targets_locus_tag and
                item[1]["feature_type"] == "gene" and
                item[1]["entry_id"] == entry_id):
                target_gene_id = item[0]
                target_get_part_of_ids_list = (
                    self._get_part_of_targets_ids_list(
                        target_gene_id))
                real_transc = self._check_if_targets_are_transcripts(
                    target_get_part_of_ids_list)
            elif item[1]["feature_type"] == "ncRNA":
                if (item[1]["item_alias"] == sRNA_name and
                    (str(item[1]["genomic_start"]) ==
                     str(sRNA_start_position)) and
                    str(item[1]["genomic_end"]) == str(sRNA_end_position)):
                    interacting_ncRNA_id = item[0]
        if (interacting_ncRNA_id is None or
            real_transc is None):
            return(None)
        else:
            return(interacting_ncRNA_id, real_transc)

    def _get_part_of_targets_ids_list(self, target_gene_id):
        target_get_part_of_ids_list = []
        property_dict = self._get_property_dict()
        item = pywikibot.ItemPage(self.repo, target_gene_id)
        item_dict = item.get()
        clm_dict = item_dict["claims"]
        clm_list = clm_dict[property_dict["part of"]]
        for clm in clm_list:
            clm_target = clm.getTarget()
            target_ID = clm_target.getID()
            target_get_part_of_ids_list.append(target_ID)
        return(target_get_part_of_ids_list)
            
    def _check_if_targets_are_transcripts(self, target_transcripts_ids_list):
        ids_of_transcripts = []
        property_dict = self._get_property_dict()
        for ID in target_transcripts_ids_list:
            item = pywikibot.ItemPage(self.repo, ID)
            item_dict = item.get()
            clm_dict = item_dict["claims"]
            clm_list = clm_dict[property_dict["instance_of"]]
            for clm in clm_list:
                clm_target = clm.getTarget()
                instance_of_target_id = clm_target.getID()
                if instance_of_target_id == property_dict["transcript"]:
                    ids_of_transcripts.append(ID)                   
        return(ids_of_transcripts)                       
                                           
    def _get_property_dict(self):
        if self.databank == "TillsWiki":
            property_dict = {"instance_of": "P6", "subclass_of": "P7",
                             "found_in_taxon": "P8", "genomic_start": "P9",
                             "genomic_end": "P10",
                             "strand_orientation": "P11", "RNA": "Q9",
                             "Forward_Strand": "Q11",
                             "Reverse_Strand": "Q12", "small_RNA": "Q194",
                             "DNA": "Q227", "gene": "Q254",
                             "protein": "Q267", "NCBI Locus tag": "P24",
                             "encodes": "P25","encoded by": "P26",
                             "NCBI Locus tag of associated gene": "P27",
                             "physically interacts with": "P15",
                             "has part": "P28",
                             "determination method": "P29",
                             "RNAplex": "Q357", "RNAup": "Q465",
                             "ribosomal RNA": "Q417", "transfer RNA": "Q424",
                             "TSS": "Q441", "part of": "P4",
                             "RefSeq": "Q36583", "ANNOgesic": "Q36584",
                             "transcript": "Q37105", "stated in": "P30"}
            return(property_dict)        
        elif self.databank == "Wikidata":
            property_dict = {"instance_of": "P31", "subclass_of": "P279",
                             "found_in_taxon": "P703",
                             "genomic_start": "P644",
                             "genomic_end": "P645",
                             "strand_orientation": "P2548", "RNA": "Q11053",
                             "Forward_Strand": "Q22809680",
                             "Reverse_Strand": "Q22809711",
                             "small_RNA": "Q24287527", "DNA": "Q7430",
                             "gene": "Q7187", "protein": "Q8054",
                             "NCBI Locus tag": "P2393", "encodes": "P688",
                             "encoded by": "P702",
                             "physically interacts with": "P129",
                             "has part": "P527",
                             "determination method": "P459",
                             "ribosomal RNA": "Q215980",
                             "transfer RNA": "Q201448",
                             "TSS": "Q2449354", "part of": "P361",
                             "RefSeq": "Q7307074",
                             "transcript": "Q26944990", "stated in": "P248",
                             "RNAplex": "Q27907827", "RNAup": "Q27907828"}
            return(property_dict)
                       
    def _add_claim_string(self, item_id, claim, target):
        new_item = pywikibot.ItemPage(self.repo, item_id) 
        claim = pywikibot.Claim(self.repo, claim)
        claim.setTarget(target)
        new_item.addClaim(claim)
        
    def _add_claim_item(self, item_id, claim, target):
        new_item = pywikibot.ItemPage(self.repo, item_id) 
        claim = pywikibot.Claim(self.repo, claim)
        target = pywikibot.ItemPage(self.repo, target)
        claim.setTarget(target)
        new_item.addClaim(claim)

    def _add_claim_item_qualifier(self, item_id, claim, target, qualifier,
                                  qualifier_target, qualifier2,
                                  qualifier_target2, qualifier3,
                                  qualifier_target3):
        new_item = pywikibot.ItemPage(self.repo, item_id) 
        claim = pywikibot.Claim(self.repo, claim)
        target = pywikibot.ItemPage(self.repo, target)
        claim.setTarget(target)
        new_item.addClaim(claim)
        qualifier = pywikibot.Claim(self.repo, qualifier)
        qualifier.setTarget(qualifier_target)
        claim.addQualifier(qualifier)
        qualifier2 = pywikibot.Claim(self.repo, qualifier2)
        qualifier2.setTarget(qualifier_target2)
        claim.addQualifier(qualifier2)
        qualifier3 = pywikibot.Claim(self.repo, qualifier3)
        qualifier_target3 = pywikibot.ItemPage(self.repo, qualifier_target3)
        qualifier3.setTarget(qualifier_target3)
        claim.addQualifier(qualifier3)
   
    def _create_new_item(self, label, item_description):
        data = self._get_data_for_new_item(label, item_description)
        item = pywikibot.ItemPage(self.repo)
        item.editEntity(data)
        item_id = item.getID()
        new_item = pywikibot.ItemPage(self.repo, item_id) 
        return(item_id)

    def _get_data_for_new_item(self, label, item_description):
        return {
             'labels': {
                'en': {
                    'language': 'en',
                    'value': label
                }
            },
            'descriptions': {
                'en': {
                    'language': 'en',
                    'value': item_description
                }
            }
        }
         


