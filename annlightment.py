import argparse
import pywikibot
from annlightmentlib.bacterialannotationbot import BacterialAnnotationBot
from annlightmentlib.delete_items import DeleteItems

def main():
    parser = argparse.ArgumentParser()    
    parser.add_argument("--version", "-v", default=False, action="store_true",
                        help="show version")

    subparsers = parser.add_subparsers(help="commands")
 
    upload_parser = subparsers.add_parser("upload", help="subcommand to upload "
                                          "genes, products (proteins, rRNAs, "
                                          "tRNAs, sRNAs), transcripts and sRNA "
                                          "interactions")    
    upload_parser.add_argument("ANNOgesic_merge_gff", help="the path to the "
                               "ANNOgesic ...merge_features.gff file")
    upload_parser.add_argument("ANNOgesic_merge_csv", help="the path to "
                               "the ANNOgesic ..._merge.csv file that "
                               "contains the sRNA interactions")
    upload_parser.add_argument("strain_id", help="the ID (Q-number) of the "
                               "item that describes the strain you are working "
                               "with")
    upload_parser.add_argument("--databank", default = "TillsWiki", help = ""
                               "the databank you want to use. Chose from Tills"
                               "Wiki or Wikidata. Default is TillsWiki") 
    upload_parser.set_defaults(func=upload_items)

    test_parser = subparsers.add_parser("test", help="subcommand to upload "
                                          "genes, products (proteins, rRNAs, "
                                          "tRNAs, sRNAs), transcripts and sRNA "
                                          "interactions") 
    test_parser.add_argument("ANNOgesic_merge_gff", help="the path to the "
                               "ANNOgesic ...merge_features.gff file")
    test_parser.add_argument("ANNOgesic_merge_csv", help="the path to "
                               "the ANNOgesic ..._merge.csv file that "
                               "contains the sRNA interactions")
    test_parser.add_argument("strain_id", help="the ID (Q-number) of the "
                               "item that describes the strain you are working "
                               "with")
    test_parser.add_argument("--databank", default = "TillsWiki", help = ""
                               "the databank you want to use. Chose from Tills"
                               "Wiki or Wikidata. Default is TillsWiki") 
    test_parser.set_defaults(func=test)

    
    
    delete_parser = subparsers.add_parser("delete", help="subcommand to delete "
                                          "items")
    delete_parser.add_argument("path_to_logfile", help="the path to the log "
                               "file that contains the IDs of the items which "
                               "should be deleted")
    delete_parser.set_defaults(func=delete_items)

    args = parser.parse_args()
       
    if args.version is True:
        print("ANNlightment version 0")
    elif "func" in dir(args):
        args.func(args)
    else:
        parser.print_help()
        
def upload_items(args):
    site = _return_database_site(args)
    BAB = BacterialAnnotationBot(
        site, args.ANNOgesic_merge_gff, args.ANNOgesic_merge_csv,
        args.strain_id, args.databank)
    BAB.all_features()
    
def _return_database_site(args):
    if args.databank == "Wikidata":
        print("working with Wikidata")
        site = pywikibot.Site("wikidata", "wikidata")
    elif args.databank == "TillsWiki":
        print("working with TillsWiki")
        site = pywikibot.Site("en", "TillsWiki")
    return(site)

def test(args):
    site = _return_database_site(args)
    BAB = BacterialAnnotationBot(
        site, args.ANNOgesic_merge_gff, args.ANNOgesic_merge_csv,
        args.strain_id, args.databank)
    BAB.test()

def delete_items(args):
    deletion = DeleteItems(args.path_to_logfile)
    id_list_for_deletion = deletion._return_ID_list_for_deletion()
    deletion._delete(id_list_for_deletion)

main()
