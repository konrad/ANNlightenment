import pywikibot

class DeleteItems():

    def __init__(self, path_to_logfile):
        self.path_to_logfile = path_to_logfile

    def _return_ID_list_for_deletion(self):
        ID_for_deletion_list = []
        placeholder_item_type_list = [
            "transcript", "sRNA", "ncRNA", "TSS", "gene", "protein",
            "sRNA gene", "rRNA", "tRNA"]
        logfile = open(self.path_to_logfile, 'r')
        for line in logfile.readlines():
            for item_type in placeholder_item_type_list:
                if line.startswith("INFO:root:created {} item with "
                                   "ID".format(item_type)):
                    line_list = line.split()
                    ID_for_deletion_list.append(line_list[-1])
        return(ID_for_deletion_list)

    def _delete(self, item_id_list):
        for item_id in item_id_list:
            site = pywikibot.Site("en", "TillsWiki")
            repo = site.data_repository()
            repoitem = pywikibot.ItemPage(repo, item_id)
            reason = "deleted by ANNlightenment"
            try:
                repoitem.delete(reason, prompt=False)
                print("Deleted item {}".format(item_id))
            except:
                continue
