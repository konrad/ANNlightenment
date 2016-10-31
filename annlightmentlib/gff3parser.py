import csv
from annlightmentlib.gff3entry import Gff3Entry

class Gff3Parser(object):
    """
    A format description can be found at:
    http://genome.ucsc.edu/FAQ/FAQformat.html#format3
    http://www.sequenceontology.org/gff3.shtml
    a validator can be found here:
    http://modencode.oicr.on.ca/cgi-bin/validate_gff3_online
    """

    def entries(self, input_gff_fh):
        """
        """
        for entry_dict in csv.DictReader(
            input_gff_fh, delimiter="\t",
            fieldnames=["seq_id", "source", "feature", "start",
                        "end", "score", "strand", "phase", "attributes"]):
            if entry_dict["seq_id"].startswith("#"):
                continue
            yield(self._dict_to_entry(entry_dict))
            # try:
            #     yield(self._dict_to_entry(entry_dict))
            # except:
            #     sys.stderr.write(
            #         "Error! Please make sure that you use GFF3 formated "
            #         "annotation files. GTF2/GTF is not valid and the usage "
            #         "of that format is not recommended anymore (see "
            #         "http://www.sequenceontology.org/gff3.shtml for more "
            #         "information).\n")
            #     sys.exit(0)
    
    def _dict_to_entry(self, entry_dict):
        return Gff3Entry(entry_dict)
