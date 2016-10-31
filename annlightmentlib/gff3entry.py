class Gff3Entry(object):

    def __init__(self, entry_dict):
        self.seq_id = entry_dict["seq_id"]
        self.source = entry_dict["source"]
        self.feature = entry_dict["feature"]
        # 1-based coordinates
        # Make sure that start <= end
        start, end = sorted([int(entry_dict["start"]), int(entry_dict["end"])])
        self.start = start
        self.end = end
        self.score = entry_dict["score"]
        self.strand = entry_dict["strand"]
        self.phase = entry_dict["phase"]
        self.attributes = self._attributes(entry_dict["attributes"])
        self.attribute_string = entry_dict["attributes"]
    
    def _attributes(self, attributes_string):
        """Translate the attribute string to dictionary"""
        if attributes_string is None:
            return({})
        if attributes_string.endswith(";"):
            attributes_string = attributes_string[:-1]
        return dict(
            [key_value_pair.split("=")
             for key_value_pair in attributes_string.split(";")])

    def __str__(self):
        return "\t".join([str(field) for field in [
            self.seq_id, self.source, self.feature, self.start,
            self.end, self.score, self.strand, self.phase,
            self.attribute_string]])
