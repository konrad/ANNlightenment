import unittest
from annlightenmentlib.bacterialannotationbot import BacterialAnnotationBot 

class TestLabelDescriptionData(unittest.TestCase):

    def test_label_description_data(self):
        item_name = ("Staphylococcus aureus subsp. aureus NCTC "
                     "8325 tran814 826489 831124")
        item_description = ("bacterial transcript found in Staphylococcus "
                            "aureus subsp. aureus NCTC 8325")
        data = {
             'labels': {
                'en': {
                    'language': 'en',
                    'value': ("Staphylococcus aureus subsp. aureus NCTC "
                              "8325 tran814 826489 831124")
                }
            },
            'descriptions': {
                'en': {
                    'language': 'en',
                    'value': ("bacterial transcript found in Staphylococcus "
                              "aureus subsp. aureus NCTC 8325")
                }
            }
        }
        self.assertEqual(
            BacterialAnnotationBot._get_data_for_new_item(
                self, item_name, item_description), data)
