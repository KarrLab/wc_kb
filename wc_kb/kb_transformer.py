""" Parsing of the core.xlsx unstructured data file and insering fields to kb_core.xlsx

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2019-02-23
:Copyright: 2019, Karr Lab
:License: MIT
"""

import wc_kb
import json
import openpyxl
import pdb
import re

class KbTransformer:

    def __init__(self, core=None, kb=None):
        if core is None:
            self.core = openpyxl.load_workbook('kbs/core.original.xlsx')
        if kb is None:
            self.kb = openpyxl.load_workbook('kbs/kb_core.xlsx')

    @staticmethod
    def parseJsonStrs(fieldStr):
        """ Detects and return individual JSON strings within fieldStr """

        jsons = []
        for match in re.findall(r'{.*?}', fieldStr):

            try:
                jsons.append(json.loads(match))
            except:
                pdb.set_trace()
                raise Exception('Can not parse JSONs from line: {}'.format(fieldStr))

        return jsons

    def transcodeReferences(self):
        """ Parses information from the References's 'Cross references' column and inserts it to the appropiate field(s) in the new kb structure. """

        coreRefs = self.core['References']
        kbRefs   = self.kb['References']
        cell2str = lambda myCell: myCell or ""

        for row_idx in range(3, coreRefs.max_row+1):

            # Make sure IDs are matching;
            assert(kbRefs.cell(column=1, row=row_idx-1).value == coreRefs.cell(column=1, row=row_idx).value)

            # Skip if there is no database ref
            coreDbref = coreRefs.cell(column=4, row=row_idx).value
            if coreDbref is None:
                continue

            # Get list of JSONs in string
            jsons = self.parseJsonStrs(coreDbref)
            for json in jsons:
                print(json)
                if json['source'] == 'URL':
                    myCell = kbRefs.cell(column=12, row=row_idx-1).value
                    kbRefs.cell(column=12, row=row_idx-1).value = cell2str(myCell) + str(json['xid']) + ', '
                else:
                    myCell = kbRefs.cell(column=11, row=row_idx-1).value
                    kbRefs.cell(column=11, row=row_idx-1).value = cell2str(myCell) + "{}:{}".format(json['source'], json['xid']) + ', '

            # Remove " ," fromt end of list
            if kbRefs.cell(column=11, row=row_idx-1).value is not None:
                kbRefs.cell(column=11, row=row_idx-1).value  = kbRefs.cell(column=11, row=row_idx-1).value[:-2]

            if kbRefs.cell(column=12, row=row_idx-1).value is not None:
                kbRefs.cell(column=12, row=row_idx-1).value  = kbRefs.cell(column=12, row=row_idx-1).value[:-2]

        self.kb.save('kbs/kb_core_ParseAdded2.xlsx')


        """
        if json['source'] == 'URL':
            # Either create new str or add to existing string
            if kbRefs.cell(column=12, row=row_idx-1).value is None:
                kbRefs.cell(column=12, row=row_idx-1).value = "{}:{}".format(json['source'], json['xid'])
            else:
                assert isinstance(kbRefs.cell(column=12, row=row_idx-1).value, str)
                kbRefs.cell(column=12, row=row_idx-1).value = "{}:{}".format(json['source'], json['xid'])
        else:
            # Either create new str or add to existing string
            if kbRefs.cell(column=11, row=row_idx-1).value is None:
                kbRefs.cell(column=11, row=row_idx-1).value = "{}:{}".format(json['source'], json['xid'])
            else:
                assert isinstance(kbRefs.cell(column=11, row=row_idx-1).value, str)
                kbRefs.cell(column=11, row=row_idx-1).value = "{}:{}".format(json['source'], json['xid'])
        """
