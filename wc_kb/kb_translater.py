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

class KbTranslater:

    def __init__(self, core_path=None, kb_path=None):

        if core_path is None:
            core_path = 'wc_kb/kbs/original_core.xlsx'
        if kb_path is None:
            kb_path = 'wc_kb/kbs/kb_core_empty.xlsx'

        self.core = openpyxl.load_workbook(core_path)
        self.kb = openpyxl.load_workbook(kb_path)

    def translateChromosomeFeatures(self):
        """ Parses information from the core's 'Chromosome Features' sheet and inserts them to the appropiate field(s) in the KB structure. """

        # Names and types will need cleanup

        coreFeats = self.core['Chromosome features']
        kbFeats   = self.kb['Chromosome features']

        for rowIdx in range(3, coreFeats.max_row+1):

            kbFeats.cell(column=1, row=rowIdx-1).value = coreFeats.cell(column=1, row=rowIdx).value  # copy ID
            kbFeats.cell(column=2, row=rowIdx-1).value = coreFeats.cell(column=2, row=rowIdx).value  # copy name
            kbFeats.cell(column=3, row=rowIdx-1).value = coreFeats.cell(column=5, row=rowIdx).value  # copy type
            kbFeats.cell(column=4, row=rowIdx-1).value = (coreFeats.cell(column=6, row=rowIdx).value).lower()  # copy polymer
            kbFeats.cell(column=5, row=rowIdx-1).value = coreFeats.cell(column=7, row=rowIdx).value  # copy coordinate
            kbFeats.cell(column=6, row=rowIdx-1).value = coreFeats.cell(column=8, row=rowIdx).value  # copy length
            kbFeats.cell(column=7, row=rowIdx-1).value = self.decodeDirection(coreFeats.cell(column=9, row=rowIdx).value) # copy direction
            kbFeats.cell(column=8, row=rowIdx-1).value = coreFeats.cell(column=10, row=rowIdx).value # copy intensity
            kbFeats.cell(column=9, row=rowIdx-1).value = coreFeats.cell(column=11, row=rowIdx).value # copy unit
            kbFeats.cell(column=12, row=rowIdx-1).value = coreFeats.cell(column=14, row=rowIdx).value # copy references
            kbFeats.cell(column=13, row=rowIdx-1).value = coreFeats.cell(column=13, row=rowIdx).value # copy comments

            # Add evidence
            if coreFeats.cell(column=12, row=rowIdx).value is None:
                continue

            eviDict = json.loads(coreFeats.cell(column=12, row=rowIdx).value)
            eviId   = 'EVI({}:intensity)'.format(coreFeats.cell(column=1, row=rowIdx).value)
            eviDict['temperature_unit'] = 'C' # Looked up from paper
            self.addEvidence(eviId, eviDict)

            # Add experiment - only need to do once since all values have same origin
            if rowIdx ==3:
                self.addExperiment(eviDict)

    def translateGenes(self):
        """ Parses information from the core's Genes sheet and inserts them to the appropiate field(s) in the KB structure. """

        coreGenes = self.core['Genes']
        kbGenes   = self.kb['Genes']

        for rowIdx in range(3, coreGenes.max_row+1):

            kbGenes.cell(column=1, row=rowIdx-1).value = coreGenes.cell(column=1, row=rowIdx).value  # copy ID
            kbGenes.cell(column=2, row=rowIdx-1).value = coreGenes.cell(column=2, row=rowIdx).value  # copy names
            kbGenes.cell(column=4, row=rowIdx-1).value = coreGenes.cell(column=3, row=rowIdx).value  # copy symbol
            kbGenes.cell(column=7, row=rowIdx-1).value = coreGenes.cell(column=7, row=rowIdx).value  # copy type
            kbGenes.cell(column=8, row=rowIdx-1).value = (coreGenes.cell(column=8, row=rowIdx).value).lower()  # copy polymer
            kbGenes.cell(column=10, row=rowIdx-1).value = coreGenes.cell(column=9, row=rowIdx).value  # copy start position
            kbGenes.cell(column=11, row=rowIdx-1).value = (coreGenes.cell(column=9, row=rowIdx).value + coreGenes.cell(column=10, row=rowIdx).value)-1  # copy end position
            kbGenes.cell(column=12, row=rowIdx-1).value = self.decodeIsEssential(coreGenes.cell(column=15, row=rowIdx).value) # copy IsEssential
            kbGenes.cell(column=16, row=rowIdx-1).value = coreGenes.cell(column=18, row=rowIdx).value  # copy comments
            kbGenes.cell(column=15, row=rowIdx-1).value = coreGenes.cell(column=19, row=rowIdx).value

            # Add synonyms
            coreSynonyms = coreGenes.cell(column=4, row=rowIdx).value
            if coreSynonyms is None:
                continue

            synonymsDict = json.loads(coreSynonyms)
            if synonymsDict
            kbGenes.cell(column=3, row=rowIdx-1).value = synonymsDict['name']





    def translateReferences(self):
        """ Parses information from the References's 'Cross references' column and inserts it to the appropiate field(s) in the new kb structure. """

        coreRefs = self.core['References']
        kbRefs   = self.kb['References']
        cell2str = lambda myCell: myCell or "" # Convert None to string or return original string

        for rowIdx in range(3, coreRefs.max_row+1):

            # Make sure IDs are matching; skip if there is no database ref
            assert(kbRefs.cell(column=1, row=rowIdx-1).value == coreRefs.cell(column=1, row=rowIdx).value)
            coreDbref = coreRefs.cell(column=4, row=rowIdx).value
            if coreDbref is None:
                continue

            # Get list of JSONs in string add them to DBrefs / comments
            jsons = self.parseJsonStrs(coreDbref)
            for json in jsons: # Each JSON is a dict with a single key (DB name) - value pair (xid)
                if json['source'] == 'URL':
                    myCell = kbRefs.cell(column=12, row=rowIdx-1).value
                    kbRefs.cell(column=12, row=rowIdx-1).value = cell2str(myCell) + str(json['xid']) + ', '
                else:
                    myCell = kbRefs.cell(column=11, row=rowIdx-1).value
                    kbRefs.cell(column=11, row=rowIdx-1).value = cell2str(myCell) + "{}:{}".format(json['source'], json['xid']) + ', '
                    self.addDatabaseReference(json)

            # Remove " ," fromt end of DBrefs / comments
            if kbRefs.cell(column=11, row=rowIdx-1).value is not None:
                kbRefs.cell(column=11, row=rowIdx-1).value  = kbRefs.cell(column=11, row=rowIdx-1).value[:-2]

            if kbRefs.cell(column=12, row=rowIdx-1).value is not None:
                kbRefs.cell(column=12, row=rowIdx-1).value  = kbRefs.cell(column=12, row=rowIdx-1).value[:-2]


    """ Auxiliary functions """
    def addEvidence(self, eviId, eviDict):
        wsName = 'Evidence'
        nextRow = self.kb[wsName].max_row+1

        if not self.idExists(wsName, eviId):
            self.kb[wsName].cell(column=1, row=nextRow).value = eviId
            self.kb[wsName].cell(column=2, row=nextRow).value = '' # Need to pass in object ID
            self.kb[wsName].cell(column=3, row=nextRow).value = 'intensity'
            self.kb[wsName].cell(column=5, row=nextRow).value = float(eviDict['value'])
            self.kb[wsName].cell(column=7, row=nextRow).value = eviDict['units']
            self.kb[wsName].cell(column=10, row=nextRow).value = eviDict['comments']

    def addExperiment(self, expDict):
        wsName = 'Experiment'
        nextRow = self.kb[wsName].max_row+1
        expId = 'EXP_{}'.format(str(nextRow-1).zfill(4))

        if not self.idExists(wsName, expId):
            self.kb[wsName].cell(column=1, row=nextRow).value  = expId
            self.kb[wsName].cell(column=5, row=nextRow).value  = expDict['species']
            self.kb[wsName].cell(column=7, row=nextRow).value  = expDict['media']
            self.kb[wsName].cell(column=8, row=nextRow).value  = expDict['temperature']
            self.kb[wsName].cell(column=9, row=nextRow).value  = expDict['temperature_unit']
            self.kb[wsName].cell(column=13, row=nextRow).value = expDict['references'][0] # Will not work for mulitple references

    def addDatabaseReference(self, json):
        wsName = 'Database references'
        nextRow  = self.kb[wsName].max_row+1 #self.nextRow(wsName)
        dbRefId = "{}:{}".format(json['source'], json['xid'])

        if not self.idExists(wsName, dbRefId):
            self.kb[wsName].cell(column=1, row=nextRow).value = "{}:{}".format(json['source'], json['xid'])
            self.kb[wsName].cell(column=2, row=nextRow).value = json['source']
            self.kb[wsName].cell(column=3, row=nextRow).value = str(json['xid'])

    def idExists(self, wsName, id):
        """ Check if ID exists within worksheet """

        for row in range(1, self.kb[wsName].max_row+1):
            if self.kb[wsName].cell(column=1, row=row).value == id:
                return True

        return False

    @staticmethod
    def parseJsonStrs(fieldStr):
        """ Detects and return individual JSON strings within fieldStr """

        jsons = []
        for match in re.findall(r'{.*?}', fieldStr):
            try:
                jsons.append(json.loads(match))
            except:
                raise Exception('Can not parse JSONs from line: {}'.format(fieldStr))

        return jsons

    @staticmethod
    def decodeIsEssential(fieldStr):
        if fieldStr == 'E':
            return True
        elif fieldStr == 'NE' or fieldStr == 'F' or fieldStr == 'NE*':
            return False
        elif fieldStr is None:
            return None
        else:
            raise ValueError('Unrecognized IsEssential value: {}.'.format(fieldStr))

    @staticmethod
    def decodeDirection(fieldStr):
        if fieldStr == 'f':
            return 'forward'
        elif fieldStr == 'r':
            return 'backward'
        else:
            raise ValueError('Unrecognized direction value.')
