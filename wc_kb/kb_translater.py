""" Parsing of the core.xlsx unstructured data file and insering fields to kb_core.xlsx

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2019-02-23
:Copyright: 2019, Karr Lab
:License: MIT
"""

import wc_kb
import json
import openpyxl
import warnings
import pdb
import re

class KbTranslater:

    def __init__(self, core_path=None, kb_path=None):

        if core_path is None:
            core_path = 'wc_kb/kbs/original_core.xlsx'
        if kb_path is None:
            kb_path = 'wc_kb/kbs/kb_core_empty.xlsx'

        self.core = openpyxl.load_workbook(core_path, read_only=True)
        self.kb = openpyxl.load_workbook(kb_path)

    def translateChromosomeFeatures(self):
        """ Parses information from the core's 'Chromosome Features' sheet and inserts them to the appropiate field(s) in the KB structure. """

        coreChromFeats = self.core['Chromosome features']
        kbChromFeats   = self.kb['Chromosome features']

        assert kbChromFeats.cell(column=1, row=1).value  == 'Id'
        assert kbChromFeats.cell(column=2, row=1).value  == 'Name'
        assert kbChromFeats.cell(column=3, row=1).value  == 'Type'
        assert kbChromFeats.cell(column=4, row=1).value  == 'Polymer'
        assert kbChromFeats.cell(column=5, row=1).value  == 'Coordinate'
        assert kbChromFeats.cell(column=6, row=1).value  == 'Length'
        assert kbChromFeats.cell(column=7, row=1).value  == 'Direction'
        assert kbChromFeats.cell(column=8, row=1).value  == 'Intensity'
        assert kbChromFeats.cell(column=9, row=1).value  == 'Unit'
        assert kbChromFeats.cell(column=10, row=1).value == 'Evidence'
        assert kbChromFeats.cell(column=11, row=1).value == 'Database references' #NO DB REFS IN CORE
        assert kbChromFeats.cell(column=12, row=1).value == 'References'
        assert kbChromFeats.cell(column=13, row=1).value == 'Comments'
        eviProperty = 'intensity'

        for rowIdx in range(3, 300): # coreChromFeats.max_row+1):
            print(rowIdx)

            kbChromFeats.cell(column=1, row=rowIdx-1).value = coreChromFeats.cell(column=1, row=rowIdx).value.replace('-','_') # ID
            kbChromFeats.cell(column=2, row=rowIdx-1).value = coreChromFeats.cell(column=2, row=rowIdx).value # Name
            kbChromFeats.cell(column=3, row=rowIdx-1).value = coreChromFeats.cell(column=5, row=rowIdx).value # Type
            kbChromFeats.cell(column=4, row=rowIdx-1).value = (coreChromFeats.cell(column=6, row=rowIdx).value).lower() #Polymer
            kbChromFeats.cell(column=5, row=rowIdx-1).value = coreChromFeats.cell(column=7, row=rowIdx).value #Coordinate
            kbChromFeats.cell(column=6, row=rowIdx-1).value = coreChromFeats.cell(column=8, row=rowIdx).value #Length
            kbChromFeats.cell(column=7, row=rowIdx-1).value = self.decodeDirection(coreChromFeats.cell(column=9, row=rowIdx).value) #Direction
            kbChromFeats.cell(column=8, row=rowIdx-1).value = coreChromFeats.cell(column=10, row=rowIdx).value #Intensity
            # NO 'unit' entry in core: kbChromFeats.cell(column=9, row=rowIdx-1).value = coreChromFeats.cell(column=11, row=rowIdx).value #Unit

            # Add evidence
            eviDictStr = coreChromFeats.cell(column=12, row=rowIdx).value
            if eviDictStr is not None:
                eviDict = json.loads(eviDictStr)

                if rowIdx==3:   # Add experiment once - all data from same experiment
                    expId = self.addExperiment(eviDict)

                self.addEvidence(
                    objectId = kbChromFeats.cell(column=1, row=rowIdx-1).value,
                    property = eviProperty,
                    eviDict = eviDict,
                    eviColumn = 10,
                    expId = expId)

            # NO 'Cross references' entry in core: kbChromFeats.cell(column=11, row=rowIdx-1).value = coreChromFeats.cell(column=11, row=rowIdx).value
            # kbChromFeats.cell(column=12, row=rowIdx-1).value = coreChromFeats.cell(column=14, row=rowIdx).value # copy references
            kbChromFeats.cell(column=13, row=rowIdx-1).value = coreChromFeats.cell(column=13, row=rowIdx).value # copy comments

    def translateTranscriptionUnits(self):
        """ Parses information from the core's Transcription unit's sheet and inserts them to the appropiate field(s) in the KB structure. """

        coreTUs = self.core['Transcription units']
        kbTUs   = self.kb['Transcription units']

        assert kbTUs.cell(column=1, row=1).value == 'Id'
        assert kbTUs.cell(column=2, row=1).value == 'Name'
        assert kbTUs.cell(column=3, row=1).value == 'Type'
        assert kbTUs.cell(column=4, row=1).value == 'Polymer'
        assert kbTUs.cell(column=5, row=1).value == 'Strand'
        assert kbTUs.cell(column=6, row=1).value == 'Pribnow_start'
        assert kbTUs.cell(column=7, row=1).value == 'Pribnow_end'
        assert kbTUs.cell(column=8, row=1).value == 'Start'
        assert kbTUs.cell(column=9, row=1).value == 'End'
        assert kbTUs.cell(column=10, row=1).value == 'Genes'
        assert kbTUs.cell(column=11, row=1).value == 'Database references'
        assert kbTUs.cell(column=12, row=1).value == 'References'
        assert kbTUs.cell(column=13, row=1).value == 'Comments'

        for rowIdx in range(3, 300): #coreTUs.max_row+1):
            print(rowIdx)

            kbTUs.cell(column=1, row=rowIdx-1).value = coreTUs.cell(column=1, row=rowIdx).value #ID
            kbTUs.cell(column=2, row=rowIdx-1).value = coreTUs.cell(column=2, row=rowIdx).value #NAME
            kbTUs.cell(column=3, row=rowIdx-1).value = self.decodeRnaType(coreTUs.cell(column=9, row=rowIdx).value) #TYPE
            kbTUs.cell(column=4, row=rowIdx-1).value = coreTUs.cell(column=5, row=rowIdx).value.lower() #POLYMET
            kbTUs.cell(column=5, row=rowIdx-1).value = coreTUs.cell(column=8, row=rowIdx).value #STRAND
            kbTUs.cell(column=6, row=rowIdx-1).value = coreTUs.cell(column=11, row=rowIdx).value #PRIBNOW START
            kbTUs.cell(column=7, row=rowIdx-1).value = coreTUs.cell(column=12, row=rowIdx).value #PRIBNOW END
            kbTUs.cell(column=8, row=rowIdx-1).value = coreTUs.cell(column=6, row=rowIdx).value # START
            kbTUs.cell(column=9, row=rowIdx-1).value = coreTUs.cell(column=7, row=rowIdx).value # END
            kbTUs.cell(column=10, row=rowIdx-1).value = coreTUs.cell(column=10, row=rowIdx).value #GENES
            #NO 'Cross references' entry in core: kbTUs.cell(column=11, row=rowIdx-1).value = coreTUs.cell(column=4, row=rowIdx).value # DB REFS
            kbTUs.cell(column=12, row=rowIdx-1).value = coreTUs.cell(column=14, row=rowIdx).value #REFS
            kbTUs.cell(column=13, row=rowIdx-1).value = coreTUs.cell(column=13, row=rowIdx).value #COMMENTS

    def translateGenes(self):
        """ Parses information from the core's Genes sheet and inserts them to the appropiate field(s) in the KB structure. """

        coreGenes = self.core['Genes']
        kbGenes   = self.kb['Genes']

        assert kbGenes.cell(column=1, row=1).value == 'Id'
        assert kbGenes.cell(column=2, row=1).value == 'Name'
        assert kbGenes.cell(column=3, row=1).value == 'Synonyms'
        assert kbGenes.cell(column=4, row=1).value == 'Symbol'
        assert kbGenes.cell(column=5, row=1).value == 'Homologs'
        assert kbGenes.cell(column=7, row=1).value == 'Type'
        assert kbGenes.cell(column=8, row=1).value == 'Polymer'
        assert kbGenes.cell(column=9, row=1).value == 'Direction'
        assert kbGenes.cell(column=10, row=1).value == 'Start'
        assert kbGenes.cell(column=11, row=1).value == 'End'
        assert kbGenes.cell(column=12, row=1).value == 'Is essential'
        assert kbGenes.cell(column=14, row=1).value == 'Database references'
        assert kbGenes.cell(column=15, row=1).value == 'References'
        assert kbGenes.cell(column=16, row=1).value == 'Comments'
        eviProperty = 'is_essential'

        for rowIdx in range(3, 300): #coreGenes.max_row+1):
            print(rowIdx)

            kbGenes.cell(column=1, row=rowIdx-1).value = coreGenes.cell(column=1, row=rowIdx).value # ID
            kbGenes.cell(column=2, row=rowIdx-1).value = coreGenes.cell(column=2, row=rowIdx).value #NAME
            synonymsDicts = self.delinateJSONs(coreGenes.cell(column=4, row=rowIdx).value) #SYNONYMS
            kbGenes.cell(column=3, row=rowIdx-1).value = self.concatenateValues(synonymsDicts, 'name')
            kbGenes.cell(column=4, row=rowIdx-1).value = coreGenes.cell(column=3, row=rowIdx).value #SYMBOL
            homologsDicts = self.delinateJSONs(coreGenes.cell(column=6, row=rowIdx).value) #HOMOLOGS
            kbGenes.cell(column=5, row=rowIdx-1).value = self.concatenateValues(homologsDicts, ('species', 'xid'))

            # COLUMN 6: add COG categories

            kbGenes.cell(column=7, row=rowIdx-1).value = coreGenes.cell(column=7, row=rowIdx).value # TYPE
            kbGenes.cell(column=8, row=rowIdx-1).value = (coreGenes.cell(column=8, row=rowIdx).value).lower() # POLYMER
            kbGenes.cell(column=9, row=rowIdx-1).value = self.decodeDirection(coreGenes.cell(column=11, row=rowIdx).value) # DIRECTION
            kbGenes.cell(column=10, row=rowIdx-1).value = coreGenes.cell(column=9, row=rowIdx).value # START
            kbGenes.cell(column=11, row=rowIdx-1).value = (coreGenes.cell(column=9, row=rowIdx).value + coreGenes.cell(column=10, row=rowIdx).value)-1 # END
            kbGenes.cell(column=12, row=rowIdx-1).value = self.decodeIsEssential(coreGenes.cell(column=15, row=rowIdx).value) # IS ESSENTIAL

            # Add evidence
            eviDictStr = coreGenes.cell(column=17, row=rowIdx).value
            if eviDictStr is not None:
                eviDict = json.loads(eviDictStr)
                if rowIdx==3:   # Add experiment once - all data from same experiment
                    expId = self.addExperiment(eviDict)

                self.addEvidence(
                    objectId = kbGenes.cell(column=1, row=rowIdx-1).value,
                    property = eviProperty,
                    eviDict = eviDict,
                    eviColumn = 13,
                    expId = expId)

            if rowIdx==3:   # Add experiment once - all data from same experiment
                self.addExperiment(eviDict)


            dbRefsDicts = self.delinateJSONs(coreGenes.cell(column=5, row=rowIdx).value)
            #kbGenes.cell(column=14, row=rowIdx-1).value = self.addDatabaseReference(dbRefsDicts, returnIds=True) # Add DB refs
            #kbGenes.cell(column=15, row=rowIdx-1).value = coreGenes.cell(column=19, row=rowIdx).value # References
            #kbGenes.cell(column=16, row=rowIdx-1).value = coreGenes.cell(column=18, row=rowIdx).value # Comments

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
            jsons = self.delinateJSONs(coreDbref)
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


    def addEvidence(self, objectId, property, eviDict, eviColumn, expId):
        if eviDict is None:
            return None

        wsName = 'Evidence'
        nextRow = self.kb[wsName].max_row+1
        eviId   = 'EVI({}:{})'.format(objectId, property)

        #if not self.idExists(wsName, eviId):
        self.kb[wsName].cell(column=1, row=nextRow).value = eviId
        self.kb[wsName].cell(column=2, row=nextRow).value = 'cell' # Need to pass in object ID
        self.kb[wsName].cell(column=3, row=nextRow).value = objectId
        self.kb[wsName].cell(column=4, row=nextRow).value = property
        self.kb[wsName].cell(column=9, row=nextRow).value = expId

        try: self.kb[wsName].cell(column=5, row=nextRow).value = eviDict['value']
        except: pass

        try: self.kb[wsName].cell(column=6, row=nextRow).value = eviDict['means']
        except: pass

        try: self.kb[wsName].cell(column=7, row=nextRow).value = eviDict['STD']
        except: pass

        try: self.kb[wsName].cell(column=8, row=nextRow).value = eviDict['units']
        except: pass

        try: self.kb[wsName].cell(column=11, row=nextRow).value = eviDict['comments']
        except: pass

        #else:
        #    warnings.warn('\nTried to create evidence object twice: \n\t {}'.format(eviId))

        kbChromFeats = self.kb['Chromosome features']
        kbChromFeats.cell(column=eviColumn, row=kbChromFeats.max_row).value = self.addComma(
            kbChromFeats.cell(column=eviColumn, row=kbChromFeats.max_row).value) + eviId

    def addExperiment(self, expDict):
        wsName = 'Experiment'
        nextRow = self.kb[wsName].max_row+1
        expId = 'EXP_{}'.format(str(nextRow-1).zfill(4))

        if not self.idExists(wsName, expId):
            self.kb[wsName].cell(column=1, row=nextRow).value  = expId
            self.kb[wsName].cell(column=5, row=nextRow).value  = expDict['species']
            self.kb[wsName].cell(column=7, row=nextRow).value  = expDict['media']
            self.kb[wsName].cell(column=8, row=nextRow).value  = expDict['temperature']
            #self.kb[wsName].cell(column=9, row=nextRow).value  = expDict['temperature_unit']
            self.kb[wsName].cell(column=13, row=nextRow).value = expDict['references'][0] # Will not work for mulitple references

        return expId

    def addDatabaseReference(self, jsons, returnIds=False):

        if jsons is None:
            return None

        assert isinstance(jsons, list)
        assert len(jsons) > 0

        wsName = 'Database references'
        dbRefIds = ''

        for json in jsons:
            nextRow  = self.kb[wsName].max_row+1
            dbRefId = "{}:{}".format(json['source'], str(json['xid']))

            # Clean up formatting; can condense code?
            dbRefId = dbRefId.replace('-','_')

            if returnIds is True:
                dbRefIds += dbRefId + ', '

            if not self.idExists(wsName, dbRefId):
                self.kb[wsName].cell(column=1, row=nextRow).value = dbRefId
                self.kb[wsName].cell(column=2, row=nextRow).value = json['source']
                self.kb[wsName].cell(column=3, row=nextRow).value = str(json['xid'])

        if returnIds is True:
            return dbRefIds[:-2]

    """ Auxiliary functions """
    def idExists(self, wsName, id):
        """ Check if ID exists within worksheet """

        for row in range(1, self.kb[wsName].max_row+1):
            if self.kb[wsName].cell(column=1, row=row).value == id:
                return True

        return False

    @staticmethod
    def delinateJSONs(fieldStr):
        """ Detects and return individual JSON strings within fieldStr """

        if fieldStr is None or fieldStr=='':
            return None

        jsons = []
        jsonStrLen = 0
        matches = re.findall(r'{.*?}', fieldStr)
        for match in matches:
            jsonStrLen += len(match)
            try:
                jsons.append(json.loads(match))
            except:
                raise Exception('\nCan not parse JSONs from line: \n\t {}'.format(fieldStr))

        # Check if whole fieldStr has been parsed to JSONs; assumes JSONS are separated as "}, {"
        if jsonStrLen+(len(matches)-1)*2 != len(fieldStr):
            warnings.warn('\nString contains non-JSON substrings: \n\t {}'.format(fieldStr))

        # Return None if only non-JSON substrings are found
        if jsons ==[]:
            return None

        return jsons

    @staticmethod
    def concatenateValues(jsons, keys):
        """ Concatenate values from jsons """

        assert isinstance(jsons, list) or jsons is None
        assert isinstance(keys, tuple) or isinstance(keys, str)
        if jsons is None or len(jsons)==0:
            return None
        if isinstance(keys, tuple):
            assert 0<len(keys)<=2

        values = ''
        for json in jsons:
            if len(keys)==2:
                entry = str(json[keys[0]]) + ':' + str(json[keys[1]])
            else:
                entry = str(json[keys])
            values += entry + ', '

        # Return string without final ', '
        return values[:-2]

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

    @staticmethod
    def decodeCOGcats(fieldStr):
        pass

    @staticmethod
    def decodeRnaType(fieldStr):
        if fieldStr.startswith("mixed_"):
            return 'mixed'
        else:
            return fieldStr

    @staticmethod
    def addComma(fieldStr):
        if fieldStr is None:
            return ''
        else:
            return fieldStr + ', '
