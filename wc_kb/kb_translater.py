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

    def __init__(self, core_path=None, kb_path=None, column=None):

        if core_path is None:
            core_path = 'wc_kb/kbs/original_core.xlsx'
        if kb_path is None:
            kb_path = 'wc_kb/kbs/kb_core_empty.xlsx'

        self.core = openpyxl.load_workbook(core_path, read_only=True)
        self.kb = openpyxl.load_workbook(kb_path)

        # Construct header dictionary
        self.header={}
        for ws in self.kb.worksheets:
            self.header[ws.title.lower()]={}
            self.header[ws.title]={}
            for columnIdx in range(1,50):
                if ws.cell(row=1, column=columnIdx).value is not None:
                    self.header[ws.title.lower()][ws.cell(row=1, column=columnIdx).value.lower()] = columnIdx
                    self.header[ws.title][ws.cell(row=1, column=columnIdx).value] = columnIdx

    def translateMetabolites(self):
        """ Parses information from the core's Metabolite sheet and inserts them to the appropiate field(s) in the KB structure. """

        coreMetas = self.core['Metabolites']
        kbMetas   = self.kb['Metabolites']

        for rowIdx in range(3, coreMetas.max_row+1):
            print(rowIdx)
            speciesTypeId = kbMetas.cell(column=1, row=rowIdx-1).value
            assert speciesTypeId == coreMetas.cell(column=1, row=rowIdx).value

            #MANUALLY COPIED: kbMetas.cell(column=1, row=rowIdx-1).value = coreMetas.cell(column=1, row=rowIdx).value # ID
            #MANUALLY COPIED: kbMetas.cell(column=2, row=rowIdx-1).value = coreMetas.cell(column=2, row=rowIdx).value # NAME
            #MANUALLY COPIED: kbMetas.cell(column=3, row=rowIdx-1).value = coreMetas.cell(column=5, row=rowIdx).value # Synonyms

            # Type
            metaType = coreMetas.cell(column=7, row=rowIdx).value # TYPE
            if metaType is None:
                kbMetas.cell(column=4, row=rowIdx-1).value = 'unknown'
            elif metaType not in wc_kb.core.MetaboliteSpeciesTypeType:
                kbMetas.cell(column=4, row=rowIdx-1).value = 'misc'
            elif metaType not in wc_kb.core.MetaboliteSpeciesTypeType:
                kbMetas.cell(column=4, row=rowIdx-1).value = metaType

            # Concentration
            if coreMetas.cell(column=18, row=rowIdx).value is not None:

                concDict = json.loads(coreMetas.cell(column=18, row=rowIdx).value)

                concentrationEntry = self.createConcentration(
                    speciesTypeId = speciesTypeId ,
                    values = concDict['concentration'],
                    units = 'mM')

                #evidenceEntry = self.createEvidence(
                #    objectId = speciesTypeId,
                #    property = 'concentration',
                #    values = coreMetas.cell(column=22, row=rowIdx).value,
                #    units = 'molecules',
                #    experiment = concExpEntry['id'],
                #    comments=coreMetas.cell(column=23, row=rowIdx).value)

                #if concentrationEntry is not None:
                #    kbMetas.cell(column=self.header['Metabolites']['Concentration'], row=rowIdx-1).value = concentrationEntry['id']
                #    if evidenceEntry is not None:
                #        self.insertCell(insertTo=concentrationEntry, header='Evidence', value=evidenceEntry['id'])

            # Species properties
            speciesPropsEntry = self.createSpeciesTypeProperties(
                speciesTypeId = speciesTypeId,
                structure = coreMetas.cell(column=9, row=rowIdx).value)
            if speciesPropsEntry is not None:
                kbMetas.cell(column=self.header['Metabolites']['Species properties'], row=rowIdx-1).value = speciesPropsEntry['id']

            #MANUALLY COPIED: kbMetas.cell(column=8, row=rowIdx-1).value = coreMetas.cell(column=refsColumn, row=rowIdx).value # REFERENCES
            #MANUALLY COPIED: kbMetas.cell(column=9, row=rowIdx-1).value = coreMetas.cell(column=commentsColumn, row=rowIdx).value #comments
            kbMetas.cell(column=self.header['Metabolites']['Database references'], row=rowIdx-1).value = \
                self.parseDatabaseReferences(cell = coreMetas.cell(column=6, row=rowIdx)) # Add DB refs

    def translateProteins(self):
        """ Parses information from the core's 'protein monomers' sheet and inserts them to the appropiate field(s) in the KB structure. """

        header=self.header
        coreProteins = self.core['Protein monomers']
        kbProteins = self.kb['Proteins']
        kbSpeciesProperties = self.kb['Species properties']
        kbConcentrations = self.kb['Concentrations']

        # Create experiment entries
        concExpEntry = self.createExperiment(
                              species='Mycoplasma Pneumoniae',
                              geneticVariant='M129',
                              externalMedia='Hayflick',
                              temperature=37,
                              refs='PUB_0932',
                              comments='Concentrations available at multiple time points, see the Comments in the evidence field')
        hlExpEntry = self.createExperiment(
                              species='Mycoplasma Pneumoniae',
                              geneticVariant='M129',
                              externalMedia='Hayflick',
                              temperature=37,
                              refs='PUB_0956')
        methionineExpEntry = self.createExperiment(
                              species='Mycoplasma Pneumoniae',
                              geneticVariant='M129',
                              refs='PUB_0937')

        for rowIdx in range(3, coreProteins.max_row+1):
            print(rowIdx)

            # Check if row is correct
            speciesTypeId = kbProteins.cell(column=1, row=rowIdx-1).value
            assert speciesTypeId == coreProteins.cell(column=1, row=rowIdx).value

            #MANUALLY COPIED: kbProteins.cell(column=1, row=rowIdx-1).value = coreProteins.cell(column=1, row=rowIdx).value # ID
            #MANUALLY COPIED: kbProteins.cell(column=2, row=rowIdx-1).value = coreProteins.cell(column=2, row=rowIdx).value # Name
            #MANUALLY COPIED: kbProteins.cell(column=3, row=rowIdx-1).value = coreProteins.cell(column=3, row=rowIdx).value # Synonyms

            # Types
            proteinType = coreProteins.cell(column=5, row=rowIdx).value
            if proteinType is None:
                kbProteins.cell(column=4, row=rowIdx-1).value = 'uncategorized'
            else:
                kbProteins.cell(column=4, row=rowIdx-1).value = coreProteins.cell(column=5, row=rowIdx).value

            # Species properties
            speciesPropsEntry = self.createSpeciesTypeProperties(
                speciesTypeId = speciesTypeId,
                halfLife = coreProteins.cell(column=27, row=rowIdx).value,
                domains= coreProteins.cell(column=14, row=rowIdx).value,
                prostheticGroups = coreProteins.cell(column=15, row=rowIdx).value)
            evidenceEntry = self.createEvidence(
                objectId = speciesTypeId,
                property = 'half life',
                values = coreProteins.cell(column=27, row=rowIdx).value,
                units = 'min',
                experiment = hlExpEntry['id'])
            if speciesPropsEntry is not None:
                kbProteins.cell(column=self.header['Proteins']['Species properties'], row=rowIdx-1).value = speciesPropsEntry['id']
                if evidenceEntry is not None:
                    self.insertCell(insertTo=speciesPropsEntry, header='Evidence', value=evidenceEntry['id'])

            # Concentration
            concentrationEntry = self.createConcentration(
                speciesTypeId = speciesTypeId ,
                values = coreProteins.cell(column=22, row=rowIdx).value,
                units = 'molecules')
            evidenceEntry = self.createEvidence(
                objectId = speciesTypeId,
                property = 'concentration',
                values = coreProteins.cell(column=22, row=rowIdx).value,
                units = 'molecules',
                experiment = concExpEntry['id'],
                comments=coreProteins.cell(column=23, row=rowIdx).value)
            if concentrationEntry is not None:
                kbProteins.cell(column=self.header['Proteins']['Concentration'], row=rowIdx-1).value = concentrationEntry['id']
                if evidenceEntry is not None:
                    self.insertCell(insertTo=concentrationEntry, header='Evidence', value=evidenceEntry['id'])

            #MANUALLY COPIED: kbProteins.cell(column=5, row=rowIdx-1).value = coreProteins.cell(column=6, row=rowIdx).value # Gene
            kbProteins.cell(column=self.header['Proteins']['Localization'], row=rowIdx-1).value = coreProteins.cell(column=9, row=rowIdx).value.replace('t', '') # localization: treat trans- membrane
            #MANUALLY COPIED: kbProteins.cell(column=9, row=rowIdx-1).value = coreProteins.cell(column=6, row=rowIdx).value # Translation rates
            #MANUALLY COPIED: kbProteins.cell(column=10, row=rowIdx-1).value = coreProteins.cell(column=6, row=rowIdx).value # Translation rate units

            evidenceEntry = self.createEvidence(
                objectId = speciesTypeId,
                property = 'translationRate',
                values = coreProteins.cell(column=24, row=rowIdx).value,
                units = '1/s/mRNA',
                experiment = hlExpEntry['id']) #measured in the same experiment as half lifes
            if evidenceEntry is not None:
                proteinEntry ={'ws':'Proteins', 'id':speciesTypeId, 'row':rowIdx-1}
                self.insertCell(insertTo=proteinEntry, header='Evidence', value=evidenceEntry['id'])

            #MANUALLY COPIED: kbProteins.cell(column=11, row=rowIdx-1).value = coreProteins.cell(column=6, row=rowIdx).value # Is methionine cleaved?

            evidenceEntry = self.createEvidence(
                objectId = speciesTypeId,
                property = 'Is methionine cleaved',
                values = coreProteins.cell(column=7, row=rowIdx).value,
                units = 'dimensionless',
                experiment = methionineExpEntry['id'])
            if evidenceEntry is not None:
                proteinEntry ={'ws':'Proteins', 'id':speciesTypeId, 'row':rowIdx-1}
                self.insertCell(insertTo=proteinEntry, header='Evidence', value=evidenceEntry['id'])

            #MANUALLY COPIED: kbProteins.cell(column=7, row=rowIdx-1).value = coreProteins.cell(column=10, row=rowIdx).value # signal sequence type
            #MANUALLY COPIED: kbProteins.cell(column=8, row=rowIdx-1).value = coreProteins.cell(column=11, row=rowIdx).value # signal sequence localization
            #MANUALLY COPIED: kbProteins.cell(column=9, row=rowIdx-1).value = coreProteins.cell(column=12, row=rowIdx).value # signal sequence length
            #MANUALLY COPIED: kbProteins.cell(column=10, row=rowIdx-1).value = coreProteins.cell(column=16, row=rowIdx).value # DNA footprint length
            #MANUALLY COPIED: kbProteins.cell(column=11, row=rowIdx-1).value = coreProteins.cell(column=17, row=rowIdx).value # DNA footprint biding
            #MANUALLY COPIED: kbProteins.cell(column=16, row=rowIdx-1).value = coreProteins.cell(column=33, row=rowIdx).value # References
            #MANUALLY COPIED: kbProteins.cell(column=17, row=rowIdx-1).value = coreProteins.cell(column=34, row=rowIdx).value # Comments

            kbProteins.cell(column=self.header['Proteins']['Database references'], row=rowIdx-1).value = \
                self.parseDatabaseReferences(cell = coreProteins.cell(column=4, row=rowIdx)) # Add DB refs

    def createSpeciesTypeProperties(self, speciesTypeId, structure=None, halfLife=None, halfLifeUnits='min', domains=None, prostheticGroups=None, evidence=None, dbRefs=None, refs=None, comments=None):
        inputs = [structure, halfLife, domains, prostheticGroups]
        if all(input is None for input in inputs):
            return None

        wsName = 'Species properties'
        nextRow = self.kb[wsName].max_row+1
        speciesPropId = 'props_{}'.format(speciesTypeId)

        self.kb[wsName].cell(column=self.header['species properties']['id'], row=nextRow).value = speciesPropId
        self.kb[wsName].cell(column=self.header['species properties']['structure'], row=nextRow).value = structure
        self.kb[wsName].cell(column=self.header['species properties']['half life'], row=nextRow).value = halfLife
        self.kb[wsName].cell(column=self.header['species properties']['half life units'], row=nextRow).value = halfLifeUnits
        self.kb[wsName].cell(column=self.header['species properties']['domains'], row=nextRow).value = domains
        self.kb[wsName].cell(column=self.header['species properties']['prosthetic groups'], row=nextRow).value = prostheticGroups
        self.kb[wsName].cell(column=self.header['species properties']['evidence'], row=nextRow).value = evidence
        self.kb[wsName].cell(column=self.header['species properties']['database references'], row=nextRow).value = dbRefs
        self.kb[wsName].cell(column=self.header['species properties']['references'], row=nextRow).value = refs
        self.kb[wsName].cell(column=self.header['species properties']['comments'], row=nextRow).value = comments

        entry ={'ws':wsName, 'id':speciesPropId, 'row':nextRow}
        return entry

    def createConcentration(self, speciesTypeId, compartment='c', values=None, units=None, evidence=None, dbRefs=None, refs=None, comments=None):
        assert compartment in ['c', 'e', 'm']
        if values is None:
            return None

        wsName = 'Concentrations'
        nextRow = self.kb[wsName].max_row+1
        concentrationId = 'conc_{}_{}'.format(speciesTypeId, compartment)

        self.kb[wsName].cell(column=self.header['concentrations']['id'], row=nextRow).value = concentrationId
        self.kb[wsName].cell(column=self.header['concentrations']['species'], row=nextRow).value = '{}[{}]'.format(speciesTypeId, compartment)
        self.kb[wsName].cell(column=self.header['concentrations']['values'], row=nextRow).value = values
        self.kb[wsName].cell(column=self.header['concentrations']['units'], row=nextRow).value = units
        self.kb[wsName].cell(column=self.header['concentrations']['evidence'], row=nextRow).value = evidence
        self.kb[wsName].cell(column=self.header['concentrations']['database references'], row=nextRow).value = dbRefs
        self.kb[wsName].cell(column=self.header['concentrations']['references'], row=nextRow).value = refs
        self.kb[wsName].cell(column=self.header['concentrations']['comments'], row=nextRow).value = comments

        entry ={'ws':'Concentrations', 'id':concentrationId, 'row':nextRow}
        return entry

    def createEvidence(self, objectId, property=None, values=None, units=None, experiment=None, dbRefs=None, refs=None, comments=None):
        inputs = [property, values, units, experiment, dbRefs, refs, comments]
        #if all(input is None for input in inputs):
        if values is None:
            return None

        wsName = 'Evidence'
        nextRow = self.kb[wsName].max_row+1
        eviId   = 'EVI{}'.format(str(nextRow-1).zfill(4))

        self.kb[wsName].cell(column=self.header['evidence']['id'], row=nextRow).value = eviId
        self.kb[wsName].cell(column=self.header['evidence']['cell'], row=nextRow).value = 'cell'
        self.kb[wsName].cell(column=self.header['evidence']['object'], row=nextRow).value = objectId
        self.kb[wsName].cell(column=self.header['evidence']['property'], row=nextRow).value = property
        self.kb[wsName].cell(column=self.header['evidence']['values'], row=nextRow).value = values
        self.kb[wsName].cell(column=self.header['evidence']['units'], row=nextRow).value = units
        self.kb[wsName].cell(column=self.header['evidence']['experiment'], row=nextRow).value = experiment
        self.kb[wsName].cell(column=self.header['evidence']['database references'], row=nextRow).value = dbRefs
        self.kb[wsName].cell(column=self.header['evidence']['comments'], row=nextRow).value = comments

        entry ={'ws':'Evidence', 'id':eviId, 'row':nextRow}
        return entry

    def createExperiment(self, experimentDesign=None, measurmentTechnology=None, analysisType=None, species=None, geneticVariant=None, externalMedia=None, temperature=None, temperatureUnits='C', ph=None, dbRefs=None, refs=None, comments=None):
        wsName = 'Experiment'
        nextRow = self.kb[wsName].max_row+1
        expId   = 'EXP{}'.format(str(nextRow-1).zfill(4))

        self.kb[wsName].cell(column=self.header['experiment']['id'], row=nextRow).value = expId
        self.kb[wsName].cell(column=self.header['experiment']['experiment design'], row=nextRow).value = experimentDesign
        self.kb[wsName].cell(column=self.header['experiment']['measurment technology'], row=nextRow).value = measurmentTechnology
        self.kb[wsName].cell(column=self.header['experiment']['analysis type'], row=nextRow).value = analysisType
        self.kb[wsName].cell(column=self.header['experiment']['species'], row=nextRow).value = species
        self.kb[wsName].cell(column=self.header['experiment']['genetic variant'], row=nextRow).value = geneticVariant
        self.kb[wsName].cell(column=self.header['experiment']['external media'], row=nextRow).value = externalMedia
        self.kb[wsName].cell(column=self.header['experiment']['temperature'], row=nextRow).value = temperature
        self.kb[wsName].cell(column=self.header['experiment']['temperature units'], row=nextRow).value = temperatureUnits
        self.kb[wsName].cell(column=self.header['experiment']['ph'], row=nextRow).value = ph
        self.kb[wsName].cell(column=self.header['experiment']['database references'], row=nextRow).value = dbRefs
        self.kb[wsName].cell(column=self.header['experiment']['references'], row=nextRow).value = refs
        self.kb[wsName].cell(column=self.header['experiment']['comments'], row=nextRow).value = comments

        entry ={'ws':'Experiment', 'id':expId, 'row':nextRow}
        return entry

    def insertCell(self, insertTo, header, value, concatenate=True):
        wsName = insertTo['ws']
        rowIdx = insertTo['row']

        if  concatenate is True:
            self.concatenateId(self.kb[wsName].cell(column=self.header[wsName][header], row=rowIdx), value)
        else:
            self.kb[wsName].cell(column=self.header[wsName][header], row=rowIdx).value =value

    @staticmethod
    def concatenateId(cell, entry):
        if entry is None:
            return

        if cell.value is None:
            cell.value = '' +  entry
        else:
            cell.value = cell.value + ', ' +  entry

    def parseDatabaseReferences(self, cell):
        if cell.value is None:
            return None

        jsons = self.delinateJSONs(cell.value)

        if jsons is None:
            return None
        assert isinstance(jsons, list)
        assert len(jsons) > 0

        dbRefIds = ''
        for json in jsons:

            if json['source']=='URL':
                dbRefId = "{}:{}".format(json['source'].lower(), str(nextRow).zfill(4))
            else:
                dbRefId = "{}:{}".format(str(json['source']).lower(), str(json['xid']))

            dbRefId = dbRefId.replace('-','_')
            dbRefIds += dbRefId + ', '

        return dbRefIds[:-2]





    """ Old Translate methods """
    def translateChromosomeFeatures(self):
        """ Parses information from the core's 'Chromosome Features' sheet and inserts them to the appropiate field(s) in the KB structure. """

        coreChromFeats = self.core['Chromosome features']
        kbChromFeats   = self.kb['Chromosome features']
        kbEviColumn = 10

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

        for rowIdx in range(3, 500): # coreChromFeats.max_row+1):
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

            # Add evidence and experiment(s)
            eviDictStr = coreChromFeats.cell(column=12, row=rowIdx).value
            if eviDictStr is not None:

                eviDict = json.loads(eviDictStr)
                objectId = kbChromFeats.cell(column=1, row=rowIdx-1).value
                if rowIdx==3:   # Add experiment once - all chrom feature data from same experiment
                    expId = self.addExperiment(eviDict)

                eviId = self.addEvidence(
                            objectId = objectId,
                            property = 'intensity',
                            eviDict = eviDict,
                            expId = expId)

                self.concatenateId(kbChromFeats.cell(column=kbEviColumn, row=kbChromFeats.max_row), eviId)

            # NO 'Cross references' entry in core: kbChromFeats.cell(column=11, row=rowIdx-1).value = coreChromFeats.cell(column=11, row=rowIdx).value
            kbChromFeats.cell(column=12, row=rowIdx-1).value = coreChromFeats.cell(column=14, row=rowIdx).value # copy references
            kbChromFeats.cell(column=13, row=rowIdx-1).value = coreChromFeats.cell(column=13, row=rowIdx).value

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

        for rowIdx in range(3, coreTUs.max_row+1):
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
            kbTUs.cell(column=13, row=rowIdx-1).value = coreTUs.cell(column=13, row=rowIdx).value

    def translateGenes(self):
        """ Parses information from the core's Genes sheet and inserts them to the appropiate field(s) in the KB structure. """

        coreGenes = self.core['Genes']
        kbGenes   = self.kb['Genes']
        kbEviColumn = 13

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

        for rowIdx in range(3, coreGenes.max_row+1):
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

            # Add evidence and experiment(s)
            eviDictStr = coreGenes.cell(column=17, row=rowIdx).value

            if eviDictStr is not None and eviDictStr!='':
                eviDict = json.loads(eviDictStr)
                objectId = kbGenes.cell(column=1, row=rowIdx-1).value

                if rowIdx==3: # Add experiment once - all essentiality data from same experiment
                    expId = self.addExperiment(eviDict)

                eviId = self.addEvidence(
                    objectId = objectId,
                    property = 'is_essential',
                    eviDict = eviDict,
                    #kbEviColumn = kbEviColumn,
                    expId = expId)

                self.concatenateId(kbGenes.cell(column=kbEviColumn, row=kbGenes.max_row),
                                   eviId)

            dbRefsDicts = self.delinateJSONs(coreGenes.cell(column=5, row=rowIdx).value)
            kbGenes.cell(column=14, row=rowIdx-1).value = self.addDatabaseReference(dbRefsDicts, returnIds=True) # Add DB refs
            kbGenes.cell(column=15, row=rowIdx-1).value = coreGenes.cell(column=19, row=rowIdx).value # References
            kbGenes.cell(column=16, row=rowIdx-1).value = coreGenes.cell(column=18, row=rowIdx).value

    def translateReactions(self):
        """ Parses information from the core's 'Reactions' sheet and inserts them to the appropiate field(s) in the KB structure. """

        coreReactions = self.core['Reactions']
        kbReactions   = self.kb['Reactions']

        for rowIdx in range(3, coreReactions.max_row+1):
            print(rowIdx)
            """ missing
            modification
            """

            # Check if row is correct
            assert kbReactions.cell(column=1, row=rowIdx-1).value == coreReactions.cell(column=1, row=rowIdx).value

            #kbReactions.cell(column=1, row=rowIdx-1).value = coreReactions.cell(column=1, row=rowIdx).value # ID
            #kbReactions.cell(column=2, row=rowIdx-1).value = coreReactions.cell(column=2, row=rowIdx).value # Name
            #kbReactions.cell(column=3, row=rowIdx-1).value = coreReactions.cell(column=3, row=rowIdx).value # Synonyms

            # Copy types
            reactionType = coreReactions.cell(column=5, row=rowIdx).value
            if reactionType is None:
                kbReactions.cell(column=4, row=rowIdx-1).value = 'Uncategorized'
            else:
                kbReactions.cell(column=4, row=rowIdx-1).value = coreReactions.cell(column=5, row=rowIdx).value

    def translateReferences(self):
        """ Parses information from the core's 'References' sheet and inserts it to the appropiate field(s) in the new kb structure. """

        coreRefs = self.core['References']
        kbRefs   = self.kb['References']

        assert kbRefs.cell(column=1, row=1).value  == 'Id'
        assert kbRefs.cell(column=2, row=1).value  == 'Name'
        assert kbRefs.cell(column=3, row=1).value  == 'Type'
        assert kbRefs.cell(column=4, row=1).value  == 'Title'
        assert kbRefs.cell(column=5, row=1).value  == 'Authors'
        assert kbRefs.cell(column=6, row=1).value  == 'Journal'
        assert kbRefs.cell(column=7, row=1).value  == 'Volume'
        assert kbRefs.cell(column=8, row=1).value  == 'Issue'
        assert kbRefs.cell(column=9, row=1).value  == 'Pages'
        assert kbRefs.cell(column=10, row=1).value == 'Year'
        assert kbRefs.cell(column=11, row=1).value == 'Database references' #NO DB REFS IN CORE
        assert kbRefs.cell(column=12, row=1).value == 'Comments'

        for rowIdx in range(3, coreRefs.max_row+1):
            print(rowIdx)
            kbRefs.cell(column=1, row=rowIdx-1).value = coreRefs.cell(column=1, row=rowIdx).value # ID
            kbRefs.cell(column=2, row=rowIdx-1).value = coreRefs.cell(column=2, row=rowIdx).value # Name
            kbRefs.cell(column=3, row=rowIdx-1).value = coreRefs.cell(column=5, row=rowIdx).value # Type
            kbRefs.cell(column=4, row=rowIdx-1).value = coreRefs.cell(column=9, row=rowIdx).value # Title
            kbRefs.cell(column=5, row=rowIdx-1).value = coreRefs.cell(column=6, row=rowIdx).value # Authors
            kbRefs.cell(column=6, row=rowIdx-1).value = coreRefs.cell(column=10, row=rowIdx).value # Journal
            kbRefs.cell(column=7, row=rowIdx-1).value = coreRefs.cell(column=12, row=rowIdx).value # Volume
            kbRefs.cell(column=8, row=rowIdx-1).value = coreRefs.cell(column=13, row=rowIdx).value # Issue
            kbRefs.cell(column=9, row=rowIdx-1).value = coreRefs.cell(column=14, row=rowIdx).value # Pages
            kbRefs.cell(column=10, row=rowIdx-1).value = coreRefs.cell(column=8, row=rowIdx).value # Year

            dbRefsDicts = self.delinateJSONs(coreRefs.cell(column=4, row=rowIdx).value)
            kbRefs.cell(column=11, row=rowIdx-1).value = self.addDatabaseReference(dbRefsDicts, returnIds=True) # Add DB refs
            kbRefs.cell(column=12, row=rowIdx-1).value = coreRefs.cell(column=15, row=rowIdx).value

    """ Nested entries """
    def addDatabaseReferenceOLD(self, cell, returnIds=False):

        jsons = self.delinateJSONs(cell.value)

        if jsons is None:
            return None

        assert isinstance(jsons, list)
        assert len(jsons) > 0

        wsName = 'Database references'
        dbRefIds = ''

        for json in jsons:
            nextRow  = self.kb[wsName].max_row+1

            if json['source']=='URL':
                dbRefId = "{}:{}".format(json['source'].lower(), str(nextRow).zfill(4))
            else:
                dbRefId = "{}:{}".format(str(json['source']).lower(), str(json['xid']))

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

    def addSpeciesProperties(self, speciesId, structure=None, half_life=None, half_life_units=None, domains=None, evidence=None, refs=None):
        wsName = 'Species properties'
        nextRow = self.kb[wsName].max_row+1
        speciePropId = 'prop_{}'.format(speciesId)

        if structure==None and half_life==None and half_life_units=='min' and domains==None and evidence==None and refs==None:
            return None

        self.kb[wsName].cell(column=1, row=nextRow).value = speciePropId

        structureColumn = self.getColumnId(wbName='kb', sheetName=wsName, header = 'Structure')
        self.kb[wsName].cell(column=structureColumn, row=nextRow).value = structure

        hlColumn = self.getColumnId(wbName='kb', sheetName=wsName, header = 'Half life')
        self.kb[wsName].cell(column=hlColumn, row=nextRow).value = half_life

        hluColumn = self.getColumnId(wbName='kb', sheetName=wsName, header = 'Half life units')
        self.kb[wsName].cell(column=hluColumn, row=nextRow).value = half_life_units

        domainsColumn = self.getColumnId(wbName='kb', sheetName=wsName, header = 'Domains')
        self.kb[wsName].cell(column=domainsColumn, row=nextRow).value = domains

        eviColumn = self.getColumnId(wbName='kb', sheetName=wsName, header = 'Evidence')
        self.kb[wsName].cell(column=eviColumn, row=nextRow).value = evidence

        refsColumn = self.getColumnId(wbName='kb', sheetName=wsName, header = 'References')
        self.kb[wsName].cell(column=refsColumn, row=nextRow).value = refs

        return speciePropId

    def addConcentrationOLD(self, objectId, fieldStr):
        objectId = objectId.value
        fieldStr = fieldStr.value
        if fieldStr is None:
            return None

        wsName = 'Concentrations'
        nextRow = self.kb[wsName].max_row+1
        coreConc  = json.loads(fieldStr)
        #concId = 'CONC({}:{})'.format(objectId, coreConc['compartment'])
        concId = 'conc_{}_{}'.format(objectId, coreConc['compartment'])
        assert coreConc['compartment'] in ['c', 'e', 'm']

        self.kb[wsName].cell(column=1, row=nextRow).value = concId

        meanColumn = self.getColumnId(wbName='kb', sheetName='Concentrations', header = 'Mean')
        self.kb[wsName].cell(column=meanColumn, row=nextRow).value = coreConc['concentration']

        unitsColumn = self.getColumnId(wbName='kb', sheetName='Concentrations', header = 'Units')
        if 'evidence' in coreConc.keys():
            self.kb[wsName].cell(column=unitsColumn, row=nextRow).value = coreConc['evidence'][0]['units'] #check if all evidence has same unit

        speciesColumn = self.getColumnId(wbName='kb', sheetName='Concentrations', header = 'Species')
        self.kb[wsName].cell(column=speciesColumn, row=nextRow).value = '{}[{}]'.format(objectId, coreConc['compartment'])

        evidenceColumn = self.getColumnId(wbName='kb', sheetName='Concentrations', header = 'Evidence')
        if 'evidence' not in coreConc.keys():
            return concId

        for evidence in coreConc['evidence']:
            evidenceId = self.addEvidence(objectId, 'conc', evidence, evidenceColumn)
            self.concatenateId(self.kb[wsName].cell(column=evidenceColumn, row=nextRow), evidenceId)

        return concId

    def addEvidenceOLD(self, objectId, property, eviDict, expId = None):
        if eviDict is None:
            return None

        wsName = 'Evidence'
        nextRow = self.kb[wsName].max_row+1
        eviId   = self.get_eviId(nextRow, objectId, property)

        #if not self.idExists(wsName, eviId):
        self.kb[wsName].cell(column=1, row=nextRow).value = eviId
        self.kb[wsName].cell(column=2, row=nextRow).value = 'cell' # Need to pass in object ID
        self.kb[wsName].cell(column=3, row=nextRow).value = objectId
        self.kb[wsName].cell(column=4, row=nextRow).value = property

        if expId is None:
            self.kb[wsName].cell(column=9, row=nextRow).value = self.addExperiment(eviDict)
        else:
            self.kb[wsName].cell(column=9, row=nextRow).value = expId

        # Not all evidence objects in core contans the following keys, thus use try
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

        return eviId

    def addExperiment(self, eviDict):

        # References and species uniquely identify experiments in coreT
        refs = None
        specie = None

        try:
            eviDict = json.loads(eviDict)
        except:
            print('Experiment could nto be JSON parsed!')

        if eviDict['references'] is not None:
            refs = self.cleanJsonDump(json.dumps(eviDict['references']))
        if eviDict['species'] is not None:
            specie = eviDict['species']

        wsName = 'Experiment'
        nextRow = self.kb[wsName].max_row+1

        speciesColumn = self.getColumnId(wbName='kb', sheetName='Experiment', header = 'Species')
        refsColumn   = self.getColumnId(wbName='kb', sheetName='Experiment', header = 'References')

        for rowIdx in range(2,nextRow):
             if self.kb[wsName].cell(row=rowIdx, column=speciesColumn).value == specie and \
                self.kb[wsName].cell(row=rowIdx, column=refsColumn).value == refs:
                experimentId = self.kb[wsName].cell(row=rowIdx, column=1).value
                return experimentId

        experimentId = 'EXP_{}'.format(str(nextRow-1).zfill(4))
        try: self.kb[wsName].cell(column=1, row=nextRow).value  = experimentId
        except: pass
        try: self.kb[wsName].cell(column=5, row=nextRow).value  = eviDict['species']
        except: pass
        try: self.kb[wsName].cell(column=7, row=nextRow).value  = eviDict['media']
        except: pass
        try: self.kb[wsName].cell(column=8, row=nextRow).value  = eviDict['temperature']
        except: pass
        try: self.kb[wsName].cell(column=13, row=nextRow).value = self.cleanJsonDump(json.dumps(eviDict['references']))
        except: pass

        return experimentId

    def getColumnId(self, wbName, sheetName, header):

        assert isinstance(wbName, str)
        assert isinstance(sheetName, str)
        assert isinstance(header, str)

        if wbName == 'kb':
            wb = self.kb
            rowIdx = 1
        elif wbName == 'core':
            wb = self.core
            rowIdx = 2
        else:
            raise Exception('Unrecognised workbook name: \n\t {}'.format(wb))

        sheet = wb[sheetName]
        id=[]

        for columnIdx in range(1, sheet.max_column+1):
            #if sheet.cell(column = columnIdx, row = rowIdx).value is None:
            #    continue
            if header.lower() == sheet.cell(column = columnIdx, row = rowIdx).value.lower():
                id.append(columnIdx)

        if len(id)>1:
            print("Mulitple '{}' headers are matched in sheet {} of workbook {}!".format(
                header, sheetName, wbName))
            return id
        elif len(id)==0:
            print("Header '{}', was not found in sheet '{}' of workbook '{}'.".format(
                header, sheetName, wbName))
            return None
        elif len(id)==1:
            return id[0]
        else:
            raise Exception('Unexpected match number in getColumnId.')

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
    def cleanJsonDump(fieldStr):
        fieldStr = fieldStr.replace('[', '')
        fieldStr = fieldStr.replace(']', '')
        fieldStr = fieldStr.replace('"', '')
        return fieldStr
