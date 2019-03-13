""" Core schema to represent a knowledge base to build models

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Jonathan Karr <karr@mssm.edu>
:Author: Bilal Shaikh  <bilal.shaikh@columbia.edu>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2018-02-07
:Copyright: 2018, Karr Lab
:License: MIT
"""

from natsort import natsorted, ns
from math import ceil, floor, exp, log, log10, isnan
from pyfaidx import Fasta
from wc_utils.util import chem
from wc_utils.util.list import det_dedupe
from wc_utils.util.units import unit_registry
import abc
import pdb
import Bio.Alphabet
import Bio.Seq
import enum
import math
import obj_model.abstract
import obj_model
import obj_model.units
import openbabel
import pkg_resources
import re
import six
import token
from obj_model import (BooleanAttribute, EnumAttribute, FloatAttribute, IntegerAttribute, PositiveIntegerAttribute,
                       RegexAttribute, SlugAttribute, StringAttribute, LongStringAttribute, UrlAttribute,
                       OneToOneAttribute, ManyToOneAttribute, ManyToManyAttribute,
                       InvalidModel, InvalidObject, InvalidAttribute, TabularOrientation)
from obj_model.expression import (ExpressionOneToOneAttribute, ExpressionManyToOneAttribute,
                                  ExpressionStaticTermMeta, ExpressionDynamicTermMeta,
                                  ExpressionExpressionTermMeta, Expression,
                                  ParsedExpression, ParsedExpressionError)
from wc_utils.util.enumerate import CaseInsensitiveEnum
from wc_utils.util.types import get_subclasses

with open(pkg_resources.resource_filename('wc_kb', 'VERSION'), 'r') as file:
    wc_kb_version = file.read().strip()

#####################
#####################
# Enumeration classes

PolymerStrand = enum.Enum(value='PolymerStrand', names=[
    ('positive', 1),
    ('+', 1),
    ('negative', -1),
    ('-', -1), ])


class RnaType(enum.Enum):
    """ Type of RNA """
    mRna = 0
    rRna = 1
    sRna = 2
    tRna = 3
    mixed = 4


class GeneType(enum.Enum):
    """ Type of gene """
    mRna = 0
    rRna = 1
    sRna = 2
    tRna = 3
    asRna = 4
    pseudogene = 5


class DirectionType(enum.Enum):
    """ Type of direction """
    forward = 1
    backward = -1


class CogCategoryType(enum.Enum):
    """ Type of Clusters of Orthologous Groups (COGs)
        List of categories obtained from: http://ecoliwiki.net/colipedia/index.php/Clusters_of_Orthologous_Groups_(COGs)
    """

    RNA_processing = 0
    chromatin_structure_dynamics = 1
    energy_production_conversion = 2
    cell_cycle_control_mitosis = 3
    amino_acid_metabolism_transport = 4
    nucleotide_metabolism_transport = 5
    carbohydrate_metabolism_transport = 6
    coenzyme_metabolis = 7
    lipid_metabolism = 8
    tranlsation = 9
    transcription = 10
    replication_repair = 11
    cell_wall_membrane_envelop_biogenesis = 12
    cell_motility = 13
    post_translational_modification_protein_turnover_chaperone_functions = 14
    inorganic_ion_transport_metabolism = 15
    secondary_structure = 16
    signal_transduction = 17
    intracellular_trafficing_secretion = 18
    nuclear_structure = 19
    cytoskeleton = 20
    general_functional_prediction = 21
    unknown = 22


class ComplexType(enum.Enum):
    """ Type of complex """
    tRnaSynthClassII = 0
    FattyAcylAcp = 1


class ComplexFormationType(enum.Enum):
    """ Type of complex formation"""
    process_ChromosomeCondensation = 0
    process_FtsZPolymerization = 1
    process_MacromolecularComplexation = 2
    process_Metabolism = 3
    process_ProteinModification = 4
    process_Replication = 5
    process_ReplicationInitiation = 6
    process_RibosomeAssembly = 7
    process_Transcription = 8
    process_Translation = 9


class ChromosomeFeatureType(enum.Enum):
    """ Type of complex formation"""
    LongStructuralRegion = 0
    DnaMethylation = 1
    GeneWizPrediction = 2
    DnaBindingSite = 3
    DnaBindingSite_Lon = 4
    FunctionalDnaABox = 5
    DnaABox = 6


class MetaboliteSpeciesTypeType(enum.Enum):
    """ Types of metabolites """

    vitamin = 0
    amino_acid = 1
    dipeptide = 2
    nucleobase = 3
    modified_nucleobase = 4
    carbohydrate_sugar = 5
    carbohydrate_sugar_phosphate = 6
    ribonucleotide_monophosphate = 7
    ribonucleotide_biphosphate = 8
    ribonucleotide_triphosphate = 9
    carboxy_acid=10
    misc=12
    unknown = 11

class SignalSequenceType(enum.Enum):
    """ Types of signal sequences """
    secretory = 0
    lipoprotein = 1

class ProteinType(enum.Enum):
    """ Types of signal sequences """
    TrnaSynthClassI = 0
    TrnaSynthClassIB = 1
    TrnaSynthClassII = 2
    uncategorized = 3

class DnaBindingType(enum.Enum):
    """ Types of DNA binding """
    ssDNA = 0
    dsDNA = 1

class ReactionType(enum.Enum):
    """ Types of DNA binding """

    Uncategorized = 0
    DnaDamageBaseAlkylationReaction = 1
    DnaDamageBaseEthylationReaction =2
    DnaDamageRadiationInducedBaseOxidation =3
    DnaDamageBaseMethylationReaction =4
    DnaDamageBaseAminationReaction =5
    DnaDamageUvBPhotodimerization =6
    DnaDamagePhotooxidationReaction =7
    DnaDamageStrandBreakReaction=8
    DnaDamageBaseGlucosylTransferReaction=9
    DnaDamageSpontaneousBaseDeaminationReaction=10
    DnaDamageSpontaneousBaseLossReaction=11
    DnaDamageBaseReductionReaction=12
    DnaRepairBaseExcisionRepairReaction=13
    DnaRepairDnaLigationReaction=14
    DnaRepairDnaPolymerizationReaction=15
    DnaRepairDnaRestrictionModificationReaction=16
    DnaRepairDnaCleavageReaction=17
    DnaRepairHomologousRecombinationReaction=18
    DnaRepairBaseExcisionRepairBaseExcisionReaction=19
    DnaRepairNucleotideExcisionRepairReaction=20
    ChemicalReaction=21
    TransportReaction=22
    ModifiedBaseTransportReaction=23
    IonTransportReaction=24
    OxidationInactivatingProteinModificationReaction=25
    ProteinModificationAdductionReaction=26
    ProteinModificationLigationReaction=27
    TrnaTransferReaction=28
    TrnaAminoacylationReaction=29
    GlycationInactivatingProteinModificationReaction=30
    PhosphorylationInactivatingProteinModificationReaction=31
    DephosphorylationActivatingProteinModificationReaction=32

class ReferenceType(enum.Enum):
    """ Types of references """

    article = 0
    preprint = 1
    supplementary_material = 2
    book = 3
    thesis = 4
    misc = 5


#####################
#####################
# Attributes

class SubunitAttribute(ManyToManyAttribute):
    """ Subunits """

    def __init__(self, related_name='', verbose_name='', verbose_related_name='', help=''):
        """
        Args:
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            help (:obj:`str`, optional): help message
        """

        super(SubunitAttribute, self).__init__('SpeciesTypeCoefficient',
                                               related_name=related_name,
                                               verbose_name=verbose_name,
                                               verbose_related_name=verbose_related_name,
                                               help=help)

    def serialize(self, subunits, encoded=None):
        """ Serialize related object

        Args:
            subunits (:obj:`list` of :obj:`SpeciesTypeCoefficient`): Python representation of subunits
            encoded (:obj:`dict`, optional): dictionary of objects that have already been encoded

        Returns:
            :obj:`str`: simple Python representation
        """
        if not subunits:
            return ''

        subunits = natsorted(subunits, lambda unit: (
            unit.species_type.id), alg=ns.IGNORECASE)

        lhs = []
        for unit in subunits:
            lhs.append(unit.serialize())

        return '{}'.format(' + '.join(lhs))

    def deserialize(self, value, objects, decoded=None):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            decoded (:obj:`dict`, optional): dictionary of objects that have already been decoded

        Returns:
            :obj:`tuple` of `object`, `InvalidAttribute` or `None`: tuple of cleaned value and cleaning error
        """
        return SpeciesTypeCoefficient.deserialize(self, value, objects)


class OneToOneSpeciesAttribute(OneToOneAttribute):
    """ Species attribute """

    def __init__(self, related_name='', verbose_name='', verbose_related_name='', help=''):
        """
        Args:
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            help (:obj:`str`, optional): help message
        """
        super(OneToOneSpeciesAttribute, self).__init__('Species',
                                                       related_name=related_name, min_related=1, min_related_rev=0,
                                                       verbose_name=verbose_name, verbose_related_name=verbose_related_name, help=help)

    def serialize(self, value, encoded=None):
        """ Serialize related object
        Args:
            value (:obj:`Model`): Python representation
            encoded (:obj:`dict`, optional): dictionary of objects that have already been encoded
        Returns:
            :obj:`str`: simple Python representation
        """
        return value.serialize()

    def deserialize(self, value, objects, decoded=None):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            decoded (:obj:`dict`, optional): dictionary of objects that have already been decoded

        Returns:
            :obj:`tuple` of :obj:`list` of :obj:`Species`, :obj:`InvalidAttribute` or :obj:`None`: :obj:`tuple` of cleaned value
                and cleaning error
        """
        return Species.deserialize(self, value, objects)


class DatabaseReferenceAttribute(ManyToManyAttribute):
    """ Database reference attribute """

    def __init__(self, related_name='', verbose_name='', verbose_related_name='', help=''):
        """
        Args:
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            help (:obj:`str`, optional): help message
        """
        super(DatabaseReferenceAttribute, self).__init__(DatabaseReference,
                                                         related_name=related_name, min_related=0, min_related_rev=0,
                                                         verbose_name=verbose_name, verbose_related_name=verbose_related_name, help=help)

    def serialize(self, database_references, encoded=None):
        """ Serialize related object
        Args:
            database_references (:obj:`list` of :obj:`Model`): a list of instances of DatabaseReference Python representation
            encoded (:obj:`dict`, optional): dictionary of objects that have already been encoded
        Returns:
            :obj:`str`: simple Python representation
        """
        if not database_references:
            return ''

        return ', '.join(obj_model.serialize() for obj_model in database_references)

    def deserialize(self, value, objects, decoded=None):
        """ Deserialize value
        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            decoded (:obj:`dict`, optional): dictionary of objects that have already been decoded
        Returns:
            :obj:`tuple` of :obj:`list` of :obj:`DatabaseReference`, :obj:`InvalidAttribute` or :obj:`None`: :obj:`tuple` of cleaned value
                and cleaning error
        """
        if not value:
            return ([], None)

        pattern = r'([a-z][a-z0-9_\-]*)\:([a-z0-9_\-]*)'
        #pattern = r'DBREF[a-z0-9_\-]*\([a-zA-Z0-9]*:[a-zA-Z0-9]*\)'

        if not re.match(pattern, value, flags=re.I):
            return (None, InvalidAttribute(self, ['Incorrectly formatted list of database references: {}'.format(value)]))

        objs = []
        for pat_match in re.findall(pattern, value, flags=re.I):
            match = re.match(pattern, value, flags=re.I)
            database_name = match.group(1)
            data_id = match.group(2)
            if self.related_class not in objects:
                objects[self.related_class] = {}
            serialized_value = value
            if serialized_value in objects[self.related_class]:
                obj = objects[self.related_class][serialized_value]
            else:
                obj = self.related_class(database=database_name, id=data_id)
                objects[self.related_class][serialized_value] = obj
            objs.append(obj)
        return (objs, None)


class ReactionParticipantAttribute(ManyToManyAttribute):
    """ Reaction participants """

    def __init__(self, related_name='', verbose_name='', verbose_related_name='', help=''):
        """
        Args:
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            help (:obj:`str`, optional): help message
        """
        super(ReactionParticipantAttribute, self).__init__('SpeciesCoefficient', related_name=related_name,
                                                           verbose_name=verbose_name,
                                                           verbose_related_name=verbose_related_name,
                                                           help=help)

    def serialize(self, participants, encoded=None):
        """ Serialize related object

        Args:
            participants (:obj:`list` of :obj:`SpeciesCoefficient`): Python representation of reaction participants
            encoded (:obj:`dict`, optional): dictionary of objects that have already been encoded

        Returns:
            :obj:`str`: simple Python representation
        """
        if not participants:
            return ''

        comps = set([part.species.compartment for part in participants])
        if len(comps) == 1:
            global_comp = comps.pop()
        else:
            global_comp = None

        if global_comp:
            participants = natsorted(
                participants, lambda part: part.species.species_type.id, alg=ns.IGNORECASE)
        else:
            participants = natsorted(participants, lambda part: (
                part.species.species_type.id, part.species.compartment.id), alg=ns.IGNORECASE)

        lhs = []
        rhs = []
        for part in participants:
            if part.coefficient < 0:
                lhs.append(part.serialize(
                    show_compartment=global_comp is None, show_coefficient_sign=False))
            elif part.coefficient > 0:
                rhs.append(part.serialize(
                    show_compartment=global_comp is None, show_coefficient_sign=False))

        if global_comp:
            return '[{}]: {} ==> {}'.format(global_comp.get_primary_attribute(), ' + '.join(lhs), ' + '.join(rhs))
        else:
            return '{} ==> {}'.format(' + '.join(lhs), ' + '.join(rhs))

    def deserialize(self, value, objects, decoded=None):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            decoded (:obj:`dict`, optional): dictionary of objects that have already been decoded

        Returns:
            :obj:`tuple` of `list` of `SpeciesCoefficient`, `InvalidAttribute` or `None`: tuple of cleaned value
                and cleaning error
        """
        errors = []

        id = r'[a-z][a-z0-9_]*'
        stoch = r'\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\)'
        gbl_part = r'({} )*({})'.format(stoch, id)
        lcl_part = r'({} )*({}\[{}\])'.format(stoch, id, id)
        gbl_side = r'{}( \+ {})*'.format(gbl_part, gbl_part)
        lcl_side = r'{}( \+ {})*'.format(lcl_part, lcl_part)
        gbl_pattern = r'^\[({})\]: ({}) ==> ({})$'.format(
            id, gbl_side, gbl_side)
        lcl_pattern = r'^({}) ==> ({})$'.format(lcl_side, lcl_side)

        global_match = re.match(gbl_pattern, value, flags=re.I)
        local_match = re.match(lcl_pattern, value, flags=re.I)

        if global_match:
            if global_match.group(1) in objects[Compartment]:
                global_comp = objects[Compartment][global_match.group(1)]
            else:
                global_comp = None
                errors.append('Undefined compartment "{}"'.format(
                    global_match.group(1)))
            lhs = global_match.group(2)
            rhs = global_match.group(14)

        elif local_match:
            global_comp = None
            lhs = local_match.group(1)
            rhs = local_match.group(13)

        else:
            return (None, InvalidAttribute(self, ['Incorrectly formatted participants: {}'.format(value)]))

        lhs_parts, lhs_errors = self.deserialize_side(
            -1., lhs, objects, global_comp)
        rhs_parts, rhs_errors = self.deserialize_side(
            1., rhs, objects, global_comp)

        parts = lhs_parts + rhs_parts
        errors.extend(lhs_errors)
        errors.extend(rhs_errors)

        if errors:
            return (None, InvalidAttribute(self, errors))
        return (parts, None)

    def deserialize_side(self, direction, value, objects, global_comp):
        """ Deserialize the LHS or RHS of a reaction equation
        Args:
            direction (:obj:`float`): -1. indicates LHS, +1. indicates RHS
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            global_comp (:obj:`Compartment`): global compartment of the reaction

        Returns:
            :obj:`tuple`:
                * :obj:`list` of :obj:`SpeciesCoefficient`: list of species coefficients
                * :obj:`list` of :obj:`Exception`: list of errors
        """
        parts = []
        errors = []

        for part in re.findall(r'(\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\) )*([a-z][a-z0-9_]*)(\[([a-z][a-z0-9_]*)\])*', value, flags=re.I):
            part_errors = []

            species_type = None
            for species_type_cls in get_subclasses(SpeciesType):
                if species_type_cls in objects and part[4] in objects[species_type_cls]:
                    species_type = objects[species_type_cls][part[4]]
                    break
            if not species_type:
                part_errors.append(
                    'Undefined species type "{}"'.format(part[4]))

            if global_comp:
                compartment = global_comp
            elif part[6] in objects[Compartment]:
                compartment = objects[Compartment][part[6]]
            else:
                part_errors.append(
                    'Undefined compartment "{}"'.format(part[6]))

            coefficient = direction * float(part[1] or 1.)

            if part_errors:
                errors += part_errors
            else:
                spec_primary_attribute = Species.gen_id(species_type.get_primary_attribute(),
                                                        compartment.get_primary_attribute())
                species, error = Species.deserialize(
                    self, spec_primary_attribute, objects)
                if error:
                    raise ValueError('Invalid species "{}"'.format(
                        spec_primary_attribute))
                    # pragma: no cover # unreachable due to error checking above

                if coefficient != 0:
                    if SpeciesCoefficient not in objects:
                        objects[SpeciesCoefficient] = {}
                    serialized_value = SpeciesCoefficient._serialize(
                        species, coefficient)
                    if serialized_value in objects[SpeciesCoefficient]:
                        rxn_part = objects[SpeciesCoefficient][serialized_value]
                    else:
                        rxn_part = SpeciesCoefficient(
                            species=species, coefficient=coefficient)
                        objects[SpeciesCoefficient][serialized_value] = rxn_part
                    parts.append(rxn_part)

        return (parts, errors)


#####################
#####################
# Base classes

class KnowledgeBaseObject(obj_model.Model):
    """ Knowledge of a biological entity

    Attributes:
        id (:obj:`str`): identifier
        name (:obj:`str`): name
        synonyms (:obj:`str`): synonyms
        comments (:obj:`str`): comments
    """

    id = obj_model.SlugAttribute(primary=True, unique=True)
    name = obj_model.StringAttribute()
    synonyms = obj_model.StringAttribute()
    comments = obj_model.LongStringAttribute()


class KnowledgeBase(KnowledgeBaseObject):
    """ A knowledge base

    Attributes:
        version (:obj:`str`): version
        translation_table (:obj:`int`): translation table
        version (:obj:`str`): version of the KB
        url (:obj:`str`): url of the KB Git repository
        branch (:obj:`str`): branch of the KB Git repository
        revision (:obj:`str`): revision of the KB Git repository
        wc_kb_version (:obj:`str`): version of ``wc_kb``

    Related attributes:
        cell (:obj:`Cell`): cell
    """
    translation_table = obj_model.IntegerAttribute()
    version = RegexAttribute(
        min_length=1, pattern=r'^[0-9]+\.[0-9+]\.[0-9]+', flags=re.I)
    url = obj_model.StringAttribute(verbose_name='URL')
    branch = obj_model.StringAttribute()
    revision = obj_model.StringAttribute()
    wc_kb_version = RegexAttribute(min_length=1, pattern=r'^[0-9]+\.[0-9+]\.[0-9]+', flags=re.I,
                                   default=wc_kb_version, verbose_name='wc_kb version')

    class Meta(obj_model.Model.Meta):
        verbose_name = 'KB'
        attribute_order = ('id', 'name', 'translation_table', 'version',
                           'url', 'branch', 'revision', 'wc_kb_version', 'comments')
        tabular_orientation = obj_model.TabularOrientation.column


class Cell(KnowledgeBaseObject):
    """ Knowledge of a cell

    Attributes:
        knowledge_base (:obj:`KnowledgeBase`): knowledge base
        taxon (:obj:`int`): NCBI taxon identifier

    Related attributes:
        references (:obj:`list` of :obj:`Reference`): references
        compartments (:obj:`list` of :obj:`Compartment`): compartments
        species_types (:obj:`list` of :obj:`SpeciesType`): species types
        concentrations (:obj:`list` of :obj:`Concentration`): concentrations
        observables (:obj:`list` or :obj:`Observable`) : observables
        loci (:obj:`list` of :obj:`PolymerLocus`): locus
        reactions (:obj:`list` of :obj:`Reaction`): reactions
    """
    knowledge_base = obj_model.OneToOneAttribute(
        KnowledgeBase, related_name='cell')
    taxon = obj_model.IntegerAttribute()

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'taxon', 'comments')
        tabular_orientation = obj_model.TabularOrientation.column


class DatabaseReference(obj_model.Model):
    """ Reference to an entity in an external database

    Attributes:
        database (:obj:`str`): name of the external database
        id (:obj:`str`): identifier within the database

    Related attributes:
        compartments (:obj:`list` of :obj:`Compartment`): compartments
        species_types (:obj:`list` of :obj:`SpeciesType`): species_types
        concentrations (:obj:`list` of :obj:`Concentration`): concentrations
        loci (:obj:`list` of :obj:`PolymerLocus`): loci
        properties (:obj:`list` of :obj:`Property`): properties
        reactions (:obj:`list` of :obj:`Reaction`): reactions
        rate_laws (:obj:`list` of :obj:`RateLaw`): rate_laws
        observables (:obj:`list` of :obj:`Observable`): observables
    """
    database = obj_model.StringAttribute()
    id = obj_model.StringAttribute()
    #id = obj_model.SlugAttribute(primary=True, unique=True)
    entry_id = obj_model.StringAttribute()
    comments = obj_model.LongStringAttribute()

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'database', 'entry_id', 'comments')
        #tabular_orientation = TabularOrientation.inline
        #unique_together = (('database', 'id'), )
        #ordering = ('database', 'id')

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: value of primary attribute
        """
        return '{}:{}'.format(self.database, self.id)


class Reference(obj_model.Model):
    """ Reference to the literature

    Attributes:
        id (:obj:`str`): identifier
        standard_id (:obj:`str`): standard identifier such as DOI or PubMed ID
        cell (:obj:`Cell`): cell

    Related attributes:
        compartments (:obj:`list` of :obj:`Compartment`): compartments
        species_types (:obj:`list` of :obj:`SpeciesType`): species_types
        concentrations (:obj:`list` of :obj:`Concentration`): concentrations
        loci (:obj:`list` of :obj:`PolymerLocus`): loci
        properties (:obj:`list` of :obj:`Property`): properties
        reactions (:obj:`list` of :obj:`Reaction`): reactions
        rate_laws (:obj:`list` of :obj:`RateLaw`): rate_laws
        observables (:obj:`list` of :obj:`Observable`): observables
    """

    id = obj_model.SlugAttribute(primary=True, unique=True)
    name = obj_model.StringAttribute()
    type = obj_model.EnumAttribute(ReferenceType)
    authors = obj_model.LongStringAttribute()
    title = obj_model.LongStringAttribute()
    volume = obj_model.StringAttribute()
    issue = obj_model.StringAttribute()
    journal = obj_model.StringAttribute()
    pages = obj_model.StringAttribute()
    year = obj_model.IntegerAttribute()
    cell = obj_model.ManyToOneAttribute(Cell, related_name='references')
    database_references = DatabaseReferenceAttribute(related_name='references')
    comments = obj_model.LongStringAttribute()

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'type', 'title', 'authors', 'journal', 'volume', 'issue', 'pages', 'year', 'database_references', 'comments')


class Compartment(KnowledgeBaseObject):
    """ Knowledge of a subcellular compartment

    Attributes:
        cell (:obj:`Cell`): cell
        volumetric_fraction (:obj:`float`): average volumetric fraction relative to the cell volume
        references (:obj:`list` of :obj:`Reference`): references
        database_references (:obj:`list` of :obj:`DatabaseReference`): database references

    Related attributes:
        reaction_participants (:obj:`list` of :obj:`ReactionParticipant`): reaction participants
    """
    cell = obj_model.ManyToOneAttribute(Cell, related_name='compartments')
    volumetric_fraction = obj_model.FloatAttribute(min=0., max=1.)
    references = obj_model.ManyToManyAttribute(Reference, related_name='compartments')
    database_references = DatabaseReferenceAttribute(related_name='compartments')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'volumetric_fraction', 'database_references', 'references', 'comments')


class SpeciesType(six.with_metaclass(obj_model.abstract.AbstractModelMeta, KnowledgeBaseObject)):
    """ Knowledge of a molecular species

    Attributes:
        cell (:obj:`Cell`): cell
        half_life  (:obj:`float`): half life (s)
        references (:obj:`list` of :obj:`Reference`): references
        database_references (:obj:`list` of :obj:`DatabaseReference`): database references

    Related attributes:
        reaction_participants (:obj:`list` of :obj:`ReactionParticipant`): reaction participants
    """

    cell = obj_model.ManyToOneAttribute(Cell, related_name='species_types')
    half_life = obj_model.FloatAttribute(min=0)
    references = obj_model.ManyToManyAttribute(Reference, related_name='species_types')
    database_references = DatabaseReferenceAttribute(related_name='species_types')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'half_life', 'comments', 'references', 'database_references')

    @abc.abstractmethod
    def get_empirical_formula(self):
        """ Get the empirical formula

        Returns:
            :obj:`chem.EmpiricalFormula`: empirical formula
        """
        pass  # pragma: no cover

    @abc.abstractmethod
    def get_charge(self):
        """ Get the charge

        Returns:
            :obj:`int`: charge
        """
        pass  # pragma: no cover

    @abc.abstractmethod
    def get_mol_wt(self):
        """ Get the molecular weight

        Returns:
            :obj:`float`: molecular weight
        """
        pass  # pragma: no cover


class Species(obj_model.Model):
    """ Species (tuple of species type, compartment)

    Attributes:
        species_type (:obj:`SpeciesType`): species type
        compartment (:obj:`Compartment`): compartment

    Related attributes:
        concentration (:obj:`Concentration`): concentration
        species_coefficients (:obj:`list` of :obj:`SpeciesCoefficient`): participations in reactions
        rate_law_expressions (:obj:`list` of :obj:`RateLawExpression`): participations in the evaluation of rates
        observable_expressions (:obj:`list` of :obj:`ObservableExpression`): participations in observables
    """
    species_type = ManyToOneAttribute(
        SpeciesType, related_name='species', min_related=1)
    compartment = ManyToOneAttribute(
        Compartment, related_name='species', min_related=1)

    class Meta(obj_model.Model.Meta):
        attribute_order = ('species_type', 'compartment')
        frozen_columns = 1
        tabular_orientation = TabularOrientation.inline
        unique_together = (('species_type', 'compartment', ), )
        ordering = ('species_type', 'compartment')
        expression_term_token_pattern = (token.NAME, token.LSQB, token.NAME, token.RSQB)

    @staticmethod
    def gen_id(species_type, compartment):
        """ Generate a Species' primary identifier

        Args:
            species_type (:obj:`object`): a `SpeciesType`, or its id
            compartment (:obj:`object`): a `Compartment`, or its id

        Returns:
            :obj:`str`: canonical identifier for a specie in a compartment, 'species_type_id[compartment_id]'
        """
        if isinstance(species_type, SpeciesType):
            species_type_id = species_type.get_primary_attribute()
        elif isinstance(species_type, six.string_types):
            species_type_id = species_type
        else:
            raise ValueError(
                "gen_id: incorrect species type: {}".format(species_type))

        if isinstance(compartment, Compartment):
            compartment_id = compartment.get_primary_attribute()
        elif isinstance(compartment, six.string_types):
            compartment_id = compartment
        else:
            raise ValueError(
                "gen_id: incorrect compartment type: {}".format(compartment))

        return '{}[{}]'.format(species_type_id, compartment_id)

    def id(self):
        """ Provide a Species' primary identifier

        Returns:
            :obj:`str`: canonical identifier for a specie in a compartment, 'specie_id[compartment_id]'
        """
        return self.serialize()

    def serialize(self):
        """ Provide a Species' primary identifier

        Returns:
            :obj:`str`: canonical identifier for a specie in a compartment, 'specie_id[compartment_id]'
        """
        return self.gen_id(self.species_type, self.compartment)

    @classmethod
    def deserialize(cls, attribute, value, objects):
        """ Deserialize value

        Args:
            attribute (:obj:`Attribute`): attribute
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model

        Returns:
            :obj:`tuple` of `object`, `InvalidAttribute` or `None`: tuple of cleaned value and cleaning error
        """
        if cls in objects and value in objects[cls]:
            return (objects[cls][value], None)

        match = re.match(
            r'^([a-z][a-z0-9_]*)\[([a-z][a-z0-9_]*)\]$', value, flags=re.I)
        if match:
            errors = []

            species_type = None
            for species_type_cls in get_subclasses(SpeciesType):
                if species_type_cls in objects and match.group(1) in objects[species_type_cls]:
                    species_type = objects[species_type_cls][match.group(1)]
                    break
            if not species_type:
                errors.append(
                    'Species type "{}" is not defined'.format(match.group(1)))

            if Compartment in objects and match.group(2) in objects[Compartment]:
                compartment = objects[Compartment][match.group(2)]
            else:
                errors.append(
                    'Compartment "{}" is not defined'.format(match.group(2)))

            if errors:
                return (None, InvalidAttribute(attribute, errors))
            else:
                obj = cls(species_type=species_type, compartment=compartment)
                if cls not in objects:
                    objects[cls] = {}
                objects[cls][obj.serialize()] = obj
                return (obj, None)

        return (None, InvalidAttribute(attribute, ['Invalid species']))


class Concentration(KnowledgeBaseObject):
    """ Species concentration

    Attributes:
        cell (:obj:`Cell`): cell
        species (:obj:`Species`): species
        value (:obj:`float`): value
        units (:obj:`unit_registry.Unit`): units; default units is 'M'
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
        database_references (:obj:`list` of :obj:`DatabaseReference`): database references
    """
    cell = obj_model.ManyToOneAttribute(Cell, related_name='concentrations')
    species = OneToOneSpeciesAttribute(related_name='concentration')
    medium = obj_model.StringAttribute()
    mean = FloatAttribute(min=0)
    std_dev = FloatAttribute(min=0)
    database_references = DatabaseReferenceAttribute(related_name='concentrations')
    references = ManyToManyAttribute(Reference, related_name='concentrations')
    units = obj_model.units.UnitAttribute(unit_registry,
                          choices=(
                              unit_registry.parse_units('molecule'),
                              unit_registry.parse_units('mM'),
                              unit_registry.parse_units('uM'),
                              unit_registry.parse_units('nM'),
                              unit_registry.parse_units('pM'),
                              unit_registry.parse_units('fM'),
                              unit_registry.parse_units('aM'),
                          ),
                          default=unit_registry.parse_units('M'))

    evidence = obj_model.OneToManyAttribute('Evidence', related_name='concentrations')
    database_references = DatabaseReferenceAttribute(related_name='concentrations')
    references = ManyToManyAttribute(Reference, related_name='concentrations')
    comments = LongStringAttribute()


    class Meta(obj_model.Model.Meta):
        unique_together = (('species', 'medium' ), )
        attribute_order = ('id', 'name', 'species', 'medium', 'mean', 'std_dev', 'units', 'evidence', 'database_references', 'references', 'comments')
        frozen_columns = 1
        ordering = ('species',)

    def serialize(self):
        """ Generate string representation
        Returns:
            :obj:`str`: value of primary attribute
        """
        return self.species.serialize()


class SpeciesTypeCoefficient(obj_model.Model):
    """ A tuple of a species type and a coefficient

    Attributes:
        species_type (:obj:`SpeciesType`): species_type
        coefficient (:obj:`float`): coefficient

    Related attributes:
        complex (:obj:`ComplexSpeciesType`): complex
    """

    species_type = ManyToOneAttribute(SpeciesType, related_name='species_type_coefficients')
    coefficient = FloatAttribute(min=0., nan=False)

    class Meta(obj_model.Model.Meta):
        attribute_order = ('species_type', 'coefficient')
        frozen_columns = 1
        tabular_orientation = TabularOrientation.inline
        ordering = ('species_type',)

    def serialize(self):
        """ Serialize related object

        Returns:
            :obj:`str`: string representation of a species type and a coefficient
        """
        return self._serialize(self.species_type, self.coefficient)

    @staticmethod
    def _serialize(species_type, coefficient):
        """ Serialize values

        Args:
            species_type (:obj:`SpeciesType`): species_type
            coefficient (:obj:`float`): coefficient

        Returns:
            :obj:`str`: string representation of a species type and a coefficient
        """
        coefficient = float(coefficient)

        if coefficient == 1:
            coefficient_str = ''
        elif coefficient % 1 == 0 and abs(coefficient) < 1000:
            coefficient_str = '({:.0f}) '.format(coefficient)
        else:
            coefficient_str = '({:e}) '.format(coefficient)

        return '{}{}'.format(coefficient_str, species_type.get_primary_attribute())

    @classmethod
    def deserialize(cls, attribute, value, objects):
        """ Deserialize value

        Args:
            attribute (:obj:`Attribute`): attribute
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model

        Returns:
            :obj:`tuple` of `list` of `SpeciesTypeCoefficient`, `InvalidAttribute` or `None`: tuple of cleaned value
                and cleaning error
        """
        parts = []
        errors = []
        id = r'[a-z][a-z0-9_]*'
        stoch = r'\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\)'
        gbl_part = r'({} )*({})'.format(stoch, id)
        gbl_side = r'{}( \+ {})*'.format(gbl_part, gbl_part)
        gbl_pattern = r'^({})$'.format(gbl_side)

        global_match = re.match(gbl_pattern, value, flags=re.I)

        if global_match:
            subunits_str = global_match.group(1)
        else:
            attr = cls.Meta.attributes['species_type']
            return (None, InvalidAttribute(attr, ['Incorrectly formatted participants: {}'.format(value)]))

        for part in re.findall(r'(\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\) )*([a-z][a-z0-9_]*)',
                               subunits_str, flags=re.I):

            species_type = None
            for species_type_cls in get_subclasses(SpeciesType):
                if species_type_cls in objects and part[4] in objects[species_type_cls]:
                    species_type = objects[species_type_cls][part[4]]
                    break

            if not species_type:
                errors.append('Undefined species type "{}"'.format(part[4]))

            coefficient = float(part[1] or 1.)

            if not errors:

                if coefficient != 0:
                    if cls not in objects:
                        objects[cls] = {}
                    serialized_value = cls._serialize(species_type, coefficient)
                    if serialized_value in objects[cls]:
                        subunit_part = objects[cls][serialized_value]
                    else:
                        subunit_part = cls(species_type=species_type, coefficient=coefficient)
                        objects[cls][serialized_value] = subunit_part
                    parts.append(subunit_part)

        if errors:
            return (None, InvalidAttribute(cls, errors))
        return (parts, None)


class SpeciesCoefficient(obj_model.Model):
    """ A tuple of a species and a coefficient

    Attributes:
        species (:obj:`Species`): species
        coefficient (:obj:`float`): coefficient

    Related attributes:
        reaction (:obj:`Reaction`): reaction
    """

    species = ManyToOneAttribute(Species, related_name='species_coefficients')
    coefficient = FloatAttribute(nan=False)

    class Meta(obj_model.Model.Meta):
        attribute_order = ('species', 'coefficient')
        frozen_columns = 1
        tabular_orientation = TabularOrientation.inline
        ordering = ('species',)

    def serialize(self, show_compartment=True, show_coefficient_sign=True):
        """ Serialize related object

        Args:
            show_compartment (:obj:`bool`, optional): if true, show compartment
            show_coefficient_sign (:obj:`bool`, optional): if true, show coefficient sign

        Returns:
            :obj:`str`: string representation of a species and a coefficient
        """
        return self._serialize(self.species, self.coefficient,
                               show_compartment=show_compartment,
                               show_coefficient_sign=show_coefficient_sign)

    @staticmethod
    def _serialize(species, coefficient, show_compartment=True, show_coefficient_sign=True):
        """ Serialize values

        Args:
            species (:obj:`Species`): species
            coefficient (:obj:`float`): coefficient
            show_compartment (:obj:`bool`, optional): if true, show compartment
            show_coefficient_sign (:obj:`bool`, optional): if true, show coefficient sign

        Returns:
            :obj:`str`: string representation of a species and a coefficient
        """
        coefficient = float(coefficient)

        if not show_coefficient_sign:
            coefficient = abs(coefficient)

        if coefficient == 1:
            coefficient_str = ''
        elif coefficient % 1 == 0 and abs(coefficient) < 1000:
            coefficient_str = '({:.0f}) '.format(coefficient)
        else:
            coefficient_str = '({:e}) '.format(coefficient)

        if show_compartment:
            return '{}{}'.format(coefficient_str, species.serialize())
        else:
            return '{}{}'.format(coefficient_str, species.species_type.get_primary_attribute())

    @classmethod
    def deserialize(cls, attribute, value, objects, compartment=None):
        """ Deserialize value

        Args:
            attribute (:obj:`Attribute`): attribute
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            compartment (:obj:`Compartment`, optional): compartment

        Returns:
            :obj:`tuple` of `list` of `SpeciesCoefficient`, `InvalidAttribute` or `None`: tuple of cleaned value
                and cleaning error
        """
        errors = []

        if compartment:
            pattern = r'^(\(((\-?\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\) )*([a-z][a-z0-9_]*)$'
        else:
            pattern = r'^(\(((\-?\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\) )*([a-z][a-z0-9_]*\[[a-z][a-z0-9_]*\])$'

        match = re.match(pattern, value, flags=re.I)
        if match:
            errors = []

            coefficient = float(match.group(2) or 1.)

            if compartment:
                species_id = Species.gen_id(match.group(
                    5), compartment.get_primary_attribute())
            else:
                species_id = match.group(5)

            species, error = Species.deserialize(
                attribute, species_id, objects)
            if error:
                return (None, error)

            serial_val = cls._serialize(species, coefficient)
            if cls in objects and serial_val in objects[cls]:
                return (objects[cls][serial_val], None)

            obj = cls(species=species, coefficient=coefficient)
            if cls not in objects:
                objects[cls] = {}
            objects[cls][obj.serialize()] = obj
            return (obj, None)

        else:
            attr = cls.Meta.attributes['species']
            return (None, InvalidAttribute(attr, ['Invalid species coefficient']))


class PolymerSpeciesType(SpeciesType):
    """ Knowledge of a polymer

    Attributes:
        circular (:obj:`bool`): is the polymer circular
        double_stranded (:obj:`bool`): is the polymer double stranded

    Related attributes:
        loci (:obj:`list` of :obj:`PolymerLocus`): loci
    """
    circular = obj_model.BooleanAttribute()
    double_stranded = obj_model.BooleanAttribute()

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'circular', 'double_stranded',
                           'half_life', 'comments', 'references', 'database_references')

    @abc.abstractmethod
    def get_seq(self):
        """ Get the polymer sequence

        Returns:
            :obj:`Bio.Seq.Seq`: sequence
        """
        pass  # pragma: no cover

    def get_len(self):
        """ Get the polymer length

        Returns:
            :obj:`int`: length
        """
        return len(self.get_seq())

    def get_subseq(self, start, end, strand=PolymerStrand.positive):
        """ Get a subsequence

        Args:
            start (:obj:`int`): start coordinate (1-indexed)
            end (:obj:`int`): end coordinate (1-indexed)
            strand (:obj:`PolymerStrand`, optional): strand

        Returns:
            :obj:`Bio.Seq.Seq`: sequence

        Raises:
            :obj:`ValueError`: if the polymer is linear and the start or end coordinates
                are less than 1 or greater than the length of the sequence
        """
        seq = self.get_seq()
        seq_len = len(seq)

        # convert to zero-based indexing
        start -= 1

        if self.circular:
            n_wrap = int(math.floor(start / seq_len))
            start = start - seq_len * n_wrap
            end = end - seq_len * n_wrap
        elif start < 0 or end > seq_len:
            raise ValueError('Start and end coordinates for linear polymers must be at '
                             'least 1 and less than the length of the sequence')

        if end <= seq_len:
            pos_seq = seq[start:end]
        else:
            pos_seq = seq[start:] + \
                str(seq) * (int(math.floor(end / seq_len)) - 1) + \
                seq[0:end % seq_len]

        if strand == PolymerStrand.positive:
            return pos_seq
        else:
            return pos_seq.reverse_complement()


class PolymerLocus(KnowledgeBaseObject):
    """ Knowledge about a locus of a polymer

    Attributes:
        polymer (:obj:`PolymerSpeciesType`): polymer
        start (:obj:`int`): start position
        end (:obj:`int`): end position
        strand (:obj:`PolymerStrand`): strand
        references (:obj:`list` of :obj:`Reference`): references
        database_references (:obj:`list` of :obj:`DatabaseReference`): database references
    """

    cell = obj_model.ManyToOneAttribute(Cell, related_name='loci')
    polymer = obj_model.ManyToOneAttribute(
        PolymerSpeciesType, related_name='loci')
    strand = obj_model.EnumAttribute(
        PolymerStrand, default=PolymerStrand.positive)
    start = obj_model.IntegerAttribute()
    end = obj_model.IntegerAttribute()
    references = obj_model.ManyToManyAttribute(Reference, related_name='loci')
    database_references = DatabaseReferenceAttribute(related_name='loci')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'polymer', 'strand',
                           'start', 'end', 'comments', 'references', 'database_references')

    def get_seq(self):
        """ Get the sequence

        Returns:
            :obj:`Bio.Seq.Seq`: sequence
        """
        return self.polymer.get_subseq(self.start, self.end, strand=self.strand)

    def get_len(self):
        """ Get the length

        Returns:
            :obj:`int`: length
        """
        return abs(self.start - self.end) + 1


class ObservableExpression(obj_model.Model, Expression):
    """ A mathematical expression of Observables and Species

    The expression used by a `Observable`.

    Attributes:
        expression (:obj:`str`): mathematical expression for an Observable
        species (:obj:`list` of :obj:`Species`): Species used by this Observable expression
        observables (:obj:`list` of :obj:`Observable`): other Observables used by this Observable expression

    Related attributes:
        observable (:obj:`Observable`): observable
    """

    expression = LongStringAttribute(primary=True, unique=True, default='')
    species = ManyToManyAttribute(Species, related_name='observable_expressions')
    observables = ManyToManyAttribute('Observable', related_name='observable_expressions')

    class Meta(obj_model.Model.Meta, Expression.Meta):
        tabular_orientation = TabularOrientation.inline
        expression_term_models = ('Species', 'Observable')
        expression_is_linear = True
        expression_unit_registry = unit_registry

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: string representation
        """
        return Expression.serialize(self)

    @classmethod
    def deserialize(cls, value, objects):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model

        Returns:
            :obj:`tuple` of :obj:`ObservableExpression`, `InvalidAttribute` or `None`:
                tuple of cleaned value and cleaning error
        """
        return Expression.deserialize(cls, value, objects)


class Observable(KnowledgeBaseObject):
    """ Observable: a linear function of other Observables and Species

    Attributes:
        cell (:obj:`Cell`): cell
        expression (:obj:`ObservableExpression`): mathematical expression for an Observable
        units (:obj:`unit_registry.Unit`): units of expression
        references (:obj:`list` of :obj:`Reference`): references
        database_references (:obj:`list` of :obj:`DatabaseReference`): database references

    Related attributes:
        observable_expressions (:obj:`list` of :obj:`ObservableExpression`): observable expressions
        rate_law_expressions (:obj:`list` of :obj:`RateLawExpression`): rate law expressions
    """

    cell = ManyToOneAttribute(Cell, related_name='observables')
    expression = ExpressionManyToOneAttribute(ObservableExpression, related_name='observable',
                                              min_related=1, min_related_rev=1)
    units = obj_model.units.UnitAttribute(unit_registry,
                          choices=(unit_registry.parse_units('molecule'),),
                          default=unit_registry.parse_units('molecule'))
    references = obj_model.ManyToManyAttribute(Reference, related_name='observables')
    database_references = DatabaseReferenceAttribute(related_name='observables')

    class Meta(obj_model.Model.Meta, ExpressionExpressionTermMeta):
        attribute_order = ('id', 'name', 'expression', 'units',
                           'comments', 'references', 'database_references')
        expression_term_model = ObservableExpression
        expression_term_units = 'units'

    def deserialize(self, value, objects, decoded=None):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            decoded (:obj:`dict`, optional): dictionary of objects that have already been decoded

        Returns:
            :obj:`tuple` of :obj:`ObservableExpression`, `InvalidAttribute` or `None`:
                tuple of cleaned value and cleaning error
        """
        return expression.deserialize()


class Parameter(KnowledgeBaseObject):
    """ Knowledge of parameters

    Attributes:
        cell (:obj:`Cell`): cell
        value (:obj:`float`): value
        error (:obj:`float`): measurement error
        units (:obj:`unit_registry.Unit`): units of value
        references (:obj:`list` of :obj:`Reference`): references
        database_references (:obj:`list` of :obj:`DatabaseReference`): database references

    Related attributes:
        rate_law_expressions (:obj:`list` of :obj:`RateLawExpression`): rate law expressions that use a Parameter
    """

    cell = obj_model.ManyToOneAttribute(Cell, related_name='parameters')
    value = FloatAttribute(min=0)
    error = FloatAttribute(min=0)
    units = obj_model.units.UnitAttribute(unit_registry, none=True)
    references = obj_model.ManyToManyAttribute(Reference, related_name='parameters')
    evidence = obj_model. OneToManyAttribute('Evidence', related_name='parameters')
    database_references = DatabaseReferenceAttribute(related_name='parameters')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'synonyms', 'value', 'units', 'evidence', 'database_references', 'references', 'comments')
        expression_term_token_pattern = (token.NAME, )


class Validator(obj_model.Validator):
    def run(self, knowledge_base, get_related=True):
        """ Validate a knowledge_base and return its errors

        Args:
            knowledge_base (:obj:`KnowledgeBase`): knowledge base
            get_related (:obj:`bool`, optional): if true, get all related objects

        Returns:
            :obj:`InvalidObjectSet` or `None`: list of invalid objects/models and their errors
        """
        return super(Validator, self).run(knowledge_base, get_related=get_related)


#####################
#####################
# Species types


class MetaboliteSpeciesType(SpeciesType):
    """ Knowledge of a metabolite

    Attributes:
        structure (:obj:`str`): InChI-encoded structure
    """
    type = obj_model.EnumAttribute(MetaboliteSpeciesTypeType)
    concentration = obj_model.OneToManyAttribute('Concentration', related_name='metabolites')
    species_properties = obj_model.OneToOneAttribute('SpeciesTypeProperty', related_name='metabolites')
    #evidence = obj_model. OneToManyAttribute('Evidence', related_name='metabolites')

    class Meta(obj_model.Model.Meta):
        verbose_name = 'Metabolite'
        attribute_order = ('id', 'name', 'synonyms', 'type', 'concentration', 'species_properties',
                           'database_references', 'references', 'comments')

    def get_structure(self, ph=7.95):
        """ Get the structure

        Returns:
            :obj:`str`: structure
        """

        # return self.structure

        mol = openbabel.OBMol()
        conversion = openbabel.OBConversion()
        conversion.SetInFormat('inchi')
        conversion.ReadString(mol, self.structure)
        mol.CorrectForPH(ph)
        conversion.SetOutFormat('inchi')
        protonated_inchi = conversion.WriteString(mol)

        return protonated_inchi

    def to_openbabel_mol(self):
        """ Convert species type to an Open Babel molecule

        Returns:
            :obj:`openbabel.OBMol`: Open Babel molecule
        """
        mol = openbabel.OBMol()
        obConversion = openbabel.OBConversion()
        obConversion.SetInFormat('inchi')
        obConversion.ReadString(mol, self.structure)

        return mol

    def get_empirical_formula(self, ph=7.95):
        """ Get the empirical formula

        Returns:
            :obj:`chem.EmpiricalFormula`: empirical formula
        """

        #mol = self.to_openbabel_mol()
        # return chem.EmpiricalFormula(mol.GetFormula().rstrip('+-'))

        mol = self.to_openbabel_mol()
        conversion = openbabel.OBConversion()
        conversion.SetInFormat('inchi')
        conversion.ReadString(mol, self.structure)
        mol.CorrectForPH(ph)
        conversion.SetOutFormat('inchi')
        protontated_inchi = conversion.WriteString(mol)
        protonated_formula = mol.GetFormula().rstrip('+-')

        return chem.EmpiricalFormula(protonated_formula)

    def get_charge(self, ph=7.95):
        """ Get the charge

        Returns:
            :obj:`int`: charge
        """

        #mol = self.to_openbabel_mol()
        # return mol.GetTotalCharge()

        mol = self.to_openbabel_mol()
        conversion = openbabel.OBConversion()
        conversion.SetInFormat('inchi')
        conversion.ReadString(mol, self.structure)
        mol.CorrectForPH(ph)
        conversion.SetOutFormat('inchi')

        return mol.GetTotalCharge()

    def get_mol_wt(self):
        """ Get the molecular weight

        Returns:
            :obj:`float`: molecular weight
        """
        mol = self.to_openbabel_mol()

        return mol.GetMolWt()


class DnaSpeciesType(PolymerSpeciesType):
    """ Knowledge of a DNA species

    Attributes:
        seq_path (:obj:`str`): path to sequence fasta file
        ploidy (:obj:`int`): ploidy
    """

    sequence_path = obj_model.StringAttribute()
    ploidy = obj_model.IntegerAttribute(min=0)

    class Meta(obj_model.Model.Meta):
        verbose_name = 'DNA'
        attribute_order = ('id', 'name', 'sequence_path', 'circular', 'double_stranded',
                           'ploidy', 'database_references', 'references', 'comments')

    def get_seq(self, start=None, end=None):
        """ Get the sequence

        Args:
            start (:obj:`int`, optional): start coordinate of the queried subsequence,
                default is the start of the full sequence
            end (:obj:`int`, optional): end coordinate of the queried subsequence,
                default is the end of the full sequence

        Returns:
            :obj:`Bio.Seq.Seq`: structure
        """

        seq_idx = Fasta(self.sequence_path, as_raw=True)
        start = start or 1
        end = end or len(seq_idx[self.id][:])

        seq = seq_idx[self.id][start-1:end]

        return Bio.Seq.Seq(seq, alphabet=Bio.Alphabet.DNAAlphabet())

    def get_empirical_formula(self):
        """ Get the empirical formula for a DNA molecule with

        * 5' monophosphate (for linear molecules)
        * Deprotonated phosphate oxygens

        * Linear DNA

            :math:`N_A * dAMP + N_C * dCMP + N_G * dGMP + N_T * dTMP - (L - 1) * OH`

        * Circular DNA

            :math:`N_A * dAMP + N_C * dCMP + N_G * dGMP + N_T * dTMP - L * OH`

        N's in the sequence will be distributed into the four bases by preserving the original ratio

        Returns:
           :obj:`chem.EmpiricalFormula`: empirical formula
        """
        seq = self.get_seq()
        n_a = seq.upper().count('A')
        n_c = seq.upper().count('C')
        n_g = seq.upper().count('G')
        n_t = seq.upper().count('T')
        n_n = seq.upper().count('N')

        l = len(seq)
        known_bases = n_a + n_c + n_g + n_t
        n_a += round(n_a / known_bases * n_n)
        n_c += round(n_c / known_bases * n_n)
        n_g += round(n_g / known_bases * n_n)
        n_t = l - (n_a + n_c + n_g)

        if self.double_stranded:
            n_a = n_a + n_t
            n_t = n_a
            n_c = n_c + n_g
            n_g = n_c

        formula = chem.EmpiricalFormula()
        formula.C = 10 * n_a + 9 * n_c + 10 * n_g + 10 * n_t
        formula.H = 12 * n_a + 12 * n_c + 12 * n_g + 13 * n_t - \
            (l - 1 + self.circular) * (1 + self.double_stranded)
        formula.N = 5 * n_a + 3 * n_c + 5 * n_g + 2 * n_t
        formula.O = 6 * n_a + 7 * n_c + 7 * n_g + 8 * n_t - \
            (l - 1 + self.circular) * (1 + self.double_stranded)
        formula.P = n_a + n_c + n_g + n_t

        return formula

    def get_charge(self):
        """ Get the charge for a DNA molecule with

        * 5' monophosphate (for linear molecules)
        * Deprotonated phosphate oxygens

        * Linear DNA

            :math:`-L - 1`

        * Circular DNA

            :math:`-L`

        Returns:
            :obj:`int`: charge
        """
        return (-self.get_len() - 1 + self.circular) * (1 + self.double_stranded)

    def get_mol_wt(self):
        """ Get the molecular weight for a DNA molecule with

        * 5' monophosphate (for linear molecules)
        * Deprotonated phosphate oxygens

        * Linear DNA

            :math:`N_A * MW_{dAMP} + N_C * MW_{dCMP} + N_G * MW_{dGMP} + N_T * MW_{dTMP} - (L - 1) * MW_{OH}`

        * Circular DNA

            :math:`N_A * MW_{dAMP} + N_C * MW_{dCMP} + N_G * MW_{dGMP} + N_T * MW_{dTMP} - L * MW_{OH}`

        Returns:
            :obj:`float`: molecular weight
        """
        return self.get_empirical_formula().get_molecular_weight()


class ComplexSpeciesType(SpeciesType):
    """ Knowledge of a protein complex

    Attributes:
        formation_process (:obj:`ComplexFormationType`): type of formation process
        subunits (:obj:`list` of :obj:`SpeciesTypeCoefficient`): subunits
        composition_in_uniprot (:obj:`str`): protein subunit composition in uniprot IDs
        complex_type (:obj:`ComplexType`): type of complex
        binding (:obj:`str`): strand of DNA bound if involved
        region (:obj:`str`): region where DNA is bound if involved

    """

    formation_process = obj_model.EnumAttribute(ComplexFormationType)
    subunits = SubunitAttribute(related_name='complex')
    composition_in_uniprot = obj_model.StringAttribute()
    complex_type = obj_model.StringAttribute()  # EnumAttribute(ComplexType)
    binding = obj_model.StringAttribute()
    region = obj_model.StringAttribute()

    class Meta(obj_model.Model.Meta):
        verbose_name = 'Complexes'
        attribute_order = ('id', 'name', 'formation_process', 'subunits',
                           'composition_in_uniprot', 'complex_type', 'binding', 'region',
                           'half_life', 'comments', 'references', 'database_references')

    def get_empirical_formula(self):
        """ Get the empirical formula

        Returns:
            :obj:`chem.EmpiricalFormula`: empirical formula
        """
        # Formula addition
        formula = chem.EmpiricalFormula()
        for subunit in self.subunits:
            for coeff in range(0, abs(int(subunit.coefficient))):
                formula = formula + subunit.species_type.get_empirical_formula()

        return formula

    def get_charge(self):
        """ Get the charge at physiological pH

        Returns:
            :obj:`int`: charge
        """
        charge = 0
        for subunit in self.subunits:
            charge += abs(subunit.coefficient)*subunit.species_type.get_charge()

        return charge

    def get_mol_wt(self):
        """ Get the molecular weight

        Returns:
            :obj:`float`: molecular weight
        """
        weight = 0
        for subunit in self.subunits:
            weight += abs(subunit.coefficient)*subunit.species_type.get_mol_wt()

        return weight


#####################
#####################
# Reactions and related classes

class RateLawDirection(int, CaseInsensitiveEnum):
    """ Rate law directions """
    backward = -1
    forward = 1


class RateLawExpression(obj_model.Model, Expression):
    """ Rate law expression

    Attributes:
        expression (:obj:`str`): mathematical expression of the rate law
        parameters (:obj:`list` of :obj:`Parameter`): parameters whose values are used in the rate law
        species (:obj:`list` of :obj:`Species`): species whose dynamic concentrations are used in the rate law
        observables (:obj:`list` of :obj:`Observable`): observables whose values are used in the rate law

    Related attributes:
        rate_law (:obj:`RateLaw`): the `RateLaw` which uses this `RateLawExpression`
    """
    expression = LongStringAttribute(primary=True, unique=True, default='')
    parameters = ManyToManyAttribute('Parameter', related_name='rate_law_expressions')
    species = ManyToManyAttribute(Species, related_name='rate_law_expressions')
    observables = ManyToManyAttribute(Observable, related_name='rate_law_expressions')

    class Meta(obj_model.Model.Meta, Expression.Meta):
        attribute_order = ('expression', 'parameters', 'species', 'observables')
        tabular_orientation = TabularOrientation.inline
        ordering = ('expression',)
        expression_term_models = ('Parameter', 'Species', 'Observable')
        expression_unit_registry = unit_registry

    def serialize(self):
        """ Generate string representation
        Returns:
            :obj:`str`: value of primary attribute
        """
        return Expression.serialize(self)

    @classmethod
    def deserialize(cls, value, objects):
        """ Deserialize value
        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
        Returns:
            :obj:`tuple` of :obj:`RateLawExpression`, `InvalidAttribute` or `None`: tuple of cleaned value
                and cleaning error
        """
        return Expression.deserialize(cls, value, objects)


class RateLaw(KnowledgeBaseObject):
    """ Rate law

    Attributes:
        reaction (:obj:`Reaction`): reaction
        direction (:obj:`RateLawDirection`): direction
        expression (:obj:`RateLawExpression`): expression
        units (:obj:`unit_registry.Unit`): units
        references (:obj:`list` of :obj:`Reference`): references
        database_references (:obj:`list` of :obj:`DatabaseReference`): database references
    """
    reaction = ManyToOneAttribute('Reaction', related_name='rate_laws')
    direction = EnumAttribute(RateLawDirection, default=RateLawDirection.forward)
    expression = ExpressionManyToOneAttribute(RateLawExpression, min_related=1, min_related_rev=1, related_name='rate_laws')
    units = obj_model.units.UnitAttribute(unit_registry,
                          choices=(unit_registry.parse_units('s^-1'),),
                          default=unit_registry.parse_units('s^-1'))
    references = obj_model.ManyToManyAttribute(Reference, related_name='rate_laws')
    database_references = DatabaseReferenceAttribute(related_name='rate_laws')

    class Meta(obj_model.Model.Meta, ExpressionExpressionTermMeta):
        attribute_order = ('id', 'reaction', 'direction', 'expression', 'units',
                           'comments', 'references', 'database_references')
        expression_term_model = RateLawExpression
        expression_term_units = 'units'

    def gen_id(self):
        """ Generate identifier
        Returns:
            :obj:`str`: identifier
        """
        return '{}_{}'.format(self.reaction.id, self.direction.name)

    def deserialize(self, value, objects, decoded=None):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            decoded (:obj:`dict`, optional): dictionary of objects that have already been decoded

        Returns:
            :obj:`tuple` of :obj:`ObservableExpression`, `InvalidAttribute` or `None`:
                tuple of cleaned value and cleaning error
        """
        return expression.deserialize()


class Reaction(KnowledgeBaseObject):
    """ Knowledge of reactions

    Attributes:
        cell (:obj:`Cell`): cell
        submodel (:obj:`str`): submodel where reaction belongs to
        participants (:obj:`list` of :obj:`SpeciesCoefficient`): participants
        reversible (:obj:`boolean`): denotes whether reaction is reversible
        references (:obj:`list` of :obj:`Reference`): references
        database_references (:obj:`list` of :obj:`DatabaseReference`): database references

    Related attributes:
        rate_laws (:obj:`list` of :obj:`RateLaw`): rate laws; if present, rate_laws[0] is the forward
            rate law, and rate_laws[1] is the backward rate law
    """

    cell = obj_model.ManyToOneAttribute(Cell, related_name='reactions')
    participants = ReactionParticipantAttribute(related_name='reactions')
    submodel = obj_model.StringAttribute()
    reversible = obj_model.BooleanAttribute()
    references = obj_model.ManyToManyAttribute(Reference, related_name='reactions')
    database_references = DatabaseReferenceAttribute(related_name='reactions')
    evidence = obj_model.OneToManyAttribute('Evidence', related_name='reactions')
    type = obj_model.EnumAttribute(ReactionType)
    enzyme = obj_model.ManyToManyAttribute(SpeciesType, related_name='reactions')
    coenzymes = obj_model.ManyToManyAttribute(SpeciesType, related_name='reactions')
    spontenaeous =obj_model.BooleanAttribute()
    deltaG = obj_model.FloatAttribute()
    parameters = obj_model.OneToManyAttribute('Parameter', related_name='reactions')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'synonyms', 'type', 'submodel', 'participants', 'reversible', 'spontenaeous', 'enzyme', 'coenzymes',
                           'deltaG', 'parameters', 'evidence', 'database_references', 'references', 'comments')

#####################
#####################
# Expansion classes

class ChromosomeFeature(KnowledgeBaseObject):
    """ Knowledge of chromosoe features

    Attributes:
        seq_path (:obj:`str`): path to sequence fasta file
        ploidy (:obj:`int`): ploidy

    Related sttributes:
        seq_path (:obj:`str`): path to sequence fasta file
        ploidy (:obj:`int`): ploidy
    """

    coordinate = obj_model.IntegerAttribute(min=0)
    length = obj_model.IntegerAttribute(min=0)
    direction = obj_model.EnumAttribute(DirectionType)
    type = obj_model.EnumAttribute(ChromosomeFeatureType)
    intensity = obj_model.FloatAttribute(min=0)
    unit = obj_model.units.UnitAttribute(unit_registry, none=True)
    polymer = obj_model.ManyToOneAttribute('DnaSpeciesType', related_name='chromosome_features')
    evidence   = obj_model.OneToManyAttribute('Evidence', related_name='chromosome_features')
    database_references = DatabaseReferenceAttribute(related_name='chromosome_features')
    references = obj_model.ManyToManyAttribute('Reference', related_name='chromosome_features')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'type', 'polymer', 'coordinate', 'length',
                           'direction', 'intensity', 'unit', 'evidence', 'database_references', 'references', 'comments')


class Evidence(obj_model.Model):
    """ Represents the measurement / observation of a property
        Attributes:
        Related attributes:
    """

    id       = obj_model.SlugAttribute(primary=True, unique=True)
    cell = obj_model.ManyToOneAttribute('Cell', related_name='evidence')
    object   =  obj_model.StringAttribute() #obj_model.ManyToOneAttribute(obj_model.Model, related_name='evidences')
    property = obj_model.StringAttribute()
    values = obj_model.FloatAttribute()
    mean = obj_model.IntegerAttribute()
    standard_error = obj_model.IntegerAttribute()
    units = obj_model.units.UnitAttribute(unit_registry, none=True) # False allows None units
    database_references = DatabaseReferenceAttribute(related_name='evidence')
    experiment = obj_model.ManyToOneAttribute('Experiment', related_name ='evidence')
    comments = obj_model.LongStringAttribute()

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'cell', 'object', 'property', 'values', 'mean', 'standard_error', 'units', 'experiment', 'database_references', 'comments')


class Experiment(obj_model.Model):
    """ Represents an experiment in which a property was measured
        Attributes:
        Related attributes:
    """

    id = obj_model.SlugAttribute(primary=True, unique=True)
    species = obj_model.StringAttribute()
    genetic_variant = obj_model.StringAttribute()
    external_media  = obj_model.StringAttribute()
    temperature	= obj_model.FloatAttribute()
    temperature_units = obj_model.units.UnitAttribute(unit_registry,
                        choices=(unit_registry.parse_units('F'),
                                 unit_registry.parse_units('C'),
                                 unit_registry.parse_units('K')),
                        default= unit_registry.parse_units('C'))
    ph = obj_model.FloatAttribute()
    ph_units = obj_model.units.UnitAttribute(unit_registry, none=True)
    experiment_design = obj_model.StringAttribute()
    measurment_technology = obj_model.StringAttribute()
    analysis_type = obj_model.StringAttribute()
    database_references = DatabaseReferenceAttribute(related_name='experiment')
    references = obj_model.ManyToManyAttribute('Reference', related_name='experiment')
    comments = obj_model.LongStringAttribute()

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'experiment_design', 'measurment_technology', 'analysis_type', 'species', 'genetic_variant', 'external_media',
                           'temperature', 'temperature_units', 'ph', 'ph_units', 'database_references',	'references', 'comments')


class Interaction(KnowledgeBaseObject):
    """ Knowledge of interactions
        Attributes:
        Related attributes:
    """

    type = obj_model.StringAttribute() # Need to convert to enumeration / ontology
    participants = obj_model.ManyToManyAttribute('SpeciesTypeCoefficient', related_name='interactions')
    binding_site_coordinate = obj_model.IntegerAttribute()
    binding_site_length = obj_model.IntegerAttribute()
    binding_site_direction = obj_model.EnumAttribute(DirectionType)
    affinity = obj_model.IntegerAttribute()
    units = obj_model.units.UnitAttribute(unit_registry, none=False) # False allows None units
    database_references = DatabaseReferenceAttribute(related_name='interactions')
    references = obj_model.ManyToManyAttribute('Reference', related_name='interactions')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'synonyms', 'type', 'participants', 'binding_site_coordinate', 'binding_site_length',
                           'binding_site_direction', 'affinity', 'units', 'database_references', 'references', 'comments')


class SpeciesTypeProperty(KnowledgeBaseObject):
    """ Knowledge of interactions
        Attributes:
        Related attributes:
    """

    structure = obj_model.LongStringAttribute()
    half_life =  obj_model.FloatAttribute()
    half_life_units = obj_model.units.UnitAttribute(unit_registry,
                        choices=(unit_registry.parse_units('s'),
                                 unit_registry.parse_units('min'),
                                 unit_registry.parse_units('hr')),
                        default= unit_registry.parse_units('min'))

    domains =  obj_model.LongStringAttribute()
    prosthetic_groups =  obj_model.LongStringAttribute()
    evidence = obj_model.OneToManyAttribute('Evidence', related_name='species_type_properties')
    database_references = DatabaseReferenceAttribute(related_name='species_type_properties')
    references = obj_model.ManyToManyAttribute('Reference', related_name='species_type_properties')

    class Meta(obj_model.Model.Meta):
        verbose_name = 'Species properties'
        attribute_order = ('id', 'name', 'synonyms', 'structure', 'half_life', 'half_life_units', 'domains',
                           'prosthetic_groups', 'evidence', 'database_references', 'references', 'comments')
