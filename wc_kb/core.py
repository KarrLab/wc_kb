""" Core schema to represent a knowledge base to build models

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Bilal Shaikh  <bilal.shaikh@columbia.edu>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2018-02-07
:Copyright: 2018, Karr Lab
:License: MIT
"""

from natsort import natsorted, ns
from math import ceil, floor, exp, log, log10, isnan
from wc_utils.util import chem
from wc_utils.util.list import det_dedupe
import abc
import Bio.SeqUtils
import enum
import math
import obj_model.abstract
import obj_model
import obj_model.extra_attributes
import openbabel
import pkg_resources
import re
import six
from obj_model import (BooleanAttribute, EnumAttribute, FloatAttribute, IntegerAttribute, PositiveIntegerAttribute,
                       RegexAttribute, SlugAttribute, StringAttribute, LongStringAttribute, UrlAttribute,
                       OneToOneAttribute, ManyToOneAttribute, ManyToManyAttribute,
                       InvalidModel, InvalidObject, InvalidAttribute, TabularOrientation)
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


class RegulatoryElementType(enum.Enum):
    """ Type of regulatory element """
    promoter = 1
    promoter_flanking_region = 2
    enhancer = 3
    CTCF_binding_site = 4
    TF_binding_site = 5
    open_chromatin_region = 6


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


ConcentrationUnit = enum.Enum('ConcentrationUnit', type=int, names=[
    ('molecules', 1),
    ('M', 2),
    ('mM', 3),
    ('uM', 4),
    ('nM', 5),
    ('pM', 6),
    ('fM', 7),
    ('aM', 8),
    ('mol dm^-2', 9),
])


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


class ObservableSpeciesParticipantAttribute(ManyToManyAttribute):
    """ Inline separated list of species and their weights of an observable

    Attributes:
        separator (:obj:`str`): list separator
    """

    def __init__(self, related_class, separator=' + ', related_name='', verbose_name='', verbose_related_name='', help=''):
        """
        Args:
            related_class (:obj:`class`): related class
            separator (:obj:`str`, optional): list separator
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            help (:obj:`str`, optional): help message
        """
        super(ObservableSpeciesParticipantAttribute, self).__init__(related_class, related_name=related_name,
                                                                    verbose_name=verbose_name,
                                                                    verbose_related_name=verbose_related_name,
                                                                    help=help)
        self.separator = separator

    def serialize(self, spec_coeffs, encoded=None):
        """ Serialize related object

        Args:
            spec_coeffs (:obj:`list` of :obj:`Model`): Python representation of species and their coefficients
            encoded (:obj:`dict`, optional): dictionary of objects that have already been encoded

        Returns:
            :obj:`str`: simple Python representation
        """
        if not spec_coeffs:
            return ''

        spec_coeff_strs = []
        for spec_coeff_obj in spec_coeffs:
            spec_coeff_str = spec_coeff_obj.serialize(
                show_compartment=True, show_coefficient_sign=True)
            spec_coeff_strs.append(spec_coeff_str)

        return self.separator.join(spec_coeff_strs)

    def deserialize(self, value, objects, decoded=None):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            decoded (:obj:`dict`, optional): dictionary of objects that have already been decoded

        Returns:
            :obj:`tuple` of `list` of `related_class`, `InvalidAttribute` or `None`: tuple of cleaned value
                and cleaning error
        """
        if not value:
            return ([], None)

        pat_id = r'([a-z][a-z0-9_]*)'
        pat_coeff = r'\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\)'
        pat_spec_coeff = r'({} )*({}\[{}\])'.format(pat_coeff, pat_id, pat_id)
        pat_observable = r'^{}( \+ {})*$'.format(pat_spec_coeff, pat_spec_coeff)
        if not re.match(pat_observable, value, flags=re.I):
            return (None, InvalidAttribute(self, ['Incorrectly formatted observable: {}'.format(value)]))

        spec_coeff_objs = []
        errors = []
        for spec_coeff_match in re.findall(pat_spec_coeff, value, flags=re.I):
            spec_type_errors = []
            
            spec_type_id = spec_coeff_match[5]
            
            spec_type = None
            for species_type_cls in get_subclasses(SpeciesType):
                if species_type_cls in objects and spec_type_id in objects[species_type_cls]:
                    spec_type = objects[species_type_cls][spec_type_id]
                    break
            if not spec_type:
                spec_type_errors.append(
                    'Undefined species type "{}"'.format(spec_type_id))            
            
            compartment_id = spec_coeff_match[6]
            if compartment_id in objects[Compartment]:
                compartment = objects[Compartment][compartment_id]
            else:
                spec_type_errors.append(
                    'Undefined compartment "{}"'.format(compartment_id))

            coefficient = float(spec_coeff_match[1] or 1.)

            if spec_type_errors:
                errors += spec_type_errors
            elif coefficient != 0:
                spec_id = Species.gen_id(
                    spec_type.get_primary_attribute(), compartment.get_primary_attribute())
                obj, error = Species.deserialize(self, spec_id, objects)

                if error:
                    raise ValueError('Invalid object "{}"'.format(spec_primary_attribute)
                                     )  # pragma: no cover # unreachable due to error checking above

                if self.related_class not in objects:
                    objects[self.related_class] = {}
                serialized_value = self.related_class._serialize(
                    obj, coefficient)
                if serialized_value in objects[self.related_class]:
                    spec_coeff_obj = objects[self.related_class][serialized_value]
                else:
                    spec_coeff_obj = self.related_class(
                        species=obj, coefficient=coefficient)
                    objects[self.related_class][serialized_value] = spec_coeff_obj
                spec_coeff_objs.append(spec_coeff_obj)

        if errors:
            return (None, InvalidAttribute(self, errors))
        return (spec_coeff_objs, None)


class ObservableObservableParticipantAttribute(ManyToManyAttribute):
    """ Inline separated list of observables and their weights of an observable

    Attributes:
        separator (:obj:`str`): list separator
    """

    def __init__(self, related_class, separator=' + ', related_name='', verbose_name='', verbose_related_name='', help=''):
        """
        Args:
            related_class (:obj:`class`): related class
            separator (:obj:`str`, optional): list separator
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            help (:obj:`str`, optional): help message
        """
        super(ObservableObservableParticipantAttribute, self).__init__(related_class, related_name=related_name,
                                                                       verbose_name=verbose_name,
                                                                       verbose_related_name=verbose_related_name,
                                                                       help=help)
        self.separator = separator

    def serialize(self, obs_coeffs, encoded=None):
        """ Serialize related object

        Args:
            obs_coeffs (:obj:`list` of :obj:`Model`): Python representation of observables and their coefficients
            encoded (:obj:`dict`, optional): dictionary of objects that have already been encoded

        Returns:
            :obj:`str`: simple Python representation
        """
        if not obs_coeffs:
            return ''

        obs_coeff_strs = []
        for obs_coeff_obj in obs_coeffs:
            obs_coeff_str = obs_coeff_obj.serialize()
            obs_coeff_strs.append(obs_coeff_str)

        return self.separator.join(obs_coeff_strs)

    def deserialize(self, value, objects, decoded=None):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            decoded (:obj:`dict`, optional): dictionary of objects that have already been decoded

        Returns:
            :obj:`tuple` of `list` of `related_class`, `InvalidAttribute` or `None`: tuple of cleaned value
                and cleaning error
        """
        if not value:
            return ([], None)

        pat_id = r'([a-z][a-z0-9_]*)'
        pat_coeff = r'\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\)'
        pat_obs_coeff = r'({} )*({})'.format(pat_coeff, pat_id, pat_id)
        pat_observable = r'^{}( \+ {})*$'.format(pat_obs_coeff, pat_obs_coeff)
        if not re.match(pat_observable, value, flags=re.I):
            return (None, InvalidAttribute(self, ['Incorrectly formatted observable: {}'.format(value)]))

        obs_coeff_objs = []
        errors = []
        for obs_coeff_match in re.findall(pat_obs_coeff, value, flags=re.I):
            obs_errors = []

            obs_id = obs_coeff_match[5]
            if obs_id in objects[Observable]:
                obs = objects[Observable][obs_id]
            else:
                obs_errors.append('Undefined observable "{}"'.format(obs_id))

            coefficient = float(obs_coeff_match[1] or 1.)

            if obs_errors:
                errors += obs_errors
            elif coefficient != 0:
                if self.related_class not in objects:
                    objects[self.related_class] = {}
                serialized_value = self.related_class._serialize(
                    obs, coefficient)
                if serialized_value in objects[self.related_class]:
                    obs_coeff_obj = objects[self.related_class][serialized_value]
                else:
                    obs_coeff_obj = self.related_class(
                        observable=obs, coefficient=coefficient)
                    objects[self.related_class][serialized_value] = obs_coeff_obj
                obs_coeff_objs.append(obs_coeff_obj)

        if errors:
            return (None, InvalidAttribute(self, errors))
        return (obs_coeff_objs, None)


#####################
#####################
# Base classes


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

    class Meta(obj_model.Model.Meta):
        attribute_order = ('database', 'id')
        tabular_orientation = TabularOrientation.inline
        unique_together = (('database', 'id'), )
        ordering = ('database', 'id')

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: value of primary attribute
        """
        return '{}:{}'.format(self.database, self.id)


class KnowledgeBaseObject(obj_model.Model):
    """ Knowledge of a biological entity

    Attributes:
        id (:obj:`str`): identifier
        name (:obj:`str`): name
        comments (:obj:`str`): comments
    """
    id = obj_model.SlugAttribute(primary=True, unique=True)
    name = obj_model.StringAttribute()
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
    standard_id = obj_model.StringAttribute()
    cell = obj_model.ManyToOneAttribute(Cell, related_name='references')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'standard_id')        


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
        attribute_order = ('id', 'name', 'volumetric_fraction', 'comments', 'references', 'database_references')


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
        species_coefficients (:obj:`list` of :obj:`SpeciesCoefficient`): participations in reactions and observables
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


class Concentration(obj_model.Model):
    """ Species concentration

    Attributes:
        cell (:obj:`Cell`): cell
        species (:obj:`Species`): species
        value (:obj:`float`): value
        units (:obj:`str`): units; default units is 'M'
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
        database_references (:obj:`list` of :obj:`DatabaseReference`): database references
    """
    cell = obj_model.ManyToOneAttribute(Cell, related_name='concentrations')
    species = OneToOneSpeciesAttribute(related_name='concentration')
    value = FloatAttribute(min=0)
    units = EnumAttribute(ConcentrationUnit, default=ConcentrationUnit.M)
    comments = LongStringAttribute()
    references = ManyToManyAttribute(Reference, related_name='concentrations')
    database_references = DatabaseReferenceAttribute(related_name='concentrations')

    class Meta(obj_model.Model.Meta):
        unique_together = (('species', ), )
        attribute_order = ('species', 'value', 'units', 'comments', 'references', 'database_references')
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


class Observable(six.with_metaclass(obj_model.abstract.AbstractModelMeta, KnowledgeBaseObject)):
    """Knowledge of an observable to include in a model

    Attributes:
        cell (:obj:`Cell`): The cell that the observable is in
        species (:obj:`list` of :obj:`SpeciesCoefficient`): A list of the species and the
            coefficients to be included in the observable
        observables (:obj:`list` of :obj:`ObservableCoefficient`): list of component observables
            and their coefficients
        references (:obj:`list` of :obj:`Reference`): references
        database_references (:obj:`list` of :obj:`DatabaseReference`): database references

    Related Attributes:
        observable_coefficients (:obj:`list` of :obj:`ObservableCoefficient`): Participants in observables
"""
    cell = ManyToOneAttribute(Cell, related_name='observables')
    species = ObservableSpeciesParticipantAttribute(
        'SpeciesCoefficient', related_name='observables')
    observables = ObservableObservableParticipantAttribute(
        'ObservableCoefficient', related_name='observables')
    references = obj_model.ManyToManyAttribute(Reference, related_name='observables')
    database_references = DatabaseReferenceAttribute(related_name='observables')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'cell', 'species',
                           'observables', 'comments', 'references', 'database_references')

    
class ObservableCoefficient(obj_model.Model):
    """ A tuple of observable and coefficient

    Attributes:
        observable (:obj:`Observable`): observable
        coefficient (:obj:`float`): coefficient
    """
    observable = ManyToOneAttribute(
        Observable, related_name='observable_coefficients')
    coefficient = FloatAttribute(nan=False)

    class Meta(obj_model.Model.Meta):
        attribute_order = ('observable', 'coefficient')
        frozen_columns = 1
        tabular_orientation = TabularOrientation.inline
        ordering = ('observable',)

    def serialize(self):
        """ Serialize related object

        Returns:
            :obj:`str`: simple Python representation
        """
        return self._serialize(self.observable, self.coefficient)

    @staticmethod
    def _serialize(observable, coefficient):
        """ Serialize values

        Args:
            observable (:obj:`Observable`): observable
            coefficient (:obj:`float`): coefficient

        Returns:
            :obj:`str`: simple Python representation
        """
        coefficient = float(coefficient)

        if coefficient == 1:
            coefficient_str = ''
        elif coefficient % 1 == 0 and abs(coefficient) < 1000:
            coefficient_str = '({:.0f}) '.format(coefficient)
        else:
            coefficient_str = '({:e}) '.format(coefficient)

        return '{}{}'.format(coefficient_str, observable.serialize())

    @classmethod
    def deserialize(cls, attribute, value, objects):
        """ Deserialize value

        Args:
            attribute (:obj:`Attribute`): attribute
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model

        Returns:
            :obj:`tuple` of `list` of `ObservableCoefficient`, `InvalidAttribute` or `None`: tuple of cleaned value
                and cleaning error
        """
        errors = []

        pattern = r'^(\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\) )*([a-z][a-z0-9_]*)$'

        match = re.match(pattern, value, flags=re.I)
        if match:
            errors = []

            coefficient = float(match.group(2) or 1.)

            obs_id = match.group(5)

            observable, error = Observable.deserialize(obs_id, objects)
            if error:
                return (None, error)

            serial_val = cls._serialize(observable, coefficient)
            if cls in objects and serial_val in objects[cls]:
                return (objects[cls][serial_val], None)

            if cls not in objects:
                objects[cls] = {}
            serialized_val = cls._serialize(observable, coefficient)
            if serialized_val in objects[cls]:
                obj = objects[cls][serialized_val]
            else:
                obj = cls(observable=observable, coefficient=coefficient)
                objects[cls][serialized_val] = obj
            return (obj, None)

        else:
            attr = cls.Meta.attributes['observable']
            return (None, InvalidAttribute(attr, ['Invalid observable coefficient']))


#####################
#####################
# Species types


class MetaboliteSpeciesType(SpeciesType):
    """ Knowledge of a metabolite

    Attributes:
        structure (:obj:`str`): InChI-encoded structure
    """
    structure = obj_model.StringAttribute()

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'structure',
                           'half_life', 'comments', 'references', 'database_references')

    def get_structure(self):
        """ Get the structure

        Returns:
            :obj:`str`: structure
        """
        return self.structure

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

    def get_empirical_formula(self):
        """ Get the empirical formula

        Returns:
            :obj:`chem.EmpiricalFormula`: empirical formula
        """
        mol = self.to_openbabel_mol()
        return chem.EmpiricalFormula(mol.GetFormula().rstrip('+-'))

    def get_charge(self):
        """ Get the charge

        Returns:
            :obj:`int`: charge
        """
        mol = self.to_openbabel_mol()
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
        seq (:obj:`Bio.Seq.Seq`): sequence
        ploidy (:obj:`int`): ploidy
    """

    seq = obj_model.extra_attributes.BioSeqAttribute(verbose_name='Sequence')
    ploidy = obj_model.IntegerAttribute(min=0)

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'seq', 'circular', 'double_stranded', 
            'ploidy', 'half_life', 'comments', 'references', 'database_references')
        verbose_name = 'DNA species type'

    def get_seq(self, start=None, end=None):
        """ Get the sequence

        Returns:
            :obj:`Bio.Seq.Seq`: structure
        """

        start = start or 0
        end = end or len(self.seq)

        return self.seq[start:end]

    def get_empirical_formula(self):
        """ Get the empirical formula for a DNA molecule with

        * 5' monophosphate (for linear molecules)
        * Deprotonated phosphate oxygens

        * Linear DNA

            :math:`N_A * dAMP + N_C * dCMP + N_G * dGMP + N_T * dTMP - (L - 1) * OH`

        * Circular DNA

            :math:`N_A * dAMP + N_C * dCMP + N_G * dGMP + N_T * dTMP - L * OH`

        Returns:
           :obj:`chem.EmpiricalFormula`: empirical formula
        """
        seq = self.get_seq()
        n_a = seq.count('A')
        n_c = seq.count('C')
        n_g = seq.count('G')
        n_t = seq.count('T')
        l = len(seq)

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
# Reactions

class RateLawDirection(int, CaseInsensitiveEnum):
    """ Rate law directions """
    backward = -1
    forward = 1


class RateLawEquationAttribute(ManyToOneAttribute):
    """ Rate law equation """

    def __init__(self, related_name='', verbose_name='', verbose_related_name='', help=''):
        """
        Args:
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            help (:obj:`str`, optional): help message
        """
        super(RateLawEquationAttribute, self).__init__('RateLawEquation',
                                                       related_name=related_name, min_related=1, min_related_rev=1,
                                                       verbose_name=verbose_name, verbose_related_name=verbose_related_name, help=help)

    def serialize(self, rate_law_equation, encoded=None):
        """ Serialize related object

        Args:
            rate_law_equation (:obj:`RateLawEquation`): the related `RateLawEquation`
            encoded (:obj:`dict`, optional): dictionary of objects that have already been encoded

        Returns:
            :obj:`str`: simple Python representation of the rate law equation
        """
        return rate_law_equation.serialize()

    def deserialize(self, value, objects, decoded=None):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            decoded (:obj:`dict`, optional): dictionary of objects that have already been decoded

        Returns:
            :obj:`tuple` of `object`, `InvalidAttribute` or `None`: tuple of cleaned value and cleaning error
        """
        return RateLawEquation.deserialize(self, value, objects)


class RateLaw(obj_model.Model):
    """ Rate law

    Attributes:
        reaction (:obj:`Reaction`): reaction
        direction (:obj:`RateLawDirection`): direction
        equation (:obj:`RateLawEquation`): equation
        k_cat (:obj:`float`): k_cat for law with Michaelis–Menten kinetics (units: 1/sec)
        k_m (:obj:`float`): K_m for law with Michaelis–Menten kinetics (units: mol/L)
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
        database_references (:obj:`list` of :obj:`DatabaseReference`): database references
    """

    reaction = ManyToOneAttribute('Reaction', related_name='rate_laws')
    direction = EnumAttribute(RateLawDirection, default=RateLawDirection.forward)
    equation = RateLawEquationAttribute(related_name='rate_laws')
    k_cat = FloatAttribute(min=0, nan=True)
    k_m = FloatAttribute(min=0, nan=True)
    comments = obj_model.StringAttribute()
    references = obj_model.ManyToManyAttribute(Reference, related_name='rate_laws')
    database_references = DatabaseReferenceAttribute(related_name='rate_laws')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('reaction', 'direction', 'equation', 'k_cat', 'k_m',
                           'comments', 'references', 'database_references')
        unique_together = (('reaction', 'direction'), )
        ordering = ('reaction', 'direction',)

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: value of primary attribute
        """
        return '{}.{}'.format(self.reaction.serialize(), self.direction.name)


class RateLawEquation(obj_model.Model):
    """ Rate law equation

    Attributes:
        expression (:obj:`str`): mathematical expression of the rate law
        modifiers (:obj:`list` of :obj:`Species`): species whose concentrations are used in the rate law
        parameters (:obj:`list` of :obj:`Parameter`): parameters whose values are used in the rate law

    Related attributes:
        rate_law (:obj:`RateLaw`): the `RateLaw` which uses this `RateLawEquation`
    """
    expression = LongStringAttribute(primary=True, unique=True)
    modifiers = ManyToManyAttribute(Species, related_name='rate_law_equations')
    parameters = ManyToManyAttribute('Parameter', related_name='rate_law_equations')

    class Meta(obj_model.Model.Meta):
        """
        Attributes:
            valid_functions (:obj:`tuple` of `builtin_function_or_method`): tuple of functions that
                can be used in a `RateLawEquation`s `expression`
            valid_models (:obj:`tuple` of `str`): names of `obj_model.Model`s in this module that a
                `RateLawEquation` is allowed to reference in its `expression`    
        """
        attribute_order = ('expression', 'modifiers', 'parameters')
        tabular_orientation = TabularOrientation.inline
        ordering = ('rate_law',)
        valid_functions = (ceil, floor, exp, pow, log, log10, min, max)
        valid_models = ('Species', 'Parameter')

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: value of primary attribute
        """
        return self.expression

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
        modifiers = []
        parameters = []
        errors = []
        modifier_pattern = r'(^|[^a-z0-9_])({}\[{}\])([^a-z0-9_]|$)'.format(SpeciesType.id.pattern[1:-1],
                                                                           Compartment.id.pattern[1:-1])
        parameter_pattern = r'(^|[^a-z0-9_\[\]])({})([^a-z0-9_\[\]]|$)'.format(Parameter.id.pattern[1:-1])

        reserved_names = set([func.__name__ for func in RateLawEquation.Meta.valid_functions] + ['k_cat', 'k_m'])

        try:
            for match in re.findall(modifier_pattern, value, flags=re.I):
                species, error = Species.deserialize(attribute, match[1], objects)
                if error:
                    errors += error.messages
                else:
                    modifiers.append(species)
            for match in re.findall(parameter_pattern, value, flags=re.I):
                if match[1] not in reserved_names:
                    parameter, error = Parameter.deserialize(match[1], objects)
                    if error:
                        errors += error.messages
                    else:
                        parameters.append(parameter)        
        except Exception as e:
            errors += ["deserialize fails on '{}': {}".format(value, str(e))]

        if errors:
            attr = cls.Meta.attributes['expression']
            return (None, InvalidAttribute(attribute, errors))

        # return value
        if cls not in objects:
            objects[cls] = {}
        serialized_val = value
        if serialized_val in objects[cls]:
            obj = objects[cls][serialized_val]
        else:
            obj = cls(expression=value, modifiers=det_dedupe(modifiers), parameters=det_dedupe(parameters))
            objects[cls][serialized_val] = obj
        return (obj, None)


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
    submodel   = obj_model.StringAttribute()
    reversible = obj_model.BooleanAttribute()
    references = obj_model.ManyToManyAttribute(Reference, related_name='reactions')
    database_references = DatabaseReferenceAttribute(related_name='reactions')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'submodel', 'participants', 'reversible',
                           'comments', 'references', 'database_references')


class Parameter(KnowledgeBaseObject):
    """ Knowledge of parameters

    Attributes:
        value (:obj:`float`): value
        error (:obj:`float`): measurement error
        units (:obj:`str`): units of value
        references (:obj:`list` of :obj:`Reference`): references
        database_references (:obj:`list` of :obj:`DatabaseReference`): database references
    
    Related attributes:
        rate_law_equations (:obj:`list` of :obj:`RateLawEquation`): rate law equations that use a Parameter
    """
    
    value = FloatAttribute(min=0)
    error = FloatAttribute(min=0)
    units = StringAttribute()
    references = obj_model.ManyToManyAttribute(Reference, related_name='parameters')
    database_references = DatabaseReferenceAttribute(related_name='parameters')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'value', 'error', 'units', 'comments', 
                            'references', 'database_references')


class Property(KnowledgeBaseObject):
    """ Other properties of cells

    Attributes:
        cell (:obj:`Cell`): cell
        value (:obj:`float`): value
        units (:obj:`str`): units
        references (:obj:`list` of :obj:`Reference`): references
        database_references (:obj:`list` of :obj:`DatabaseReference`): database references

    """
    cell = obj_model.ManyToOneAttribute(Cell, related_name='properties')
    value = obj_model.FloatAttribute()
    units = obj_model.StringAttribute()
    references = obj_model.ManyToManyAttribute(Reference, related_name='properties')
    database_references = DatabaseReferenceAttribute(related_name='properties')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'value', 'units', 'comments',
                           'references', 'database_references')
        verbose_name_plural = 'Properties'
