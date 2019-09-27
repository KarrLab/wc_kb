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
from wc_utils.util.chem.core import get_major_micro_species, OpenBabelUtils
from wc_utils.util.list import det_dedupe
from wc_utils.util.units import unit_registry
from wc_onto import onto as kbOnt
import abc
import Bio.Alphabet
import Bio.Seq
import enum
import math
import obj_tables.abstract
import obj_tables
import obj_tables.units
import openbabel
import pkg_resources
import re
import six
import token
from obj_tables import (BooleanAttribute, EnumAttribute, FloatAttribute, IntegerAttribute, PositiveIntegerAttribute,
                       RegexAttribute, SlugAttribute, StringAttribute, LongStringAttribute, UrlAttribute,
                       OneToOneAttribute, ManyToOneAttribute, ManyToManyAttribute,
                       InvalidModel, InvalidObject, InvalidAttribute, TableFormat)
from obj_tables.expression import (ExpressionOneToOneAttribute, ExpressionManyToOneAttribute,
                                  ExpressionStaticTermMeta, ExpressionDynamicTermMeta,
                                  ExpressionExpressionTermMeta, Expression,
                                  ParsedExpression, ParsedExpressionError)
from wc_utils.util.enumerate import CaseInsensitiveEnum
from wc_utils.util.types import get_subclasses
from obj_tables.ontology import OntologyAttribute
from wc_utils.util.ontology import are_terms_equivalent
import os


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

PolymerDirection = enum.Enum(value='PolymerDirection', names=[
    ('forward', 1),
    ('reverse', -1), ])


#####################
#####################
# Attributes

class SubunitAttribute(ManyToManyAttribute):
    """ Subunits """

    def __init__(self, related_name='', verbose_name='', verbose_related_name='', description=''):
        """
        Args:
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            description (:obj:`str`, optional): description
        """

        super(SubunitAttribute, self).__init__('SpeciesTypeCoefficient',
                                               related_name=related_name,
                                               verbose_name=verbose_name,
                                               verbose_related_name=verbose_related_name,
                                               description=description)

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

    def __init__(self, related_name='', verbose_name='', verbose_related_name='', description=''):
        """
        Args:
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            description (:obj:`str`, optional): description
        """
        super(OneToOneSpeciesAttribute, self).__init__('Species',
                                                       related_name=related_name, min_related=1, min_related_rev=0,
                                                       verbose_name=verbose_name, verbose_related_name=verbose_related_name, description=description)

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


class IdentifierAttribute(ManyToManyAttribute):
    """ Identifier attribute """

    def __init__(self, related_name='', verbose_name='', verbose_related_name='', description=''):
        """
        Args:
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            description (:obj:`str`, optional): description
        """
        super(IdentifierAttribute, self).__init__(Identifier,
                                                         related_name=related_name, min_related=0, min_related_rev=0,
                                                         verbose_name=verbose_name, verbose_related_name=verbose_related_name, description=description)

    def serialize(self, identifiers, encoded=None):
        """ Serialize related object
        Args:
            identifiers (:obj:`list` of :obj:`Model`): a list of instances of Identifier Python representation
            encoded (:obj:`dict`, optional): dictionary of objects that have already been encoded
        Returns:
            :obj:`str`: simple Python representation
        """
        if not identifiers:
            return ''

        return ', '.join(obj_tables.serialize() for obj_tables in identifiers)

    def deserialize(self, value, objects, decoded=None):
        """ Deserialize value
        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            decoded (:obj:`dict`, optional): dictionary of objects that have already been decoded
        Returns:
            :obj:`tuple` of :obj:`list` of :obj:`Identifier`, :obj:`InvalidAttribute` or :obj:`None`: :obj:`tuple` of cleaned value
                and cleaning error
        """
        if not value:
            return ([], None)

        obj_pattern = r'({}) *\: *({})'.format(Identifier.namespace.pattern[1:-1], Identifier.id.pattern[1:-1])
        lst_pattern = obj_pattern + r'( *, *{})*'.format(obj_pattern)

        if not re.match(lst_pattern, value, flags=re.I):
            return (None, InvalidAttribute(self, ['Incorrectly formatted list of identifiers: {}'.format(value)]))

        objs = []
        for pat_match in re.findall(obj_pattern, value, flags=re.I):
            namespace_name = pat_match[0]
            data_id = pat_match[1]
            if self.related_class not in objects:
                objects[self.related_class] = {}
            serialized_value = self.related_class()._serialize(namespace=namespace_name, id=data_id)
            if serialized_value in objects[self.related_class]:
                obj = objects[self.related_class][serialized_value]
            else:
                obj = self.related_class(namespace=namespace_name, id=data_id)
                objects[self.related_class][serialized_value] = obj
            objs.append(obj)
        return (objs, None)


class ReactionParticipantAttribute(ManyToManyAttribute):
    """ Reaction participants """

    def __init__(self, related_name='', verbose_name='', verbose_related_name='', description=''):
        """
        Args:
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            description (:obj:`str`, optional): description
        """
        super(ReactionParticipantAttribute, self).__init__('SpeciesCoefficient', related_name=related_name,
                                                           verbose_name=verbose_name,
                                                           verbose_related_name=verbose_related_name,
                                                           description=description)

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

        st_id = SpeciesType.id.pattern[1:-1]
        comp_id = Compartment.id.pattern[1:-1]
        stoch = r'\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\)'
        gbl_part = r'({} )*({})'.format(stoch, st_id)
        lcl_part = r'({} )*({}\[{}\])'.format(stoch, st_id, comp_id)
        gbl_side = r'{}( \+ {})*'.format(gbl_part, gbl_part)
        lcl_side = r'{}( \+ {})*'.format(lcl_part, lcl_part)
        gbl_pattern = r'^\[({})\]: ({}) ==> ({})$'.format(
            comp_id, gbl_side, gbl_side)
        lcl_pattern = r'^({}) ==> ({})$'.format(lcl_side, lcl_side)
        
        import_pattern = r'^\[({})\]: ==> ({})$'.format(comp_id, st_id)
        export_pattern = r'^\[({})\]: ({}) ==> $'.format(comp_id, st_id)

        global_match = re.match(gbl_pattern, value, flags=re.I)
        local_match = re.match(lcl_pattern, value, flags=re.I)
        import_match = re.match(import_pattern, value, flags=re.I)
        export_match = re.match(export_pattern, value, flags=re.I)

        if global_match:
            if global_match.group(1) in objects[Compartment]:
                global_comp = objects[Compartment][global_match.group(1)]
            else:
                global_comp = None
                errors.append('Undefined compartment "{}"'.format(
                    global_match.group(1)))
            lhs = global_match.group(11)
            rhs = global_match.group(41)

        elif local_match:
            global_comp = None
            lhs = local_match.group(1)
            rhs = local_match.group(49)

        elif import_match:
            if import_match.group(1) in objects[Compartment]:
                global_comp = objects[Compartment][import_match.group(1)]
            else:
                global_comp = None
                errors.append('Undefined compartment "{}"'.format(
                    import_match.group(1)))
            lhs = None
            rhs = import_match.group(11) #todo
                
        elif export_match:
            if export_match.group(1) in objects[Compartment]:
                global_comp = objects[Compartment][export_match.group(1)]
            else:
                global_comp = None
                errors.append('Undefined compartment "{}"'.format(
                    export_match.group(1)))
            lhs = export_match.group(11) #todo
            rhs = None

        else:
            return (None, InvalidAttribute(self, ['Incorrectly formatted participants: {}'.format(value)]))

        lhs_parts = []
        rhs_parts = []
        if lhs:
            lhs_parts, lhs_errors = self.deserialize_side(
                -1., lhs, objects, global_comp)
            errors.extend(lhs_errors)
        if rhs:    
            rhs_parts, rhs_errors = self.deserialize_side(
                1., rhs, objects, global_comp)
            errors.extend(rhs_errors)
        parts = lhs_parts + rhs_parts       

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

        st_id = SpeciesType.id.pattern[1:-1]
        comp_id = Compartment.id.pattern[1:-1]
        pattern = r'(\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\) )*({})(\[({})\])*'.format(st_id, comp_id)
        i_st = 4
        i_comp = 15
        for part in re.findall(pattern, value, flags=re.I):
            part_errors = []

            species_type = None
            for species_type_cls in get_subclasses(SpeciesType):
                if species_type_cls in objects and part[i_st] in objects[species_type_cls]:
                    species_type = objects[species_type_cls][part[i_st]]
                    break
            if not species_type:
                part_errors.append(
                    'Undefined species type "{}"'.format(part[i_st]))

            if global_comp:
                compartment = global_comp
            elif part[i_comp] in objects[Compartment]:
                compartment = objects[Compartment][part[i_comp]]
            else:
                part_errors.append(
                    'Undefined compartment "{}"'.format(part[i_comp]))

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


class Identifier(obj_tables.Model):
    """ Reference to an entity in an external namespace

    Attributes:
        namespace (:obj:`str`): namespace
        id (:obj:`str`): identifier within the namespace

    Related attributes:
        compartments (:obj:`list` of :obj:`Compartment`): compartments
        species_types (:obj:`list` of :obj:`SpeciesType`): species_types
        concentrations (:obj:`list` of :obj:`Concentration`): concentrations
        loci (:obj:`list` of :obj:`PolymerLocus`): loci
        properties (:obj:`list` of :obj:`SpeciesTypeProperty`): species type properties
        reactions (:obj:`list` of :obj:`Reaction`): reactions
        rate_laws (:obj:`list` of :obj:`RateLaw`): rate_laws
        observables (:obj:`list` of :obj:`Observable`): observables
    """
    namespace = obj_tables.RegexAttribute(pattern=r'^[^ \:,]+$')
    id = obj_tables.RegexAttribute(pattern=r'^[^ \:,]+$')

    class Meta(obj_tables.Model.Meta):
        attribute_order = ('namespace', 'id')
        table_format = TableFormat.cell
        ordering = ('namespace', 'id')

    @staticmethod
    def _serialize(namespace, id):
        """ Generate string representation

        Args:
            namespace (:obj:`str`): namespace
            id (:obj:`str`): identifier within the namespace

        Returns:
            :obj:`str`: value of primary attribute
        """
        return '{}:{}'.format(namespace, id)

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: value of primary attribute
        """
        return self._serialize(self.namespace, self.id)


class KnowledgeBaseObject(obj_tables.Model):
    """ Knowledge of a biological entity

    Attributes:
        id (:obj:`str`): identifier
        name (:obj:`str`): name
        synonyms (:obj:`str`): synonyms
        comments (:obj:`str`): comments

    Related Attributes:
        time_course_evidences (:obj:`list` of :obj:`TimeCourseEvidence): list of wc_core.TimeCourseEvidence
    """

    id = obj_tables.StringAttribute(primary=True, unique=True)
    name = obj_tables.StringAttribute()
    synonyms = obj_tables.StringAttribute()
    comments = obj_tables.LongStringAttribute()

    def get_nested_metadata(self):
        """ Returns a list of wc_kb.core.Reference / wc_kb.core.DatabaseReference / wc_kb.core.Comments objects that
            appear in the object's wc_kb.core.Evidence and the associated wc_kb.core.Experiment

        Returns:
            id (:obj:`list` of :obj:`Reference`): references
        """

        metadataObjs = {self.id:[], SpeciesTypeProperty:[], Evidence:[], Experiment:[]}
        metadataObjs = self._append_metadata_entries(key=self.id, metadataObjs=metadataObjs)
        metadataObjs = self._parse_EviNExperiment(metadataObjs)

        if hasattr(self,'properties') and self.properties is not None:
            for property in self.properties:
                metadataObjs = property._append_metadata_entries(key=SpeciesTypeProperty, metadataObjs=metadataObjs)
                metadataObjs = property._parse_EviNExperiment(metadataObjs)

            return metadataObjs

    def _parse_EviNExperiment(self, metadataObjs):
        if hasattr(self, 'evidence'):
            for evidence in self.evidence:
                metadataObjs = evidence._append_metadata_entries(key=Evidence, metadataObjs=metadataObjs)
                if evidence.experiment is not None:
                    metadataObjs = evidence.experiment._append_metadata_entries(key=Experiment, metadataObjs=metadataObjs)

        return metadataObjs

    def _append_metadata_entries(self, key, metadataObjs):
        """ Appends wc_kb.core.Reference / wc_kb.core.DatabaseReference / wc_kb.core.Comments objects
            to metadataObjs list

            Input:
                obj(:obj:`obj_tables.Model`): model object

            Return:
                metadataObjs (:obj:`list` of :obj:`Reference` / :obj:`Evidence` / :obj:`Comments`): list of metadata objects
        """

        if self.references!=[]:
            for reference in self.references:
                metadataObjs[key].append(reference)

        if self.identifiers!=[]:
            for database_reference in self.identifiers:
                metadataObjs[key].append(database_reference)

        if self.comments!='':
            metadataObjs[key].append(self.comments)

        return metadataObjs


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
    translation_table = obj_tables.IntegerAttribute()
    version = RegexAttribute(
        min_length=1, pattern=r'^[0-9]+\.[0-9+]\.[0-9]+[0-9a-z]*$', flags=re.I)
    url = obj_tables.StringAttribute(verbose_name='URL')
    branch = obj_tables.StringAttribute()
    revision = obj_tables.StringAttribute()
    wc_kb_version = RegexAttribute(min_length=1, pattern=r'^[0-9]+\.[0-9+]\.[0-9]+[0-9a-z]*$', flags=re.I,
                                   default=wc_kb_version, verbose_name='wc_kb version')

    class Meta(obj_tables.Model.Meta):
        verbose_name = 'KB'
        description = 'Knowledge base'
        attribute_order = ('id', 'name', 'translation_table', 'version',
                           'url', 'branch', 'revision', 'wc_kb_version', 'comments')
        table_format = obj_tables.TableFormat.column


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
    knowledge_base = obj_tables.OneToOneAttribute(
        KnowledgeBase, related_name='cell')
    taxon = obj_tables.IntegerAttribute()

    class Meta(obj_tables.Model.Meta):
        attribute_order = ('id', 'name', 'taxon', 'comments')
        table_format = obj_tables.TableFormat.column


class Reference(obj_tables.Model):
    """ Reference to the literature

    Attributes:
        id (:obj:`str`): identifier
        name (:obj:`str`): name
        authors (:obj:`str`): authors
        title (:obj:`str`): title
        volume (:obj:`str`): volume
        issue (:obj:`str`): issue
        journal (:obj:`str`): journal
        pages (:obj:``str): pages
        year (:obj:`int`): year
        cell (:obj:`Cell`) : cell
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        comments (:obj:`str`): comments
        type (:obj:`pronto`): type of reference

    Related attributes:
        compartments (:obj:`list` of :obj:`Compartment`): compartments
        species_types (:obj:`list` of :obj:`SpeciesType`): species_types
        concentrations (:obj:`list` of :obj:`Concentration`): concentrations
        loci (:obj:`list` of :obj:`PolymerLocus`): loci
        properties (:obj:`list` of :obj:`SpeciesTypeProperty`): species type properties
        reactions (:obj:`list` of :obj:`Reaction`): reactions
        rate_laws (:obj:`list` of :obj:`RateLaw`): rate_laws
        observables (:obj:`list` of :obj:`Observable`): observables
    """

    id = obj_tables.SlugAttribute(primary=True, unique=True)
    name = obj_tables.StringAttribute()
    authors = obj_tables.LongStringAttribute()
    title = obj_tables.LongStringAttribute()
    volume = obj_tables.StringAttribute()
    issue = obj_tables.StringAttribute()
    journal = obj_tables.StringAttribute()
    pages = obj_tables.StringAttribute()
    year = obj_tables.IntegerAttribute()
    cell = obj_tables.ManyToOneAttribute(Cell, related_name='references')
    identifiers = IdentifierAttribute(related_name='references')
    comments = obj_tables.LongStringAttribute()
    type = obj_tables.ontology.OntologyAttribute(kbOnt,
                                  terms = kbOnt['WC:reference'].rchildren(),
                                  default = kbOnt['WC:article'],
                                  none=True)

    class Meta(obj_tables.Model.Meta):
        attribute_order = ('id', 'name', 'type', 'title', 'authors', 'journal', 'volume', 'issue', 'pages', 'year', 'identifiers', 'comments')


class Compartment(KnowledgeBaseObject):
    """ Knowledge of a subcellular compartment

    Attributes:
        cell (:obj:`Cell`): cell
        volumetric_fraction (:obj:`float`): average volumetric fraction relative to the cell volume
        references (:obj:`list` of :obj:`Reference`): references
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers

    Related attributes:
        reaction_participants (:obj:`list` of :obj:`ReactionParticipant`): reaction participants
    """
    id = obj_tables.SlugAttribute(primary=True, unique=True)
    cell = obj_tables.ManyToOneAttribute(Cell, related_name='compartments')
    volumetric_fraction = obj_tables.FloatAttribute(min=0., max=1.)
    references = obj_tables.ManyToManyAttribute(Reference, related_name='compartments')
    identifiers = IdentifierAttribute(related_name='compartments')

    class Meta(obj_tables.Model.Meta):
        attribute_order = ('id', 'name', 'volumetric_fraction', 'identifiers', 'references', 'comments')


class SpeciesType(six.with_metaclass(obj_tables.abstract.AbstractModelMeta, KnowledgeBaseObject)):
    """ Knowledge of a molecular species

    Attributes:
        cell (:obj:`Cell`): cell
        references (:obj:`list` of :obj:`Reference`): references
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers

    Related attributes:
        reaction_participants (:obj:`list` of :obj:`ReactionParticipant`): reaction participants
    """

    id = obj_tables.SlugAttribute(primary=True, unique=True)
    cell = obj_tables.ManyToOneAttribute(Cell, related_name='species_types')
    references = obj_tables.ManyToManyAttribute(Reference, related_name='species_types')
    identifiers = IdentifierAttribute(related_name='species_types')

    class Meta(obj_tables.Model.Meta):
        attribute_order = ('id', 'name', 'comments', 'references', 'identifiers')

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


class Species(obj_tables.Model):
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

    id = obj_tables.StringAttribute(primary=True, unique=True)
    species_type = ManyToOneAttribute(
        SpeciesType, related_name='species', min_related=1)
    compartment = ManyToOneAttribute(
        Compartment, related_name='species', min_related=1)

    class Meta(obj_tables.Model.Meta):
        attribute_order = ('id', 'species_type', 'compartment')
        frozen_columns = 1
        table_format = TableFormat.cell
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

        pattern = r'^({})\[({})\]$'.format(SpeciesType.id.pattern[1:-1], Compartment.id.pattern[1:-1])
        match = re.match(pattern, value, flags=re.I)
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

            if Compartment in objects and match.group(11) in objects[Compartment]:
                compartment = objects[Compartment][match.group(11)]
            else:
                errors.append(
                    'Compartment "{}" is not defined'.format(match.group(11)))

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
        medium (:obj:`str`): medium
        value (:obj:`float`): value
        units (:obj:`unit_registry.Unit`): units; default units is 'M'
        evidence (:obj:`list` of :obj:`Evidence`): evidence
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
    """

    cell = obj_tables.ManyToOneAttribute(Cell, related_name='concentrations')
    species = OneToOneSpeciesAttribute(related_name='concentrations')
    medium = obj_tables.StringAttribute()
    value = FloatAttribute(min=0)
    units = obj_tables.units.UnitAttribute(unit_registry,
                          choices=(
                              unit_registry.parse_units('molecule'),
                              unit_registry.parse_units('mM'),
                              unit_registry.parse_units('uM'),
                              unit_registry.parse_units('nM'),
                              unit_registry.parse_units('pM'),
                              unit_registry.parse_units('fM'),
                              unit_registry.parse_units('aM')),
                          default=unit_registry.parse_units('M'))
    evidence = obj_tables.OneToManyAttribute('Evidence', related_name='concentrations')
    references = ManyToManyAttribute(Reference, related_name='concentrations')
    identifiers = IdentifierAttribute(related_name='concentrations')

    class Meta(obj_tables.Model.Meta):
        attribute_order = ('id', 'species', 'value', 'units', 'evidence', 'identifiers', 'references', 'comments')
        unique_together = (('species', ), )
        ordering = ('species',)
        frozen_columns = 1

    def serialize(self):
        """ Generate string representation
        Returns:
            :obj:`str`: value of primary attribute
        """
        return 'CONC[{}]'.format(self.species.serialize())


class SpeciesTypeCoefficient(obj_tables.Model):
    """ A tuple of a species type and a coefficient

    Attributes:
        species_type (:obj:`SpeciesType`): species_type
        coefficient (:obj:`float`): coefficient

    Related attributes:
        complex (:obj:`ComplexSpeciesType`): complex
    """

    species_type = ManyToOneAttribute(SpeciesType, related_name='species_type_coefficients')
    coefficient = FloatAttribute(min=0.)

    class Meta(obj_tables.Model.Meta):
        attribute_order = ('species_type', 'coefficient')
        frozen_columns = 1
        table_format = TableFormat.cell
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
        st_id = SpeciesType.id.pattern[1:-1]
        stoch = r'\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\)'
        gbl_part = r'({} )*({})'.format(stoch, st_id)
        gbl_side = r'{}( \+ {})*'.format(gbl_part, gbl_part)
        gbl_pattern = r'^({})$'.format(gbl_side)

        global_match = re.match(gbl_pattern, value, flags=re.I)

        if global_match:
            subunits_str = global_match.group(1)
        else:
            attr = cls.Meta.attributes['species_type']
            return (None, InvalidAttribute(attr, ['Incorrectly formatted participants: {}'.format(value)]))

        for part in re.findall(gbl_part, subunits_str, flags=re.I):

            species_type = None
            for species_type_cls in get_subclasses(SpeciesType):
                if species_type_cls in objects and part[4] in objects[species_type_cls]:
                    species_type = objects[species_type_cls][part[4]]
                    break

            if not species_type:
                errors.append('Undefined species type "{}"'.format(part[4]))

            coefficient = float(part[1] or 1.)

            if not errors:
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


class SpeciesCoefficient(obj_tables.Model):
    """ A tuple of a species and a coefficient

    Attributes:
        species (:obj:`Species`): species
        coefficient (:obj:`float`): coefficient

    Related attributes:
        reaction (:obj:`Reaction`): reaction
    """

    species = ManyToOneAttribute(Species, related_name='species_coefficients')
    coefficient = FloatAttribute(nan=False)

    class Meta(obj_tables.Model.Meta):
        attribute_order = ('species', 'coefficient')
        frozen_columns = 1
        table_format = TableFormat.cell
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

        st_id = SpeciesType.id.pattern[1:-1]
        comp_id = Compartment.id.pattern[1:-1]
        if compartment:
            pattern = r'^(\(((\-?\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\) )*({})$'.format(st_id)
        else:
            pattern = r'^(\(((\-?\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\) )*({}\[{}\])$'.format(st_id, comp_id)

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
    circular = obj_tables.BooleanAttribute()
    double_stranded = obj_tables.BooleanAttribute()

    class Meta(obj_tables.Model.Meta):
        attribute_order = ('id', 'name', 'circular', 'double_stranded',
                           'comments', 'references', 'identifiers')

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
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
    """

    cell = obj_tables.ManyToOneAttribute(Cell, related_name='loci')
    polymer = obj_tables.ManyToOneAttribute(PolymerSpeciesType, related_name='loci')
    start = obj_tables.IntegerAttribute()
    end = obj_tables.IntegerAttribute()
    references = obj_tables.ManyToManyAttribute(Reference, related_name='loci')
    identifiers = IdentifierAttribute(related_name='loci')
    strand = obj_tables.EnumAttribute(
        PolymerStrand, default=PolymerStrand.positive)

    class Meta(obj_tables.Model.Meta):
        attribute_order = ('id', 'name', 'polymer', 'strand', 'start', 'end', 'identifiers', 'references', 'comments')

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

    def get_direction(self):
        """ Returns the direction of the polymer feature defind by its strand and start/end coordinate

            Returns:
                :obj:`PolymerDirection`: direction (in ['forward', 'reverse'])

            Raises:
                :obj::obj:`ValueError`: start and end coordinate of chromosome feature can not be the same
                :obj::obj:`Exception`: strand is not member of PolymerStrand
        """

        if self.start < self.end:
            if self.strand == PolymerStrand.positive:
                return PolymerDirection.forward
            elif self.strand == PolymerStrand.negative:
                return PolymerDirection.reverse
            else:
                raise Exception('Unrecognized polymer strand ({}) found for {}.'.format(self.strand, self.id))

        elif self.start > self.end:
            if self.strand == PolymerStrand.positive:
                return PolymerDirection.reverse
            elif self.strand == PolymerStrand.negative:
                return PolymerDirection.forward
            else:
                raise Exception('Unrecognized polymer strand ({}) found for {}.'.format(self.strand, self.id))

        elif self.start == self.end:
            raise ValueError('Start and end position of chromosome feature can not be the same (Chrom feature id: {}).'.format(self.id))


class ObservableExpression(obj_tables.Model, Expression):
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

    class Meta(obj_tables.Model.Meta, Expression.Meta):
        table_format = TableFormat.cell
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
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers

    Related attributes:
        observable_expressions (:obj:`list` of :obj:`ObservableExpression`): observable expressions
        rate_law_expressions (:obj:`list` of :obj:`RateLawExpression`): rate law expressions
    """

    cell = ManyToOneAttribute(Cell, related_name='observables')
    expression = ExpressionManyToOneAttribute(ObservableExpression, related_name='observable',
                                              min_related=1, min_related_rev=1)
    units = obj_tables.units.UnitAttribute(unit_registry,
                          choices=(unit_registry.parse_units('molecule'),),
                          default=unit_registry.parse_units('molecule'))
    references = obj_tables.ManyToManyAttribute(Reference, related_name='observables')
    identifiers = IdentifierAttribute(related_name='observables')

    class Meta(obj_tables.Model.Meta, ExpressionExpressionTermMeta):
        attribute_order = ('id', 'name', 'expression', 'units', 'identifiers', 'references', 'comments')
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
        evidence (:obj:`list` of :obj:`Evidence`): evidence
        references (:obj:`list` of :obj:`Reference`): references
        identifierss (:obj:`list` of :obj:`DatabaseReference`): reference in external namespaces

    Related attributes:
        rate_law_expressions (:obj:`list` of :obj:`RateLawExpression`): rate law expressions that use a Parameter
    """

    cell = obj_tables.ManyToOneAttribute(Cell, related_name='parameters')
    value = FloatAttribute(min=0)
    error = FloatAttribute(min=0)
    units = obj_tables.units.UnitAttribute(unit_registry, none=True)
    references = obj_tables.ManyToManyAttribute(Reference, related_name='parameters')
    evidence = obj_tables.OneToManyAttribute('Evidence', related_name='parameters')
    identifiers = IdentifierAttribute(related_name='parameters')

    class Meta(obj_tables.Model.Meta):
        attribute_order = ('id', 'name', 'synonyms', 'value', 'units', 'evidence', 'identifiers', 'references', 'comments')
        expression_term_token_pattern = (token.NAME, )


class Validator(obj_tables.Validator):
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
        synonyms (:obj:`str`): synonyms
        type (:obj:`pronto`): type

    """
    synonyms = obj_tables.LongStringAttribute()
    type = obj_tables.ontology.OntologyAttribute(kbOnt,
                                  terms = kbOnt['WC:metabolite'].rchildren(),
                                  none = True)

    class Meta(obj_tables.Model.Meta):
        verbose_name = 'Metabolite'
        attribute_order = ('id', 'name', 'synonyms', 'type', 'identifiers', 'references', 'comments')

    def get_structure(self):
        """ Get the structure

        Returns:
            :obj:`str`: InChI or SMILES structure

        Raises:
            :obj:`ValueError`: if structure has not been provided
        """
        structure = self.properties.get_one(property='structure')
        if structure:
            return structure.get_value()            
        else:    
            raise ValueError('The structure of {} has not been provided'.format(self.id))
                
    def calc_structure(self, ph=7.4, major_tautomer=False, keep_hydrogens=False, dearomatize=False):
        """ Get the major microspecies

        Args:
            pH (:obj:`float`, optional): pH, default is 7.4
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, use the major tautomeric in the calculation
            keep_hydrogens (:obj:`bool`, optional): if :obj:`True`, keep explicity defined hydrogens
            dearomatize (:obj:`bool`, optional): if :obj:`True`, dearomatize molecule

        Returns:
            :obj:`str`: InChI-encoded structure
        """
        structure_str = self.get_structure()
        if 'InChI=' in structure_str:
            return get_major_micro_species(structure_str, 'inchi', 'inchi', 
                ph=ph, major_tautomer=major_tautomer, keep_hydrogens=keep_hydrogens, dearomatize=dearomatize)
        else:
            return get_major_micro_species(structure_str, 'smiles', 'smiles', 
                ph=ph, major_tautomer=major_tautomer, keep_hydrogens=keep_hydrogens, dearomatize=dearomatize)    

    def to_openbabel_mol(self):
        """ Convert species type to an Open Babel molecule

        Returns:
            :obj:`openbabel.OBMol`: Open Babel molecule
        """
        structure_str = self.get_structure()
        structure_type = 'inchi' if 'InChI=' in structure_str else 'smi'
        mol = openbabel.OBMol()
        obConversion = openbabel.OBConversion()
        obConversion.SetInFormat(structure_type)
        obConversion.ReadString(mol, structure_str)

        return mol

    def get_empirical_formula(self):
        """ Get the empirical formula

        Returns:
            :obj:`chem.EmpiricalFormula`: empirical formula
        """
        prop = self.properties.get_one(property='empirical_formula')
        if prop:
            return chem.EmpiricalFormula(prop.get_value())

        return self.calc_empirical_formula()

    def calc_empirical_formula(self):
        """ Calculate the empirical formula

        Returns:
            :obj:`chem.EmpiricalFormula`: empirical formula
        """
        mol = self.to_openbabel_mol()
        return OpenBabelUtils.get_formula(mol)

    def get_charge(self):
        """ Get the charge

        Returns:
            :obj:`int`: charge
        """
        prop = self.properties.get_one(property='charge')
        if prop:
            return prop.get_value()

        return self.calc_charge()

    def calc_charge(self):
        """ Calculate the charge

        Returns:
            :obj:`int`: charge
        """
        mol = self.to_openbabel_mol()
        return mol.GetTotalCharge()

    def get_mol_wt(self):
        """ Get the molecular weight

        Returns:
            :obj:`float`: molecular weight

        Raises:
            :obj:`ValueError`: if there is not enough information to calculate molecular weight
        """
        prop = self.properties.get_one(property='empirical_formula')
        if prop:
            return chem.EmpiricalFormula(prop.get_value()).get_molecular_weight()

        elif self.properties.get_one(property='structure'):
            mol = self.to_openbabel_mol()
            return mol.GetMolWt()
        
        else:
            raise ValueError('Molecular weight cannot be calculated because no structure or '
                'empirical formula has been provided for {}'.format(self.id))
      

class DnaSpeciesType(PolymerSpeciesType):
    """ Knowledge of a DNA species

    Attributes:
        seq_path (:obj:`str`): path to sequence fasta file
        ploidy (:obj:`int`): ploidy
    """

    sequence_path = obj_tables.StringAttribute()
    ploidy = obj_tables.IntegerAttribute(min=0)

    class Meta(obj_tables.Model.Meta):
        verbose_name = 'Chromosome'
        attribute_order = ('id', 'name', 'sequence_path', 'circular', 'double_stranded',
                           'ploidy', 'identifiers', 'references', 'comments')

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
        formation_process (:obj:`pronto`): type of formation process
        subunits (:obj:`list` of :obj:`SpeciesTypeCoefficient`): subunits
        type (:obj:`pronto`): type of complex formation

    """

    subunits = SubunitAttribute(related_name='complexes')
    type = obj_tables.ontology.OntologyAttribute(kbOnt,
                                  terms = kbOnt['WC:complex'].rchildren(),
                                  none=True)
    formation_process  = obj_tables.ontology.OntologyAttribute(kbOnt,
                                  terms = kbOnt['WC:complexFormation'].rchildren(),
                                  none=True)

    class Meta(obj_tables.Model.Meta):
        verbose_name = 'Complex'
        attribute_order = ('id', 'name', 'synonyms', 'type', 'formation_process', 'subunits',
                           'identifiers', 'references', 'comments')

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


class RateLawExpression(obj_tables.Model, Expression):
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

    class Meta(obj_tables.Model.Meta, Expression.Meta):
        attribute_order = ('expression', 'parameters', 'species', 'observables')
        table_format = TableFormat.cell
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
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
    """
    reaction = ManyToOneAttribute('Reaction', related_name='rate_laws')
    expression = ExpressionManyToOneAttribute(RateLawExpression, min_related=1, min_related_rev=1, related_name='rate_laws')
    units = obj_tables.units.UnitAttribute(unit_registry,
                          choices=(unit_registry.parse_units('s^-1'),),
                          default=unit_registry.parse_units('s^-1'))
    references = obj_tables.ManyToManyAttribute(Reference, related_name='rate_laws')
    identifiers = IdentifierAttribute(related_name='rate_laws')
    direction = EnumAttribute(RateLawDirection, default=RateLawDirection.forward)


    class Meta(obj_tables.Model.Meta, ExpressionExpressionTermMeta):
        attribute_order = ('id', 'reaction', 'direction', 'expression', 'units', 'identifiers', 'references', 'comments')
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
        participants (:obj:`list` of :obj:`SpeciesCoefficient`): participants
        reversible (:obj:`boolean`): denotes whether reaction is reversible
        references (:obj:`list` of :obj:`Reference`): references
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        evidence (:obj:`list` of :obj:`Evidence`): evidence
        enzymes (:obj:`list` of :obj:`SpeciesType`): enzymes
        coenzymes (:obj:`list` of :obj:`SpeciesType`): coenzymes
        spontaneous (:obj:`bool`): spontaneity
        parameters (:obj:`Parameter`): parameters
        type (:obj:`pronto`): type

    Related attributes:
        rate_laws (:obj:`list` of :obj:`RateLaw`): rate laws; if present, rate_laws[0] is the forward
            rate law, and rate_laws[1] is the backward rate law
    """

    cell = obj_tables.ManyToOneAttribute(Cell, related_name='reactions')
    participants = ReactionParticipantAttribute(related_name='reactions')
    reversible = obj_tables.BooleanAttribute()
    references = obj_tables.ManyToManyAttribute(Reference, related_name='reactions')
    identifiers = IdentifierAttribute(related_name='reactions')
    evidence = obj_tables.OneToManyAttribute('Evidence', related_name='reactions')
    enzymes = obj_tables.ManyToManyAttribute(SpeciesType, related_name='reactions')
    coenzymes = obj_tables.ManyToManyAttribute(SpeciesType, related_name='reactions')
    spontaneous = obj_tables.BooleanAttribute()
    parameters = obj_tables.OneToManyAttribute('Parameter', related_name='reactions')
    type = obj_tables.ontology.OntologyAttribute(kbOnt,
                                  terms = kbOnt['WC:reaction'].rchildren(),
                                  none=True)

    class Meta(obj_tables.Model.Meta):
        attribute_order = ('id', 'name', 'synonyms', 'type', 'participants', 'enzymes', 'coenzymes',
                           'reversible', 'spontaneous', 'parameters', 'evidence', 'identifiers', 'references', 'comments')

#####################
#####################
# Expansion classes

class ChromosomeFeature(PolymerLocus):
    """ Knowledge of chromosome features

    Attributes:
        cell (:obj:`Cell`): cell
        value (:obj:`float`): value
        error (:obj:`float`): measurement error
        units (:obj:`unit_registry.Unit`): units of value
        references (:obj:`list` of :obj:`Reference`): references
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers

    Related attributes:
        seq_path (:obj:`str`): path to sequence fasta file
        ploidy (:obj:`int`): ploidy
    """

    coordinate = obj_tables.IntegerAttribute(min=0)
    start = obj_tables.IntegerAttribute(min=0)
    end = obj_tables.IntegerAttribute(min=0)
    intensity = obj_tables.FloatAttribute(min=0)
    unit = obj_tables.units.UnitAttribute(unit_registry, none=True)
    polymer = obj_tables.ManyToOneAttribute('DnaSpeciesType', related_name='chromosome_features')
    evidence   = obj_tables.OneToManyAttribute('Evidence', related_name='chromosome_features')
    identifiers = IdentifierAttribute(related_name='chromosome_features')
    references = obj_tables.ManyToManyAttribute('Reference', related_name='chromosome_features')
    type = obj_tables.ontology.OntologyAttribute(kbOnt,
                                  terms = kbOnt['WC:chromosomeFeature'].rchildren(),
                                  none=True)

    class Meta(obj_tables.Model.Meta):
        attribute_order = ('id', 'name', 'type', 'polymer', 'start', 'end',
                            'intensity', 'unit', 'evidence', 'identifiers', 'references', 'comments')
        expression_term_token_pattern = (token.NAME, )

    def get_direction(self):
        """ Returns the direction of chromosome feature

            Returns:
                :obj:`PolymerDirection`: direction (in ['forward', 'reverse'])
        """

        if self.start < self.end:
            return PolymerDirection.forward
        elif self.start > self.end:
            return PolymerDirection.reverse
        elif self.start == self.end:
            raise ValueError('Start and end position of chromosome feature can not be the same (Chrom feature id: {}).'.format(self.id))


class Evidence(KnowledgeBaseObject):
    """ Represents the measurement / observation of a property

        Attributes:
            id (:obj:`str`): identifier
            cell (:obj:`Cell`): cell
            object (:obj:`str`): object
            property (:obj:`str`): property
            value (:obj:`float`): value
            units (:obj:`Units`): units
            identifiers(:obj:`list` of :obj:`Identifier`): identifiers
            references (:obj:`list` of :obj:`Reference`): references
            experiment (:obj:`Experiment`): experiment
            comments(:obj:`str`): comments

        Related attributes:

    """

    cell = obj_tables.ManyToOneAttribute('Cell', related_name='evidence')
    object   =  obj_tables.StringAttribute()
    property = obj_tables.StringAttribute()
    value = obj_tables.FloatAttribute()
    units = obj_tables.units.UnitAttribute(unit_registry, none=True) # False allows None units
    identifiers = IdentifierAttribute(related_name='evidence')
    references = obj_tables.ManyToManyAttribute('Reference', related_name='evidence')
    experiment = obj_tables.ManyToOneAttribute('Experiment', related_name ='evidence')
    comments = obj_tables.LongStringAttribute()

    class Meta(obj_tables.Model.Meta):
        attribute_order = ('id', 'cell', 'object', 'property', 'value', 'units', 'experiment', 'identifiers', 'references', 'comments')


class TimeCourseEvidence(Evidence):
    """ Represents the measurement / observation of a property. Supports time resolved
    and non-time resolved data.

    Attributes:
        id (:obj:`str`): identifier
        cell (:obj:`Cell`): cell
        object_tce (:obj:`KnowledgeBaseObject`): object
        property (:obj:`str`): property
        value_tce (:obj:`float`): value
        values_unit_tce (:obj:`Units`): units
        times (:obj:`ndarray`): time
        times_unit (:obj:`Units`): units
        identifiers(:obj:`list` of :obj:`Identifier`): identifiers
        references (:obj:`list` of :obj:`Reference`): references
        experiment (:obj:`Experiment`): experiment
        comments(:obj:`str`): comments

    Related attributes:

    """
    print('reading_object_tce')
    object_tce = obj_tables.ManyToOneAttribute('Observable', related_name='time_course_evidences') # Todo: this
    values_tce = obj_tables.obj_math.NumpyArrayAttribute()
    values_unit_tce = obj_tables.units.UnitAttribute(unit_registry, none=False)
    times = obj_tables.obj_math.NumpyArrayAttribute()
    times_unit = obj_tables.units.UnitAttribute(unit_registry, 
        choices=(unit_registry.parse_units('s'),), none=False)

    class Meta(obj_tables.Model.Meta):
        attribute_order = ('id', 'cell', 'object_tce', 'property', 'values_tce',
            'values_unit_tce', 'times', 'times_unit', 'experiment', 'identifiers',
            'references', 'comments')

    # This validation is partially obsololete, as I cannot set object_tce to any KnowledgeBaseObject:
    # Relationships to `KnowledgeBase` from classes other than `Cell` are prohibited: KnowledgeBaseObject.time_course_evidences to KnowledgeBase
    # obj_tables.ManyToOneAttribute has no method that returns its length or shape (see below)
    # obj_tables.ManyToOneAttribute can read in strings such as 1, 2, blah
    # In this specifice case we should allow only floats and np.nans in obj_tables.ManyToOneAttribute.
    '''def validate(self):
                    """ Determine if the object is valid
                    Returns:
                        :obj:`InvalidObject` or None: `None` if the object is valid,
                            otherwise return a list of errors as an instance of `InvalidObject`
                    """
                    errors = []
            
                    # attributes
                    for attr_name, attr in self.Meta.attributes.items():
                        error = attr.validate(self, getattr(self, attr_name))
                        if error:
                            errors.append(error)
            
                    # related attributes
                    for attr_name, attr in self.Meta.related_attributes.items():
                        if attr.related_name:
                            error = attr.related_validate(self, getattr(self, attr.related_name))
                            if error:
                                errors.append(error)
                        # Specifically checking if obcect_tce is instance of `Species` or `Observable`
                        if attr_name == 'object_tce':
                            if not (isinstance(attr, Species) or isinstance(attr, Observable)):
                                errors.append('{} must be wc_kb.core.Species or wc_kb.core.Observable'.format(attr_name))
            
                    # Checking if wc_kb.core.values_tce and wc_kb.core.times are of same length
                    print(self.Meta.attributes['times'].__dir__())
                    if not len(self.Meta.attributes['times'].deserialize) != len(self.Meta.attributes['values_tce']):
                        errors.append('wc_kb.core.values_tce and wc_kb.core.times must be of same length')
            
                    if errors:
                        return InvalidObject(self, errors)
                    return None'''


class Experiment(KnowledgeBaseObject):
    """ Represents an experiment in which a property was measured

        Attributes:
            id (:obj:`str`): identifier
            species (:obj:`str`): species
            genetic_variant (:obj:`str`): genetic_variant
            external_media (:obj:`str`): external_media
            temperature (:obj:`float`): temperature
            temperature_units (:obj:`Units`): temperature_units
            ph (:obj:`float`): pH
            experiment_design (:obj:`str`): experimental design
            measurement_technology (:obj:`str`): measurement technology
            analysis_type (:obj:`str`): analysis type
            identifiers(:obj:`list` of :obj:`Identifier`): identifiers
            references (:obj:`list` of :obj:`Reference`): references
            comments(:obj:`str`): comments

        Related attributes:
    """

    species = obj_tables.StringAttribute() # Todo: Is this a molecular species or a core.Cell?
    genetic_variant = obj_tables.StringAttribute()
    external_media  = obj_tables.StringAttribute()
    temperature	= obj_tables.FloatAttribute()
    temperature_units = obj_tables.units.UnitAttribute(unit_registry,
                        choices=(unit_registry.parse_units('F'),
                                 unit_registry.parse_units('C'),
                                 unit_registry.parse_units('K')),
                        default= unit_registry.parse_units('C'))
    ph = obj_tables.FloatAttribute()
    experiment_design = obj_tables.StringAttribute()
    measurement_technology = obj_tables.StringAttribute()
    analysis_type = obj_tables.StringAttribute()
    identifiers = IdentifierAttribute(related_name='experiments')
    references = obj_tables.ManyToManyAttribute('Reference', related_name='experiment')
    comments = obj_tables.LongStringAttribute()

    class Meta(obj_tables.Model.Meta):
        attribute_order = ('id', 'experiment_design', 'measurement_technology', 'analysis_type', 'species', 'genetic_variant', 'external_media',
                           'temperature', 'temperature_units', 'ph', 'identifiers', 'references', 'comments')

class TimeCourseExperiment(Experiment):
    """ Represents an experiment in which a property was measured

        Attributes:
            id (:obj:`str`): identifier
            species (:obj:`str`): species
            genetic_variant (:obj:`str`): genetic_variant
            external_media (:obj:`str`): external_media
            temperature (:obj:`float`): temperature
            temperature_units (:obj:`Units`): temperature_units
            ph (:obj:`float`): pH
            times (:obj:`ndarray`): time
            times_unit (:obj:`Units`): units
            objects_tce (:obj:`list` of :obj:`Observable`): list of observables that were perturbed
            objects_tce_units (:obj:`list` of :obj:`Units`): list of units of observables
            experiment_design_tce (:obj:`ndarray`): objects_tce x times ndarray of perturbation courses
            measurement_technology (:obj:`str`): measurement technology
            analysis_type (:obj:`str`): analysis type
            identifiers(:obj:`list` of :obj:`Identifier`): identifiers
            references (:obj:`list` of :obj:`Reference`): references
            comments(:obj:`str`): comments

        Related attributes:
    """

    objects_tce = ManyToManyAttribute('Observable', related_name='time_course_experiments') # Todo: make sure this list is ordered
    objects_tce_units = obj_tables.units.UnitAttribute(unit_registry, 
        choices=(unit_registry.parse_units('s'),), none=False) # Todo: allow this to be a list or ndarray attribute
    times = obj_tables.obj_math.NumpyArrayAttribute()
    times_unit = obj_tables.units.UnitAttribute(unit_registry, 
        choices=(unit_registry.parse_units('s'),), none=False)
    experiment_design_tce = obj_tables.obj_math.NumpyArrayAttribute()
    

    class Meta(obj_tables.Model.Meta):
        attribute_order = ('id', 'times', 'times_unit', 'objects_tce', 'objects_tce_units',
            'experiment_design_tce', 'measurement_technology', 'analysis_type', 'species',
            'genetic_variant', 'external_media', 'temperature', 'temperature_units', 'ph',
            'identifiers', 'references', 'comments')

class SpeciesTypeProperty(KnowledgeBaseObject):
    """ Knowledge of the properties of species types

        Attributes:
            species_type (:obj:`SpeciesType`): species type
            property (:obj:`str`): name of property
            units (:obj:`unit_registry`): units
            value (:obj:`str`): value
            identifiers (:obj:`list` of :obj:`Identifier`): identifiers
            references (:obj:`list` of :obj:`Reference`): references
            evidence (:obj:`list` of :obj:`Evidence`): evidence
            value_type (:obj:`pronto`): value type
    """
    species_type = ManyToOneAttribute(SpeciesType, related_name='properties') #Do we want min_related=1?
    property = obj_tables.StringAttribute()
    units = obj_tables.units.UnitAttribute(unit_registry, none=True)
    value = obj_tables.LongStringAttribute()
    identifiers = IdentifierAttribute(related_name='properties')
    references = ManyToManyAttribute(Reference, related_name='properties')
    evidence = obj_tables.OneToManyAttribute(Evidence, related_name='properties')
    value_type = obj_tables.ontology.OntologyAttribute(kbOnt,
                                terms = kbOnt['WC:valueType'].rchildren(),
                                default = kbOnt['WC:float'],
                                none=False)

    class Meta(obj_tables.Model.Meta):
        verbose_name_plural = 'Species type properties'
        unique_together = (('species_type', 'property', ), )
        attribute_order = ('id', 'species_type', 'property', 'value', 'value_type', 'units', 'evidence',
                           'identifiers', 'references', 'comments')

    def gen_id(self):
        """ Generate id
        Returns:
            :obj:`str`: identifier
        """
        return 'PROP({}:{})'.format(self.species_type.id, self.property)

    def get_value(self):
        """ SpeciesType property values are stored as strings, this function returns the value as the correct type. """

        if self.value == '':
            return None

        if are_terms_equivalent(self.value_type, kbOnt['WC:boolean']):
            return bool(self.value)
        elif are_terms_equivalent(self.value_type, kbOnt['WC:string']):
            return self.value
        elif are_terms_equivalent(self.value_type, kbOnt['WC:integer']):
            return int(self.value)
        elif are_terms_equivalent(self.value_type, kbOnt['WC:float']):
            return float(self.value)
        else:
            raise ValueError('SpeciesTypeProperty "{}" has unexpected value type "{}".'.format(self.id, self.value_type))
