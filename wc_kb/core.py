""" Schema to represent a knowledge base to build models

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-02-07
:Copyright: 2018, Karr Lab
:License: MIT
"""

from natsort import natsorted, ns
from wc_utils.util import chem
import abc
import Bio.Seq
import Bio.SeqUtils
import enum
import math
import obj_model.abstract
import obj_model.core
import obj_model.extra_attributes

from six import with_metaclass, string_types
from obj_model import (BooleanAttribute, EnumAttribute, FloatAttribute, IntegerAttribute, PositiveIntegerAttribute,
                       RegexAttribute, SlugAttribute, StringAttribute, LongStringAttribute, UrlAttribute,
                       OneToOneAttribute, ManyToOneAttribute, ManyToManyAttribute,
                       InvalidModel, InvalidObject, InvalidAttribute, TabularOrientation)
import openbabel
import pkg_resources
import six
import re

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


class ComplexType(enum.Enum):
    """ Type of gene """
    tRnaSynthClassII = 0
    FattyAcylAcp = 1


class ComplexFormationType(enum.Enum):
    """ Type of gene """
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


#####################
#####################
# Base classes


class KnowledgeBaseObject(obj_model.core.Model):
    """ Knowlege of a biological entity

    Attributes:
        id (:obj:`str`): identifier
        name (:obj:`str`): name
        comments (:obj:`str`): comments
    """
    id = obj_model.core.StringAttribute(primary=True, unique=True)
    name = obj_model.core.StringAttribute()
    comments = obj_model.core.StringAttribute()


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
    translation_table = obj_model.core.IntegerAttribute()
    version = RegexAttribute(min_length=1, pattern='^[0-9]+\.[0-9+]\.[0-9]+', flags=re.I)
    url = obj_model.core.StringAttribute(verbose_name='URL')
    branch = obj_model.core.StringAttribute()
    revision = obj_model.core.StringAttribute()
    wc_kb_version = RegexAttribute(min_length=1, pattern='^[0-9]+\.[0-9+]\.[0-9]+', flags=re.I,
                                   default=wc_kb_version, verbose_name='wc_kb version')

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'name', 'translation_table', 'version', 'url', 'branch', 'revision', 'wc_kb_version', 'comments')
        tabular_orientation = obj_model.core.TabularOrientation.column


class Cell(KnowledgeBaseObject):
    """ Knowledge of a cell

    Attributes:
        knowledge_base (:obj:`KnowledgeBase`): knowledge base

    Related attributes:
        compartments (:obj:`list` of :obj:`Compartment`): compartments
        species_types (:obj:`list` of :obj:`SpeciesType`): species types
        loci (:obj:`list` of :obj:`PolymerLocus`): locus
        reactions (:obj:`list` of :obj:`Reaction`): reactions
    """
    knowledge_base = obj_model.core.OneToOneAttribute(KnowledgeBase, related_name='cell')

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'name', 'comments')
        tabular_orientation = obj_model.core.TabularOrientation.column


class Compartment(KnowledgeBaseObject):
    """ Knowledge of a subcellular compartment

    Attributes:
        cell (:obj:`Cell`): cell
        volume (:obj:`float`): average volume at the begining of the cell cycle (L)

    Related attributes:
        reaction_participants (:obj:`list` of :obj:`ReactionParticipant`): reaction participants
    """
    cell = obj_model.core.ManyToOneAttribute(Cell, related_name='compartments')
    volume = obj_model.core.FloatAttribute(min=0)

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'name', 'volume', 'comments')


class SpeciesType(six.with_metaclass(obj_model.abstract.AbstractModelMeta, KnowledgeBaseObject)):
    """ Knowledge of a molecular species

    Attributes:
        cell (:obj:`Cell`): cell
        concentration (:obj:`float`): concentration (M)
        half_life  (:obj:`float`): half life (s)

    Related attributes:
        reaction_participants (:obj:`list` of :obj:`ReactionParticipant`): reaction participants
    """

    cell = obj_model.core.ManyToOneAttribute(Cell, related_name='species_types')
    concentration = obj_model.core.FloatAttribute(min=0)
    half_life = obj_model.core.FloatAttribute(min=0)

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'name', 'concentration', 'half_life', 'comments')

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
        species_coefficients (:obj:`list` of `SpeciesCoefficient`): participations in reactions and observables
        rate_law_equations (:obj:`RateLawEquation`): rate law equations
    """

    species_type = ManyToOneAttribute(SpeciesType, related_name='species', min_related=1)
    compartment = ManyToOneAttribute(Compartment, related_name='species', min_related=1)

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
        if isinstance(species_type, SpeciesType) and isinstance(compartment, Compartment):
            species_type_id = species_type.get_primary_attribute()
            compartment_id = compartment.get_primary_attribute()
        elif isinstance(species_type, string_types) and isinstance(compartment, string_types):
            species_type_id = species_type
            compartment_id = compartment
        else:
            raise ValueError("gen_id: incorrect parameter types: {}, {}".format(species_type, compartment))
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

        match = re.match('^([a-z][a-z0-9_]*)\[([a-z][a-z0-9_]*)\]$', value, flags=re.I)
        if match:
            errors = []
            if match.group(1) in objects[RnaSpeciesType]:
                species_type = objects[RnaSpeciesType][match.group(1)]
            elif match.group(1) in objects[ProteinSpeciesType]:
                species_type = objects[ProteinSpeciesType][match.group(1)]
            elif match.group(1) in objects[MetaboliteSpeciesType]:
                species_type = objects[MetaboliteSpeciesType][match.group(1)]
            elif match.group(1) in objects[ComplexSpeciesType]:
                species_type = objects[ComplexSpeciesType][match.group(1)]
            else:
                errors.append('Species type "{}" is not defined'.format(match.group(1)))

            if match.group(2) in objects[Compartment]:
                compartment = objects[Compartment][match.group(2)]
            else:
                errors.append('Compartment "{}" is not defined'.format(match.group(2)))

            if errors:
                return (None, InvalidAttribute(attribute, errors))
            else:
                obj = cls(species_type=species_type, compartment=compartment)
                if cls not in objects:
                    objects[cls] = {}
                objects[cls][obj.serialize()] = obj
                return (obj, None)

        return (None, InvalidAttribute(attribute, ['Invalid species']))


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
            :obj:`str`: simple Python representation
        """
        return self._serialize(self.species, self.coefficient,
                               show_compartment=show_compartment, show_coefficient_sign=show_coefficient_sign)

    @staticmethod
    def _serialize(species, coefficient, show_compartment=True, show_coefficient_sign=True):
        """ Serialize values

        Args:
            species (:obj:`Species`): species
            coefficient (:obj:`float`): coefficient
            show_compartment (:obj:`bool`, optional): if true, show compartment
            show_coefficient_sign (:obj:`bool`, optional): if true, show coefficient sign

        Returns:
            :obj:`str`: simple Python representation
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
            pattern = '^(\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\) )*([a-z][a-z0-9_]*)$'
        else:
            pattern = '^(\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\) )*([a-z][a-z0-9_]*\[[a-z][a-z0-9_]*\])$'

        match = re.match(pattern, value, flags=re.I)
        if match:
            errors = []

            coefficient = float(match.group(2) or 1.)

            if compartment:
                species_id = Species.gen_id(match.group(5), compartment.get_primary_attribute())
            else:
                species_id = match.group(5)

            species, error = Species.deserialize(attribute, species_id, objects)
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
    circular = obj_model.core.BooleanAttribute()
    double_stranded = obj_model.core.BooleanAttribute()

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'name', 'circular', 'double_stranded', 'concentration', 'half_life', 'comments')

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
            :obj:`ValueError`: if the polymer is linear and the start or end coordinates are less than 1 or greater than the length
                of the sequence
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
            pos_seq = seq[start:] + str(seq) * (int(math.floor(end / seq_len)) - 1) + seq[0:end % seq_len]

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
    """

    cell = obj_model.core.ManyToOneAttribute(Cell, related_name='loci')
    polymer = obj_model.core.ManyToOneAttribute(PolymerSpeciesType, related_name='loci')
    strand = obj_model.core.EnumAttribute(PolymerStrand, default=PolymerStrand.positive)
    start = obj_model.core.IntegerAttribute()
    end = obj_model.core.IntegerAttribute()

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'name', 'polymer', 'strand', 'start', 'end', 'comments')

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


#####################
#####################
# Species types


class MetaboliteSpeciesType(SpeciesType):
    """ Knowledge of a metabolite

    Attributes:
        structure (:obj:`str`): InChI-encoded structure
    """
    structure = obj_model.core.StringAttribute()

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'name', 'structure', 'concentration', 'half_life', 'comments')

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
        return chem.EmpiricalFormula(mol.GetFormula().rstrip('-'))

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
    """

    seq = obj_model.extra_attributes.BioSeqAttribute(verbose_name='Sequence')

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'name', 'seq', 'circular', 'double_stranded', 'concentration', 'half_life', 'comments')
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
        formula.H = 12 * n_a + 12 * n_c + 12 * n_g + 13 * n_t - (l - 1 + self.circular) * (1 + self.double_stranded)
        formula.N = 5 * n_a + 3 * n_c + 5 * n_g + 2 * n_t
        formula.O = 6 * n_a + 7 * n_c + 7 * n_g + 8 * n_t - (l - 1 + self.circular) * (1 + self.double_stranded)
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


class RnaSpeciesType(PolymerSpeciesType):
    """ Knowledge of an RNA species

    Attributes:
        transcription_units (:obj:`list` of :obj:`TranscriptionUnitLocus`): transcription units
        type (:obj:`RnaType`): type

    Related attributes:
        proteins (:obj:`list` of :obj:`ProteinSpeciesType`): protein(s)
    """

    transcription_units = obj_model.core.ManyToManyAttribute('TranscriptionUnitLocus', related_name='rna')
    type = obj_model.core.EnumAttribute(RnaType)

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'name', 'type', 'transcription_units',
                           'circular', 'double_stranded', 'concentration', 'half_life', 'comments')

    def get_seq(self):
        """ Get the sequence

        Returns:
            :obj:`Bio.Seq.Seq`: sequence
        """
        tu_start = self.transcription_units[0].start
        tu_end = self.transcription_units[0].end
        dna_seq = self.transcription_units[0].polymer.get_subseq(start=tu_start, end=tu_end)
        return dna_seq.transcribe()

    def get_empirical_formula(self):
        """ Get the empirical formula for an RNA transcript with

        * 5' monophosphate
        * Deprotonated phosphate oxygens

        :math:`N_A * AMP + N_C * CMP + N_G * GMP + N_U * UMP - (L-1) * OH`

        Returns:
           :obj:`chem.EmpiricalFormula`: empirical formula
        """
        seq = self.get_seq()
        n_a = seq.count('A')
        n_c = seq.count('C')
        n_g = seq.count('G')
        n_u = seq.count('U')
        l = len(seq)

        formula = chem.EmpiricalFormula()
        formula.C = 10 * n_a + 9 * n_c + 10 * n_g + 9 * n_u
        formula.H = 12 * n_a + 12 * n_c + 12 * n_g + 11 * n_u - (l - 1)
        formula.N = 5 * n_a + 3 * n_c + 5 * n_g + 2 * n_u
        formula.O = 7 * n_a + 8 * n_c + 8 * n_g + 9 * n_u - (l - 1)
        formula.P = n_a + n_c + n_g + n_u

        return formula

    def get_charge(self):
        """ Get the charge for a transcript with

        * 5' monophosphate
        * Deprotonated phosphate oxygens

        :math:`-L - 1`

        Returns:
           :obj:`int`: charge
        """
        return -self.get_len() - 1

    def get_mol_wt(self):
        """ Get the molecular weight for a transcript with

        * 5' monophosphate
        * Deprotonated phosphate oxygens

        Returns:
            :obj:`float`: molecular weight (Da)
        """
        return self.get_empirical_formula().get_molecular_weight()


class ProteinSpeciesType(PolymerSpeciesType):
    """ Knowledge of a protein monomer

    Attributes:
        gene (:obj:`GeneLocus`): gene
        rna (:obj:`RnaSpeciesType`): rna
    """

    gene = obj_model.core.ManyToOneAttribute('GeneLocus', related_name='proteins')
    rna = obj_model.core.ManyToOneAttribute('RnaSpeciesType', related_name='proteins')

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'name', 'gene', 'rna', 'circular', 'double_stranded', 'concentration', 'half_life', 'comments')

    def get_seq(self):
        """ Get the sequence

        Returns:
            :obj:`Bio.Seq.Seq`: sequence
        """
        trans_table = self.gene.transcription_units[0].polymer.cell.knowledge_base.translation_table
        return self.gene.get_seq().translate(trans_table)

    def get_empirical_formula(self):
        """ Get the empirical formula

        Returns:
            :obj:`chem.EmpiricalFormula`: empirical formula
        """

        seq = self.get_seq()
        l = len(seq)

        n_a = seq.count('A')  # Ala: Alanine (C3 H7 N O2)
        n_r = seq.count('R')  # Arg: Arginine (C6 H14 N4 O2)
        n_n = seq.count('N')  # Asn: Asparagine (C4 H8 N2 O3)
        n_d = seq.count('D')  # Asp: Aspartic acid (C4 H7 N O4)
        n_c = seq.count('C')  # Cys: Cysteine (C3 H7 N O2 S)

        n_q = seq.count('Q')  # Gln: Glutamine (C5 H10 N2 O3)
        n_e = seq.count('E')  # Glu: Glutamic acid (C5 H9 N O4)
        n_g = seq.count('G')  # Gly: Glycine (C2 H5 N O2)
        n_h = seq.count('H')  # His: Histidine (C6 H9 N3 O2)
        n_i = seq.count('I')  # Ile: Isoleucine (C6 H13 N O2)

        n_l = seq.count('L')  # Leu: Leucine (C6 H13 N O2)
        n_k = seq.count('K')  # Lys: Lysine (C6 H14 N2 O2)
        n_m = seq.count('M')  # Met: Methionine (C5 H11 N O2 S)
        n_f = seq.count('F')  # Phe: Phenylalanine (C9 H11 N O2)
        n_p = seq.count('P')  # Pro: Proline (C5 H9 N O2)

        n_s = seq.count('S')  # Ser: Serine (C3 H7 N O3)
        n_t = seq.count('T')  # Thr: Threonine (C4 H9 N O3)
        n_w = seq.count('W')  # Trp: Tryptophan (C11 H12 N2 O2)
        n_y = seq.count('Y')  # Tyr: Tyrosine (C9 H11 N O3)
        n_v = seq.count('V')  # Val: Valine (C5 H11 N O2)

        formula = chem.EmpiricalFormula()

        formula.C = 3 * n_a + 6 * n_r + 4 * n_n + 4 * n_d + 3 * n_c + \
            5 * n_q + 5 * n_e + 2 * n_g + 6 * n_h + 6 * n_i + \
            6 * n_l + 6 * n_k + 5 * n_m + 9 * n_f + 5 * n_p + \
            3 * n_s + 4 * n_t + 11 * n_w + 9 * n_y + 5 * n_v

        formula.H = 7 * n_a + 14 * n_r + 8 * n_n + 7 * n_d + 7 * n_c + \
            10 * n_q + 9 * n_e + 5 * n_g + 9 * n_h + 13 * n_i + \
            13 * n_l + 14 * n_k + 11 * n_m + 11 * n_f + 9 * n_p + \
            7 * n_s + 9 * n_t + 12 * n_w + 11 * n_y + 11 * n_v - 2 * (l - 1)

        formula.N = 1 * n_a + 4 * n_r + 2 * n_n + 1 * n_d + 1 * n_c + \
            2 * n_q + 1 * n_e + 1 * n_g + 3 * n_h + 1 * n_i + \
            1 * n_l + 2 * n_k + 1 * n_m + 1 * n_f + 1 * n_p + \
            1 * n_s + 1 * n_t + 2 * n_w + 1 * n_y + 1 * n_v

        formula.O = 2 * n_a + 2 * n_r + 3 * n_n + 4 * n_d + 2 * n_c + \
            3 * n_q + 4 * n_e + 2 * n_g + 2 * n_h + 2 * n_i + \
            2 * n_l + 2 * n_k + 2 * n_m + 2 * n_f + 2 * n_p + \
            3 * n_s + 3 * n_t + 2 * n_w + 3 * n_y + 2 * n_v - (l - 1)

        formula.S = n_c + n_m
        return formula

    def get_charge(self):
        """ Get the charge at physiological pH

        Returns:
            :obj:`int`: charge
        """
        seq = self.get_seq()

        n_r = seq.count('R')
        n_h = seq.count('H')
        n_k = seq.count('K')
        n_d = seq.count('D')
        n_e = seq.count('E')

        return (n_r + n_h + n_k) - (n_d + n_e)

    def get_mol_wt(self):
        """ Get the molecular weight

        Returns:
            :obj:`float`: molecular weight
        """
        return self.get_empirical_formula().get_molecular_weight()


class ComplexSpeciesType(SpeciesType):
    """ Knowledge of a protein complexe

    Attributes:
        complex_type (:obj:`ComplexType`): type of complex
        formation_process (:obj:`ComplexFormationType`): type of formation process
        biosynthesis (:obj:`string`): eq governing the formaiton of complex
    """

    complex_type = obj_model.core.StringAttribute()  # EnumAttribute(ComplexType)
    formation_process = obj_model.core.EnumAttribute(ComplexFormationType)
    binding = obj_model.core.StringAttribute()
    region = obj_model.core.StringAttribute()

    formation_reaction = obj_model.core.OneToOneAttribute('Reaction', related_name='complex')
    #subunits = obj_model.core.ManyToManyAttribute('SpeciesCoefficient', related_name='complexes')

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'name', 'formation_process', 'formation_reaction',
                           'complex_type', 'binding', 'region', 'concentration', 'half_life', 'comments')

    def get_subunits(self):
        """ Get the subunit composition of the complex

        Returns:
            `list` of :obj:`SpeciesType`: list of Speciestype objects that compose the complex
        """
        # todo: represent composition directly and eliminate this method
        subunits = []

        for participant in self.formation_reaction.participants:
            if participant.species.species_type.id != self.id:
                for n in range(0, abs(int(participant.coefficient))):
                    subunits.append(participant.species.species_type)

        return subunits

    def get_empirical_formula(self):
        """ Get the empirical formula

        Returns:
            :obj:`chem.EmpiricalFormula`: empirical formula
        """
        # Formula addition
        formula = chem.EmpiricalFormula()
        for subunit in self.get_subunits():
            formula = formula + subunit.get_empirical_formula()

        return formula

    def get_charge(self):
        """ Get the charge at physiological pH

        Returns:
            :obj:`int`: charge
        """
        charge = 0
        for subunit in self.get_subunits():
            charge = charge + subunit.get_charge()

        return charge

    def get_mol_wt(self):
        """ Get the molecular weight

        Returns:
            :obj:`float`: molecular weight
        """
        weight = 0
        for subunit in self.get_subunits():
            weight = weight + subunit.get_mol_wt()

        return weight


#####################
#####################
# Locus types


class PromoterLocus(PolymerLocus):
    """ Knowledge of a promoter for a transcription unit

    Attributes:
        pribnow_start (:obj:`int`): Pribnow box start coordinate
        pribnow_end (:obj:`int`): Pribnow box end coordinate

    Related attributes:
        transcription_units (:obj:`list` of :obj:`TranscriptionUnitLocus`)

    """
    pribnow_start = obj_model.core.IntegerAttribute()
    pribnow_end = obj_model.core.IntegerAttribute()

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'polymer', 'name', 'pribnow_start', 'pribnow_end', 'strand', 'start', 'end', 'comments')


class TranscriptionUnitLocus(PolymerLocus):
    """ Knowledge about an open reading frame

    Attributes:
        promoter (:obj:`PromoterLocus`): promoter controlling the TU
        genes (:obj:`list` of :obj:`GeneLocus`): genes
    """

    promoter = obj_model.core.ManyToOneAttribute('PromoterLocus', related_name='transcription_units')
    genes = obj_model.core.ManyToManyAttribute('GeneLocus', related_name='transcription_units')

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'polymer', 'name', 'strand', 'promoter', 'start', 'end', 'genes', 'comments')

    def get_3_prime(self):
        """ Get the 3' coordinate

        Returns:
            :obj:`int`: 3' coordinate
        """
        if self.strand == PolymerStrand.positive:
            return self.end
        else:
            return self.start

    def get_5_prime(self):
        """ Get the 5' coordinate

        Returns:
            :obj:`int`: 5' coordinate
        """
        if self.strand == PolymerStrand.positive:
            return self.start
        else:
            return self.end


class GeneLocus(PolymerLocus):
    """ Knowledge of a gene

    Attributes:
        symbol (:obj:`str`): symbol

    Related attributes:
        proteins (:obj:`list` of :obj:`ProteinSpeciesType`): protein
    """

    symbol = obj_model.core.StringAttribute()

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'polymer', 'name', 'symbol', 'strand', 'start', 'end', 'comments')


#####################
#####################
# Reactions


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
            participants (:obj:`list` of `SpeciesCoefficient`): Python representation of reaction participants
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
            participants = natsorted(participants, lambda part: part.species.species_type.id, alg=ns.IGNORECASE)
        else:
            participants = natsorted(participants, lambda part: (
                part.species.species_type.id, part.species.compartment.id), alg=ns.IGNORECASE)

        lhs = []
        rhs = []
        for part in participants:
            if part.coefficient < 0:
                lhs.append(part.serialize(show_compartment=global_comp is None, show_coefficient_sign=False))
            elif part.coefficient > 0:
                rhs.append(part.serialize(show_compartment=global_comp is None, show_coefficient_sign=False))

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

        id = '[a-z][a-z0-9_]*'
        stoch = '\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\)'
        gbl_part = '({} )*({})'.format(stoch, id)
        lcl_part = '({} )*({}\[{}\])'.format(stoch, id, id)
        gbl_side = '{}( \+ {})*'.format(gbl_part, gbl_part)
        lcl_side = '{}( \+ {})*'.format(lcl_part, lcl_part)
        gbl_pattern = '^\[({})\]: ({}) ==> ({})$'.format(id, gbl_side, gbl_side)
        lcl_pattern = '^({}) ==> ({})$'.format(lcl_side, lcl_side)

        global_match = re.match(gbl_pattern, value, flags=re.I)
        local_match = re.match(lcl_pattern, value, flags=re.I)

        if global_match:
            if global_match.group(1) in objects[Compartment]:
                global_comp = objects[Compartment][global_match.group(1)]
            else:
                global_comp = None
                errors.append('Undefined compartment "{}"'.format(global_match.group(1)))
            lhs = global_match.group(2)
            rhs = global_match.group(14)

        elif local_match:
            global_comp = None
            lhs = local_match.group(1)
            rhs = local_match.group(13)

        else:
            return (None, InvalidAttribute(self, ['Incorrectly formatted participants: {}'.format(value)]))

        lhs_parts, lhs_errors = self.deserialize_side(-1., lhs, objects, global_comp)
        rhs_parts, rhs_errors = self.deserialize_side(1., rhs, objects, global_comp)

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

        for part in re.findall('(\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\) )*([a-z][a-z0-9_]*)(\[([a-z][a-z0-9_]*)\])*', value, flags=re.I):
            part_errors = []

            if part[4] in objects[RnaSpeciesType]:
                species_type = objects[RnaSpeciesType][part[4]]
            elif part[4] in objects[ProteinSpeciesType]:
                species_type = objects[ProteinSpeciesType][part[4]]
            elif part[4] in objects[MetaboliteSpeciesType]:
                species_type = objects[MetaboliteSpeciesType][part[4]]
            elif part[4] in objects[ComplexSpeciesType]:
                species_type = objects[ComplexSpeciesType][part[4]]
            else:
                part_errors.append('Undefined species type "{}"'.format(part[4]))

            if global_comp:
                compartment = global_comp
            elif part[6] in objects[Compartment]:
                compartment = objects[Compartment][part[6]]
            else:
                part_errors.append('Undefined compartment "{}"'.format(part[6]))

            coefficient = direction * float(part[1] or 1.)

            if part_errors:
                errors += part_errors
            else:
                spec_primary_attribute = Species.gen_id(species_type.get_primary_attribute(),
                                                        compartment.get_primary_attribute())
                species, error = Species.deserialize(self, spec_primary_attribute, objects)
                if error:
                    raise ValueError('Invalid species "{}"'.format(spec_primary_attribute))
                    # pragma: no cover # unreachable due to error checking above

                if coefficient != 0:
                    if SpeciesCoefficient not in objects:
                        objects[SpeciesCoefficient] = {}
                    serialized_value = SpeciesCoefficient._serialize(species, coefficient)
                    if serialized_value in objects[SpeciesCoefficient]:
                        rxn_part = objects[SpeciesCoefficient][serialized_value]
                    else:
                        rxn_part = SpeciesCoefficient(species=species, coefficient=coefficient)
                        objects[SpeciesCoefficient][serialized_value] = rxn_part
                    parts.append(rxn_part)

        return (parts, errors)


class Reaction(KnowledgeBaseObject):
    """ Knowledge of reactions

    Attributes:
        cell (:obj:`Cell`): cell
        participants (:obj:`list` of :obj:`ReactionParticipant`): participants
        v_max (:obj:`float`):V_max value of reaction (unit: mol/L/min)
        k_m (:obj:`float`): K_m value of reaction (unit: mol/L)
        reversible (:obj:`boolean`): denotes whether reaction is reversible
    """

    cell = obj_model.core.ManyToOneAttribute(Cell, related_name='reactions')
    participants = ReactionParticipantAttribute(related_name='reactions')
    v_max = obj_model.core.FloatAttribute(min=0, verbose_name='Vmax')
    k_m = obj_model.core.FloatAttribute(min=0, verbose_name='Km')
    reversible = obj_model.core.BooleanAttribute()

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'name', 'participants', 'v_max', 'k_m', 'reversible', 'comments')
