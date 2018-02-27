""" Read and represent a knowledge base

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-02-07
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_utils.util import chem
import abc
import Bio.Seq
import Bio.SeqUtils
import enum
import math
import obj_model.abstract
import obj_model.core
import obj_model.extra_attributes
import six

PolymerStrand = enum.Enum(value='PolymerStrand', names=[
    ('positive', 1),
    ('+', 1),
    ('negative', -1),
    ('-', -1),
])


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

    Related attributes:
        cell (:obj:`Cell`): cell
    """

    version = obj_model.core.StringAttribute()
    translation_table = obj_model.core.IntegerAttribute()

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'name', 'version', 'translation_table')
        tabular_orientation = obj_model.core.TabularOrientation.column


class Cell(KnowledgeBaseObject):
    """ Knowledge of a cell

    Related attributes:
        species_types (:obj:`list` of :obj:`SpeciesType`): species types
        chromosomes (:obj:`list` of :obj:`Chromosome`): chromosomes
        reactions (:obj:`list` of :obj:`Reaction`): reactions
    """

    knowledge_base = obj_model.core.OneToOneAttribute(KnowledgeBase, related_name='cell')

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'knowledge_base', 'name')
        tabular_orientation = obj_model.core.TabularOrientation.column


class Compartment(KnowledgeBaseObject):
    """ Knowledge of a subcellular compartment

    Attributes:
        cell (:obj:`Cell`): cell
        volume (:obj:`float`): average volume at the begining of the cell cycle (L)
    """
    cell = obj_model.core.ManyToOneAttribute(Cell, related_name='species_types')
    volume = obj_model.core.FloatAttribute(min=0)

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'cell', 'name', 'volume')


class SpeciesType(six.with_metaclass(obj_model.abstract.AbstractModelMeta, KnowledgeBaseObject)):
    """ Knowledge of a molecular species

    Attributes:
        cell (:obj:`Cell`): cell
        concentration (:obj:`float`): concentration (M)
        half_life  (:obj:`float`): half life (s)

    Related attributes:
        reactions (:obj:`Reaction`): reactions
    """

    cell = obj_model.core.ManyToOneAttribute(Cell, related_name='species_types')
    concentration = obj_model.core.FloatAttribute(min=0)
    half_life = obj_model.core.FloatAttribute(min=0)

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'cell', 'name', 'structure', 'concentration', 'half_life')

    @abc.abstractmethod
    def get_structure(self):
        """ Get the structure

        Returns:
            :obj:`str`: structure
        """
        pass

    @abc.abstractmethod
    def get_empirical_formula(self):
        """ Get the empirical formula

        Returns:
            :obj:`str`: empirical formula
        """
        pass

    @abc.abstractmethod
    def get_molecular_weight(self):
        """ Get the molecular weight

        Returns:
            :obj:`float`: molecular weight
        """
        pass

    @abc.abstractmethod
    def get_charge(self):
        """ Get the charge

        Returns:
            :obj:`int`: charge
        """
        pass


class MetaboliteSpeciesType(SpeciesType):
    """ Knowledge of a metabolite

    Attributes:
        structure (:obj:`str`): InChI-encoded structure
    """
    structure = obj_model.core.StringAttribute()

    def get_structure(self):
        """ Get the structure

        Returns:
            :obj:`str`: structure
        """
        return self.structure

    def get_empirical_formula(self):
        """ Get the empirical formula

        Returns:
            :obj:`str`: empirical formula
        """
        pass  # todo calculate the empirical formula

    def get_molecular_weight(self):
        """ Get the molecular weight

        Returns:
            :obj:`float`: molecular weight
        """
        pass  # todo calculate the molecular weight from the structure

    def get_charge(self):
        """ Get the charge

        Returns:
            :obj:`int`: charge
        """
        pass  # todo calculate the charge from the structure


class PolymerSpeciesType(SpeciesType):
    """ Knowledge of a polymer

    Attributes:        
        is_circular (:obj:`bool`): is the polymer circular
        is_double_stranded (:obj:`bool`): is the polymer double stranded

    Related attributes:
        loci (:obj:`list` of :obj:`PolymerLocus`): loci
    """
    is_circular = obj_model.core.BooleanAttribute()
    is_double_stranded = obj_model.core.BooleanAttribute()

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'cell', 'name', 'is_circular', 'is_double_stranded')

    def get_len(self):
        """ Get the polymer length

        Returns:
            :obj:`int`: length
        """
        return len(self.seq)

    def get_subseq(self, start, end, strand=PolymerStrand.positive):
        """ Get a subsequence

        Args:
            start (:obj:`int`): start coordinate (1-indexed)
            end (:obj:`int`): end coordinate (1-indexed)
            strand (:obj:`PolymerStrand`, optional): strand

        Returns:
            :obj:`str`: sequence
        """
        seq_len = self.get_len()

        # convert to zero-based indexing
        start -= 1

        n_wrap = int(math.floor(start / seq_len))
        start = start - seq_len * n_wrap
        end = end - seq_len * n_wrap

        if end <= seq_len:
            pos_seq = self.seq[start:end]
        else:
            pos_seq = self.seq[start:] + str(self.seq) * (int(math.floor(end / seq_len)) - 1) + self.seq[0:end % seq_len]

        if strand == PolymerStrand.positive:
            return pos_seq
        else:
            return pos_seq.reverse_complement()


class DnaSpeciesType(PolymerSpeciesType):
    """ Knowledge of a DNA species 

    Attributes:
        seq (:obj:`Bio.Seq.Seq`): sequence
    """
    seq = obj_model.extra_attributes.BioSeqAttribute()

    def get_structure(self):
        """ Get the structure

        Returns:
            :obj:`Bio.Seq.Seq`: structure
        """
        return self.seq

    def get_empirical_formula(self):
        """ Get the empirical formula

        Returns:
            :obj:`str`: empirical formula
        """
        pass  # todo calculate the empirical formula from the sequence

    def get_molecular_weight(self):
        """ Get the molecular weight

        Returns:
            :obj:`float`: molecular weight
        """
        pass  # todo calculate the molecular weight from the sequence

    def get_charge(self):
        """ Get the charge

        Returns:
            :obj:`int`: charge
        """
        pass  # todo calculate the charge from the sequence


class ChromosomeSpeciesType(DnaSpeciesType):
    """ Knowledge of a chromosome

    Related attributes:
        transcription_units (:obj:`list` of :obj:`TranscriptionUnit`): transcription units
    """
    pass


class RnaType(enum.Enum):
    """ Type of RNA """
    mRna = 0
    rRna = 1
    sRna = 2
    tRna = 3


class RnaSpeciesType(PolymerSpeciesType):
    """ Knowledge of an RNA species

    Attributes:
        transcription_unit (:obj:`TranscriptionUnit`): transcription unit
        type (:obj:`RnaType`): type

    Related attributes:
        genes (:obj:`list` of :obj:`Gene`): genes
    """
    transcription_unit = obj_model.core.ManyToOneAttribute(TranscriptionUnit, related_name='rna')
    type = obj_model.core.EnumAttribute(RnaType)

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'transcription_unit', 'name', 'type',
                           'copy_number', 'half_life')

    def get_seq(self):
        """ Get the sequence 

        Returns:
            :obj:`Bio.Seq.Seq`: sequence
        """
        return self.transcription_unit.get_seq()

    def get_structure(self):
        """ Get the structure

        Returns:
            :obj:`Bio.Seq.Seq`: structure
        """
        return self.get_seq()

    def get_empirical_formula(self):
        """ Get the empirical formula for a transcript with

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

    def get_molecular_weight(self):
        """ Get the molecular weight for a transcript with

        * 5' monophosphate
        * Deprotonated phosphate oxygens

        Returns:
            :obj:`float`: molecular weight (Da)
        """
        return Bio.SeqUtils.molecular_weight(self.get_seq())


class NascentRnaSpeciesType(RnaSpeciesType):
    pass


class MatureRnaSpeciesType(RnaSpeciesType):
    pass


class ProteinSpeciesType(PolymerSpeciesType):
    """ Knowledge of a protein monomer """

    def get_structure(self):
        """ Get the structure

        Returns:
            :obj:`Bio.Seq.Seq`: structure
        """
        return self.get_seq()

    def get_empirical_formula(self):
        """ Get the empirical formula

        Returns:
            :obj:`str`: empirical formula
        """
        pass  # todo calculate the empirical formula from the sequence

    def get_molecular_weight(self):
        """ Get the molecular weight

        Returns:
            :obj:`float`: molecular weight
        """
        pass  # todo calculate the molecular weight from the sequence

    def get_charge(self):
        """ Get the charge

        Returns:
            :obj:`int`: charge
        """
        pass  # todo calculate the charge from the sequence


class PolymerLocus(KnowledgeBaseObject):
    """ Represents knowledge about a locus of a polymer

    Attributes:
        polymer (:obj:`PolymerSpeciesType`): polymer
        start (:obj:`int`): start position
        end (:obj:`int`): end position
        strand (:obj:`PolymerStrand`): strand
    """
    polymer = obj_model.core.ManyToOneAttribute(PolymerSpeciesType, related_name='loci')
    start = obj_model.core.IntegerAttribute()
    end = obj_model.core.IntegerAttribute()
    strand = obj_model.core.EnumAttribute(PolymerStrand, default=PolymerStrand.positive)

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'polymer', 'name', 'start', 'end', 'strand')


class TranscriptionUnitLocus(PolymerLocus):
    """ Knowledge of a transcription unit

    Attributes:
        pribnow_start (:obj:`int`): start position of the Pribnow box (promoter) relative to the 5' coordinate
        pribnow_end (:obj:`int`): end position of the Pribnow box (promoter) relative to the 5' coordinate

    Related attributes:
        rnas (:obj:`list` of :obj:`Rna`): RNAs
    """

    pribnow_start = obj_model.core.IntegerAttribute()
    pribnow_end = obj_model.core.IntegerAttribute()

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'polymer', 'name', 'start', 'end', 'strand', 'pribnow_start', 'pribnow_end')

    @property
    def chromosome(self):
        """ Get the chromosome

        Returns:
            :obj:`Chromosome`: chromosome
        """
        return self.polymer

    @chromosome.setter
    def chromosome(self, chromosome):
        """ Set the chromosome

        Args:
            chromosome (:obj:`Chromosome`): chromosome
        """
        self.polymer = chromosome

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

    def get_len(self):
        """ Get the transcription unit length

        Returns:
            :obj:`int`: length
        """
        return self.end - self.start + 1

    def get_seq(self):
        """ Get the transcription unit sequence

        Returns:
            :obj:`Bio.Seq.Seq`: sequence
        """
        dna_seq = self.chromosome.get_subseq(self.start, self.end, self.strand)
        return dna_seq.transcribe()

    def get_pribnow_len(self):
        """ Get the length of the Pribnow box

        Returns:
            :obj:`int`: length
        """
        return self.pribnow_start - self.pribnow_end + 1

    def get_pribnow_seq(self):
        """ Get the sequence of the Pribnow box

        Returns:
            :obj:`Bio.Seq.Seq`: sequence
        """
        if self.strand == PolymerStrand.positive:
            return self.chromosome.get_subseq(
                self.get_5_prime() + self.pribnow_end,
                self.get_5_prime() + self.pribnow_start,
                self.strand)[::-1]
        else:
            return self.chromosome.get_subseq(
                self.get_5_prime() - self.pribnow_start,
                self.get_5_prime() - self.pribnow_end).complement()


class GeneType(enum.Enum):
    """ Type of gene """
    mRna = 0
    rRna = 1
    sRna = 2
    tRna = 3


class GeneLocus(PolymerLocus):
    """ Knowledge of a gene

    Attributes:
        rna (:obj:`Rna`): RNA
        symbol (:obj:`str`): symbol
        type (:obj:`GeneType`): type
    """

    rnas = obj_model.core.ManyToManyAttribute(Rna, related_name='genes')
    symbol = obj_model.core.StringAttribute()
    type = obj_model.core.EnumAttribute(GeneType)

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'rnas', 'name', 'symbol', 'type')


class ReactionParticipant(obj_model.core.Model):
    """ Knowledge of a participant in a reaction

    Attributes:
        species_type (:obj:`SpeciesType`): species type
        compartment (:obj:`Compartment): compartment
        coefficient (:obj:`float`): coefficient

    Related attributes:
        reactions (:obj:`list` of :obj:`Reaction`): reactions
    """
    species_type = obj_model.core.ManyToManyAttribute(SpeciesType, related_name='reaction_participants')
    compartment = obj_model.core.ManyToManyAttribute(Compartment, related_name='reaction_participants')
    coefficient = obj_model.core.FloatAttribute()


class Reaction(KnowledgeBaseObject):
    """ Knowledge of reactions

    Attributes:
        id (:obj:`str`): identifier
        cell (:obj:`Cell`): cell
        name (:obj:`str`): name
        participants (:obj:`list` of :obj:`ReactionParticipant`): participants
    """

    cell = obj_model.core.ManyToOneAttribute(Cell, related_name='reactions')
    participants = obj_model.core.ManyToManyAttribute(ReactionParticipant, related_name='reactions')

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'cell', 'name', 'participants')
