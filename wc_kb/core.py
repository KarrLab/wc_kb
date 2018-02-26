""" Read and represent a knowledge base

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-02-07
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_utils.util import chem
import Bio.Seq
import Bio.SeqUtils
import enum
import math
import obj_model.core
import obj_model.extra_attributes

ChromosomeStrand = enum.Enum(value='ChromosomeStrand', names=[
    ('positive', 1),
    ('+', 1),
    ('negative', -1),
    ('-', -1),
])


class KnowledgeBase(obj_model.core.Model):
    """ A knowledge base

    Attributes:
        id (:obj:`str`): identifier
        translation_table (:obj:`int`): translation table

    Related attributes:
        cell (:obj:`Cell`): cell
    """

    id = obj_model.core.StringAttribute(primary=True, unique=True)
    name = obj_model.core.StringAttribute()
    version = obj_model.core.StringAttribute()
    translation_table = obj_model.core.IntegerAttribute()

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'name', 'version', 'translation_table')
        tabular_orientation = obj_model.core.TabularOrientation.column


class Cell(obj_model.core.Model):
    """ Knowledge of a cell

    Attributes:
        id (:obj:`str`): identifier
        knowledge_base (:obj:`KnowledgeBase`): knowledge base

    Related attributes:
        chromosomes (:obj:`list` of :obj:`Chromosome`): chromosomes
    """

    knowledge_base = obj_model.core.OneToOneAttribute(KnowledgeBase, related_name='cell')
    id = obj_model.core.StringAttribute(primary=True, unique=True)

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'knowledge_base')
        tabular_orientation = obj_model.core.TabularOrientation.column


class SpeciesTypeType(enum.Enum):
    """ Type of species type """
    metabolite = 0
    dna = 1
    rna = 2
    protein = 3
    pseudospecies = 4


class SpeciesType(obj_model.core.Model):
    """ Knowledge of a molecular species

    Attributes:
        id (:obj:`str`): identifier
        name (:obj:`str`): name
        molecular_weight (:obj:`float`): molecular weight
        type (:obj:`str`): DNA/RNA/metabolite/pseudospecies
    """

    id = obj_model.core.StringAttribute(primary=True, unique=True)
    cell = obj_model.core.ManyToOneAttribute(Cell, related_name='species_types')
    name = obj_model.core.StringAttribute()
    molecular_weight = obj_model.core.FloatAttribute()
    type = obj_model.core.EnumAttribute(SpeciesTypeType)

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'name', 'cell', 'molecular_weight', 'type')


class Chromosome(obj_model.core.Model):
    """ Knowledge of a chromosome

    Attributes:
        id (:obj:`str`): identifier
        cell (:obj:`Cell`): cell
        seq (:obj:`str`): sequence

    Related attributes:
        transcription_units (:obj:`list` of :obj:`TranscriptionUnit`): transcription units
    """

    id = obj_model.core.StringAttribute(primary=True, unique=True)
    cell = obj_model.core.ManyToOneAttribute(Cell, related_name='chromosomes')
    seq = obj_model.extra_attributes.BioSeqAttribute()

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'cell', 'seq')

    def get_len(self):
        """ Get the chromosome length

        Returns:
            :obj:`int`: length
        """
        return len(self.seq)

    def get_subseq(self, start, end, strand=ChromosomeStrand.positive):
        """ Get a subsequence

        Args:
            start (:obj:`int`): start coordinate (1-indexed)
            end (:obj:`int`): end coordinate (1-indexed)
            strand (:obj:`ChromosomeStrand`, optional): strand

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

        if strand == ChromosomeStrand.positive:
            return pos_seq
        else:
            return pos_seq.reverse_complement()


class TranscriptionUnit(obj_model.core.Model):
    """ Knowledge of a transcription unit

    Attributes:
        id (:obj:`str`): identifier
        name (:obj:`str`): identifier
        chromosome (:obj:`Chromosome`): chromosome
        start (:obj:`int`): start position
        end (:obj:`int`): end position
        strand (:obj:`ChromosomeStrand`): strand
        pribnow_start (:obj:`int`): start position of the Pribnow box (promoter) relative to the 5' coordinate
        pribnow_end (:obj:`int`): end position of the Pribnow box (promoter) relative to the 5' coordinate
    """

    id = obj_model.core.StringAttribute(primary=True, unique=True)
    chromosome = obj_model.core.ManyToOneAttribute(Chromosome, related_name='transcription_units')

    name = obj_model.core.StringAttribute()
    start = obj_model.core.IntegerAttribute()
    end = obj_model.core.IntegerAttribute()
    strand = obj_model.core.EnumAttribute(ChromosomeStrand, default=ChromosomeStrand.positive)
    pribnow_start = obj_model.core.IntegerAttribute()
    pribnow_end = obj_model.core.IntegerAttribute()
    genes = obj_model.core.StringAttribute()

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'name', 'chromosome', 'start', 'end', 'strand', 'pribnow_start', 'pribnow_end', 'genes')

    def get_3_prime(self):
        """ Get the 3' coordinate

        Returns:
            :obj:`int`: 3' coordinate
        """
        if self.strand == ChromosomeStrand.positive:
            return self.end
        else:
            return self.start

    def get_5_prime(self):
        """ Get the 5' coordinate

        Returns:
            :obj:`int`: 5' coordinate
        """
        if self.strand == ChromosomeStrand.positive:
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
        if self.strand == ChromosomeStrand.positive:
            return self.chromosome.get_subseq(
                self.get_5_prime() + self.pribnow_end,
                self.get_5_prime() + self.pribnow_start,
                self.strand)[::-1]
        else:
            return self.chromosome.get_subseq(
                self.get_5_prime() - self.pribnow_start,
                self.get_5_prime() - self.pribnow_end).complement()


class Reaction(obj_model.core.Model):
    """ Knowledge of reactions

    Attributes:
        id (:obj:`str`): identifier
        name (:obj:`str`): name
        participants (:obj:`list` of :obj:`SpeciesType`): participants
    """

    id = obj_model.core.StringAttribute(primary=True, unique=True)
    cell = obj_model.core.ManyToOneAttribute(Cell, related_name='reactions')
    name = obj_model.core.StringAttribute()
    participants = obj_model.core.ManyToManyAttribute(SpeciesType, related_name='reactions')

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'name', 'cell', 'participants')
