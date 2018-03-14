""" Schema to represent knowledge bases

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
import openbabel
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

    Attributes:
        knowledge_base (:obj:`KnowledgeBase`): knowledge base

    Related attributes:
        compartments (:obj:`list` of :obj:`Compartment`): compartments
        species_types (:obj:`list` of :obj:`SpeciesType`): species types
        reactions (:obj:`list` of :obj:`Reaction`): reactions
    """
    knowledge_base = obj_model.core.OneToOneAttribute(KnowledgeBase, related_name='cell')

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'knowledge_base', 'name')
        tabular_orientation = obj_model.core.TabularOrientation.column

    def get_dna_species_types(self):
        """ Get the DNA species types

        Returns:
            :obj:`list` of :obj:`DnaSpeciesType`: DNA species types
        """
        return filter(lambda species_type: isinstance(species_type, DnaSpeciesType), self.species_types)


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
        attribute_order = ('id', 'cell', 'name', 'volume')


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
        attribute_order = ('id', 'cell', 'name', 'concentration', 'half_life')

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
            :obj:`chem.EmpiricalFormula`: empirical formula
        """
        pass

    @abc.abstractmethod
    def get_charge(self):
        """ Get the charge

        Returns:
            :obj:`int`: charge
        """
        pass

    @abc.abstractmethod
    def get_mol_wt(self):
        """ Get the molecular weight

        Returns:
            :obj:`float`: molecular weight
        """
        pass


class MetaboliteSpeciesType(SpeciesType):
    """ Knowledge of a metabolite

    Attributes:
        structure (:obj:`str`): InChI-encoded structure
    """
    structure = obj_model.core.StringAttribute()

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'cell', 'name', 'structure', 'concentration', 'half_life')

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
        attribute_order = ('id', 'cell', 'name', 'circular', 'double_stranded', 'concentration', 'half_life')

    def get_structure(self):
        """ Get the polymer sequence

        Returns:
            :obj:`Bio.Seq.Seq`: sequence
        """
        self.get_seq()

    @abc.abstractmethod
    def get_seq(self):
        """ Get the polymer sequence

        Returns:
            :obj:`Bio.Seq.Seq`: sequence
        """
        pass

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


class DnaSpeciesType(PolymerSpeciesType):
    """ Knowledge of a DNA species

    Attributes:
        seq (:obj:`Bio.Seq.Seq`): sequence
    """
    seq = obj_model.extra_attributes.BioSeqAttribute()

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'cell', 'name', 'seq', 'circular', 'double_stranded', 'concentration', 'half_life')

    def get_seq(self):
        """ Get the sequence

        Returns:
            :obj:`Bio.Seq.Seq`: structure
        """
        return self.seq

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


class RnaType(enum.Enum):
    """ Type of RNA """
    mRna = 0
    rRna = 1
    sRna = 2
    tRna = 3
    mixed = 4


class RnaSpeciesType(PolymerSpeciesType):
    """ Knowledge of an RNA species

    Attributes:
        dna (:obj:`DnaSpeciesType`): polymer
        start (:obj:`int`): start position
        end (:obj:`int`): end position
        strand (:obj:`PolymerStrand`): strand
        type (:obj:`RnaType`): type

    Related attributes:
        genes (:obj:`list` of :obj:`GeneLocus`): genes
        promoters (:obj:`list` of :obj:`PromoterLocus`): promoters
    """
    dna = obj_model.core.ManyToOneAttribute(DnaSpeciesType, related_name='rnas')
    start = obj_model.core.IntegerAttribute()
    end = obj_model.core.IntegerAttribute()
    strand = obj_model.core.EnumAttribute(PolymerStrand, default=PolymerStrand.positive)
    type = obj_model.core.EnumAttribute(RnaType)

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'name', 'cell', 'dna', 'strand', 'start', 'end',
                           'type', 'concentration', 'half_life', 'comments')

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
        """ Get the length

        Returns:
            :obj:`int`: length
        """
        return self.end - self.start + 1

    def get_seq(self):
        """ Get the sequence

        Returns:
            :obj:`Bio.Seq.Seq`: sequence
        """
        dna_seq = self.dna.get_subseq(self.start, self.end, self.strand)
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

    Related attributes:
        orfs (:obj:`list` of :obj:`OpenReadingFrameLocus`): open reading frames
    """

    def get_seq(self):
        """ Get the sequence

        Returns:
            :obj:`Bio.Seq.Seq`: sequence
        """
        orf = self.orfs[0]
        return orf.get_seq().translate(orf.polymer.dna.cell.knowledge_base.translation_table)

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


class PolymerLocus(KnowledgeBaseObject):
    """ Knowledge about a locus of a polymer

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
        return self.start - self.end + 1


class GeneType(enum.Enum):
    """ Type of gene """
    mRna = 0
    rRna = 1
    sRna = 2
    tRna = 3


class GeneLocus(PolymerLocus):
    """ Knowledge of a gene

    Attributes:
        rnas (:obj:`Rna`): RNA
        symbol (:obj:`str`): symbol
        type (:obj:`GeneType`): type
    """
    rnas = obj_model.core.ManyToManyAttribute(RnaSpeciesType, related_name='genes')
    symbol = obj_model.core.StringAttribute()
    type = obj_model.core.EnumAttribute(GeneType)

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'polymer', 'rnas', 'name', 'symbol', 'start', 'end', 'strand', 'type')


class PromoterLocus(PolymerLocus):
    """ Knowledge of a promoter for a transcription unit

    Attributes:
        rnas (:obj:`list` of :obj:`Rna`): RNAs produced from the promoter
        pribnow_start (:obj:`int`): Pribnow box start coordinate, relative to the start site of the RNA (TSS)
        pribnow_end (:obj:`int`): Pribnow box end coordinate, relative to the start site of the RNA (TSS)
    """
    rnas = obj_model.core.OneToManyAttribute(RnaSpeciesType, related_name='promoters')
    pribnow_start = obj_model.core.IntegerAttribute()
    pribnow_end = obj_model.core.IntegerAttribute()

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'name', 'polymer', 'strand', 'rnas', 'start', 'end','pribnow_start', 'pribnow_end')

    def get_pribnow_seq(self):
        """ Get the Pribnow sequence

        Returns:
            :obj:`Bio.Seq.Seq`: sequence
        """
        rna = self.rnas[0]
        five_prime = rna.get_5_prime()

        if rna.strand == PolymerStrand.positive:
            return self.polymer.get_subseq(
                five_prime + self.pribnow_end,
                five_prime + self.pribnow_start)[::-1]
        else:
            return self.polymer.get_subseq(
                five_prime - self.pribnow_start,
                five_prime - self.pribnow_end).complement()


class OpenReadingFrameLocus(PolymerLocus):
    """ Knowledge about an open reading frame

    Attributes:
        protein (:obj:`ProteinSpeciesType`): protein
    """
    protein = obj_model.core.ManyToOneAttribute(ProteinSpeciesType, related_name='orfs')

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'polymer', 'protein', 'name', 'start', 'end', 'strand')


class ReactionParticipant(KnowledgeBaseObject):
    """ Knowledge of a participant in a reaction

    Attributes:
        species_type (:obj:`SpeciesType`): species type
        compartment (:obj:`Compartment`): compartment
        coefficient (:obj:`float`): coefficient

    Related attributes:
        reactions (:obj:`list` of :obj:`Reaction`): reactions
    """
    species_type = obj_model.core.ManyToManyAttribute(SpeciesType, related_name='reaction_participants')
    compartment = obj_model.core.ManyToManyAttribute(Compartment, related_name='reaction_participants')
    coefficient = obj_model.core.FloatAttribute()

    class Meta(obj_model.core.Model.Meta):
        unique_together = (('species_type', 'compartment', 'coefficient'), )
        tabular_orientation = obj_model.core.TabularOrientation.inline

    def serialize(self):
        """ Get value of primary attribute

        Returns:
            :obj:`str`: value of primary attribute
        """
        if self.coefficient == 1:
            return '{}[{}]'.format(self.species_type.serialize(), self.compartment.serialize())
        else:
            return '{} {}[{}]'.format(self.coefficient, self.species_type.serialize(), self.compartment.serialize())

    @classmethod
    def deserialize(cls, value, objects):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model

        Returns:
            :obj:`tuple` of `object`, `InvalidAttribute` or `None`: tuple of cleaned value and cleaning error
        """
        if value in objects[cls]:
            return (objects[cls][value], None)

        matches = re.match('^(.*?) (.*?)\[(.*?)\]$', value)
        if matches:
            coefficient = float(matches.group(1))
            species_type_id = matches.group(2)
            compartment_id = matches.group(3)
        else:
            matches = re.match('^(.*?)\[(.*?)\]$', value)
            coefficient = 1.
            species_type_id = matches.group(1)
            compartment_id = matches.group(2)

        if species_type_id in objects[SpeciesType] and compartment_id in objects[Compartment]:
            species_type = objects[SpeciesType][species_type_id]
            compartment = objects[Compartment][compartment_id]
            return cls(species_type=species_type, compartment=compartment, coefficient=coefficient)

        attr = cls.Meta.primary_attribute
        return (None, InvalidAttribute(attr, [
            'No species type and compartment with primary attribute values "{}" and "{}"'.format(
                species_type_id, compartment_id)]))


class Reaction(KnowledgeBaseObject):
    """ Knowledge of reactions
    Attributes:
        cell (:obj:`Cell`): cell
        participants (:obj:`list` of :obj:`ReactionParticipant`): participants
        k_m (:obj:`float`): K_m value of reaction (unit: mol/L)
        v_max (:obj:`float`):V_max value of reaction (unit: mol/L/min)
        reversible (:obj:`boolean`): denotes whether reaction is reversible
        todo: Handle submodel here or during model generation?
    """

    cell = obj_model.core.ManyToOneAttribute(Cell, related_name='reactions')
    participants = obj_model.core.ManyToManyAttribute(ReactionParticipant, related_name='reactions')
    k_m = obj_model.core.FloatAttribute(min=0)
    v_max = obj_model.core.FloatAttribute(min=0)
    reversible = obj_model.core.BooleanAttribute()

    class Meta(obj_model.core.Model.Meta):
        attribute_order = ('id', 'cell', 'name', 'participants', 'k_m', 'v_max', 'reversible')
