""" Schema to represent a knowledge base to build models of prokaryotes

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Bilal Shaikh  <bilal.shaikh@columbia.edu>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2018-09-10
:Copyright: 2018, Karr Lab
:License: MIT

TODO:
ProteinSpeciesType.get_seq() => cds=True (complete coding sequence) causes errors, talk to J
"""

from wc_utils.util import chem
import obj_model
import wc_kb.core as schema_core


#####################
#####################
# Species types


class RnaSpeciesType(schema_core.PolymerSpeciesType):
    """ Knowledge of an RNA species

    Attributes:
        transcription_units (:obj:`list` of :obj:`TranscriptionUnitLocus`): transcription units
        type (:obj:`RnaType`): type

    Related attributes:
        proteins (:obj:`list` of :obj:`ProteinSpeciesType`): protein(s)
    """

    transcription_units = obj_model.ManyToManyAttribute(
        'TranscriptionUnitLocus', related_name='rna')
    type = obj_model.EnumAttribute(schema_core.RnaType)

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'type', 'transcription_units',
                           'circular', 'double_stranded', 'concentration', 'half_life', 'comments')

    def get_seq(self):
        """ Get the sequence

        Returns:
            :obj:`Bio.Seq.Seq`: sequence
        """
        tu_start = self.transcription_units[0].start
        tu_end = self.transcription_units[0].end
        dna_seq = self.transcription_units[0].polymer.get_subseq(
            start=tu_start, end=tu_end)
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


class ProteinSpeciesType(schema_core.PolymerSpeciesType):
    """ Knowledge of a protein monomer

    Attributes:
        gene (:obj:`GeneLocus`): gene
        rna (:obj:`RnaSpeciesType`): rna
    """

    gene = obj_model.ManyToOneAttribute('GeneLocus', related_name='proteins')
    rna = obj_model.ManyToOneAttribute(
        'RnaSpeciesType', related_name='proteins')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'gene', 'rna', 'circular',
                           'double_stranded', 'concentration', 'half_life', 'comments')

    def get_seq(self, cds=True):
        """ Get the sequence

        Returns:
            :obj:`Bio.Seq.Seq`: sequence
        """
        return self.gene.get_seq().translate(table=self.cell.knowledge_base.translation_table, cds=cds)

    def get_empirical_formula(self, cds=True):
        """ Get the empirical formula

        Returns:
            :obj:`chem.EmpiricalFormula`: empirical formula
        """

        seq = self.get_seq(cds)
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

    def get_charge(self, cds=True):
        """ Get the charge at physiological pH

        Returns:
            :obj:`int`: charge
        """
        seq = self.get_seq(cds)

        n_r = seq.count('R')
        n_h = seq.count('H')
        n_k = seq.count('K')
        n_d = seq.count('D')
        n_e = seq.count('E')

        return (n_r + n_h + n_k) - (n_d + n_e)

    def get_mol_wt(self, cds=True):
        """ Get the molecular weight

        Returns:
            :obj:`float`: molecular weight
        """
        return self.get_empirical_formula(cds).get_molecular_weight()


#####################
#####################
# Locus types


class PromoterLocus(schema_core.PolymerLocus):
    """ Knowledge of a promoter for a transcription unit

    Attributes:
        pribnow_start (:obj:`int`): Pribnow box start coordinate
        pribnow_end (:obj:`int`): Pribnow box end coordinate

    Related attributes:
        transcription_units (:obj:`list` of :obj:`TranscriptionUnitLocus`)

    """
    pribnow_start = obj_model.IntegerAttribute()
    pribnow_end = obj_model.IntegerAttribute()

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'polymer', 'name', 'pribnow_start',
                           'pribnow_end', 'strand', 'start', 'end', 'comments')


class TranscriptionUnitLocus(schema_core.PolymerLocus):
    """ Knowledge about an open reading frame

    Attributes:
        promoter (:obj:`PromoterLocus`): promoter controlling the TU
        genes (:obj:`list` of :obj:`GeneLocus`): genes
    """

    promoter = obj_model.ManyToOneAttribute(
        'PromoterLocus', related_name='transcription_units')
    genes = obj_model.ManyToManyAttribute(
        'GeneLocus', related_name='transcription_units')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'polymer', 'name', 'strand',
                           'promoter', 'start', 'end', 'genes', 'comments')

    def get_3_prime(self):
        """ Get the 3' coordinate

        Returns:
            :obj:`int`: 3' coordinate
        """
        if self.strand == schema_core.PolymerStrand.positive:
            return self.end
        else:
            return self.start

    def get_5_prime(self):
        """ Get the 5' coordinate

        Returns:
            :obj:`int`: 5' coordinate
        """
        if self.strand == schema_core.PolymerStrand.positive:
            return self.start
        else:
            return self.end


class GeneLocus(schema_core.PolymerLocus):
    """ Knowledge of a gene

    Attributes:
        symbol (:obj:`str`): symbol

    Related attributes:
        proteins (:obj:`list` of :obj:`ProteinSpeciesType`): protein
    """

    symbol = obj_model.StringAttribute()
    type = obj_model.EnumAttribute(schema_core.GeneType)

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'polymer', 'name', 'symbol',
                           'type', 'strand', 'start', 'end', 'comments')

