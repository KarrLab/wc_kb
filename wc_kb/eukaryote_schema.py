""" Schema to represent a knowledge base to build models of eukaryotes

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2018-09-10
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_utils.util import chem
import enum
import obj_model
from wc_kb import core


#####################
#####################
# Enumeration classes

class RegulatoryElementType(enum.Enum):
    """ Type of regulatory element """
    promoter = 1
    promoter_flanking_region = 2
    enhancer = 3
    CTCF_binding_site = 4
    TF_binding_site = 5
    open_chromatin_region = 6


class ActivityLevel(enum.Enum):
    active = 1
    poised = 2
    repressed = 3
    inactive = 4
    na = 5


#####################
#####################
# Locus types


class GeneLocus(core.PolymerLocus):
    """ Knowledge of a gene

    Attributes:
        symbol (:obj:`str`): symbol
        type (:obj:`GeneType`): type of gene

    Related attributes:
        rna (:obj:`list` of :obj:`PreRnaSpeciesType`): rna
        regulatory_modules (:obj:`RegulatoryModule`): regulatory_modules
    """
    symbol = obj_model.StringAttribute()
    type = obj_model.EnumAttribute(core.GeneType)

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'polymer', 'name', 'symbol', 'type', 'strand', 'start', 
                           'end', 'comments', 'references', 'database_references')


class ExonLocus(core.PolymerLocus):
    """ Knowledge of an exon

    Related attributes:
        transcript (:obj:`TranscriptSpeciesType`): transcript
    """        
    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'polymer', 'name', 'start', 'end', 
                           'comments', 'references', 'database_references')


class CdsLocus(core.PolymerLocus):
    """ Knowledge of an coding region

    Attributes:
        exon (:obj:`ExonLocus`): exon

    Related attributes:
        protein (:obj:`ProteinSpeciesType`): protein
    """        
    exon = obj_model.OneToOneAttribute(ExonLocus, related_name='cds_loci')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'polymer', 'name', 'exon', 'start', 'end', 
                           'comments', 'references', 'database_references')
        verbose_name = 'CDS loci'                           


class RegulatoryElementLocus(core.PolymerLocus):
    """ Knowledge of a regulatory element of a gene

    Attributes:
        type (:obj:`RegulatoryElementType`): type of regulatory element
        activity (:obj:`ActivityLevel`): cell-type specific activity level
        bound_start (:obj:`int`): start coordinate of binding
        bound_end (:obj:`int`): end coordinate of binding
        motif_features (:obj:`list` of :obj:`ProteinSpeciesType`): proteins that bind to the site
        
    Related attributes:
        regulatory_modules (:obj:`list` of :obj:`RegulatoryModule`): regulatory_modules
    """
    type = obj_model.EnumAttribute(RegulatoryElementType)
    activity = obj_model.EnumAttribute(ActivityLevel)
    bound_start = obj_model.PositiveIntegerAttribute()
    bound_end = obj_model.PositiveIntegerAttribute()
    motif_features = obj_model.ManyToManyAttribute('ProteinSpeciesType', related_name='regulatory_elements')
    
    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'polymer', 'name', 'type', 'activity', 
                           'start', 'end', 'bound_start', 'bound_end', 'motif_features', 
                           'comments', 'references', 'database_references')


class RegulatoryModule(obj_model.Model):
    """ Knowledge about regulatory modules

    Attributes:        
        gene (:obj:`GeneLocus`): gene
        regulatory_elements (:obj:`list` of :obj:`RegulatoryElementLocus`): regulatory elements
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
        database_references (:obj:`list` of :obj:`DatabaseReference`): database references
    """
    gene = obj_model.OneToOneAttribute(GeneLocus, related_name='regulatory_modules')
    regulatory_elements = obj_model.ManyToManyAttribute(
        RegulatoryElementLocus, related_name='regulatory_modules')
    comments = obj_model.LongStringAttribute()
    references = obj_model.ManyToManyAttribute(core.Reference, related_name='regulatory_modules')
    database_references = core.DatabaseReferenceAttribute(related_name='regulatory_modules')    

    class Meta(obj_model.Model.Meta):
        attribute_order = ('gene', 'regulatory_elements', 'comments', 'references', 'database_references')


#####################
#####################
# Species types


class PreRnaSpeciesType(core.PolymerSpeciesType):
    """ Knowledge of a transcribed RNA species before splicing

    Attributes:
        gene (:obj:`GeneLocus`): gene
        type (:obj:`RnaType`): type

    Related attributes:
        transcripts (:obj:`list` of :obj:`TranscriptSpeciesType`): transcript(s)
    """
    gene = obj_model.OneToOneAttribute(
        'GeneLocus', related_name='rna')
    type = obj_model.EnumAttribute(core.RnaType)

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'type', 'gene', 'circular', 'double_stranded', 
                           'half_life', 'comments', 'references', 'database_references')
        verbose_name = 'pre-RNA species type'

    def get_seq(self):
        """ Get the 5' to 3' sequence

        Returns:
            :obj:`Bio.Seq.Seq`: sequence
        """
        
        dna_seq = self.gene.polymer.get_subseq(
            start=self.gene.start, end=self.gene.end, strand=self.gene.strand)
        
        return dna_seq.transcribe()    

    def get_empirical_formula(self):
        """ Get the empirical formula for a transcribed RNA species before splicing with

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
        """ Get the charge for a transcribed RNA species before splicing with

        * 5' monophosphate
        * Deprotonated phosphate oxygens

        :math:`-L - 1`

        Returns:
           :obj:`int`: charge
        """
        return -self.get_len() - 1

    def get_mol_wt(self):
        """ Get the molecular weight for a transcribed RNA species before splicing with

        * 5' monophosphate
        * Deprotonated phosphate oxygens

        Returns:
            :obj:`float`: molecular weight (Da)
        """
        return self.get_empirical_formula().get_molecular_weight()


class TranscriptSpeciesType(core.PolymerSpeciesType):
    """ Knowledge of a transcript (spliced RNA) species

    Attributes:
        rna (:obj:`RnaSpeciesType`): transcribed RNA species before splicing        
        exons (:obj:`list` of :obj:`ExonLocus`): exons

    Related attributes:
        protein (:obj:`ProteinSpeciesType`): protein    
    """
    rna = obj_model.ManyToOneAttribute(PreRnaSpeciesType, related_name='transcripts')
    exons = obj_model.ManyToManyAttribute(ExonLocus, related_name='transcripts')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'rna', 'exons', 'half_life', 
                           'comments', 'references', 'database_references')

    def get_seq(self):
        """ Get the 5' to 3' sequence

        Returns:
            :obj:`Bio.Seq.Seq`: sequence
        """
        dna_seq = ''

        for exon in self.exons:                
            dna_seq += self.rna.gene.polymer.get_subseq(
                    start=exon.start, end=exon.end)

        if self.rna.gene.strand==core.PolymerStrand.negative:
            dna_seq = dna_seq.reverse_complement()    
        
        return dna_seq.transcribe()    

    def get_empirical_formula(self):
        """ Get the empirical formula for a transcript (spliced RNA) species with

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
        """ Get the charge for a transcript (spliced RNA) species with

        * 5' monophosphate
        * Deprotonated phosphate oxygens

        :math:`-L - 1`

        Returns:
           :obj:`int`: charge
        """
        return -self.get_len() - 1

    def get_mol_wt(self):
        """ Get the molecular weight for a transcript (spliced RNA) species with

        * 5' monophosphate
        * Deprotonated phosphate oxygens

        Returns:
            :obj:`float`: molecular weight (Da)
        """
        return self.get_empirical_formula().get_molecular_weight()    


class ProteinSpeciesType(core.PolymerSpeciesType):
    """ Knowledge of a protein monomer

    Attributes:
        uniprot (:obj:`str`): uniprot id
        transcript (:obj:`TranscriptSpeciesType`): transcript
        coding_regions (:obj:1list` of :obj:`CdsLocus`): CDS Loci

    Related attributes:
        regulatory_elements (:obj:`RegulatoryElementLocus`): protein binding sites    
    """

    uniprot = obj_model.StringAttribute()
    transcript = obj_model.OneToOneAttribute(TranscriptSpeciesType, related_name='protein')
    coding_regions = obj_model.OneToManyAttribute(CdsLocus, related_name='protein')
    
    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'uniprot', 'transcript', 'coding_regions', 'half_life', 
                           'comments', 'references', 'database_references')

    def get_seq(self, table=1, cds=True):
        """ Get the 5' to 3' sequence

        Args:
            table (:obj:`int`, optional): NCBI identifier for translation table 
                                        (default = standard table)
            cds (:obj:`bool`, optional): True indicates the sequence is a complete CDS     

        Returns:
            :obj:`Bio.Seq.Seq`: sequence
        """
        dna_seq = ''

        for cds in sorted(self.coding_regions, key=lambda x: x.start):                
            dna_seq += self.transcript.rna.gene.polymer.get_subseq(
                    start=cds.start, end=cds.end)
        
        if self.transcript.rna.gene.strand == core.PolymerStrand.negative:
            dna_seq = dna_seq.reverse_complement()
            
        return dna_seq.transcribe().translate(table=table, cds=cds)

    def get_empirical_formula(self, table=1, cds=True):
        """ Get the empirical formula

        Args:
            table (:obj:`int`, optional): NCBI identifier for translation table 
                                        (default = standard table)
            cds (:obj:`bool`, optional): True indicates the sequence is a complete CDS

        Returns:
            :obj:`chem.EmpiricalFormula`: empirical formula
        """
        seq = self.get_seq(table=table, cds=cds)
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

    def get_charge(self, table=1, cds=True):
        """ Get the charge at physiological pH

        Args:
            table (:obj:`int`, optional): NCBI identifier for translation table 
                                        (default = standard table)
            cds (:obj:`bool`, optional): True indicates the sequence is a complete CDS

        Returns:
            :obj:`int`: charge
        """
        seq = self.get_seq(table=table, cds=cds)

        n_r = seq.count('R')
        n_h = seq.count('H')
        n_k = seq.count('K')
        n_d = seq.count('D')
        n_e = seq.count('E')

        return (n_r + n_h + n_k) - (n_d + n_e)

    def get_mol_wt(self, table=1, cds=True):
        """ Get the molecular weight

        Args:
            table (:obj:`int`, optional): NCBI identifier for translation table 
                                        (default = standard table)
            cds (:obj:`bool`, optional): True indicates the sequence is a complete CDS

        Returns:
            :obj:`float`: molecular weight
        """
        return self.get_empirical_formula(table=table, cds=cds).get_molecular_weight()
