
from pysam import VariantFile, FastaFile, VariantHeader
import numpy as np
from pybedtools import BedTool, Interval
from math import ceil, floor
import os

            
class VcfAltParser():
    
    '''
    Class to iterate over a vcf in batches, returns strings of DNA-sequences.
    
    Loosely inspired by janggu.dna.VarianStreamer
    
    :ivar pysam.VariantFile vcf: VariantFile, the variant calls
    :ivar pysam.FastaFile ref: FastaFile, the reference sequence
    :ivar str idx_path: Path to the exported (compatible) variants in bed-format
    :ivar int bin_size: size of the DNA-sequences
    :ivar int n_variants: number of exported (compatible) variants
    '''
    
    def __init__(self,ref_fa_path=None, vcf_path=None, idx_path=None, batch_size=32, bin_size=100, tie='r'):
        '''
        :param str ref_fa_path: Path to indexed reference fasta
        :param str vcf_path: Path to indexed vcf
        :param str idx_path: Path to bed-file which will contain the names and locations of compatible variants
        :param int batch_size: Batch size
        :param int bin_size: Length of the DNA-sequences (centered on the start position of the variant)
        '''
        self.vcf = VariantFile(vcf_path)
        self.ref = FastaFile(ref_fa_path)
        assert os.path.isfile(ref_fa_path + '.fai'), 'Error: no index found for Fasta-file: {}'.format(ref_fa_path)
        self.idx_path = idx_path
        self.batch_size = batch_size
        self.bin_size = bin_size 
        assert tie in ['l','r']
        self.tie = tie
        if not bin_size % 2:
            self.offset = 0 if tie == 'r' else 1
        else:
            self.offset = 0
        self.n_variants = self._initialize_index()
        self._verify_refmatch()
                            
    def get_flanking_centered(self, variant):
        '''
        get flanking sequence, variant will be centered
        '''
        # centers the alt variant (note: ref centering not implemented)
        # flank
        lf, rf = ceil(self.bin_size/2), floor(self.bin_size/2)
        lenref, lenalt = len(variant.ref), len(variant.alts[0])
        d = lenalt - lenref
        # len diff
        ld, rd = ceil(d/2), floor(d/2)
        pos = variant.pos - 1 # 0-based
        left_seq = self.ref.fetch(variant.chrom, pos - lf + ld + self.offset, pos)
        right_seq = self.ref.fetch(variant.chrom, pos + lenref, pos + rf - rd + self.offset)
        return left_seq, right_seq
    
    def get_flanking_right(self, variant):
        '''
        get flanking sequence, variant will be aligned to the right of the center
        '''
        # aligns the varaint to the right of the center 
        # flank
        if self.bin_size % 2:
            rf = floor(self.bin_size/2)
            lf = rf
        else:
            lf, rf = ceil(self.bin_size/2), floor(self.bin_size/2)
        lenref, lenalt = len(variant.ref), len(variant.alts[0])
        d = lenalt - lenref
        # len diff
        pos = variant.pos - 1 # 0-based
        left_seq = self.ref.fetch(variant.chrom, pos - lf + self.offset, pos)
        right_seq = self.ref.fetch(variant.chrom, pos + lenref, pos + rf - d + self.offset)
        return left_seq, right_seq
    
    def get_alt(self, variant):
        '''
        get alternative sequence for a variant
        '''
        l, r = self.get_flanking_right(variant)
        return (l + variant.alts[0] + r).upper()
    
    def get_ref(self, variant):
        '''
        get reference sequence for a variant
        '''
        if self.bin_size % 2:
            rf = floor(self.bin_size/2)
            lf = rf
        else:
            lf, rf = ceil(self.bin_size/2), floor(self.bin_size/2)
        return self.ref.fetch(variant.chrom, variant.pos - lf - 1 + self.offset, variant.pos + rf - 1 + self.offset).upper()
    
    def is_compatible(self, variant):
        '''
        simple test for compatibility 
        '''
        if len(variant.alts[0]) >= (self.bin_size/2):
            return False
        if len(variant.alts) > 1:
            return False
        return True
    
    def _verify_refmatch(self):
        variants = self.vcf.fetch()
        err = 0
        proc = 0
        varid = []
        ref = []
        true_ref = []
        for i in range(50000):
            try:
                variant = next(variants)
                proc += 1
            except StopIteration:
                break
            if variant.ref != self.ref.fetch(variant.chrom, variant.pos - 1, variant.pos - 1 + len(variant.ref)):
                err += 1
                ref.append(variant.ref)
                true_ref.append(self.ref.fetch(variant.chrom, variant.pos - 1, variant.pos - 1 + len(variant.ref)))
                varid.append(variant.id)

        if err:
            print('Warning: {} mismatches with reference based on the first {} variants.'.format(err,min(proc, 50000)))
            for i in range(min(err, 10)):
                print('variant: {}, vcf ref : {}, actual ref: {}'.format(varid[i], ref[i], true_ref[i]))
            if err > 10:
                print('...')
    
    def _initialize_index(self):
        '''
        create a bed-file containing the variant locations and ids
        '''
        bedtool = BedTool((Interval(record.chrom, record.pos - 1, record.pos - 1 + len(record.ref), name='{}_{}>{}'.format(record.id, record.ref, record.alts[0])) for record in self.vcf.fetch() if self.is_compatible(record)))
        bedtool.saveas(self.idx_path)
        with open(self.idx_path, 'r') as infile:
            for n, _ in enumerate(infile):
                pass
        try:
            import subprocess
            subprocess.call(['gzip', '-f', self.idx_path])
            self.idx_path += '.gz'
        except:
            pass
        return n + 1
        
    def batch_generator(self):
        '''
        returns a generator that iterates over pairs of reference and alternative sequences
        '''
        variants = self.vcf.fetch()
        ibatch = 0
        try:
            while True:
                br = []
                ba = []
                ids = []
                while ibatch < self.batch_size:
                    variant = next(variants)
                    if not self.is_compatible(variant):
                        continue
                    ids.append(variant.id)
                    br.append(self.get_ref(variant))
                    ba.append(self.get_alt(variant))
                    ibatch += 1
                yield np.array(ids), (np.array(br), np.array(ba))
                ibatch = 0
        except StopIteration:
            yield np.array(ids), (np.array(br), np.array(ba))
     