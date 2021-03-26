from argparse import ArgumentParser

import tensorflow as tf
from tensorflow.keras import layers
from tensorflow.keras import Model
from tensorflow.keras.models import load_model
from keras import backend as K

from progress.bar import Bar
from math import ceil
import h5py
import numpy as np

from util.dl import DnaOneHot
from util.variant import VcfAltParser

import os
import sys


def get_args(*argv):
    p = ArgumentParser()
    p.add_argument('-gpu', default=None)
    p.add_argument('-strand', required=True)
    p.add_argument('-vcf', required=True)
    p.add_argument('-ref', required=True)
    p.add_argument('-out', required=True, help='output prefix')
    p.add_argument('-batchsize', required=False, default=1024, type=int)

    args = p.parse_args(*argv)

    if args.gpu is not None:
        _ = int(args.gpu)
        os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
        os.environ["CUDA_VISIBLE_DEVICES"] = args.gpu

    if args.strand in ['minus', 'plus']:
        return args.vcf, args.ref, args.strand, args.out, args.batchsize
    else:
        raise Exception('strand must be "plus" or "minus"')


def get_deepripe_models(strand):
    
    basedir = 'data/deepripe_models/'
    seqlen = 200

    # dummy functions to avoid import errors
    def precision(y_true, y_pred):
        return K.mean(y_true)

    def recall(y_true, y_pred):
        return K.mean(y_true)

    # ENCODE HepG2, K562
    eclip_dirs = [basedir + 'eclip_model_encodeHepG2_high1_seq.h5',
                  basedir + 'eclip_model_encodeHepG2_high2_seq.h5',
                  basedir + 'eclip_model_encodeHepG2_mid1_seq.h5',
                  basedir + 'eclip_model_encodeHepG2_mid2_seq.h5',
                  basedir + 'eclip_model_encodeHepG2_low_seq.h5',
                  basedir + 'eclip_model_encodeK562_high1_seq.h5',
                  basedir + 'eclip_model_encodeK562_high2_seq.h5',
                  basedir + 'eclip_model_encodeK562_mid1_seq.h5',
                  basedir + 'eclip_model_encodeK562_mid2_seq.h5',
                  basedir + 'eclip_model_encodeK562_low_seq.h5']

    # parclip
    parclip_dirs = [basedir + 'parclip_model_high_seq.h5',
                    basedir + 'parclip_model_med_seq.h5',
                    basedir + 'parclip_model_low_seq.h5']

    eclip_models = [load_model(m, custom_objects={'precision': precision, 'recall': recall}) for m in eclip_dirs]
    parclip_models = [load_model(m, custom_objects={'precision': precision, 'recall': recall}) for m in parclip_dirs]

    # rename to avoid conflicts later
    for i, _ in enumerate(eclip_models):
        eclip_models[i]._name = 'eclip_{}'.format(i)
    for i, _ in enumerate(parclip_models):
        parclip_models[i]._name = 'parclip_{}'.format(i)

    in_string = layers.Input((1,), name='dna_string', dtype=tf.string)

    if strand == 'plus':
        '''
        the forward model (plus strand):
        '''
        print('using plus-stand')

        in_seq = DnaOneHot(seqlen + 4)(in_string)

        # eclip input
        seq_squeeze = layers.Reshape((204, 4), name='rna_reshape')(in_seq)
        # parclip input
        seq_squeeze_crp = layers.Cropping1D(cropping=(25, 25))(seq_squeeze)

        in_sq = [] # shifted
        in_sq_crp = [] # shifted, cropped

        for i in range(4):
            in_sq.append(layers.Cropping1D(cropping=(i,4-i), name='in_sq_{}'.format(i))(seq_squeeze))
            in_sq_crp.append(layers.Cropping1D(cropping=(i,4-i), name='in_sq_cr_{}'.format(i))(seq_squeeze_crp))

        out = [[],[],[],[]]
        for m in eclip_models:
            for i in range(4):
                out[i].append(m(in_sq[i]))

        for m in parclip_models:
            for i in range(4):
                out[i].append(m(in_sq_crp[i]))

        out_concat = []
        for i in range(4):
            out_concat.append(layers.Concatenate(name='concat_{}'.format(i))(out[i]))

        out_avg= layers.Average(name='avg_pred')(out_concat)

        model = Model([in_string], [out_avg], name='forward_model')

    else:
        '''
        the reverse-complement model (minus strand):
        '''
        print('using minus-stand')

        in_seq = DnaOneHot(seqlen + 4, reverse=True, complement=True)(in_string)

        # eclip input
        seq_squeeze = layers.Reshape((204, 4), name='rna_reshape')(in_seq)
        # parclip input
        seq_squeeze_crp = layers.Cropping1D(cropping=(25, 25))(seq_squeeze)

        in_sq = []  # shifted
        in_sq_crp = []  # shifted, cropped

        for i in range(4):
            in_sq.append(layers.Cropping1D(cropping=(i, 4 - i), name='in_sq_{}'.format(i))(seq_squeeze))
            in_sq_crp.append(layers.Cropping1D(cropping=(i, 4 - i), name='in_sq_cr_{}'.format(i))(seq_squeeze_crp))

        out = [[], [], [], []]
        for m in eclip_models:
            for i in range(4):
                out[i].append(m(in_sq[i]))

        for m in parclip_models:
            for i in range(4):
                out[i].append(m(in_sq_crp[i]))

        out_concat = []
        for i in range(4):
            out_concat.append(layers.Concatenate(name='concat_{}'.format(i))(out[i]))

        out_avg = layers.Average(name='avg_pred')(out_concat)

        model = Model([in_string], [out_avg], name='forward_model')

    # labels
    hepg_names = ['DDX3X', 'PCBP2', 'FAM120A', 'HNRNPL', 'RBFOX2', 'PTBP1', 'MATR3', 'EFTUD2', 'PRPF4', 'UPF1',
                  'GRWD1', 'PRPF8', 'PPIG', 'CSTF2T', 'QKI', 'U2AF2', 'SUGP2', 'HNRNPM', 'AQR', 'BCLAF1',
                  'LSM11', 'NKRF', 'SUB1', 'NCBP2', 'UCHL5', 'LIN28B', 'IGF2BP3', 'SF3A3', 'AGGF1', 'DROSHA', 'DDX59',
                  'CSTF2', 'DKC1', 'EIF3H', 'FUBP3',
                  'SFPQ', 'HNRNPC', 'ILF3', 'TIAL1', 'HLTF', 'ZNF800', 'PABPN1', 'YBX3', 'FXR2',
                  'GTF2F1', 'IGF2BP1', 'HNRNPK', 'XPO5', 'RPS3', 'SF3B4', 'LARP4', 'BUD13', 'SND1', 'G3BP1', 'AKAP1',
                  'KHSRP',
                  'RBM22', 'GRSF1', 'CDC40', 'NOLC1', 'FKBP4', 'DGCR8', 'ZC3H11A', 'XRN2', 'SLTM', 'DDX55', 'TIA1',
                  'SRSF1', 'U2AF1', 'RBM15']
    hepg_names = [n + '_hepg2' for n in hepg_names]

    k562_names = ['BUD13', 'PTBP1', 'DDX24', 'EWSR1', 'RBM15',
                  'SF3B4', 'YBX3', 'UCHL5', 'KHSRP', 'ZNF622', 'NONO', 'EXOSC5', 'PRPF8', 'CSTF2T', 'AQR', 'UPF1',
                  'U2AF2', 'AKAP8L', 'METAP2', 'SMNDC1', 'GEMIN5', 'HNRNPK', 'SLTM', 'SRSF1', 'FMR1', 'SAFB2', 'DROSHA',
                  'RPS3', 'IGF2BP2', 'ILF3',
                  'RBFOX2', 'QKI', 'PCBP1', 'ZNF800', 'PUM1',
                  'EFTUD2', 'LIN28B', 'AGGF1', 'HNRNPL', 'SND1', 'GTF2F1', 'EIF4G2', 'TIA1', 'TARDBP', 'FXR2', 'HNRNPM',
                  'IGF2BP1', 'PUM2', 'FAM120A',
                  'DDX3X', 'MATR3', 'FUS', 'GRWD1', 'PABPC4',
                  'MTPAP', 'RBM22', 'DHX30', 'DDX6', 'DDX55', 'TRA2A', 'XRN2', 'U2AF1', 'LSM11', 'ZC3H11A', 'NOLC1',
                  'KHDRBS1', 'GPKOW', 'DGCR8', 'AKAP1',
                  'FXR1', 'DDX52', 'AATF']
    k562_names = [n + '_k562' for n in k562_names]

    parclip_names = ['DND1', 'CPSF7', 'CPSF6', 'CPSF1', 'CSTF2', 'CSTF2T', 'ZC3H7B', 'FMR1iso1', 'RBM10', 'MOV10',
                     'ELAVL1',
                     'TARDBP', 'ELAVL2', 'ELAVL3', 'ELAVL4', 'RBM20', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3', 'EWSR1',
                     'HNRNPD',
                     'RBPMS', 'SRRM4', 'AGO2', 'NUDT21', 'FIP1L1', 'CAPRIN1', 'FMR1iso7', 'FXR2', 'AGO1', 'L1RE1',
                     'ORF1',
                     'MBNL1', 'P53_NONO', 'PUM2', 'QKI', 'AGO3', 'FUS', 'TAF15', 'ZFP36', 'DICER1', 'EIF3A', 'EIF3D',
                     'EIF3G', 'SSB', 'PAPD5', 'CPSF4', 'CPSF3', 'RTCB', 'FXR1', 'NOP58', 'NOP56', 'FBL', 'LIN28A',
                     'LIN28B', 'UPF1', 'G35', 'G45', 'XPO5']

    clabels = hepg_names + k562_names + parclip_names

    print(model.summary())
    sys.stdout.flush()

    return model, clabels


def predict_variant_effect(model, parser, labels):
    assert len(labels) == model.output_shape[-1], 'Error: got {} labels, but there are {} outputs!'.format(len(labels),
                                                                                                           model.output_shape[
                                                                                                               -1])
    # output file path
    if parser.idx_path.endswith('gz'):
        out_scores_path = parser.idx_path.replace('.bed.gz', '_scores.h5')
    else:
        out_scores_path = parser.idx_path.replace('.bed', '_scores.h5')

    # number of batches
    n_batch = ceil(parser.n_variants / parser.batch_size)

    print('number of batches: {}'.format(n_batch))
    sys.stdout.flush()
    
    # sequences we want to predict
    variants = parser.batch_generator()

    h5file = h5py.File(out_scores_path, 'w')

    refscores = h5file.create_dataset('refscore', (parser.n_variants, len(labels)), dtype='float16')
    altscores = h5file.create_dataset('altscore', (parser.n_variants, len(labels)), dtype='float16')
    diffscores = h5file.create_dataset('diffscore', (parser.n_variants, len(labels)), dtype='float16')

    h5file.create_dataset('labels', (len(labels),), dtype=h5py.special_dtype(vlen=str),
                          data=np.array(labels).astype(h5py.special_dtype(vlen=str)))

    bar = Bar('Parsing {} batches: '.format(n_batch), max=n_batch)

    predicted = []

    idx = 0
    b = 0
    for ids, (refbatch, altbatch) in variants:

        if not len(ids):
            bar.next()
            break

        refsc = model.predict_on_batch(refbatch)
        altsc = model.predict_on_batch(altbatch)

        refscores[idx:(idx + len(refbatch)), :] = refsc.astype('float16')
        altscores[idx:(idx + len(refbatch)), :] = altsc.astype('float16')
        diffscores[idx:(idx + len(refbatch)), :] = (altsc - refsc).astype('float16')
        idx += len(refbatch)

        predicted.append(ids)
        bar.next()
        #b += 1
        #if (b % 10) == 0:
        #    sys.stdout.flush()
        
    bar.finish()

    return np.concatenate(predicted)


def main():
    vcf_file, ref_file, strand, out, batchsize = get_args()

    if not os.path.isdir(os.path.dirname(out)):
        os.makedirs(os.path.dirname(out), exist_ok=True)

    out_bed = out + '.bed'

    # we add 2 to each side of binsize so we can perform average predictions over small shifts
    vcfparser = VcfAltParser(ref_file, vcf_file, out_bed, batch_size=batchsize, bin_size=204, tie='l')
    model, labels = get_deepripe_models(strand)

    ids = predict_variant_effect(model, vcfparser, labels)

    out_ids = out + '.ids.txt'
    with open(out_ids, 'w') as outfile:
        for rec in ids:
            outfile.write('{}\n'.format(rec))


if __name__ == '__main__':
    main()