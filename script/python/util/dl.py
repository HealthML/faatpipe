import tensorflow as tf
from tensorflow.keras import layers
import numpy as np


class DnaOneHot(layers.Layer):
    
    '''
    Layer that turns DNA-strings into one-hot encoded sequences
    '''
    
    def __init__(self, seqlen, complement=False, reverse=False, **kwargs):
        super(DnaOneHot, self).__init__(**kwargs)
        # inspired by, and extended:
        # https://medium.com/@h_76213/efficient-dna-embedding-with-tensorflow-ffce5f499083
        self.reverse = reverse
        self.complement = complement
        self.seqlen = seqlen
        _embedding_values = np.zeros([85, 4], np.float32)
        if not complement:
            _embedding_values[ord('A')] = np.array([1, 0, 0, 0])
            _embedding_values[ord('C')] = np.array([0, 1, 0, 0])
            _embedding_values[ord('G')] = np.array([0, 0, 1, 0])
            _embedding_values[ord('T')] = np.array([0, 0, 0, 1])
        else:
            _embedding_values[ord('A')] = np.array([0, 0, 0, 1])
            _embedding_values[ord('C')] = np.array([0, 0, 1, 0])
            _embedding_values[ord('G')] = np.array([0, 1, 0, 0])
            _embedding_values[ord('T')] = np.array([1, 0, 0, 0])            
        
        self.embedding_table = self.add_weight(
            shape=_embedding_values.shape,
            initializer=tf.constant_initializer(_embedding_values),
            trainable=False,
            name='dna_lookup_table'
        )
        
    def call(self, inputs):
        dna_32 = tf.cast(tf.io.decode_raw(inputs, tf.uint8), tf.int32)
        encoded = tf.nn.embedding_lookup(self.embedding_table, dna_32)
        encoded = tf.reshape(encoded, (-1,self.seqlen, 1, 4))
        if self.reverse:
            return tf.reverse(encoded, axis=[1])
        return encoded
