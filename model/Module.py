# This file is part of Lingo3DMol
#
# Lingo3DMol is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# Lingo3DMol is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with Lingo3DMol. If not, see <https://www.gnu.org/licenses/>.



import torch.nn as nn
import torch
import torch.nn.functional as F
import copy
import math
import numpy as np


def clones(module, N):
    "Produce N identical layers."
    return nn.ModuleList([copy.deepcopy(module) for _ in range(N)])


class Encoder(nn.Module):
    "Core encoder is a stack of N layers"

    def __init__(self, layer, N):
        super(Encoder, self).__init__()
        self.layers = clones(layer, N)
        self.norm = LayerNorm(layer.size)

    def forward(self, x, mask):
        "Pass the input (and mask) through each layer in turn."

        for layer in self.layers:
            x = layer(x, mask)

        return self.norm(x)


class LayerNorm(nn.Module):
    "Construct a layernorm module (See citation for details)."

    def __init__(self, features, eps=1e-6):
        super(LayerNorm, self).__init__()

        self.a_2 = nn.Parameter(torch.ones(features))
        self.b_2 = nn.Parameter(torch.zeros(features))
        self.eps = eps

    def forward(self, x):

        mean = x.mean(-1, keepdim=True)
        std = x.std(-1, keepdim=True)
        return self.a_2 * (x - mean) / (std + self.eps) + self.b_2


class SublayerConnection(nn.Module):
    """
    A residual connection followed by a layer norm.
    Note for code simplicity the norm is first as opposed to last.
    """

    def __init__(self, size, dropout):
        super(SublayerConnection, self).__init__()

        self.norm = LayerNorm(size)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x, sublayer):
        "Apply residual connection to any sublayer with the same size."
        return x + self.dropout(sublayer(self.norm(x)))


class EncoderLayer(nn.Module):
    "Encoder is made up of self-attn and feed forward (defined below)"

    def __init__(self, size, self_attn, feed_forward, dropout):
        super(EncoderLayer, self).__init__()
        self.self_attn = self_attn
        self.feed_forward = feed_forward
        self.sublayer = clones(SublayerConnection(size, dropout), 2)
        self.size = size

    def forward(self, x, mask):
        "Follow Figure 1 (left) for connections."
        x = self.sublayer[0](x, lambda x: self.self_attn(x, x, x, mask))
        return self.sublayer[1](x, self.feed_forward)


class Decoder(nn.Module):
    "Generic N layer decoder with masking."

    def __init__(self, layer, N):
        super(Decoder, self).__init__()
        self.layers = clones(layer, N)
        self.norm = LayerNorm(layer.size)

    def forward(self, x, memory, src_mask, tgt_mask):
        for layer in self.layers:
            x = layer(x, memory, src_mask, tgt_mask)
        return self.norm(x)

class DecoderBias(nn.Module):
    "Generic N layer decoder with masking."

    def __init__(self, layer, N):
        super(DecoderBias, self).__init__()
        self.layers = clones(layer, N)
        self.norm = LayerNorm(layer.size)

    def forward(self, x, memory, src_mask, tgt_mask,bias,crossbias=None):
        for layer in self.layers:
            x,bias,att = layer(x, memory, src_mask, tgt_mask,bias,crossbias)
        return self.norm(x),att

class DecoderGPT(nn.Module):
    "Generic N layer decoder with masking."

    def __init__(self, layer, N):
        super(DecoderGPT, self).__init__()
        self.norm = LayerNorm(layer.size)
        self.gpts = clones(layer, N)

    def forward(self, x, tgt_mask):
        for layer in self.gpts:
            x = layer(x, tgt_mask)
        return self.norm(x)


class DecoderLayer(nn.Module):
    "Decoder is made of self-attn, src-attn, and feed forward (defined below)"

    def __init__(self, size, self_attn, src_attn, feed_forward, dropout):
        super(DecoderLayer, self).__init__()
        self.size = size
        self.self_attn = self_attn
        self.src_attn = src_attn
        self.feed_forward = feed_forward
        self.sublayer = clones(SublayerConnection(size, dropout), 3)

    def forward(self, x, memory, src_mask, tgt_mask):
        "Follow Figure 1 (right) for connections."
        m = memory
        x = self.sublayer[0](x, lambda x: self.self_attn(x, x, x, tgt_mask))
        x = self.sublayer[1](x, lambda x: self.src_attn(x, m, m, src_mask))
        return self.sublayer[2](x, self.feed_forward)


class DecoderLayerBias(nn.Module):
    "Decoder is made of self-attn, src-attn, and feed forward (defined below)"

    def __init__(self, size, self_attn, src_attn, feed_forward, dropout):
        super(DecoderLayerBias, self).__init__()
        self.size = size
        self.self_attn = self_attn
        self.src_attn = src_attn
        self.feed_forward = feed_forward
        self.sublayer = clones(SublayerConnection(size, dropout), 1)
        self.norm1 = LayerNorm(size)
        self.norm2 = LayerNorm(size)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x, memory, src_mask, tgt_mask,bias,crossbias=None):
        "Follow Figure 1 (right) for connections."
        m = memory
        x1,new_bias =  self.self_attn(x, x, x, bias,tgt_mask)
        x2 = x+self.dropout(self.norm1(x1))

        # if crossbias is not None:
        #     x3, cross_bias = self.src_attn(x2, m, m, crossbias, src_mask)
        # else:
        x3,att = self.src_attn(x2, m, m, src_mask)
        cross_bias = att
        x4 = x2+self.dropout(self.norm2(x3))
        return self.sublayer[0](x4, self.feed_forward),new_bias,cross_bias

class DecoderLayerGPT(nn.Module):
    "Decoder is made of self-attn, src-attn, and feed forward (defined below)"

    def __init__(self, size, self_attn, feed_forward, dropout):
        super(DecoderLayerGPT, self).__init__()
        self.size = size
        self.self_attn = self_attn
        self.feed_forward = feed_forward
        self.sublayer = clones(SublayerConnection(size, dropout), 3)

    def forward(self, x, tgt_mask):
        "Follow Figure 1 (right) for connections."
        x = self.sublayer[0](x, lambda x: self.self_attn(x, x, x, tgt_mask))
        return self.sublayer[2](x, self.feed_forward)


def subsequent_mask(size):
    "Mask out subsequent positions."
    attn_shape = (1, size, size)
    subsequent_mask = np.triu(np.ones(attn_shape), k=1).astype('uint8')
    return torch.from_numpy(subsequent_mask) == 0

def subsequent_mask1(size):
    "Mask out subsequent positions."
    attn_shape = (1, size, size+1)
    subsequent_mask = np.triu(np.ones(attn_shape), k=1).astype('uint8')
    return torch.from_numpy(subsequent_mask) == 0

def attention(query, key, value, mask=None, dropout=None):
    "Compute 'Scaled Dot Product Attention'"
    d_k = query.size(-1)
    scores = torch.matmul(query, key.transpose(-2, -1)) \
             / math.sqrt(d_k)
    # if mask is not None:
    scores1 = scores.masked_fill(mask == 0, -1e9)
    p_attn = F.softmax(scores1, dim=-1)
    if dropout is not None:
        p_attn = dropout(p_attn)
    return torch.matmul(p_attn, value), scores

def attentionbias( value, bias,mask=None, dropout=None):
    "Compute 'Scaled Dot Product Attention'"
    # scores = torch.matmul(query, key.transpose(-2, -1)) \
    #          / math.sqrt(d_k)
    if mask is not None:
        # scores = scores.masked_fill(mask == 0, -1e9)
        bias1  = bias.masked_fill(mask==0,-1e9)
    weight = bias1
    p_attn = F.softmax(weight, dim=-1)
    if dropout is not None:
        p_attn = dropout(p_attn)
    return torch.matmul(p_attn, value),bias

class MultiHeadedAttention(nn.Module):
    def __init__(self, h, d_model, dropout=0.1):
        "Take in model size and number of heads."
        super(MultiHeadedAttention, self).__init__()
        assert d_model % h == 0
        # We assume d_v always equals d_k
        self.d_k = d_model // h
        self.h = h
        self.linears = clones(nn.Linear(d_model, d_model), 4)
        self.attn = None
        self.dropout = nn.Dropout(p=dropout)

    def forward(self, query, key, value, mask=None):
        "Implements Figure 2"
        if mask is not None:
            # Same mask applied to all h heads.
            mask = mask.unsqueeze(1)
        nbatches = query.size(0)


        # 1) Do all the linear projections in batch from d_model => h x d_k
        query, key, value = \
            [l(x).view(nbatches, -1, self.h, self.d_k).transpose(1, 2)
             for l, x in zip(self.linears, (query, key, value))]



        # 2) Apply attention on all the projected vectors in batch.
        x, self.attn = attention(query, key, value, mask=mask, dropout=self.dropout)

        # 3) "Concat" using a view and apply a final linear.
        x = x.transpose(1, 2).contiguous() \
            .view(nbatches, -1, self.h * self.d_k)
        return self.linears[-1](x)

class MultiHeadedAttention_att(nn.Module):
    def __init__(self, h, d_model, dropout=0.1):
        "Take in model size and number of heads."
        super(MultiHeadedAttention_att, self).__init__()
        assert d_model % h == 0
        # We assume d_v always equals d_k
        self.d_k = d_model // h
        self.h = h
        self.linears = clones(nn.Linear(d_model, d_model), 4)
        self.attn = None
        self.dropout = nn.Dropout(p=dropout)

    def forward(self, query, key, value, mask=None):
        "Implements Figure 2"
        if mask is not None:
            # Same mask applied to all h heads.
            mask = mask.unsqueeze(1)
        nbatches = query.size(0)


        # 1) Do all the linear projections in batch from d_model => h x d_k
        query, key, value = \
            [l(x).view(nbatches, -1, self.h, self.d_k).transpose(1, 2)
             for l, x in zip(self.linears, (query, key, value))]



        # 2) Apply attention on all the projected vectors in batch.
        x, self.attn = attention(query, key, value, mask=mask, dropout=self.dropout)

        # 3) "Concat" using a view and apply a final linear.
        x = x.transpose(1, 2).contiguous() \
            .view(nbatches, -1, self.h * self.d_k)
        return self.linears[-1](x),self.attn


class MultiHeadedAttentionBias(nn.Module):
    def __init__(self, h, d_model, dropout=0.1):
        "Take in model size and number of heads."
        super(MultiHeadedAttentionBias, self).__init__()
        assert d_model % h == 0

        self.d_k = d_model // h
        self.h = h
        self.linears = clones(nn.Linear(d_model, d_model), 2)
        self.attn          = None
        self.dropout       = nn.Dropout(p=dropout)

        # self.w             = nn.parameter.Parameter(torch.ones(self.h))
        # self.w1            = nn.parameter.Parameter(torch.FloatTensor([0.05]))
        # self.proj          = MLP(in_features=self.h,hidden_layer_sizes=[256,256],out_features=self.h,dropout_p=0.1)
    def forward(self, query, key, value,bias, mask=None):

        "Implements Figure 2"

        # w = self.w.unsqueeze(0).unsqueeze(-1).unsqueeze(-1)
        if mask is not None:
            # Same mask applied to all h heads.
            mask = mask.unsqueeze(1)

        nbatches = query.size(0)

        # 1) Do all the linear projections in batch from d_model => h x d_k
        value = self.linears[0](value).view(nbatches, -1, self.h, self.d_k).transpose(1, 2)
        #     [l(x).view(nbatches, -1, self.h, self.d_k).transpose(1, 2)
        #      for l, x in zip(self.linears, ( value))]


        # edge_fea1 = self.w1*bias#bias.permute(0,3,1,2)

        # 2) Apply attention on all the projected vectors in batch.
        x, self.attn = attentionbias(value,bias=bias, mask=mask, dropout=self.dropout)

        # edge_fea2 = bias+w*self.attn

        # 3) "Concat" using a view and apply a final linear.
        x = x.transpose(1, 2).contiguous() \
            .view(nbatches, -1, self.h * self.d_k)




        return self.linears[-1](x),bias



class PositionwiseFeedForward(nn.Module):
    "Implements FFN equation."

    def __init__(self, d_model, d_ff, dropout=0.1):
        super(PositionwiseFeedForward, self).__init__()

        self.w_1 = nn.Linear(d_model, d_ff)
        self.w_2 = nn.Linear(d_ff, d_model)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x):
        return self.w_2(self.dropout(F.relu(self.w_1(x))))


class Embeddings(nn.Module):
    def __init__(self, d_model, vocab):
        super(Embeddings, self).__init__()
        self.lut = nn.Embedding(vocab, d_model)
        self.d_model = d_model

    def forward(self, x):
        return self.lut(x) * math.sqrt(self.d_model)


class PositionalEncoding(nn.Module):
    "Implement the PE function."

    def __init__(self, d_model, dropout, max_len=5000):
        super(PositionalEncoding, self).__init__()
        self.dropout = nn.Dropout(p=dropout)

        # Compute the positional encodings once in log space.
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len).unsqueeze(1).float()
        div_term = torch.exp(torch.arange(0, d_model, 2).float() *
                             -(math.log(10000.0) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0)
        self.register_buffer('pe', pe)

    def forward(self, x):
        x1 = x + self.pe[:, :x.size(1)]
        return self.dropout(x1)

class MLP(torch.nn.Module):
    """
    Multi-layer perceptron. Applies SELU after every linear layer.

    Args:
    ----
        in_features (int)         : Size of each input sample.
        hidden_layer_sizes (list) : Hidden layer sizes.
        out_features (int)        : Size of each output sample.
        dropout_p (float)         : Probability of dropping a weight.
    """

    def __init__(self, in_features : int, hidden_layer_sizes : list, out_features : int,
                 dropout_p : float) -> None:
        super().__init__()

        activation_function = torch.nn.SELU

        # create list of all layer feature sizes
        fs = [in_features, *hidden_layer_sizes, out_features]

        # create list of linear_blocks
        layers = [self._linear_block(in_f, out_f,
                                     activation_function,
                                     dropout_p)
                  for in_f, out_f in zip(fs, fs[1:])]

        # concatenate modules in all sequentials in layers list
        layers = [module for sq in layers for module in sq.children()]

        # add modules to sequential container
        self.seq = torch.nn.Sequential(*layers)

    def _linear_block(self, in_f : int, out_f : int, activation : torch.nn.Module,
                      dropout_p : float) -> torch.nn.Sequential:
        """
        Returns a linear block consisting of a linear layer, an activation function
        (SELU), and dropout (optional) stack.

        Args:
        ----
            in_f (int)                   : Size of each input sample.
            out_f (int)                  : Size of each output sample.
            activation (torch.nn.Module) : Activation function.
            dropout_p (float)            : Probability of dropping a weight.

        Returns:
        -------
            torch.nn.Sequential : The linear block.
        """
        # bias must be used in most MLPs in our models to learn from empty graphs
        linear = torch.nn.Linear(in_f, out_f, bias=True)
        torch.nn.init.xavier_uniform_(linear.weight)
        return torch.nn.Sequential(linear, activation()) #, torch.nn.AlphaDropout(dropout_p)

    def forward(self, layers_input : torch.nn.Sequential) -> torch.nn.Sequential:
        """
        Defines forward pass.
        """
        return self.seq(layers_input)

class cdist(nn.Module):

    def __init__(self,dist_num=42,embedding_dim=16):
        super(cdist,self).__init__()

        '''
        :param dist_num: 0.1A
        '''
        self.dist = nn.Embedding(dist_num,embedding_dim)

    def forward(self,input1,input2,emb_dim=3):
        input_shape1 = input1.shape
        input_shape2 = input2.shape

        # Flatten input
        flat_input1 = input1.reshape(input1.shape[0],-1, emb_dim)
        flat_input2 = input2.reshape(input2.shape[0],-1,emb_dim)
        flat_input3 = flat_input2.permute(0,2,1)


        # Calculate distances
        a= torch.sum(flat_input1 ** 2, dim=-1, keepdim=True)
        b = torch.sum(flat_input2 ** 2, dim=-1).unsqueeze(dim=1)
        c = torch.matmul(flat_input1, flat_input3)
        distances = torch.sqrt(a+b-2*c)

        dist_map = torch.clamp(distances.reshape(-1,input_shape1[1],input_shape2[1]),min=0,max=4)*10

        dist_map1 = torch.round(dist_map).long()

        dist      = self.dist(dist_map1)

        return dist


class edge_vector(nn.Module):
    def __init__(self,dist_num=85,embedding_dim=16):
        super(edge_vector,self).__init__()


        self.relative_coords = nn.Embedding(dist_num,embedding_dim)
    def forward(self,input1,input2):

        input1_new = input1.unsqueeze(dim=-2)

        input2_new = input2.unsqueeze(dim=1)

        new_vec = torch.round((torch.clamp(input1_new-input2_new,min=-4.0,max=4.0)+4.0)*10.0).long()

        new_vec1 = self.relative_coords(new_vec)

        new_vec2 = new_vec1.reshape(input1.shape[0],input1.shape[1],input2.shape[1],-1)

        return new_vec2
