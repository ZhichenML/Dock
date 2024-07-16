# This file is part of Lingo3DMol
#
# Lingo3DMol is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# Lingo3DMol is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with Lingo3DMol. If not, see <https://www.gnu.org/licenses/>.

import os
os.sys.path.append(os.pardir)
import rdkit
from rdkit import Chem
from copy import deepcopy
import numpy as np
from rdkit.Chem import AllChem, rdGeometry
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
import traceback
import pandas as pd


class FragmolUtil():
    def __init__(self): #TODO  code 还可以改

        self.vocab_list = ["pad", "start", "end", 'sep',
                           "C", "c", "N", "n", "S", "s", "O", "o",
                           "F",
                           "X", "Y", "Z",
                           "/", "\\", "@", "V", "H",  # W for @， V for @@
                           # "Cl", "[nH]", "Br", # "X", "Y", "Z", for origin
                           "1", "2", "3", "4", "5", "6",
                           "#", "=", "-", "(", ")", "[", "]", "M" ,"L","+" # Misc
                           ]

        self.vocab_i2c_v1 = {i: x for i, x in enumerate(self.vocab_list)}

        self.vocab_c2i_v1 = {self.vocab_i2c_v1[i]: i for i in self.vocab_i2c_v1}

        self.vocab_list_decode = ["pad", "start", "end", "sep",
                                  "C", "c", "N", "n", "S", "s", "O", "o", # self.ele_token
                                  "F",
                                  "Cl", "[nH]", "Br",
                                  "/", "\\", "@", "@@", "H",  # W for @， V for @@
                                  # "Cl", "[nH]", "Br", # "X", "Y", "Z", for origin
                                  "1", "2", "3", "4", "5", "6",
                                  "#", "=", "-", "(", ")", "[", "]", "[*]","([*])","+"  # Misc
                                  ]

        self.vocab_i2c_v1_decode = {i: x for i, x in enumerate(self.vocab_list_decode)}

        self.vocab_c2i_v1_decode = {self.vocab_i2c_v1_decode[i]: i for i in self.vocab_i2c_v1_decode}

        self.encode_length = 100

        self.ele_token = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]

        self.resolution = 0.1

        self.dist_grid = {0: 0, 10: 1, 11: 2, 12: 3, 13: 4, 14: 5, 15: 6, 16: 7, 17: 8, 18: 9, 19: 10, 20: 11, 21: 12}
        self.dist_grid_inv = {v: k for k, v in self.dist_grid.items()}

        self.vocab_list_decode_new = [
            "pad_0", "start_0", "end_0", 'sep_0', # 3
            "C_0", "C_5", "C_6", "C_10", "C_11", "C_12", # 9
            "c_0", "c_5", "c_6", "c_10", "c_11", "c_12", # 15
            "N_0", "N_5", "N_6", "N_10", "N_11", "N_12", # 21
            "n_0", "n_5", "n_6", "n_10", "n_11", "n_12", # 27
            "S_0",
            "s_0", "s_5", "s_6", "s_10", "s_11", "s_12", # 34
            "O_0", "O_5", "O_6", "O_10", "O_11", "O_12", # 40
            "o_0", "o_5", "o_6", "+_0", "o_11", "o_12", # 46
            "F_0",
            "Cl_0",
            "[nH]_0", "[nH]_5", "[nH]_6", "[nH]_10", "[nH]_11", "[nH]_12", # 54
            "Br_0",
            "/_0", "\\_0", "@_0", "@@_0", "H_0", # 60
            "1_0", "2_0", "3_0", "4_0", "5_0", "6_0", # 66
            "#_0", "=_0", "-_0", "(_0", ")_0", "[_0", "]_0", "[*]_0","([*])_0" # 75
        ]




        self.vocab_i2c_v1_decode_new = {i: x for i, x in enumerate(self.vocab_list_decode_new)}

        self.vocab_c2i_v1_decode_new = {self.vocab_i2c_v1_decode_new[i]: i for i in self.vocab_i2c_v1_decode_new}


    def mergeSmi(self, rwmol, smi, uuid):
        # merge smi to rwmol, delete * position, and connect smi with single bond
        try:
            
            mol = Chem.MolFromSmiles(smi, sanitize=False)
            # label conn atom
            nbr_idx = -1
            orig_idx = -1
            orig_symbol = None

            for i in range(rwmol.GetNumAtoms()-1, -1, -1):
                atom = rwmol.GetAtomWithIdx(i)
                s = atom.GetSymbol()
                if s == '*':
                    nbr_idx = atom.GetIdx()
                    n = atom.GetNeighbors()
                    nbr = n[0]
                    orig_idx = nbr.GetIdx()
                    # orig_symbol = nbr.GetSymbol()
                    atom1 = rwmol.GetAtomWithIdx(orig_idx)
                    atom1.SetProp('delete', str(uuid))
                    break

            if nbr_idx == -1:
                return rwmol, True  # return modified mol, isFinish true

            # remove fake atom *
            rwmol.RemoveAtom(nbr_idx)
            node_to_idx = {}
            # add mol
            #
            pre = rwmol.GetNumAtoms()
            if mol is None:
                return None,True
            combo = Chem.CombineMols(rwmol, mol)

            # connect two atom
            # 1 find pre orig atom
            pre_index = -1
            for i in range(rwmol.GetNumAtoms()):
                atom = rwmol.GetAtomWithIdx(i)
                labels = atom.GetPropsAsDict()
                #         print(labels.items())
                if 'delete' not in labels.keys():
                    continue
                label = labels['delete']
                if str(label) == str(uuid):
                    pre_index = atom.GetIdx()
                    break

            # 2 find cur conn atom (fake)
            # cur_idx = node_to_idx[1 + pre]

            # 3 connect
            ecombo = Chem.RWMol(combo)
            ecombo.AddBond(pre_index, pre, Chem.BondType.SINGLE)

            return ecombo, False
        except Exception as e:
            print(f'{e}')
            return None, True

    def mergeSmiles3D(self, res, position):
        try:

            m = rdkit.Chem.RWMol()
            node_to_idx = {}
            conn = []

            mol = Chem.MolFromSmiles(res[0])

            m = Chem.RWMol(mol)
            for i, r in enumerate(res[1:]):
                # if i == 0:
                #     continue
                m, isFinish = self.mergeSmi(m, r, i)
                # print(Chem.MolToSmiles(m,rootedAtAtom=0))
                if m is None:
                    print('invalid m is None 2')
                    return None, None
                # print(r,Chem.MolToSmiles(m))
                if isFinish:
                    break

            # delete unnecessary fake atom
            # fakes = []
            # for i in range(m.GetNumAtoms()):
            #     atom = m.GetAtomWithIdx(i)
            #     s = atom.GetSymbol()
            #     if s == '*':
            #         fakes.append(atom.GetIdx())

            # if len(fakes)>0:
            #     print('exist fake atom')
            #     return None,None
            # fakes = sorted(fakes, reverse=True)
            # for f in fakes:
            #     # m.ReplaceAtom(f, Chem.Atom('H'))
            #     m.RemoveAtom(f)
            # Chem.SanitizeMol(m)
            # m = Chem.RemoveHs(m)
            # import pdb;pdb.set_trace()

            print('finish topo partial product')
            # if m.GetNumAtoms() != len(position):
            #     print(len(position), m.GetNumAtoms())
            #     return None, None

            conf = Chem.Conformer(m.GetNumAtoms()) # all zeros conf
            
            # conf = m.GetConformer()
            pos_num=0
            for i in range(m.GetNumAtoms()):
                atom_symbol = m.GetAtomWithIdx(i).GetSymbol()
                if atom_symbol=='*':
                    # atom_symbol='H'
                    continue
                pos = position[pos_num]
                p = rdGeometry.Point3D(pos[0], pos[1], pos[2])
                conf.SetAtomPosition(i, p)
                pos_num+=1

            m.AddConformer(conf)
            m.SetProp('_Name', Chem.MolToSmiles(m))
            Chem.SanitizeMol(m)
            # import pdb;pdb.set_trace()

            # print('geo')
            if m is None:
                return None, None

            return Chem.MolToSmiles(m), m

        except Exception as e:
            print(f'invalid {e} {traceback.format_exc()}')

            return None,None

    def decode3d(self, batch_codes, positions):
        
        gen_smiles = []

        tokens = []

        pos_res = []

        positions = positions.astype(np.float16)

        positions = positions.tolist()

        for i, sample in enumerate(batch_codes):
            try:
                res = [] # 去掉数字的 token of fsmiles fragment 列表
                csmile = ""
                token = [] # fsmiles token
                pos = positions[i]
                pos_temp = []

                for j, xchar in enumerate(sample[:]):
                    t = self.vocab_i2c_v1_decode_new[xchar]
                    new_t = t.split('_')[0]

                    if xchar == 2 or xchar== 0:
                        # res.append(deepcopy(csmile))
                        # csmile = ""
                        break

                    if xchar == 1:
                        continue

                    token.append(t)
                    # res save fragment fsmiles
                    if xchar == 3:
                        res.append(deepcopy(csmile))
                        csmile = ""
                        continue

                    if xchar == 1:
                        continue

                    csmile += new_t
                    if self.vocab_c2i_v1_decode[new_t] in self.ele_token:
                        pos_temp.append(pos[j])
                
                tokens.append(deepcopy(token))
                gen_smiles.append(deepcopy(res))
                pos_res.append(deepcopy(pos_temp))
            except Exception as e:
                print(f'decode {e}')
                continue
        reses = []
        moleculars = []
        print(len(gen_smiles))
        
        for i, res in enumerate(gen_smiles):

            try:
                if len(res) > 0:
                    smi, m = self.mergeSmiles3D(res, pos_res[i])
                    reses.append(smi)
                    moleculars.append(m)

                else:

                    reses.append(None)

                    moleculars.append(None)

            except:

                reses.append(None)

                moleculars.append(None)

        return reses, tokens, moleculars
    


    def mergeSmiles3D_local(self, res, position, dist, theta, degree, smi_map, smi_map_n1, smi_map_n2, root_dist, root_theta, root_dehe):
        try:

            m = rdkit.Chem.RWMol()
            node_to_idx = {}
            conn = []

            mol = Chem.MolFromSmiles(res[0])

            m = Chem.RWMol(mol)
            for i, r in enumerate(res[1:]):
                # if i == 0:
                #     continue
                m, isFinish = self.mergeSmi(m, r, i)
                # print(Chem.MolToSmiles(m,rootedAtAtom=0))
                if m is None:
                    print('invalid m is None 2')
                    return None, None
                # print(r,Chem.MolToSmiles(m))
                if isFinish:
                    break

            # delete unnecessary fake atom
            

            print('finish topo partial product')
            # if m.GetNumAtoms() != len(position):
            #     print(len(position), m.GetNumAtoms())
            #     return None, None

            conf = Chem.Conformer(m.GetNumAtoms()) # all zeros conf
            
            # conf = m.GetConformer()
            pos_num=0
            for i in range(m.GetNumAtoms()):
                atom_symbol = m.GetAtomWithIdx(i).GetSymbol()
                if atom_symbol=='*':
                    # atom_symbol='H'
                    continue
                pos = position[pos_num]
                p = rdGeometry.Point3D(pos[0], pos[1], pos[2])
                conf.SetAtomPosition(i, p)
                pos_num+=1

            m.AddConformer(conf)
            m.SetProp('_Name', Chem.MolToSmiles(m))
            Chem.SanitizeMol(m)
            # import pdb;pdb.set_trace()

            # print('geo')
            if m is None:
                return None, None

            return Chem.MolToSmiles(m), m

        except Exception as e:
            print(f'invalid {e} {traceback.format_exc()}')

            return None,None
        

    def decode3d_local(self, batch_codes, positions, local_coords, mask):
        # local_coords = {'dist_list': dist_list.cpu().data.numpy()*mask, 'theta_list': theta_list.cpu().data.numpy()*mask, 'degree_list': degree_list.cpu().data.numpy()*mask, 
                            # 'root_symbol_dist': root_symbol_dist.cpu().data.numpy()*mask, 'root_symbol_theta': root_symbol_theta.cpu().data.numpy()*mask, 'root_symbol_dehe': root_symbol_dehe.cpu().data.numpy()*mask}
        
        dist = local_coords['dist_list']
        theta = local_coords['theta_list']
        degree = local_coords['degree_list']
        root_symbol_dist = local_coords['root_symbol_dist']
        root_symbol_theta = local_coords['root_symbol_theta']
        root_symbol_dehe = local_coords['root_symbol_dehe']
        smi_map = local_coords['smi_map']
        smi_map_n1 = local_coords['smi_map_n1']
        smi_map_n2 = local_coords['smi_map_n2']

        
        gen_smiles = []

        tokens = []

        pos_res = []

        positions = positions.astype(np.float16)

        positions = positions.tolist()

        from collections import defaultdict
        result = {'fsmiles': [], 'smiles':[], 'step':[], 'current index':[], 'current token':[], 'dist root index':[], 'dist root token':[], 'dist predicted':[], 
                  'theta root index':[], 'theta root token':[], 'theta predicted':[], 'dehe root index':[], 'dehe root token':[], 'dehe predicted':[]}
        
        gen_dist, gen_theta, gen_degree = [], [], []
        gen_smi_map, gen_smi_map_n1, gen_smi_map_n2 = [], [], []
        gen_root_symbol_dist, gen_root_symbol_theta, gen_root_symbol_dehe = [], [], []

        for i, sample in enumerate(batch_codes):
            fsmiles = []
            cur_index = []
            cur_token = []
            dist_root_index = []
            dist_root_token = []
            dist_pred= []
            theta_root_index = []
            theta_root_token = []
            theta_pred = []
            degree_root_index = []
            degree_root_token = []
            degree_pred = []
            step = []
            gen_dist_tmp, gen_theta_tmp, gen_degree_tmp = [], [], []
            gen_smi_map_tmp, gen_smi_map_n1_tmp, gen_smi_map_n2_tmp = [], [], []
            root_symbol_dist_tmp, root_symbol_theta_tmp, root_symbol_dehe_tmp = [], [], []

            try:
                res = [] # 去掉数字的 token of fsmiles fragment 列表
                csmile = ""
                token = [] # fsmiles token
                pos = positions[i]
                
                
                pos_temp = []

                for j, xchar in enumerate(sample[:]):
                    t = self.vocab_i2c_v1_decode_new[xchar]
                    # current fsmiles token, root index, root fsmiles
                   
                    fsmiles.append(t)
                    cur_index.append(j)
                    cur_token.append(t)

                    dist_root_index.append(root_symbol_dist[i][j])
                    # dist_root_token.append(self.vocab_i2c_v1_decode_new[sample[root_symbol_dist[i][j]]])
                    dist_pred.append(self.dist_grid_inv[dist[i][j]])

                    theta_root_index.append(root_symbol_theta[i][j])
                    # theta_root_token.append(self.vocab_i2c_v1_decode_new[sample[root_symbol_theta[i][j]]])
                    theta_pred.append(theta[i][j])

                    degree_root_index.append(root_symbol_dehe[i][j])
                    # degree_root_token.append(self.vocab_i2c_v1_decode_new[sample[root_symbol_dehe[i][j]]])
                    degree_pred.append(degree[i][j])
                    
                    # step.append({'current (token, index)':(t,j),
                    #              'dist (root, token, pred)':(root_symbol_dist[i][j], self.vocab_i2c_v1_decode_new[sample[root_symbol_dist[i][j]]], self.dist_grid_inv[dist[i][j]]),
                    #              'theta (root, token, pred)':(root_symbol_theta[i][j], self.vocab_i2c_v1_decode_new[sample[root_symbol_theta[i][j]]], theta[i][j]),
                    #              'dehedraw (root, token, pred)':(root_symbol_dehe[i][j], self.vocab_i2c_v1_decode_new[sample[root_symbol_dehe[i][j]]], degree[i][j]),
                    #              '3D position (current atom, dist root atom, theta root atom, dehedraw root atom)':(pos[j], pos[root_symbol_dist[i][j]], pos[root_symbol_theta[i][j]], pos[root_symbol_dehe[i][j]]) 
                    #                })


                    new_t = t.split('_')[0]

                    if xchar == 2 or xchar== 0:
                        # res.append(deepcopy(csmile))
                        # csmile = ""
                        break

                    if xchar == 1:
                        continue

                    token.append(t)
                    # res save fragment fsmiles
                    if xchar == 3:
                        res.append(deepcopy(csmile))
                        csmile = ""
                        continue

                    csmile += new_t
                
                    if self.vocab_c2i_v1_decode[new_t] in self.ele_token:
                        pos_temp.append(pos[j])

                        gen_dist_tmp.append(dist[i][j])
                        gen_theta_tmp.append(theta[i][j])
                        gen_degree_tmp.append(degree[i][j])
                        gen_smi_map_tmp.append(smi_map[i][j])
                        gen_smi_map_n1_tmp.append(smi_map_n1[i][j])
                        gen_smi_map_n2_tmp.append(smi_map_n2[i][j])
                        
                        root_symbol_dist_tmp.append(root_symbol_dist[i][j])
                        root_symbol_theta_tmp.append(root_symbol_theta[i][j])
                        root_symbol_dehe_tmp.append(root_symbol_dehe[i][j])

                
                result['fsmiles'].append(deepcopy(fsmiles))
                result['smiles'].append(deepcopy(res))
                result['step'].append(deepcopy(step))
                result['current index'].append(deepcopy(cur_index))
                result['current token'].append(deepcopy(cur_token))
                result['dist root index'].append(deepcopy(dist_root_index))
                result['dist root token'].append(deepcopy(dist_root_token))
                result['dist predicted'].append(deepcopy(dist_pred))
                result['theta root index'].append(deepcopy(theta_root_index))
                result['theta root token'].append(deepcopy(theta_root_token))
                result['theta predicted'].append(deepcopy(theta_pred))
                result['dehe root index'].append(deepcopy(degree_root_index))
                result['dehe root token'].append(deepcopy(degree_root_token))
                result['dehe predicted'].append(deepcopy(degree_pred))
                
                
                tokens.append(deepcopy(token))
                gen_smiles.append(deepcopy(res))
                pos_res.append(deepcopy(pos_temp))
                
                gen_dist.append(deepcopy(gen_dist_tmp))
                gen_theta.append(deepcopy(gen_theta_tmp))
                gen_degree.append(deepcopy(gen_degree_tmp))
                gen_smi_map.append(deepcopy(gen_smi_map_tmp))
                gen_smi_map_n1.append(deepcopy(gen_smi_map_n1_tmp))
                gen_smi_map_n2.append(deepcopy(gen_smi_map_n2_tmp))
                gen_root_symbol_dist.append(deepcopy(root_symbol_dist_tmp))
                gen_root_symbol_theta.append(deepcopy(root_symbol_theta_tmp))
                gen_root_symbol_dehe.append(deepcopy(root_symbol_dehe_tmp))

            except Exception as e:
                print(f'Excpetion: decode {e}')
                continue
        
        pd.DataFrame.from_dict(result).to_excel('result.xlsx')
        reses = []
        moleculars = []
        print(len(gen_smiles))
        
        for i, res in enumerate(gen_smiles):
            
            try:
                if len(res) > 0:
                    import pdb; pdb.set_trace()
                    smi_map_symbol = [self.vocab_i2c_v1_decode_new[v] for v in sample[ gen_smi_map[i]]]
                    smi_map_n1_symbol = [self.vocab_i2c_v1_decode_new[v] for v in sample[ gen_smi_map_n1[i]]]
                    smi_map_n2_symbol = [self.vocab_i2c_v1_decode_new[v] for v in sample[ gen_smi_map_n2[i]]]
                    smi, m = self.mergeSmiles3D_local(res, pos_res[i], gen_dist[i], gen_theta[i], gen_degree[i], gen_smi_map[i], gen_smi_map_n1[i], gen_smi_map_n2[i], gen_root_symbol_dist[i], gen_root_symbol_theta[i], gen_root_symbol_dehe[i])
                    reses.append(smi)
                    moleculars.append(m)

                else:

                    reses.append(None)

                    moleculars.append(None)

            except:

                reses.append(None)

                moleculars.append(None)

        return reses, tokens, moleculars

