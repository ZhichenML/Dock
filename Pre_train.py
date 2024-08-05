import os
os.sys.path.append(os.pardir)
import torch
import warnings
import time
import torch.nn as nn
from torch.autograd import Variable
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from loss.loss_v2 import loss_function
# from util.fragmol_v2 import FragmolUtil as MolUtil
from loss.metric4 import metric_acc

from util.sas_qed import calculateScore
from rdkit.Chem import QED


import os
os.sys.path.append('./')
import torch
import warnings
import time
import torch.nn as nn
from util.fragmol_frag_zyh import FragmolUtil
from rdkit import Chem
import numpy as np
from rdkit.Chem import AllChem, rdGeometry
from model.transformer_v1_res_mp1 import TransformerModel as TransformerModel_contact
from model.transformer_v1_res_fac2 import TransformerModel
from dataloader.dataloader_case_nci_res_merge import pre_train_dataset as dataset
from torch.utils.data import DataLoader
from rdkit import RDLogger
import time
from model.transformer_v1_res_mp1 import topkp_random
import argparse

from torch.utils.tensorboard import SummaryWriter

os.environ["CUDA_DEVICE_ORDER"] =    "PCI_BUS_ID"
device_ids=[0,1,2]
warnings.filterwarnings('ignore')
torch.backends.cudnn.enabled= True
epochs = 200
vocab = ['C','N','O','S' ,'P','F','Cl','Br','I','total']
lambda_ = 10
vocab2 = {'C':1,'N':2,'O':3,'S':4,'P':5,'F':6,'Cl':7,'Br':8,'I':9}

def nonC_rate(mol):
    total = 0
    nonC = 0
    non = ['C']
    for i in range(mol.GetNumAtoms()):
        s = mol.GetAtomWithIdx(i).GetSymbol()
        if s in vocab2.keys():
            total+=1
            if s not in non:
               nonC+=1

    return nonC/total


def criterion(recon,target):
    # label = quantized.max(dim=1)[1]

    reconstruction_function = nn.CrossEntropyLoss()
    NLL = reconstruction_function(recon, target)
    return NLL
def dice_batch_single(vox_ref, vox_g):

    vox_g = vox_g
    vox_ref = vox_ref
    sum_ref = torch.sum(torch.square(vox_ref), dim=-1)
    sum_gt = torch.sum(torch.square(vox_g), dim=-1)
    intersection = vox_ref * vox_g
    intersection_score = torch.sum(intersection, dim=-1)
    total_score = sum_ref + sum_gt - intersection_score
    score = intersection_score / total_score
    score1 = score.reshape(-1)

    return score1.cpu().numpy().tolist()


def train_loop(caption, testloader, testset, savedir, rank, optimizer, lr, epochs, batch_size):
    
    caption.train()

    if rank==0:
        f = open(savedir+'/pred','w')
        ff = open(savedir+'/smiles.smi','w')

    losses = []
    
    total = 0
    count = 0
    unique = set()
    nonC = []
    scores_finger = []
    fnames = testset.fnames
    dist_metrics = []
    theta_metrics = []
    degree_metrics = []
    sas = []
    qed = []

    target = np.load('/public/home/gongzhichen/code/Lingo3DMol/examples.npz') # captions, gt_coords, smi_map, smi_map_n1, smi_map_n2
    # captions = torch.tensor(target['captions'][:,1:])
    captions = torch.tensor(target['captions']) #[40, 100]
    gt_coords = torch.tensor(target['gt_coords']) #[40, 100, 3]  [:, 1:99]
    smi_map = torch.tensor(target['smi_map'])
    smi_map_n1 = torch.tensor(target['smi_map_n1'])
    smi_map_n2 = torch.tensor(target['smi_map_n2'])
    dist = torch.tensor(target['dist_list']) # [:, 1:]
    theta = torch.tensor(target['theta_list'])
    degree = torch.tensor(target['degree_list'])
    
    code = captions
    tgt_coords = gt_coords
    batch_size=1
    epoch_losses = []
    for epoch in range(epochs):
        # for i, (code, residue, contact, coords, atom_type,mask, tgt_coords, smi_map, smi_map_n1, smi_map_n2,theta, dist,degree, finger, center, index) in enumerate(testloader):
        for i, (center, atom_types, positions, mask, 
                code, new_position, neighbor_positoins, neighbor_positoins_n1, neighbor_positoins_n2, theta, new_dist, degree, index
        ) in enumerate(testloader):
        # for i, (coords, residue, atom_type, mask, center, index, contact_prob, contact_scaffold_prob) in enumerate(testloader):
                    
                    # target output
                    
                    code = code.cuda()
                    captions = Variable(code) #[40, 100]
                    gt_coords = new_position.cuda() #[40, 100, 3]  [:, 1:99]
                    smi_map = neighbor_positoins.cuda()
                    smi_map_n1 = neighbor_positoins_n1.cuda()
                    smi_map_n2 = neighbor_positoins_n2.cuda()
                    dist = new_dist.cuda() # [:, 1:]
                    theta = theta.cuda()
                    degree = degree.cuda()
                    
                    # model input
                    atom_type = atom_types.cuda()
                    coords = positions.cuda()
                    mask = mask.cuda()
                    
                    # finger = finger.repeat(2,1)
                    index = index.repeat(batch_size, 1)
                    center = center.repeat(batch_size, 1)
                    
                    targets = code[:,1:].reshape(-1)
                    label_coords = gt_coords[:, 1:, :]
                    
                    src_mask_repeat = mask.squeeze(1).repeat(batch_size,1)
                    

                    
                    residue = None
                    contact_prob1 = None
                    contact_idx = None
                    recon, pred_pos, dist_pred, dist_pred_aux, theta_pred, theta_pred_aux, degree_pred, degree_pred_aux = caption(coords=coords, residue=residue, 
                                                    type=atom_type, critical_anchor=contact_prob1, contact_idx=contact_idx, sample_num=batch_size,
                                                    src_mask= mask, 
                                                    captions=code,gt_coords=gt_coords,smi_map=smi_map, smi_map_n1=smi_map_n1, smi_map_n2=smi_map_n2,isTrain=True)
                    # recon, pred_pos,cap_mask, Value_prob, local_coords = caption(coords=coords, type=atom_type,
                    # residue=residue, src_mask=mask,
                    # critical_anchor=contact_prob1,
                    # contact_idx=contact_idx,
                    # sample_num=sample_num, isTrain=args.isTrain, 
                    # isMultiSample=args.isMultiSample,
                    # USE_THRESHOLD=args.USE_THRESHOLD,
                    # isGuideSample=args.isGuideSample,guidepath=[np.tile(gudiecap,(args.gen_set_size,1)),np.tile(gudiepos,(args.gen_set_size,1,1))],
                    # start_this_step=(gudiecap>0).sum(),
                    # OnceMolGen=args.OnceMolGen,frag_len=args.frag_len_add,tempture=args.tempture)

                    # caption(coords, type, residue,critical_anchor,src_mask, 
                    #             captions=None,gt_coords=None, smi_map=None, smi_map_n1=None, 
                    #             smi_map_n2=None,isTrain=True, contact_idx=None, sample_num=1, max_sample_step=100, 
                    #             isMultiSample=False,USE_THRESHOLD=False,isGuideSample=False,guidepath=None,
                    #             start_this_step=0,OnceMolGen=False,frag_len=0,tempture=1.0)

                    label_ind = torch.where(code[:, 1:] != 0)
                    targets = code[:, 1:][label_ind].reshape(-1)
                    
                    dist = dist[:, 1:][label_ind]
                    theta = theta[:, 1:][label_ind]
                    degree = degree[:, 1:][label_ind]
                    label_coords = label_coords[label_ind]

                    recon = recon[:, :-1]
                    recon = recon[label_ind]
                    recon = recon.reshape(-1, 76)
                    pred_pos = pred_pos[label_ind]
                    dist_pred = dist_pred[label_ind]
                    dist_pred_aux = dist_pred_aux[label_ind]
                    theta_pred = theta_pred[label_ind]
                    theta_pred_aux = theta_pred_aux[label_ind]
                    degree_pred = degree_pred[label_ind]
                    degree_pred_aux = degree_pred_aux[label_ind]

                    
                    loss,type_loss,coords_nll,dist_loss1,dist_loss1_aux,theta_loss,theta_loss_aux,degree_loss,degree_loss_aux= loss_function(recon,targets,
                                                label_coords,pred_pos,dist,dist_pred,dist_pred_aux,theta,theta_pred,theta_pred_aux,degree,degree_pred,degree_pred_aux)

                    epoch_losses.append(loss.item())
                    writer.add_scalar('loss', loss.item(), epoch)
                    writer.add_scalar('type_loss', type_loss.item(), epoch)
                    writer.add_scalar('coords_nll', coords_nll.item(), epoch)
                    writer.add_scalar('dist_loss1', dist_loss1.item(), epoch)
                    writer.add_scalar('dist_loss1_aux', dist_loss1_aux.item(), epoch)
                    writer.add_scalar('theta_loss', theta_loss.item(), epoch)
                    writer.add_scalar('theta_loss_aux', theta_loss_aux.item(), epoch)
                    writer.add_scalar('degree_loss', degree_loss.item(), epoch)
                    writer.add_scalar('degree_loss_aux', degree_loss_aux.item(), epoch)

                    if i % 100 == 0:
                        print(f'Epoch: {epoch}, Batch: {i}, Loss: {loss.item()}')

                    loss.backward()
                    optimizer.step()
                    optimizer.zero_grad()
                        
def get_pdb_files(ints):
    pdbs = []
    with open(ints,'r') as f:
        lines = f.readlines()
        for li in lines:
            pdbs.append(li.strip())
    return pdbs

def main(args):
    os.environ["CUDA_VISIBLE_DEVICES"] = args.cuda

    cases_dude = get_pdb_files(args.input_list)
    caption_contact = TransformerModel_contact()
    path_model = args.contact_path
    
    dict_ = torch.load(path_model, map_location='cpu')

    caption_contact.load_state_dict(dict_,strict=False)
    caption_contact = nn.DataParallel(caption_contact)
    caption_contact.cuda()
    caption_contact.eval()

    caption = TransformerModel()
    caption.load_state_dict(torch.load(args.caption_path, map_location='cpu'))
    caption = nn.DataParallel(caption)
    caption.cuda()
    caption.eval()

    num_workers = 4
    epochs = 200
    lr = 0.0001
    optimizer = torch.optim.Adam(caption.parameters(), lr=lr)
    batch_size = 40
    
    for i,case in enumerate(cases_dude):
        print(case,'-'*100)

        cases = [case]
        name = case.split('/')[-1]
        print(f'{i} {case} .........................................................')
        savedir = f'{args.savedir}{args.cuda}/{name}'
        print(savedir)
        if args.saveMol:
            os.system(f'rm -rf {savedir}')
            os.makedirs(savedir, exist_ok=True)

        testset = dataset(cases)

        testloader = DataLoader(dataset=testset,
                                batch_size=1,
                                shuffle=False,pin_memory=True, num_workers=0)

        print("Prep data done", len(testloader))
        
        train_loop(caption, testloader, testset, savedir, rank=0, optimizer=optimizer, lr=lr, epochs=epochs, batch_size=batch_size)

if __name__ == '__main__':
    import pytz
    from datetime import datetime

    start = time.time()
    parser = argparse.ArgumentParser(description='inference')
    parser.add_argument('--savedir',type=str, help='savepath')
    parser.add_argument('--contact_path', type=str,default='checkpoint/contact.pkl')
    parser.add_argument('--caption_path', type=str,default='checkpoint/gen_mol.pkl')
    parser.add_argument('--cuda', type=str)
    parser.add_argument('--coc_dis', type=float, default=2.5)
    parser.add_argument('--nci_thrs', type=float, default=0.7)
    parser.add_argument('--topk', type=int, default=5)
    parser.add_argument('--max_run_hours', type=int)
    parser.add_argument('--gennums', type=int)
    parser.add_argument('--cuda_list', type=int,nargs='+')
    parser.add_argument('--input_list', type=str)
    parser.add_argument('--saveMol', action='store_true', default=True)
    parser.add_argument('--isTrain', action='store_true')
    parser.add_argument('--USE_THRESHOLD', action='store_true', default=True)
    parser.add_argument('--isMultiSample', action='store_true', default=True)
    parser.add_argument('--isGuideSample', action='store_true', default=True)
    parser.add_argument('--OnceMolGen', action='store_true')
    parser.add_argument('--gen_set_size', type=int, default=1)
    parser.add_argument('--prod_time', type=int, default=1)
    parser.add_argument('--tempture', type=float, default=1.0)
    parser.add_argument('--frag_len_add', type=int, default=0)
    args = parser.parse_args()
    import logging
    if args.cuda == str(args.cuda):
        logging.basicConfig(filename='logs.log', level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')
        tz = pytz.timezone('Asia/Shanghai')    
        now = datetime.now(tz)
        for arg, value in args.__dict__.items():
            logging.info('%s: %s', arg, value)

    writer = SummaryWriter(log_dir=f'logs/{args.cuda}')
    main(args)