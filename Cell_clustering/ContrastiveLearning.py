#!/usr/bin/env python
# coding: utf-8


import torch as tc
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
from torch.autograd import Function
from torch.utils.data import Dataset, DataLoader
import torch.nn as nn 
from torch.optim import Adam
import torch.nn.functional as F
from tqdm import tqdm


dat_raw = pd.read_parquet('../data/MERGED_normalized_5000genesTIERS.parquet')#.iloc[:500,:]

dat = dat_raw.iloc[:,6:]

class DS(Dataset):
    def __init__(self, df, patients):
        self.df = df 

        self.patients = np.array(patients)
        self.unique_patients = patients.drop_duplicates().to_numpy()

        self.data_tensor = tc.FloatTensor((self.df.values))

        self.data_dict = {patient: self.data_tensor[self.patients == patient,:] for patient in tqdm(self.unique_patients[:])}

        self.patient_ids_dict = {patient: self.make_one_hot(self.unique_patients.shape[0], i) for i,patient in enumerate(self.unique_patients)}


    def __len__(self):
        return len(self.data_dict)


    def __getitem__(self, idx):
        current_data_tensor =  self.data_dict[self.unique_patients[idx]]
        nsamples = current_data_tensor.shape[0]


        current_patient_id = self.patient_ids_dict[self.unique_patients[idx]]
        
        if current_data_tensor.shape[0]>500:
            return current_data_tensor[tc.randperm(nsamples)[:500],:], current_patient_id * tc.ones(500,1)

        return current_data_tensor[tc.randperm(nsamples),:], current_patient_id * tc.ones(nsamples,1)


    def make_one_hot(self, length, idx):
        x = tc.zeros(length)
        x[idx]=1
        return x


class ReverseLayerF(Function):

    @staticmethod
    def forward(ctx, x, alpha):
        ctx.alpha = alpha

        return x.view_as(x)

    @staticmethod
    def backward(ctx, grad_output):
        output = grad_output.neg() * ctx.alpha

        return output, None
    
    

class Backbone(nn.Module):
    def __init__(self, inp, hidden, outp):
        super().__init__()

        self.layers = nn.Sequential(
                    #nn.Dropout(0.8),          #0.5
                    nn.Linear(inp, hidden), 
                    nn.BatchNorm1d(hidden),
                    nn.LeakyReLU(), 
                    nn.Dropout(0.0),      
                    nn.Linear(hidden, hidden),
                    nn.BatchNorm1d(hidden),
                    nn.LeakyReLU(),
                    nn.Dropout(0.0),      
                    nn.Linear(hidden, hidden),
                    nn.BatchNorm1d(hidden),
                    nn.LeakyReLU(),
                    nn.Dropout(0.0),      
                    nn.Linear(hidden, 1000),
                    nn.LeakyReLU(),
                    nn.Dropout(0.0),      
                    nn.Linear(1000, outp))

    def forward(self,x):
        x = self.layers(x)

        return x
        


class Model(nn.Module):
    def __init__(self, inp, hidden, npatients):
        super().__init__()

        self.npatients = npatients
        self.backbone = Backbone(inp, hidden, 64)
        
        self.projector = nn.Sequential(
            nn.Linear(64, 512), 
            nn.LeakyReLU(),
            nn.Linear(512,512))

        self.classifier = nn.Sequential(
            nn.Linear(64,512),
            nn.LeakyReLU(),
            nn.Linear(512, self.npatients)
        )
    
    def forward(self, x, y):


        def off_diagonal(arr):
            n = len(arr)
            return arr.flatten()[:-1].view(n - 1, n + 1)[:, 1:].flatten()

        x = self.projector(self.backbone(x))
        y = self.projector(self.backbone(y))


        repr_loss = F.mse_loss(x, y)  # invariance (2)
        x = x - x.mean(dim=0)
        y = y - y.mean(dim=0)
        

        std_x = tc.sqrt(x.var(dim=0) + 0.0001)  # variance (1)
        std_y = tc.sqrt(y.var(dim=0) + 0.0001)
        std_loss = tc.mean(F.relu(1 - std_x)) / 2 + tc.mean(F.relu(1 - std_y)) / 2

        cov_x = (x.T @ x) / (len(x) - 1) # covariance (3)
        cov_y = (y.T @ y) / (len(y) - 1)
        cov_loss = off_diagonal(cov_x).pow_(2).sum().div(
            256) + off_diagonal(cov_y).pow_(2).sum().div(256)

        
        loss = (25. * repr_loss + 25. * std_loss + 1. * cov_loss)
        return loss

    def dal(self, x,patient_y):

        x = self.backbone(x)
        reverse_x =  ReverseLayerF.apply(x, 1.0)

        y_hat = self.classifier(reverse_x)
                
        criterion = nn.CrossEntropyLoss()
        loss = criterion(y_hat, patient_y,)


        return loss




ds = DS(dat, dat_raw['Pseudo'])
model = Model(dat.shape[1], 5000, npatients = ds.unique_patients.shape[0])
device = tc.device('cuda:1')

dl = DataLoader(ds, batch_size = 1, shuffle = True)


optimizer = Adam(model.parameters(), lr=1e-4)
model.to(device).train()
tc.set_num_threads(5)
for epoch in range(10001):
    do_dal = False
    cum_loss = 0
    cum_dal_loss = 0
    cum_contrastive_loss = 0
    
    for X, patient in tqdm(dl):
        optimizer.zero_grad()
        
        X = X.to(device).squeeze()
        patient = patient.to(device).squeeze()

        noise1 = tc.cat((X[1:,:].clone(),X[0,:].clone().unsqueeze(0)), axis= 0)
        noise2 = tc.cat((X[-1,:].clone().unsqueeze(0) , X[:-1,:].clone()), axis=0)

        ps = tc.rand(2) * 0.9 #dropout between 0.0 and 0.9
        dropout1 = nn.Dropout(ps[0])
        dropout2 = nn.Dropout(ps[1])

        X1 = dropout1.forward(X.clone()) #* 0.9 + 0.1 * noise1
        X2 = dropout2.forward(X.clone()) #* 0.9 + 0.1 * noise2

        contrastive_loss = model(X1, X2)

        if do_dal:
            dal_loss = model.dal(X, patient)*1

        loss = contrastive_loss + dal_loss if do_dal else contrastive_loss

        loss.backward()
        optimizer.step()
        cum_loss += loss.item()
        if do_dal:
            cum_dal_loss += dal_loss.item()
        cum_contrastive_loss += contrastive_loss.item()
    

    if epoch%500==0:
        tc.save(model.state_dict(), './save/vicreg/model_per_patient_batchnorm_highinputdropout' + str(epoch) + '.pt' )
    print(epoch, cum_loss, cum_contrastive_loss, cum_dal_loss)
    
        

