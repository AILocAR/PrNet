# Standard API
import torch
import pandas as pd
import time
import numpy as np
from torch import nn
from d2l import torch as d2l
import matplotlib.pyplot as plt

# User defined API
from ReadingRawGnssDataset import readingRawGnssDataset
from PrNet_parallel import MlpFeatureExtractor
from PrNet_parallel import PrNet
from DataLoader_SingleFile_NoTime import GNSSSingleDataFileLoader
from Evaluate_PrNet_parallel import evaluate_gnss_net
# from Evaluate_PrNet_noMask import evaluate_gnss_net

# Set the global data type
torch.set_default_tensor_type(torch.DoubleTensor)

# A minibatch is organized in a tensor with size [batch_size * PRN_size * Input_size]
input_size = 39
# input_size = 55
PRN_size = 32
res_size = 1
label_size = 3

# %% *********************************** Reading testing data  ***********************************
data_file_eval = "../data/Dynamic/Data4QE/RouteR/Testing/SvPVT3D_Error_label_dynamic_2020-05-14-US-MTV-1.csv"
data_eval = pd.read_csv(data_file_eval)

inputs_eval, outputs_eval = torch.tensor(data_eval.iloc[:, 0:input_size].values), torch.tensor(data_eval.iloc[:, 14:input_size].values)

# A minibatch is organized in a tensor with size [batch_size * MovingWindowSize * PRN_size * Input_size]
batch_size_eval = 1

# features, pseudorange residuals, labels
x_features_eval, res_labels_eval, y_labels_eval = readingRawGnssDataset(inputs_eval, outputs_eval, input_size, res_size,
                                                                        label_size, PRN_size)  # a list of tensors
# plt.figure()
# # plt.plot(torch.stack([y_labels_eval[i][1, 0] for i in range(len(y_labels_eval))]),
# #          torch.stack([y_labels_eval[i][1, 1] for i in range(len(y_labels_eval))]), 'v')
# plt.plot(torch.stack([sum(y[y[:,0]!=0, 0])/sum(y[:,0]!=0) for y in y_labels_eval]),
#          torch.stack([sum(y[y[:,0]!=0, 1])/sum(y[:,0]!=0) for y in y_labels_eval]), 'v')         
# plt.show()

# Define the data iterator for evaluation
data_iter_eval = GNSSSingleDataFileLoader(x_features_eval, res_labels_eval, batch_size_eval)

# %%  *********************************** PrNet ***********************************
# Size of input features of Encoder's MLP
# CN0, sinE, cosE, PRN, Wls_lon*3, Wls_lat*3, Unit_geometry_vector*3, heading*3
input_size_debiasing = 16 

# Number of hidden neurons on the hidden layer of the Encoder's MLP
num_hiddens_debiasing = 40 

# Number of layers of the Encoder's MLP
num_debiasing_layers = 20 

extractor = MlpFeatureExtractor(input_size_debiasing, num_hiddens_debiasing, num_debiasing_layers, dropout = 0)
model_eval = PrNet(extractor)

# The number of trainable parameters
pytorch_total_params = sum(p.numel() for p in model_eval.parameters() if p.requires_grad)
print('Tunable parameters: ', pytorch_total_params)

# Load the trained GnssNet.tar
checkpoint = torch.load('Weights/RouteR/PrNet_Layer20_H40_heading_RouteR_500.tar')
model_eval.load_state_dict(checkpoint['model_state_dict'])

# Evaluation
model_eval.eval()
# start_time = time.time()
prm_bias = evaluate_gnss_net(model_eval, data_iter_eval, batch_size_eval, d2l.try_gpu())
# end_time = time.time()

# # Estimate the inference time per sample
# elapsed_time = (end_time - start_time)/len(y_labels_eval)
# print('Inference time per sample 2: ', elapsed_time)

# # Write denoised results to a .csv file
# prm_bias_np = prm_bias.cpu().detach().numpy()
# prm_bias_np_df = pd.DataFrame(prm_bias_np)
# prm_bias_np_df.to_csv('PrM_Bias_2020-06-05-US-MTV-1-comparision.csv')