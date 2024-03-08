import torch
import statistics
from torch import nn
from d2l import torch as d2l
torch.set_default_tensor_type(torch.DoubleTensor)


def data_preprocessing(enc_x, device):
    # Shape of enc_x: ('batch_size', 'PRN_size', 'input_size')
    # If a satellite is visible at the current epoch, this satellite will be counted.
    # Shape of valid_prn_index: (batch_size, PRN_size)
    valid_prn_index = enc_x[:, :, 1] != 0
    
    # Features of inputs: 
    # Sequence of post-preprocessing inputs
    post_enc_x_seq = []

    # 0. CN0 is scaled by 50 dBHz
    # Shape of CN0: ('batch_size', 'PRN_size', 1)
    post_enc_x_seq.append(enc_x[:, :, 8].unsqueeze(-1) / 50)

    # 1. sin and 2. cos of elevation and azimuth are calculated
    # Shape of sinE: ('batch_size', 'PRN_size', 1)
    post_enc_x_seq.append(torch.sin(enc_x[:, :, 6]).unsqueeze(-1))
    post_enc_x_seq.append(torch.cos(enc_x[:, :, 6]).unsqueeze(-1))
    # post_enc_x_seq.append(torch.sin(enc_x[:, :, :, 7]).unsqueeze(-1))
    # post_enc_x_seq.append(torch.cos(enc_x[:, :, :, 7]).unsqueeze(-1))

    # 3. PRN is normalized by 32
    # Shape of PRN: ('batch_size', 'PRN_size', 1)
    post_enc_x_seq.append(enc_x[:, :, 1:2]/32)

    # 4, 5, 6. Satellite positions
    post_enc_x_seq.append(enc_x[:, :, 2:5])

    # 7, 8, 9. WLS-based estimation of longitude
    post_enc_x_seq.append(enc_x[:, :, 22:25]/torch.tensor([180, 60, 60]).to(device))

    # 10, 11, 12. WLS-based estimation of latitude
    post_enc_x_seq.append(enc_x[:, :, 25:28]/torch.tensor([90, 60, 60]).to(device))

    # 13, 14, 15. Unit geometry vector computed with WLS-based estimation of positions
    post_enc_x_seq.append(enc_x[:, :, 28:31])

    # 16. Smoothed Pseudorange errors
    post_enc_x_seq.append(enc_x[:, :, 34:35])
    # 16. Unsmoothed Pseudorange errors
    # post_enc_x_seq.append(enc_x[:, :, 31:32])

    # 17, 18, 19. Heading of smartphones: Unit vector in NED coordinate system
    post_enc_x_seq.append(enc_x[:, :, 35:38])

    # 20. The last row of matrix H
    post_enc_x_seq.append(enc_x[:, :, 38:39])

    # 21. Unsmoothed Pseudorange errors
    post_enc_x_seq.append(enc_x[:, :, 31:32])
        
    # Concatenate all input features
    # Shape of post_enc_x: (batch_size, PRN_size, input_feature_size)
    # Shape of valid_prn_index: (batch_size, PRN_size)
    return torch.cat(post_enc_x_seq, dim=-1), valid_prn_index


