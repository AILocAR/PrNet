import torch
import time
import statistics
import matplotlib.pyplot as plt
from torch import nn
from d2l import torch as d2l
from DataPreprocessing_PrNet_parallel import data_preprocessing
torch.set_default_tensor_type(torch.DoubleTensor)
torch.autograd.set_detect_anomaly(True)


def train_gnss_net(net, data_iter, lr, num_epochs, num_iterations, device):
    """Train a model with GNSS moving horizon DAE."""

    # Delegate computation to CPU or GPU
    net.to(device)

    
    # for param in net.parameters():
    #     print(type(param), param.size())
    # print(sum([p.numel() for p in net.parameters()]))
    # assert 1==0
    # Determine the optimizer
    optimizer = torch.optim.Adam(net.parameters(), lr=lr)
    
    # Define loss function
    loss = nn.MSELoss()

    # Set the neural network to training mode
    net.train()

    # Set figure to plot training loss
    animator = d2l.Animator(xlabel='epoch', ylabel='Training Loss (m$^2$)')
    # animator1 = d2l.Animator(xlabel='epoch', ylabel='loss', xlim=[0, 200], ylim=[0, 1000])
    time_step = 0 
    validation_error_sum = 0

    # Evaluate training time
    time_sum = []

    # Training epoch by epoch
    for epoch in range(num_epochs):
        # timer = d2l.Timer()
        # metric = d2l.Accumulator(2)  # Sum of training loss, no. of tokens
        
        # Determine the optimizer
        # if epoch == 99 or epoch == 199 or epoch == 299 or epoch == 399:
        # if epoch == 999 or epoch == 1999:
        #     decay_lr = (epoch+1)//1000
        #     optimizer = torch.optim.Adam(net.parameters(), lr=lr/(10**decay_lr))        
        for batch in data_iter:           
            optimizer.zero_grad()
            start_time = time.time()
            
            time_step = time_step + 1

            # Read a batch of training data and delegate the data to our device
            # Shape of X: (batch_size, PRN_size, input_size)
            # Shape of Y: (batch_size, PRN_size, res_size)
            X, Y, num_batches = [x for x in batch]
            X = X.to(device)
            Y = Y.to(device)
            
            # Shape of enc_x: (batch_size, PRN_size, input_size)
            enc_x = X

            # Data Preprocessing
            # Shape of post_enc_x:       ('batch_size', 'PRN_size', 'input_feature_size')
            # Shape of valid_prn_index:  ('batch_size', 'PRN_size')
            post_enc_x, valid_prn_index = data_preprocessing(enc_x, device)

            # Pass input data through the neural network
            # Shape of dec_y_scaled:    ('batch_size', 'PRN_size', 1)
            # Shape of dec_y:           ('batch_size', 'PRN_size', 1)
            # Shape of total_prm_error: ('batch_size', 'PRN_size', 1)
            total_prm_error = net(post_enc_x)
            
            # Loss Computation
            # Shape of broadcast_index_valid:  ('batch_size', 'PRN_size',1)
            broadcast_index_valid = valid_prn_index.unsqueeze(-1)

            # Compute the error of clock bias estimation
            # # Training using unsmoothed pseudoranges
            # J = loss((total_prm_error-torch.bmm(post_enc_x[:, :, 20:21].permute(0,2,1), total_prm_error))[broadcast_index_valid]
            #             , post_enc_x[:, :, 21:22][broadcast_index_valid])
            
            # # Training using unsmoothed pseudoranges without considering receiver clock estimating residuals
            # J = loss(total_prm_error[broadcast_index_valid]
            #             , post_enc_x[:, :, 21:22][broadcast_index_valid])
            
            # Training using smoothed pseudoranges
            J = loss((total_prm_error-torch.bmm(post_enc_x[:, :, 20:21].permute(0,2,1), total_prm_error))[broadcast_index_valid]
                        , post_enc_x[:, :, 16:17][broadcast_index_valid])

            # # Training using smoothed pseudoranges without considering receiver clock estimating residuals
            # J = loss(total_prm_error[broadcast_index_valid]
            #             , post_enc_x[:, :, 16:17][broadcast_index_valid])

            if time_step < num_batches:
                # Backward Gradient Descent
                J.sum().backward()
                d2l.grad_clipping(net, 1)
                optimizer.step()
            
            end_time = time.time()
            elapsed_time = end_time - start_time
            time_sum.append(elapsed_time)
            

        if (epoch + 1) % 1 == 0:
            animator.add(epoch + 1, [J.cpu().detach().numpy()])
            # animator1.add(epoch + 1, [J1.cpu().detach().numpy()])
            validation_error_sum = validation_error_sum+J.cpu().detach().numpy()
            

        # if (epoch + 1) % 1 == 0:
        #     print('Epoch', epoch+1, 'is done. ', 'Loss is',J.cpu().detach().numpy())
        # if (epoch + 1) % 100 == 0:           
        #     filename = 'PrNet_2023V_' + str(epoch+1+1900) + '.tar'       
        #     torch.save({
        #         'model_state_dict': net.state_dict(),
        #         'optimizer_state_dict': optimizer.state_dict(),
        #         }, filename)
    elapsed_time_per_batch = sum(time_sum)/len(time_sum)
    print('Training time per batch: ', elapsed_time_per_batch)
    print('Validation Error: ', validation_error_sum/num_epochs)
    return optimizer
