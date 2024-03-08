import torch
import random
torch.set_default_tensor_type(torch.DoubleTensor)


def eval_data_iter(longSeqInput, longSeqLabel, batch_size):  # @save
    """Generate a minibatch of subsequences using sequential selection."""
    
    # number of input data for partitioning
    num_data = len(longSeqInput)

    # Randomly shuffle the indices of all data
    data_indices = range(0, num_data)
    # random.shuffle(data_indices)

    # Number of minibatches
    num_batches = num_data // batch_size

    """Generate a minibatch of data"""
    for i in range(0, batch_size * num_batches, batch_size):
        data_indices_per_batch = data_indices[i: i + batch_size]
        Xs_seq = []
        Ys_seq = []
        for j in data_indices_per_batch:
            # Find the jth data with shape: (PRN_size * Input_size) 
            Xs = longSeqInput[j]
            Ys = longSeqLabel[j]

            # Add the batch_size dimension: from (PRN_size * Input_size) to (Batch_size * PRN_size * Input_size)
            Xs_seq.append(Xs.unsqueeze(0))
            Ys_seq.append(Ys.unsqueeze(0))
        yield torch.cat(Xs_seq, 0), torch.cat(Ys_seq, 0)


class GNSSSingleDataFileLoader:
    """An iterator to load sequence data."""
    def __init__(self, long_seq_input, long_seq_label, batch_size):
        self.data_iter_fn = eval_data_iter
        self.long_seq_input, self.long_seq_label = long_seq_input, long_seq_label
        self.batch_size = batch_size

    def __iter__(self):
        return self.data_iter_fn(self.long_seq_input, self.long_seq_label, self.batch_size)
