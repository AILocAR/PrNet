import torch
import random
torch.set_default_tensor_type(torch.DoubleTensor)


# longSeqInput_list is a list of input files, each component is a single data file that is a list of tensors.
def seq_data_iter_in_multiple_files(longSeqInput_list, longSeqLabel_list, batch_size):  # @save
    # Create a dictionary to index each entry in all input data files
    index_all_data = []
    
    num_files = len(longSeqInput_list)

    index = 0

    for i in range(num_files):
        # Get data in the ith file
        longSeqInput = longSeqInput_list[i]
        
        # Aligh 2D index with 1D index
        for j in range(0,len(longSeqInput)):
            index_all_data.append(torch.tensor([i, j]))
            index = index + 1
    
    # Number of data in total
    total_num = index

    # Randomly shuffle the indices of all data
    data_indices = list(range(0, total_num))
    random.shuffle(data_indices)

    # Number of minibatches
    num_batches = total_num // batch_size
    
    """Generate a minibatch of data"""
    for i in range(0, batch_size * num_batches, batch_size):
        data_indices_per_batch = data_indices[i: i + batch_size]
        Xs_seq = []
        Ys_seq = []
        for j in data_indices_per_batch:
            file_index = index_all_data[j][0]
            data_in_file_index = index_all_data[j][1]

            # Find the jth data with shape: (PRN_size * Input_size) 
            Xs = longSeqInput_list[file_index][data_in_file_index]
            Ys = longSeqLabel_list[file_index][data_in_file_index]

            # Add the batch_size dimension: from (PRN_size * Input_size) to (Batch_size * PRN_size * Input_size)
            Xs_seq.append(Xs.unsqueeze(0))
            Ys_seq.append(Ys.unsqueeze(0))
        yield torch.cat(Xs_seq, 0), torch.cat(Ys_seq, 0), num_batches

 
class GnssMultipleDataFileLoader:
    """An iterator to load sequence data."""
    def __init__(self, long_seq_input, long_seq_label, batch_size):
        self.data_iter_fn = seq_data_iter_in_multiple_files
        self.long_seq_input, self.long_seq_label = long_seq_input, long_seq_label
        self.batch_size  =  batch_size

    def __iter__(self):
        return self.data_iter_fn(self.long_seq_input, self.long_seq_label, self.batch_size)
