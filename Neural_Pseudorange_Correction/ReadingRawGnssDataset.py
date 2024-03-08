import torch
torch.set_default_tensor_type(torch.DoubleTensor)

# 'inputs' is the input feature of training data
# 'outputs' is the label of training data
def readingRawGnssDataset(inputs, outputs, input_size, res_size, label_size, PRN_size):
    # Initialize the space for a input tensor
    InputPerStep = torch.zeros(PRN_size, input_size)

    # Initialize the space for a pseudorange residual tensor
    ResPerStep = torch.zeros(PRN_size, res_size)

    # Initialize the space for a label tensor
    LabelPerStep = torch.zeros(PRN_size, label_size)


    # Sequence of features and labels 
    featureSeq = []
    labelSeq = []
    resSeq = []


    # The index of epoch at the start of the data file
    index = inputs[0, 0]

    # Row index of the input data file
    j = 0

    for i in inputs[:,0]:
        if i == index:
            # The visible satellites at the current time step are less than 4 satellites
            # Delete this entry
            if inputs[j, 1] == 0:
                index = index + 1
            else:
                PRN_index = (inputs[j,1]-1).type(torch.long)
                InputPerStep[PRN_index,:] = inputs[j, 0:input_size] 
                ResPerStep[PRN_index,:] = outputs[j, -1]
                LabelPerStep[PRN_index,:] = outputs[j, 0:label_size]
            j = j+1 # j is the row index of input data 
        else:     
            # Wrap training data for the current time step
            featureSeq.append(InputPerStep)
            resSeq.append(ResPerStep)
            labelSeq.append(LabelPerStep)
            # Increment index by 1
            index = index + 1
            # Clear cache for input and output
            InputPerStep = torch.zeros(PRN_size, input_size)
            ResPerStep = torch.zeros(PRN_size, res_size)
            LabelPerStep = torch.zeros(PRN_size, label_size)

            # Read next time step
            # The visible satellites at the current time step are less than 4 satellites
            # Delete this entry
            if inputs[j, 1] == 0:
                index = index + 1
            else:
                PRN_index = (inputs[j, 1]-1).type(torch.long)
                InputPerStep[PRN_index, :] = inputs[j, 0:input_size]
                ResPerStep[PRN_index, :] = outputs[j, -1]
                LabelPerStep[PRN_index, :] = outputs[j, 0:label_size]
            j = j+1
    # Wrap training data for the last time step
    featureSeq.append(InputPerStep)
    resSeq.append(ResPerStep)
    labelSeq.append(LabelPerStep)
    return featureSeq, resSeq, labelSeq
