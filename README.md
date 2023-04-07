# <em>PrNet</em>: A Neural <em>Net</em>work for Correcting <em>P</em>seudo<em>r</em>anges to Improve Positioning with Android Raw GNSS Measurements
PrNet is a neural network 🤖 for correcting pseudoranges to improve positioning with Android 📱 raw GNSS 🛰️ measurements. This repository includes the preprocessing code, the code of PrNet, and the evaluation data set. 

## Preprocessing Android Raw GNSS Measurements (Coming soon)
The preprocessing code (MATLAB) is based on Google's open-source [gnss measurement tools](https://github.com/google/gps-measurement-tools). Through preprocessing, we extract input features as well as label the training data. Compared to Google's original code, we add more functionalities, including:
* Moving horizon estimation-based positioning
* Extended Kalman filter-based positioning
* Rauch–Tung–Striebel (RTS) smoother-based positioning

👩🏽‍💻 The code is put under PrNet/Preprocessing.
* Set the directory of your Android raw GNSS data file in ProcessGnssMeasScriptPrNet.m, e.g.:
      
      `dirName ='../Data/GSDC2021/Route1/2020-05-29-US-MTV-1'`;
* Specify the name of your Android raw GNSS data file in ProcessGnssMeasScriptPrNet.m, e.g.:

      `prFileName = 'Pixel4_GnssLog.txt'`;
* Specify the name of the ground truth data file, e.g.:

      `gtNmeaFileName = 'SPAN_Pixel4_10Hz.nmea'`;
* Run ProcessGnssMeasScriptPrNet.m to process Android raw GNSS measurements;
* More details can be found in [gnss measurement tools](https://github.com/google/gps-measurement-tools). 

Then, the processed files contain the input features and labels and can be found in the directory you just set. While one of them has a header, the other one only consists of data. For example:

      `SvPVT3D_Error_label_dynamic_2020-05-29-US-MTV-1.csv` with a header
      `SvPVT3D_Error_label_dynamic_data_2020-05-29-US-MTV-1.csv` without headers

## Data Set (Coming soon)
We use the open data set for [Google Smartphone Decimeter Challenge (GSDC) 2021](https://www.kaggle.com/competitions/google-smartphone-decimeter-challenge/overview) to evaluate our method. Most data files were collected in Moutain view. Therefore, the data captured along the following two routes are selected for evaluation:

* 🛣️ Route 1: 🚗 from San Bruno to Mountain View along Interstate 280 (I-280) highway
* 🛣️ Route 2: 🚗 from Brisbane to Mountain View along U.S. Highway 101 
<img src="AllRoutes.png" width="400" height="300">

🧑🏼‍💻 The original GSDC 2021 data files are put under:

      `PrNet/Data/GSDC2021/Route1 or Route2`
The preprocessed data files are put under:

      `PrNet/Data/Route1 or Route2` 


## PrNet Implementation (Coming soon)
PrNet is based on a simple Multilayer perceptron (MLP) structure and implemented using PyTorch and d2l libraries. The related code is included under:

      `PrNet/DNN`
And the weights we trained are stored in:

      `PrNet/DNN/Weights/Route1 or Route2`
👨🏻‍💻 Our code was developed in a `conda` environment running on Ubuntu 18.04.6 LTS. To create the `conda` environment, use: 

      `conda env create -f environment.yml`
If you wanna train it by yourself, run the cells in `PrNet/DNN/PrNet_MultipleFile_parallel.ipynb` in turn. Specifically,
* Set the directory for training data files, e.g.,

      `training_data_dir = "../Data/Route2/Training/"`
* Tune the number of training epochs and learning rate in the cell "Training Process", e.g.,

      `num_epochs, lr = 100, 0.001`
* Save the weights of PrNet in the cell "Save Weights":

      `torch.save({
            'model_state_dict': net.state_dict(),
            'optimizer_state_dict': optimizer.state_dict(),
            }, 'PrNet_Layer20_H40_heading_400.tar')`
* After training the neural network, set the directory for the testing data file in the cell "Evaluation Process", e.g.,

      `data_file_eval = "../Data/Route1/Testing/SvPVT3D_Error_label_dynamic_2020-05-14-US-MTV-1.csv"`
* Load the weight file, e.g.,

      `checkpoint = torch.load('Weights/Route2/PrNet_Layer20_H40_heading_400.tar')
       model_eval.load_state_dict(checkpoint['model_state_dict'])`
* The predicted pseudorange bias will be logged into a .csv file, e.g.,

      `prm_bias_np_df.to_csv('PrM_Bias_2020-05-14-US-MTV-1.csv')`
Put the predicted pseudorange bias file to `PrNet\Preprocessing` and modify the following lines of code.
* In GpsWlsPvtEKF.m:

      `#101 GT_data = load('PrM_Bias_2020-05-14-US-MTV-1.csv');
       #216 prM = prM - GT_data(index_GT,end);
       #434 prM = prM - GT_data(index_GT,end);`
* In MHEstimator.m:

      `#50  prM = prM - GT_data(index_GT,end);`
* Run ProcessGnssMeasScriptPrNet.m to process the corrected pseudoranges and get the positioning results.
