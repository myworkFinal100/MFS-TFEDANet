As the proposed framework utilizes the Fourier Decomposition Method (FDM), 
the corresponding FDM implementation was adopted from the original source code. 
Please cite the original FDM work if this code or repository is used in your research: 
"P. Singh, S. D. Joshi, R. K. Patney, K. Saha, The Fourier Decomposition Method for 
Nonlinear and Non-Stationary Time Series Analysis, Proceedings of the Royal Society A: 
Mathematical, Physical and Engineering Sciences 473 (2202) (2017) 20160871. 
doi:10.1098/rspa. 2016.0871. URL http://dx.doi.org/10.1098/rspa.2016.0871"

### Code Files

* ExampleCode.m
  Contains the complete end-to-end pipeline using the already trained **TFEDANet** model.
  Download the required files and run the script to reproduce the results.

* WeightedComputation_Performance_Of_OA_DetectionMethod.m
  Computes the weights of individual features and evaluates the performance of the
  **OA-contaminated EEG channel detection** method. The obtained feature weights
  are kept fixed throughout the study.

* SS_DatasetCode
  Contains the complete code for the **MFS-TFEDANet** study.
  The code can be executed by providing the required dataset details.

