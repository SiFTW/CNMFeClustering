# Generate outputs from CNMF_E analysis
Performs bespoke clustering and plotting on analysis from CNMF_e analysis.
Generates output from Samarasinghe et al.

# How to run
- Download and install CNMF_e from https://github.com/zhoupc/CNMF_E
- Save the output file from CNMF_e
- Run cnfme_setup.m if not already run (likely if generating outputs in a different session or machine than the one on which CNMFE was run).
- Download stacked plot from  MATLAB File Exchange and place in working folder or PATH www.mathworks.com/matlabcentral/fileexchange/24368-stacked-plot/
- Run plotOuputs.m or plotSpikeSynchronization_update.m to generate the clustering and synchronisation plots from the paper

Note some of the lines to produce movies are commented out in the code and need to be commented back in to run. This is done as they can take a long time to go fram by frame through a large tif stack.
