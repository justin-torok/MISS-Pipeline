# MISS-Pipeline
[Chris Mezias](https://github.com/chm2062) & [Justin Torok](https://github.com/justin-torok)

The set of code necessary to run the MISS (Matrix Inversion with Subset Selection) pipeline to extract cell counts from scRNAseq and ISH transcriptomic data. Full methodological details can be found in our [preprint](https://www.biorxiv.org/content/10.1101/833566v1).  

## 1. Setup
All code is written in MATLAB and requires version 2018b or later.

Step 1: Clone this repository into a local directory of your choice.<br>
Step 2: Download the following [Zip File](https://drive.google.com/file/d/1fznvXSCCQNr3YYYqOLi0TfOzIKJ9oAgg/view?usp=sharing) and unpack its contents in the local copy of the repository, which contains all of the data dependencies necessary to run the code.

## 2. Files
- `brainframe.m` (input_struct): Tool that plots regional or voxel-wise densities on a 3-D rendering of the brain, using MATLAB's built-in isosurface and point cloud functionalities. 
    - input_struct: MATLAB struct with pre-defined fields that are used to set the visualization parameters (Required).
- `Cell_Type_Brainframe.m` (outstruct,idx,types,*savenclose*,*view_*,*directory*): Wrapper function that calls `brainframe.m` with an input_struct that is optimal for visualizing voxel-wise cell counts in the mouse brain. 
    - outstruct: MATLAB struct that is output by either `nG_ParameterFitter.m` or `lambda_ParameterFitter.m` and contains the inferred cell counts per cell type for an array of parameters (Required)
    - idx: numeric index specifying input parameter in outstruct (Required) 
    - types: numeric array of indices indicating which cell types to plot, with the indices corresponding to the order of cell types in Preset_Inputs.mat\classkey (Required)
    - savenclose: logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure (Optional; default "0")
    - view_: 1x3 numeric array specifying an argument to MATLAB view(). If supplied with savenclose = 1, this custom view is saved along with the three on-axis views (Optional; default "[]")
    - directory: character array indicating the file path of the MatFiles folder (Optional; default "[cd filesep MatFiles]")
- `Cell_Type_Data_Extract.m`
