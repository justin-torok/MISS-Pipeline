# MISS-Pipeline
[Chris Mezias](https://github.com/chm2062) & [Justin Torok](https://github.com/justin-torok)

The set of code necessary to run the MISS (Matrix Inversion with Subset Selection) pipeline to extract cell counts from scRNAseq and ISH transcriptomic data. Full methodological details can be found in our [preprint](https://www.biorxiv.org/content/10.1101/833566v1).  

## 1. Setup
All code is written in MATLAB and requires version 2018b or later.

Step 1: Clone this repository into a local directory of your choice.<br>
Step 2: Download the following [Zip File](https://drive.google.com/file/d/1fznvXSCCQNr3YYYqOLi0TfOzIKJ9oAgg/view?usp=sharing) and unpack its contents in the local copy of the repository, which contains all of the data dependencies necessary to run the code.

## 2. Files
Below is a short description of each of the code files contained in the MISS-Pipeline folder, grouped by general functionality in alphabetical order. The "ExtraCode" folder contains code that is either outdated or incomplete, and none of the functions in the main folder require any of those scripts to run. Scripts that are also functions have their inputs and outputs described, with required inputs in boldface text and optional inputs with their default setting in parentheses.

### Data Preprocessing
- `Cell_Type_Data_Extract.m`: Data pre-processing function that outputs average expression profiles per cell type
    - __Inputs__:
        - **classstruct**: MATLAB struct that is output by `scRNAseq_Data_Extract.m` that contains raw scRNAseq data grouped by cell type
        - directory (default "[cd filesep MatFiles]"): character array indicating the file path of the MatFiles folder
    - __Outputs__:
        - meanexprmat: An n_genes x n_types numeric array of average gene expression profiles for each cell type
        - classkey: A 1 x n_types cell array of string descriptors of each cell type, in the same order as the columns of meanexprmat
        - entrez_names: A n_genes x 1 cell array of standard gene abbreviations of each gene, in the same order as the rows of meanexprmat

### Cell Count Inference and Analysis
- `CellDensityInference.m`: Function that performs the non-negative matrix inversion to determine cell density per voxel from ISH and scRNAseq data, after gene subsetting and data normalization. Called by `nG_ParameterFitter.m`.
    - __Inputs__:
        - **E_red**: A n_genes x n_voxels numeric array of gene expression per voxel
        - **C_red**: A n_genes x n_types numeric array of gene expression per cell type
    - __Outputs__:
        - B: A n_voxels x n_types numeric array of inferred cell type density per voxel, in arbitrary units
- `CorrelationsCalc.m`: Function that performs the non-negative matrix inversion to determine cell density per voxel from ISH and scRNAseq data, after gene subsetting and data normalization. Called by `nG_ParameterFitter.m`.
    - __Inputs__:
        - **E_red**: A n_genes x n_voxels numeric array of gene expression per voxel
        - **C_red**: A n_genes x n_types numeric array of gene expression per cell type
    - __Outputs__:
        - B: A n_voxels x n_types numeric array of inferred cell type density per voxel, in arbitrary units

### Visualization
- `brainframe.m`: Tool that plots regional or voxel-wise densities on a 3-D rendering of the brain, using MATLAB's built-in isosurface and point cloud functionalities.
    - __Inputs__:
        - **input_struct**: MATLAB struct with pre-defined fields that are used to set the visualization parameters (Required).
    - __Outputs__:
        - None
- `Cell_Type_Brainframe.m`: Wrapper function that calls `brainframe.m` with an input_struct that is optimal for visualizing voxel-wise cell counts in the mouse brain. 
    - __Inputs__:
        - **outstruct**: MATLAB struct that is output by either `nG_ParameterFitter.m` or `lambda_ParameterFitter.m` and contains the inferred cell counts per cell type for an array of parameters
        - **idx**: numeric index specifying input parameter in outstruct 
        - **types**: numeric array of indices indicating which cell types to plot, with the indices corresponding to the order of cell types in Preset_Inputs.mat\classkey
        - savenclose (default "0"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
        - view_ (default "[]"): 1x3 numeric array specifying an argument to MATLAB view(). If supplied with savenclose = 1, this custom view is saved along with the three on-axis views 
        - directory (default "[cd filesep MatFiles]"): character array indicating the file path of the MatFiles folder
    - __Outputs__:
        - None