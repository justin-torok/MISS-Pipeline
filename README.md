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
    - ***Inputs***:
        - **classstruct**: MATLAB struct that is output by `scRNAseq_Data_Extract.m` that contains raw scRNAseq data grouped by cell type
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - meanexprmat: An n_genes x n_types numeric array of average gene expression profiles for each cell type
        - classkey: A 1 x n_types cell array of string descriptors of each cell type, in the same order as the columns of meanexprmat
        - entrez_names: A n_genes x 1 cell array of standard gene abbreviations of each gene, in the same order as the rows of meanexprmat

### Cell Count Inference and Analysis
- `CellDensityInference.m`: Function that performs the non-negative matrix inversion to determine cell density per voxel from ISH and scRNAseq data, after gene subsetting and data normalization. Called by `nG_ParameterFitter.m`.
    - ***Inputs***:
        - **E_red**: A n_genes x n_voxels numeric array of gene expression per voxel
        - **C_red**: A n_genes x n_types numeric array of gene expression per cell type
    - ***Outputs***:
        - B: A n_voxels x n_types numeric array of inferred cell type density per voxel, in arbitrary units
- `CorrelationsCalc.m`: Function that performs Pearson correlations and Lin's concordance correlations between the inferred and empirically determined cell counts from various previously published studies for three groups of brain regions (neocortical regions only, forebrain regions only, and all brain regions). *Pvalb*+, *Sst*+, and *Vip*+ counts are compared against [Kim *et al*, 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5870827/); total cell counts are compared against [Murakami *et al*, 2018](https://www.nature.com/articles/s41593-018-0109-1); microglia cell counts are compared against [Keller *et al*, 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6205984/); total neuron counts are compared against [Herculano-Houzel *et al*, 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3800983/). Called by `nG_ParameterFitter.m`.
    - ***Inputs***:
        - **outstruct**: MATLAB struct that is output by either `nG_ParameterFitter.m` or `lambda_ParameterFitter.m` and contains the inferred cell counts per cell type for an array of parameters
        - **idx**: numeric index specifying input parameter in outstruct
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - LinRvals: A MATLAB struct with six fields ('pv', 'sst', 'vip', 'all', 'micro', 'neuron'), each of which contains a 1x3 numeric array of Lin's concordance correlation coefficients between MISS-inferred counts and data for neocortical regions only, forebrain (neocortical, subcortical, olfactory, and thalamic) regions only, and all brain regions, respectively
        - PeaRvals: A MATLAB struct with six fields ('pv', 'sst', 'vip', 'all', 'micro', 'neuron'), each of which contains a 1x3 numeric array of Pearson correlation coefficients between MISS-inferred counts and data for neocortical regions only, forebrain (neocortical, subcortical, olfactory, and thalamic) regions only, and all brain regions, respectively
- `Density_to_Counts.m`: Function that performs rescales the voxel-wise densities output by `CellDensityInference.m` into cell counts using a global rescaling factor, relying upon data from [Murakami *et al*, 2018](https://www.nature.com/articles/s41593-018-0109-1). Called by `nG_ParameterFitter.m`.
    - ***Inputs***:
        - **B**: A n_voxels x n_types numeric array of inferred cell type density per voxel, in arbitrary units
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - B_corrected: A n_voxels x n_types numeric array of inferred cell counts per voxel
        - Bfactor: The scalar multiplication factor used to correct B

### Visualization
- `brainframe.m`: Tool that plots regional or voxel-wise densities on a 3-D rendering of the brain, using MATLAB's built-in isosurface and point cloud functionalities.
    - ***Inputs***:
        - **input_struct**: MATLAB struct with pre-defined fields that are used to set the visualization parameters (Required).
    - ***Outputs***:
        - None
- `Cell_Type_Brainframe.m`: Wrapper function that calls `brainframe.m` with an input_struct that is optimal for visualizing voxel-wise cell counts in the mouse brain. 
    - ***Inputs***:
        - **outstruct**: MATLAB struct that is output by either `nG_ParameterFitter.m` or `lambda_ParameterFitter.m` and contains the inferred cell counts per cell type for an array of parameters
        - **idx**: numeric index specifying input parameter in outstruct 
        - **typeinds**: numeric array of indices indicating which cell types to plot, with the indices corresponding to the order of cell types in Preset_Inputs.mat\classkey
        - savenclose (default "0"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
        - view_ (default "[]"): 1x3 numeric array specifying an argument to MATLAB view(). If supplied with savenclose = 1, this custom view is saved along with the three on-axis views 
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - None
 - `Density_Slice_Maps.m`: Function that plots voxel-wise cell type densities on user-selected coronal slices. 
    - ***Inputs***:
        - **outstruct**: MATLAB struct that is output by either `nG_ParameterFitter.m` or `lambda_ParameterFitter.m` and contains the inferred cell counts per cell type for an array of parameters
        - **idx**: numeric index specifying input parameter in outstruct 
        - **typeinds**: numeric array of indices indicating which cell types to plot, with the indices corresponding to the order of cell types in Preset_Inputs.mat\classkey
        - slicelocs (default "[25,30,47]"): numeric array of indices indicating which coronal slices to render
        - savenclose (default "0"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - None