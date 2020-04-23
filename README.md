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
- `Cell_Type_Data_Extract.m`: Data preprocessing function that outputs average expression profiles per cell type
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
        - **input_struct**: MATLAB struct with pre-defined fields that are used to set the visualization parameters
    - ***Outputs***:
        - None
- `Cell_Type_Brainframe.m`: Wrapper function that calls `brainframe.m` with an input_struct that is optimal for visualizing voxel-wise cell counts in the mouse brain. Used to generate figure panels 3a, 4c, and 5e in the manuscript. 
    - ***Inputs***:
        - **outstruct**: MATLAB struct that is output by either `nG_ParameterFitter.m` or `lambda_ParameterFitter.m` and contains the inferred cell counts per cell type for an array of parameters
        - **idx**: numeric index specifying input parameter in outstruct 
        - **typeinds**: numeric array of indices indicating which cell types to plot, with the indices corresponding to the order of cell types in Preset_Inputs.mat\classkey
        - savenclose (default "0"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
        - view_ (default "[]"): 1x3 numeric array specifying an argument to MATLAB view(). If supplied with savenclose = 1, this custom view is saved along with the three on-axis views 
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - None
 - `Density_Slice_Maps.m`: Function that plots voxel-wise cell type densities on user-selected coronal slices. Used to generate figure panel 2d in the manuscript.
    - ***Inputs***:
        - **outstruct**: MATLAB struct that is output by either `nG_ParameterFitter.m` or `lambda_ParameterFitter.m` and contains the inferred cell counts per cell type for an array of parameters
        - **idx**: numeric index specifying input parameter in outstruct 
        - **typeinds**: numeric array of indices indicating which cell types to plot, with the indices corresponding to the order of cell types in Preset_Inputs.mat\classkey
        - slicelocs (default "[25,30,47]"): numeric array of indices indicating which coronal slices to render
        - savenclose (default "0"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - None
 - `Figure_1a_tsne.m`: Function that generates the tSNE representation of the AIBS scRNAseq data, as shown in panel 1a of the manuscript. (Note: may appear different to the panel in the manuscript due to random seeding)
    - ***Inputs***:
        - savenclose (default "0"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - None
- `Figure_1b_methodshistograms.m`: Function that generates histograms of differential expression score by subset selection method, as shown in panel 1b of the manuscript.
    - ***Inputs***:
        - **method**: character array indicating which subset selection method to use
        - **ngen_param**: scalar indicating the cutoff for gene inclusion
        - **typeinds** (required only for method = "DBSCAN"): numeric array of indices indicating which cell types to plot, with the indices corresponding to the order of cell types in Preset_Inputs.mat\classkey
        - savenclose (default "0"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - None
- `Figure_2b_methodcomp_sumfit.m`: Function that generates line plots of sum fit score with respect to number of genes for each method, as shown in panel 2b of the manuscript. (Note: although not directly generated by `nG_ParameterFitter.m`, the results in the `methodcomp_output_final.mat` data file loaded by this function could be regenerated using `nG_ParameterFitter.m` using the appropriate n*G* parameter ranges)
    - ***Inputs***:
        - savenclose (default "0"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - None      
 - `Figure_2c_bigscatter.m`: Function that creates a scatterplot of inferred vs. empirical cell counts for all studies with regional quantification used in this manuscript (see the above comments for `CorrelationsCalc.m`), as shown in panel 2c.
    - ***Inputs***:
        - **outstruct**: MATLAB struct that is output by either `nG_ParameterFitter.m` or `lambda_ParameterFitter.m` and contains the inferred cell counts per cell type for an array of parameters
        - **idx**: numeric index specifying input parameter in outstruct 
        - savenclose (default "0"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - None
- `Figure_2d_correlationmaps.m`: Function that plots voxel-wise correlations between cell type expression profiles and ISH expression on user-selected coronal slices following methodology described in [Zeisel *et al*, 2018](https://www.cell.com/cell/fulltext/S0092-8674(18)30789-X?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS009286741830789X%3Fshowall%3Dtrue), with the gene set used subject to subset selection. Used to generate figure panel 2d in the manuscript.
    - ***Inputs***:
        - **typeinds**: numeric array of indices indicating which cell types to plot, with the indices corresponding to the order of cell types in Preset_Inputs.mat\classkey
        - **method**: character array indicating which subset selection method to use
        - **ngen_param**: scalar indicating the cutoff for gene inclusion
        - lambda (default "150"): scalar indicating the lambda value for MRx3-based subset selection
        - slicelocs (default "[25,30,47]"): numeric array of indices indicating which coronal slices to render
        - savenclose (default "0"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - None
- `Figure_3b_interneuron.m`: Function that creates scatterplots of inferred vs. empirical cell counts for interneurons (using data from [Kim *et al*, 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5870827/)), as shown in panel 3b of the manuscript.
    - ***Inputs***:
        - **outstruct**: MATLAB struct that is output by either `nG_ParameterFitter.m` or `lambda_ParameterFitter.m` and contains the inferred cell counts per cell type for an array of parameters
        - **idx**: numeric index specifying input parameter in outstruct 
        - savenclose (default "0"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - None
- `Figure_3c_interneuron_ratio.m`: Function that creates the bar plot and associated heat map of relative interneuron density per region, as shown in panel 3c of the manuscript.
    - ***Inputs***:
        - **outstruct**: MATLAB struct that is output by either `nG_ParameterFitter.m` or `lambda_ParameterFitter.m` and contains the inferred cell counts per cell type for an array of parameters
        - **idx**: numeric index specifying input parameter in outstruct 
        - slicelocs (default "[25,30,47]"): numeric array of indices indicating which coronal slices to render
        - savenclose (default "0"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - None
- `Figure_3de_glia.m`: Function that creates scatterplots of inferred vs. empirical cell counts for microglia (using data from [Keller *et al*, 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6205984/)) as well as the bar plot of mean glia density per major region group, as shown in panels 3d and 3e of the manuscript.
    - ***Inputs***:
        - **outstruct**: MATLAB struct that is output by either `nG_ParameterFitter.m` or `lambda_ParameterFitter.m` and contains the inferred cell counts per cell type for an array of parameters
        - **idx**: numeric index specifying input parameter in outstruct 
        - savenclose (default "0"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - None
- `Figure_4ab_taulayerslice.m`: Function that generates line plot of corrected Kendall's tau value for layer ordering per coronal slice (refer to the description of `TauCalc.m` for more details) as well as the coronal slice panel of layer-specific glutamatergic neurons in the data set. Used to generate panels 4a and 4b of the manuscript.
    - ***Inputs***:
        - **outstruct**: MATLAB struct that is output by either `nG_ParameterFitter.m` or `lambda_ParameterFitter.m` and contains the inferred cell counts per cell type for an array of parameters
        - **idx**: numeric index specifying input parameter in outstruct 
        - savenclose (default "0"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - None - 
`Figure_5ab_exinhplots.m`: Function that creates a bar-and-whisker plot of glutamatergic fraction of neurons per major region group, based on the cell types present in the data set. Used to generate panels 5a and 5b in the manuscript.
    - ***Inputs***:
        - **outstruct**: MATLAB struct that is output by either `nG_ParameterFitter.m` or `lambda_ParameterFitter.m` and contains the inferred cell counts per cell type for an array of parameters
        - **idx**: numeric index specifying input parameter in outstruct 
        - savenclose (default "0"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
    - ***Outputs***:
        - None
`Figure_5c_homogeneityindex.m`: Function that creates a bar-and-whisker plot of the homogeneity index per major region group, based on the cell types present in the data set. Used to generate panel 5c in the manuscript.
    - ***Inputs***:
        - **outstruct**: MATLAB struct that is output by either `nG_ParameterFitter.m` or `lambda_ParameterFitter.m` and contains the inferred cell counts per cell type for an array of parameters
        - **idx**: numeric index specifying input parameter in outstruct 
        - savenclose (default "0"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
    - ***Outputs***:
        - None
- `Figure_5d_thaldens.m`: Function that creates a heatmap of relative cell type density in thalamic regions, as shown in panel 5d of the manuscript.
    - ***Inputs***:
        - **outstruct**: MATLAB struct that is output by either `nG_ParameterFitter.m` or `lambda_ParameterFitter.m` and contains the inferred cell counts per cell type for an array of parameters
        - **idx**: numeric index specifying input parameter in outstruct 
        - savenclose (default "0"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - None -    
- `Figure_6a_dendrograms.m`: Function that generates dendrograms for the regional clustering by cell type similarity and the clustering of regions by developmental ontology. Used to generate panel 6a in the manuscript.
    - ***Inputs***:
        - **outstruct**: MATLAB struct that is output by either `nG_ParameterFitter.m` or `lambda_ParameterFitter.m` and contains the inferred cell counts per cell type for an array of parameters
        - **idx**: numeric index specifying input parameter in outstruct 
        - onttype (default "fore"): character array specifying whether to use forebrain regions or mid/hindbrain regions
        - savenclose (default "0"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - None - 
`Figure_6bc_aristats.m`: Function that creates the line and bar plots displaying the agreement between clusterings for every cluster number *k* using Adjusted Rand Index, as shown in panels 6b and 6c of the manuscript
    - ***Inputs***:
        - **randstruct**: MATLAB struct that is output by `Rand_Index_Calc.m` containing the relevant Rand Index-related metrics and associated cell-type-based regional clusterings
        - savenclose (default "0"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
    - ***Outputs***:
        - None
`Figure_6d_ribrainframe.m`: Function that creates `brainframe.m` renderings of the regional clustering by developmental ontology and cell type as point clouds centered on the center-of-mass of each region, with color-coding used to indicate cluster identity. Used to generate panel 6d in the manuscript.
    - ***Inputs***:
        - **randstruct**: MATLAB struct that is output by either `nG_ParameterFitter.m` or `lambda_ParameterFitter.m` and contains the inferred cell counts per cell type for an array of parameters
        - savenclose (default "1"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - None
`Figure_6d_rihistograms.m`: Function that creates histograms for each cluster number *k* showing the null distribution of Rand Index values along with the Rand Index between developmental ontology and cell type clusterings. Used to generate panel 6d in the manuscript.
    - ***Inputs***:
        - **randstruct**: MATLAB struct that is output by either `nG_ParameterFitter.m` or `lambda_ParameterFitter.m` and contains the inferred cell counts per cell type for an array of parameters
        - savenclose (default "1"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
    - ***Outputs***:
        - None