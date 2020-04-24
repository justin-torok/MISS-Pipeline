# MISS-Pipeline
[Chris Mezias](https://github.com/chm2062) & [Justin Torok](https://github.com/justin-torok)

The set of code necessary to run the MISS (Matrix Inversion with Subset Selection) pipeline to extract cell counts from scRNAseq and ISH transcriptomic data. Full methodological details can be found in our [preprint](https://www.biorxiv.org/content/10.1101/833566v1).  

## 1. Setup
All code is written in MATLAB and requires version 2018b or later.

Step 1: Clone this repository into a local directory of your choice.<br>
Step 2: Download the following [Zip File](https://drive.google.com/file/d/1fznvXSCCQNr3YYYqOLi0TfOzIKJ9oAgg/view?usp=sharing) and unpack its contents in the local copy of the repository, which contains all of the data dependencies necessary to run the code.

## 2. Files
Below is a short description of each of the code files contained in the MISS-Pipeline folder, grouped by general functionality in alphabetical order. The "ExtraCode" folder contains code that is either outdated or incomplete, and none of the functions in the main folder require any of those scripts to run. Scripts that are also functions have their inputs and outputs described, with required inputs in boldface text and optional inputs with their default setting in parentheses.

`MISS_demo.m` walks through the basic functionality of the MISS pipeline and demonstrates several of the plotting functions.

### Data Preprocessing
- `Cell_Type_Data_Extract.m`: Data preprocessing function that outputs average expression profiles per cell type
    - ***Inputs***:
        - **classstruct**: MATLAB struct that is output by `scRNAseq_Data_Extract.m` that contains raw scRNAseq data grouped by cell type
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - meanexprmat: An n_genes x n_types numeric array of average gene expression profiles for each cell type
        - classkey: A 1 x n_types cell array of string descriptors of each cell type, in the same order as the columns of meanexprmat
        - entrez_names: A n_genes x 1 cell array of standard gene abbreviations of each gene, in the same order as the rows of meanexprmat
- `ISH_Data_Extract.m`: Data preprocessing function that outputs the voxel-wise gene expression for the consensus genes between the scRNAseq and ISH data sets, averaging duplicate probes where necessary. 
    - ***Inputs***:
        - **classstruct**: MATLAB struct that is output by `scRNAseq_Data_Extract.m` that contains raw scRNAseq data grouped by cell type
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - regvgene: An n_genes x n_voxels numeric array of ISH gene expression per voxel, normalized by gene
- `scRNAseq_Data_Extract.m`: Data preprocessing function that groups raw scRNAseq data downloaded from the AIBS by cell type according to their given taxonomy.
    - ***Inputs***:
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - classstruct: MATLAB struct that is hierarchically organized by region and cell type taxonomy containing raw scRNAseq data
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
        - **B**: n_voxels x n_types numeric array of inferred cell type density per voxel, in arbitrary units
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - B_corrected: n_voxels x n_types numeric array of inferred cell counts per voxel
        - Bfactor: scalar multiplication factor used to correct B
- `GeneSelector.m`: Function that performs gene subset selection according to the various information-theoretic methods used in the manuscript. Called by `nG_ParameterFitter.m`, among others. See the manuscript for further details.
    - ***Inputs***:
        - **genevct**: n_genes x n_types numeric array of gene expression per cell type
        - **voxvgene**: n_genes x n_voxels numeric array of gene expression per voxel
        - **gene_names**: n_genes x 1 cell array of standard gene abbreviations of each gene, in the same order as the rows of genevct and voxvgene
        - **ngen_param**: scalar indicating the cutoff for gene inclusion
        - lambda (default "150"): scalar indicating the lambda value for MRx3-based subset selection
        - method (default "MRx3"): character array specifying the gene subset selection method to use
        - preloadinds (default "[]"): 1 x n_genes numeric array of indices ranked in the order of decreasing criterion value (e.g. genes with higher MRx3 criterion values are closer to the top); only useful from a computational efficiency perspective when the method chosen is MRx3 and multiple ngen_param values are being tested for a given lambda value (see `nG_ParameterFitter.m` line 40).
    - ***Outputs***:
        - E_red: n_genes x n_voxels numeric array of gene expression per voxel, where n_genes has been constrained by ngen_param
        - C_red: n_genes x n_types numeric array of gene expression per cell type, where n_genes has been constrained by ngen_param
        - nGen: scalar value of n_genes (size of E_red/C_red row dimension)
        - reduced_gene_names: n_genes x 1 cell array of standard gene abbreviations of each gene, in the same order as the rows of E_red and C_red, where n_genes has been constrained by ngen_param
- `lambda_ParameterFitter.m`: Function that determines the optimal (n*G*,lambda) pair for given lists of n*G* parameters and lambda values in a brute-force fashion for MRx3-based subset selection. Functions as a wrapper for `nG_ParameterFitter.m` for multiple lambda values.
    - ***Inputs***:
        - **voxvgene**: n_genes x n_voxels numeric array of gene expression per voxel
        - **genevct**: n_genes x n_types numeric array of gene expression per cell type
        - **gene_names**: n_genes x 1 cell array of standard gene abbreviations of each gene, in the same order as the rows of genevct and voxvgene
        - ng_param_list (default "100:70:800"): numeric array of n*G* cutoff values
        - lambda_param_list (default "50:50:500"): numeric array of lambda values
    - ***Outputs***:
        - lambdastruct: 1 x (length(ng_param_list) * length(lambda_param_list)) MATLAB struct with the following fields:
            - Bvals: n_voxels x n_types numeric array of cell densities per voxel (output of `CellDensityInference.m`)
            - nGen: scalar value of n*G*
            - lambda: scalar value of lambda
            - corrB: n_voxels x n_types numeric array of cell counts per voxel (output of `Density_to_Counts.m`)
            - Bsums: n_regions x n_types numeric array of cell counts per CCF region (output of `Voxel_to_Region.m`)
            - Bmeans: n_regions x n_types numeric array of cell densities per CCF region in units of counts/(0.2mm)^3 (output of `Voxel_to_Region.m`)
            - LinR: 1 x 1 MATLAB struct of Lin's concordance correlation coefficients (output of `CorrelationsCalc.m`)
            - PearsonR: 1 x 1 MATLAB struct of Pearson correlation coefficients (output of `CorrelationsCalc.m`)
            - sumfit: scalar value of the sum fit metric
        - peakind: scalar index of the (n*G*,lambda) pair in lambdastruct that has the peak sum fit value
- `mRMR_Selector.m`: Function that performs mRMR-based gene subset selection; called by `GeneSelector.m`. See manuscript for more details.
    - ***Inputs***:
        - **C**: n_genes x n_types numeric array of gene expression per type, column-normalized
        - **n**: scalar cutoff of number of genes
        - **method**: character array indicating whether the quotient or difference criterion is to be used
    - ***Outputs***:
        - geneinds: 1 x n numeric array of gene indices to include from PresetInputs.mat\entrez_names
- `MRx3_Selector.m`: Function that performs MRx3-based gene subset selection; called by `GeneSelector.m`. See manuscript for more details.
    - ***Inputs***:
        - **C_raw**: n_genes x n_types numeric array of gene expression per type
        - **E**: n_genes x n_voxels numeric array of gene expression per voxel
        - **n**: scalar cutoff of number of genes
        - **lambda**: scalar lambda parameter value
    - ***Outputs***:
        - geneinds: 1 x n numeric array of gene indices to include from PresetInputs.mat\entrez_names        
- `nG_ParameterFitter.m`: Function that determines the optimal n*G* for a given list of n*G* parameters in a brute-force fashion for a user-defined gene subset selection method.
    - ***Inputs***:
        - **voxvgene**: n_genes x n_voxels numeric array of gene expression per voxel
        - **genevct**: n_genes x n_types numeric array of gene expression per cell type
        - **gene_names**: n_genes x 1 cell array of standard gene abbreviations of each gene, in the same order as the rows of genevct and voxvgene
        - **method**: character array specifying the gene subset selection method to use
        - ng_param_list (default depends on method): numeric array of ng_param cutoff values
        - lambda (default "150"): lambda value (only relevant for MRx3 method)
    - ***Outputs***:
        - outstruct: 1 x length(ng_param_list MATLAB struct with the following fields:
            - Bvals: n_voxels x n_types numeric array of cell densities per voxel (output of `CellDensityInference.m`)
            - nGen: scalar value of n*G*
            - lambda: scalar value of lambda (only present if method used is MRx3)
            - corrB: n_voxels x n_types numeric array of cell counts per voxel (output of `Density_to_Counts.m`)
            - Bsums: n_regions x n_types numeric array of cell counts per CCF region (output of `Voxel_to_Region.m`)
            - Bmeans: n_regions x n_types numeric array of cell densities per CCF region in units of counts/(0.2mm)^3 (output of `Voxel_to_Region.m`)
            - LinR: 1 x 1 MATLAB struct of Lin's concordance correlation coefficients (output of `CorrelationsCalc.m`)
            - PearsonR: 1 x 1 MATLAB struct of Pearson correlation coefficients (output of `CorrelationsCalc.m`)
            - sumfit: scalar value of the sum fit metric
        - peakind: scalar index of the (n*G*,lambda) pair in lambdastruct that has the peak sum fit value
- `Rand_Index_Calc.m`: Function that calculates the Adjusted Rand Index for a given set of regional cell type densities. See manuscript for further details.
    - ***Inputs***:
        - **outstruct**: MATLAB struct that is output by either `nG_ParameterFitter.m` or `lambda_ParameterFitter.m` and contains the inferred cell counts per cell type for an array of parameters
        - **idx**: numeric index specifying input parameter in outstruct 
        - onttype (default "fore"): character array specifying whether to use forebrain regions or mid/hindbrain regions
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - randstruct: 1 x 1 MATLAB struct with the following fields:
            - nclust: 1 x n_clusterings numeric array of cluster numbers *k* in the AIBS developmental ontology
            - RandIndexT: 1 x n_clusterings numeric array of uncorrected Rand Index values per cluster number *k* between the cell type clustering and developmental ontology
            - RandIndexRandom: n_clusterings x iters numeric array of uncorrected Rand Index values per cluster number *k* between the random clusterings and developmental ontology
            - AdjustedRand: 1 x n_clusterings numeric array of Adjusted Rand Index values per cluster number *k* between the cell type clustering and developmental ontology given the null distribution of RI values
            - StdAboveMean: 1 x n_clusterings numeric array of numbers of standard deviations above the mean per cluster number *k* of the RI between the cell type clustering and developmental ontology relative to the null distribution of RI values
            - T: n_regions x n_clusterings numeric array of cluster identities per region
- `TauCalc.m`: Function that calculates an adjusted Kendall's rank correlation coefficient (tau) between the expected ordering of neocortical layers based on the annotation of layer-specific glutamatergic neurons in the scRNAseq dataset and the empirical ordering, as determined by ranking the distance between the cortical surface and the neocortical bands of cell density for each of these cell types. The calculation is corrected per coronal slice for off-target and insufficient signal in one or more cell types. See the manuscript for more details.
    - ***Inputs***:
        - **outstruct**: MATLAB struct that is output by either `nG_ParameterFitter.m` or `lambda_ParameterFitter.m` and contains the inferred cell counts per cell type for an array of parameters
        - **idx**: numeric index specifying input parameter in outstruct 
        - cell_names (default "{'L2n3','L4','L5IT','L5PT','L6CT','L6IT','L6b'}"): 1 x n_layer_types cell array of layer-specific glutamatergic cell names
        - cell_inds (default "9:15"): 1 x n_layer_types numeric array of indices corresponding to the location of the types in cell_names within PresetInputs.mat\classkey
        - ranks (default "[1 2 3 3 4 4 4]"): 1 x n_layer_types numeric array of expected type ordering based on functional annotation; ties are allowed
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - metric_struct: 1 x 1 MATLAB struct with the following fields:
            - paramlist: scalar index for outstruct (identical to idx)
            - tau: scalar, whole-brain adjusted tau value
            - tau_perslice: scalar, mean per-slice tau value across the whole brain (in practice very similar to tau)
            - pval: scalar, assessment of significance of the tau value
            - bic_pct: scalar, BIC value (legacy)
- `Voxel_to_Region.m`: Function that takes cell counts per voxel and groups them into cell counts and cell densities in 426 CCF regions, with right hemispheric regions constituting the first 213 entries of the resulting matrices. See listB.mat for the order of regions in our convention.
    - ***Inputs***:
        - **D**: n_voxels x n_types numeric array of cell counts per voxel (output of `Density_to_Counts.m`)
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - cell_types_sum: n_regions (426) x n_types numeric array of cell counts per CCF region
        - cell_types_mean: n_regions (426) x n_types numeric array of cell densities (in counts/(0.2mm)^3)  per CCF region
- `Voxel_to_Region_Bilateral.m`: Function that takes cell counts per voxel and groups them into cell counts and cell densities in 213 bilateral regions in the CCF. See listB.mat for the order of regions in our convention.
    - ***Inputs***:
        - **D**: n_voxels x n_types numeric array of cell counts per voxel (output of `Density_to_Counts.m`)
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - cell_types_sum: n_regions (213) x n_types numeric array of cell counts per bilateral CCF region
        - cell_types_mean: n_regions (213) x n_types numeric array of cell densities (in counts/(0.2mm)^3)  per bilateral CCF region
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
        - None
- `Figure_6bc_aristats.m`: Function that creates the line and bar plots displaying the agreement between clusterings for every cluster number *k* using Adjusted Rand Index, as shown in panels 6b and 6c of the manuscript
    - ***Inputs***:
        - **randstruct**: MATLAB struct that is output by `Rand_Index_Calc.m` containing the relevant Rand Index-related metrics and associated cell-type-based regional clusterings
        - savenclose (default "0"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
    - ***Outputs***:
        - None
- `Figure_6d_ribrainframe.m`: Function that creates `brainframe.m` renderings of the regional clustering by developmental ontology and cell type as point clouds centered on the center-of-mass of each region, with color-coding used to indicate cluster identity. Used to generate panel 6d in the manuscript. (Note: does not work for mid/hindbrain clustering visualization)
    - ***Inputs***:
        - **randstruct**: MATLAB struct that is output by `Rand_Index_Calc.m` containing the relevant Rand Index-related metrics and associated cell-type-based regional clusterings
        - savenclose (default "1"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - None
- `Figure_6d_rihistograms.m`: Function that creates histograms for each cluster number *k* showing the null distribution of Rand Index values along with the Rand Index between developmental ontology and cell type clusterings. Used to generate panel 6d in the manuscript.
    - ***Inputs***:
        - **randstruct**: MATLAB struct that is output by `Rand_Index_Calc.m` containing the relevant Rand Index-related metrics and associated cell-type-based regional clusterings
        - savenclose (default "1"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
    - ***Outputs***:
        - None
- `Gene_Expression_Slice_Maps.m`: Function that plots voxel-wise gene expression energy from the AIBS ISH data set on user-selected coronal slices. Used to generate figure panel 1a in the manuscript.
    - ***Inputs***:
        - **gene_names**: cell array of character arrays of standard gene symbols corresponding to entries of PresetInputs.mat\entrez_names
        - slicelocs (default "34"): numeric array of indices indicating which coronal slices to render
        - savenclose (default "0"): logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure
        - directory (default "[cd filesep 'MatFiles']"): character array indicating the file path of the MatFiles folder
    - ***Outputs***:
        - None