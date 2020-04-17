% load PresetInputs.mat
justpresets = 0;
uniquetypes = 1;
% redinds = randperm(3855,2500);
% rednames = entrez_names(redinds);
load ExtractedHippData.mat

if justpresets
    [genevct,voxvgene,ctkey,icell_key,genes] = ProcessedData_Generator();
else
    if uniquetypes
        [genevct,voxvgene,ctkey,icell_key,genes] = ProcessedData_Generator(allcells,cell_labvec,gennames,hippkey);
    else
        [genevct,voxvgene,ctkey,icell_key,genes] = ProcessedData_Generator(allcells,cell_labvec,gennames,hippkey);
    end
end