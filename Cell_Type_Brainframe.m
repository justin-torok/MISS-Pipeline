function Cell_Type_Brainframe(outstruct,idx,typeinds,savenclose,view_,study,directory)
if nargin < 7
    directory = [cd filesep 'MatFiles'];
    if nargin < 6
        study = 'tasic';
        if nargin < 5
            view_ = [];
            if nargin < 4
                savenclose = 0;
            end
        end
    end
end
load([directory filesep 'PresetInputs.mat'],'GENGDmod','structList',...
        'structIndex','nonzerovox');
if strcmp(study,'zeisel')
    load([directory filesep 'Zeisel_extract.mat'],'classkey');
elseif strcmp(study,'tasic')
    load([directory filesep 'PresetInputs.mat'],'classkey');
end
load([directory filesep 'input_struct_voxelrender'],'input_struct');

input_struct.savenclose = savenclose;

for j = 1:length(typeinds)
    testvals = outstruct(idx).corrB(:,typeinds(j));
    datamap = zeros(size(GENGDmod));
    for i = 1:length(structList)
        [~,voxinds] = ismember(structIndex{i},nonzerovox);
        allvox = nonzerovox(voxinds);
        curvox = testvals(voxinds);
        datamap(allvox) = curvox;
        clear allvox curvox voxinds
    end
    interpmap = datamap;
    interpmap = imresize3(interpmap,[133 81 115]);
    interpmap(interpmap<0) = 0;
    datinput = interpmap;
    input_struct.data = datinput;
    input_struct.img_labels = classkey{typeinds(j)};
    if ~isempty(view_) && savenclose == 1
        input_struct.savenclose = 0;
        figure;
        brainframe(input_struct);
        view(view_);
        print([input_struct.img_labels '_customview'],'-dtiffn');
        close
        input_struct.savenclose = 1;
        figure;
        brainframe(input_struct);
    elseif ~isempty(view_) && savenclose == 0
        figure;
        brainframe(input_struct);
        view(view_)
    else
        figure;
        brainframe(input_struct);
    end
end
end