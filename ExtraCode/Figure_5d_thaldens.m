function Figure_5d_thaldens(outstruct,idx,savenclose,directory)
if nargin < 4
    directory = [cd filesep 'MatFiles'];
    if nargin < 3
        savenclose = 0;
    end
end
load([directory filesep 'PresetInputs.mat'],'classkey');  
regsums = sum(outstruct(idx).Bsums);
regsums = repmat(regsums,426,1);
ct_regprop = outstruct(idx).Bsums ./ regsums;
ct_regprop([12 12+213],:) = [];
amg = 1:11; cer = 12:22; sub = 23:25; hip = 26:36; hyp = 37:56; neo = 57:94;
med = 95:119; mid = 120:140; olf = 141:148; pal = 149:156; pns = 157:169;
str = 170:177; tha = 178:212;
sumprop(1,:) = sum(ct_regprop(amg,:));
sumprop(2,:) = sum(ct_regprop(cer,:));
sumprop(3,:) = sum(ct_regprop(sub,:));
sumprop(4,:) = sum(ct_regprop(hip,:));
sumprop(5,:) = sum(ct_regprop(hyp,:));
sumprop(6,:) = sum(ct_regprop(neo,:));
sumprop(7,:) = sum(ct_regprop(med,:));
sumprop(8,:) = sum(ct_regprop(mid,:));
sumprop(9,:) = sum(ct_regprop(olf,:));
sumprop(10,:) = sum(ct_regprop(pal,:));
sumprop(11,:) = sum(ct_regprop(pns,:));
sumprop(12,:) = sum(ct_regprop(str,:));
sumprop(13,:) = sum(ct_regprop(tha,:));
labkey = classkey;
labkey{9} = 'L2+3-IT';
labkey{11} = 'L5-IT';
labkey{12} = 'L5-PT';
labkey{13} = 'L6-CT';
labkey{14} = 'L6-IT';

figure('Position',[0 0 600 150]);
imagesc(sumprop(end,:),[0.01 0.1]); colormap copper;
xticks(1:25);
xticklabels(labkey);
xtickangle(90);
yticks([]);
set(gca,'TickLength',[0 0]);
set(gca,'FontSize',16);
title('CT Concentration in Thalamus','FontSize',18);

if savenclose
    print('Figure_5d_thaldens','-dtiffn');
    close
end
end