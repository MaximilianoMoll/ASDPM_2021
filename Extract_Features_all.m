function Feat = Extract_Features(RR,Mu,Var,pVLF,pLF,pHF,BAL,pTOT,phase_lbl)

% Copyright (C) Maximiliano Mollura and Riccardo Barbieri, 2019-2020.
% All Rights Reserved. See LICENSE.TXT for license details.
% maximiliano.mollura@polimi.it


pLFn = pLF./(pTOT-pVLF);
pHFn = pHF./(pTOT-pVLF);

featnames = {'AVNN','SDNN','NN20','NN50','pNN20','pNN50',...
    'RMSSD','SD1','SD2','SD_ratio','avgMU','avgVAR','stdVAR','avgVLF',...
    'stdVLF','avgLF','stdLF','avgHF','stdHF','avgBAL','stdBAL','avgLFn',...
    'stdLFn','avgHFn','stdHFn','class'};

Feat = [];

% Time Domain
[AVNN] = get_AVNN(RR);
[SDNN] = get_SDNN(RR);
[NN20, pNN20] = get_NN20(RR);
[NN50, pNN50] = get_NN50(RR);
[RMSSD] = get_RMSSD(RR); % Output in mseconds
[SD1, SD2, SD_ratio]=get_poincare(RR);

feats = {AVNN,SDNN,NN20,NN50,pNN20,pNN50,RMSSD,SD1,...
    SD2,SD_ratio,mean(Mu,'omitnan'),mean(Var,'omitnan'),std(Var,'omitnan'),...
    mean(pVLF,'omitnan'),std(pVLF,'omitnan'),mean(pLF,'omitnan'),std(pLF,'omitnan'),...
    mean(pHF,'omitnan'),std(pHF,'omitnan'),mean(BAL,'omitnan'),std(BAL,'omitnan'),...
    mean(pLFn,'omitnan'),std(pLFn,'omitnan'),mean(pHFn,'omitnan'),std(pHFn,'omitnan'),...
    phase_lbl};

Feat = [Feat;cell2table(feats,'VariableNames',featnames)];


end
