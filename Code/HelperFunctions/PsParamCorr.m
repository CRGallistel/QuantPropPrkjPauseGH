function Corrs = PsParamCorr(On,Off,Wdth,MxIsI)
LV  = ~isnan(On);
Corrs = corr([On(LV) Off(LV) Wdth(LV) MxIsI(LV)]);