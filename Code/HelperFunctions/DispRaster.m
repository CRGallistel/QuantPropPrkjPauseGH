function DispRaster(tsdat,CSdur,Xlm1,Xlm2)
% displays rasters for a specified CS duration with pause and then Xlm1 and
% the Xlm2
if CSdur>1
    CSdur = CSdur/1000;
end
figure
TSraster(tsdat,{[20 20] [20 inf]},[20 0;40 0],['r+';'k.'])
E = evalin('caller','sub');
title(['Electrode Index ' num2str(E)])
pause
hold on
xlim(Xlm1)
plot([.2 .2],ylim,'b',[.2+CSdur .2+CSdur],ylim,'b')
pause
xlim(Xlm2)
pause
close all