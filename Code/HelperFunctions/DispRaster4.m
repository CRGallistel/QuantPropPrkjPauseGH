function DispRaster3(tsdat,CSdur,Bs,Es,Xlm1,Ax)
% displays rasters for a specified CS duration with Xlm1 and CS delimited
% and pause beginnings and endings marked. Bs is the 'PsOn' field
% and Es the 'PsOff' field at the Session level
CSdur = CSdur/1000;
% figure
TSraster(tsdat,{[20 20] [20 inf]},[20 0;40 0],Ax,['r+';'k.'])
    E = evalin('caller','sub');
    C = evalin('caller','ses');
    title(['CS' num2str(1000*CSdur) ': Sub ' num2str(E) '; Cell ' num2str(C)])
hold on
xlim(Xlm1)
plot([.2 .2],ylim,'b',[.2+CSdur .2+CSdur],ylim,'b')
%% Plotting pause beginnings and endings on raster plot
for ll = 1:length(Bs)
    plot(.2+Bs(ll),ll,'*g')
end

for ll = 1:length(Es)
    plot(.2+Es(ll),ll,'*r')
end