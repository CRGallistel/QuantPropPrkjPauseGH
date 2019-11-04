function [r150,r200,r300,r400,r450] = CSdurElectrodes(D)
% creates logical vectors flagging electrodes in a given CS duration group
r150 = unique(D(D(:,1)==150,2))'; % col 2 for all rows with 150 in 1st col
r200 = unique(D(D(:,1)==200,2))';
r300 = unique(D(D(:,1)==300,2))';
r400 = unique(D(D(:,1)==400,2))';
r450 = unique(D(D(:,1)==450,2))';
