function [SUCCESS, exp, eIDnum, phase, box, MatlabStartDate, Duration,...
                tsdata, Notes, Weight, timeunit] = LoadFredrik(filename)
% loads Fredrik's parallel fiber stimulation data earlier version loaded
% his 112-cell paw-stimulation data. This version will probably still
% perform that function, perhaps better than the older version
global Experiment
SUCCESS=0;
exp = 101; % experiment ID #
eIDnum=[];phase=[];box=[];tsdata=[];MatlabStartDate=[];Duration=[];Notes=[];Weight=[];timeunit=[];
c1 = regexp(filename,'_e\d\d\d','start')+2;
if isempty(c1) % for the 1120 & 1131 cases
    c1 = regexp(filename,'\d\d\d\d','start');
end
c2 = regexp(filename,'\d\d\d[_p]','end')-1;
c2 = c2(find(c2>c1));
eIDnum = str2num(filename(c1:c2)); % the e# in the file name
DirName = cd;
c = regexp(DirName,'ISI');
phase = str2num(DirName(c+4:c+6));
try
    S=Experiment.Subjects==eIDnum; % flags subject index #
    MatlabStartDate = now + Experiment.Subject(S).NumSessions; % the nominal
% start date is the present time + a number of days equal to the number of
% "sessions" so far loaded for this electrode
catch ME
    keyboard
end
% CA = {'penitremA_e1113pc3a_pre.mat' 'TRAM34_e1116pc1a_pre.mat' ...
%     'TRAM34_e1118pc1a_pre.mat' 'CGP_1120pc1a_pre.mat'};
% if any(strcmp(filename,CA))
%     keyboard
% end
load(filename)
Duration = Ch4.times(end);
Notes = filename;
SpkTms = Ch4.times;
MrkTms = Ch2.times;
LV = Ch2.codes(:,1)~=8; % email from Fredrik (4/21/18) explains that only
% the 8 codes are valid marker events
MrkTms(LV)=[]; % deleting bogus marker events
tsdata =sortrows([MrkTms 20*ones(length(MrkTms),1);...
    SpkTms 40*ones(length(SpkTms),1)]); % data, consisting of spike times
% (eventcode 40) and CS onset times (flagged with eventcode 20)
SUCCESS=1;

%%
% for d = 1:length(Dr)
%     if regexp(Dr(d).name,'.mat')
% %         if regexp(Dr(d).name,'980');keyboard;end
%         c1 = regexp(Dr(d).name,'_e\d\d\d','start')+2;
%         c2 = regexp(Dr(d).name,'\d\d\d[_p]','end')-1;
%         c2 = c2(find(c2>c1));
%         disp(Dr(d).name(c1:c2));
%     end
% end