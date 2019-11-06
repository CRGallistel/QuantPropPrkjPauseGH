% script m-file Exper1acquisition
% The code for analyzing Dan-Anders single-ISI acquisition data. Each data
% file contains 2 or 4 structures. The two of interest are those that
% report the CS pulse times (Ch2 or Ch3, or CS_stim_pulses' the one with
% "CS" in the Title field) and those that report the spike times (Ch5,Ch7
% or Ch9 or "Spikes", the one with "Spikes" in the Title field). The ISIs
% (CS-US intervals) were either 300 ms or 200 ms, as indicated by the
% latter part of the file name, while the duration of the CS was either
% 300, 400, 600 or 800 ms, as indicated in the early part of the file name.
% When the sessions are loaded, the ISI is extracted from the file name and
% put in the Phase field. The duration of the CS may be determined from the
% cs pulse time data. For subjects (cells) 7-11, the CS persisted beyond
% the US. The US stimulation varied unsystematically between subjects:
%{
Extracted from Dan-Anders email on July 28, 2017
"Here are the US pulse times in [the] original data files."
436: 0, 2, 4, 6, 8 and 20, 22, 24, 26, 28
457: 0, 2, 4, 6, 8 and 20, 22, 24, 26, 28
474: 0, 2, 4, 6, 8 and 20, 22, 24, 26, 28  
479: 0, 2, 4, 6, 8 and 20, 22, 24, 26, 28
483: 0, 2, 4, 6, 8 and 20, 22, 24, 26, 28
485: 0, 2, 4, 6, 8 and 20, 22, 24, 26, 28
487: 0, 2, 4, 6, 8 and 20, 22, 24, 26, 28
505: 0, 2, 4, 6, 8 and 12, 14, 16, 18, 20
524: 0, 2, 4, 6, 8 and 12, 14, 16, 18, 20
536: 0, 2, 4, 6, 8 and 12, 14, 16, 18, 20
543: 0, 2, 4, 6, 8 and 12, 14, 16, 18, 20
%}


TSinitexperiment('DanAndersExperiment',101,...
    [436 457 474 479 483 485 487 505 524 536 543],'electrodes','Hesslow')
Experiment.Info.LoadFunction='LoadDanAnders';
Experiment.Info.InputTimeUnit = 1;
Experiment.Info.OutputTimeUnit=1;
Experiment.Info.FileExtension ='.mat';
TSsaveexperiment
%%
TSimporteventcodes('DanAndersEventCodes.txt')
% spike = 40;
% cspulse = 30;
% CSon = 50;
% CSoff = 60;
% USon = 70;
% USoff = 75;
% TrlOn = 80;
% TrlOff = 90;
% The USoff event was added after I obtained from Dan-Anders the actual
% times at which the US pulses were delivered, so that I knew when the US
% ended. The following code inserted the USoff events into the data
%%
TSlimit('Subjects','all')
TSlimit('Sessions','all')
USdurs =[[436.00         28.00];...
        [457.00         28.00];...
        [474.00         28.00];...
        [479.00         28.00];...
        [483.00         28.00];...
        [485.00         28.00];...
        [487.00         28.00];...
        [505.00         20.00];...
        [524.00         20.00];...
        [536.00         20.00];...
        [543.00         20.00]]; % sub ids & US durations
    
TSapplystat('TSDataA','TSDataA',@AddUSoff,USdurs,USon,USoff)
%{
function tsd=AddUSoff(tsd,USdurs,ONcode,OFFcode)
S = evalin('caller','Experiment.Subject(sub).SubId'); % current subject
LVon = tsd(:,2)==ONcode; % flags USons
dur = USdurs(USdurs(:,1)==S,2)/1000; % duration of US for current subject, in ms
tsd = [tsd;[tsd(LVon,1)+dur OFFcode*ones(sum(LVon),1)]]; % tsd(LVon,1)+dur is a
    % the vector of session times for the US offs; OFFcode*ones(sum(LVon),1) is
    % a vector of the same length whose elements are the event code for USoff
tsd = sortrows(tsd);
%}
%%
TSloadsessions
TSsaveexperiment
%% Combining the two parts of Subject 9's data (ID# = 524)
% This session had to be interrupted to start a new file; thus the data 
% come in 2 files, which appear as 2 sessions after loading
DataPart2 = Experiment.Subject(9).Session(2).TSData;
EndPart1 =  Experiment.Subject(9).Session(1).TSData(end,1);
DataPart2(:,1) = DataPart2(:,1) + EndPart1 + 3600; % Dan-Anders had to take
% a 1-hr break to make save initial data file and open a new one
Experiment.Subject(9).Session(1).TSData = ...
    [Experiment.Subject(9).Session(1).TSData;DataPart2];
% Above code concatenates data from the two parts   
% It was later discovered that the second "part" is a duplication of the
% "first" part
%{
Looking at the two raster plots for the cell/experiment e524, I realized
that I made an error when I exported the data for your analysis. The spikes
in the two raster plots are identical! The data from e524 contains a total
of 260 CS+US trials. These data were initially saved in two separate files
 - Part One contained 240 trials (one hour) and Part Two contained an
additional 20 trials when I did the experiment back in 2006. I then
combined the two data files into one, containing all 260 trials. But when
I exported the data for you, I had forgotten that I had already combined
the two, and I therefore mistakenly gave you two identical copies of the
same data set, called Part One and Part Two respectively. The solution is
simply to throw one of them out.
%}
% Hence, I subsequently deleted the second "session" for this subject
%
TSremovesession(524,2)
% Now there is only one session for each cell/subject
%% Adding TrlOn, TrlOff, CSon, CSoff and USon events. The CSon's and off's are
% computed from the cs pulses; the USon times are computed from the CS on
TSapplystat('TSDataA',{'TSData' 'Phase'},@AddEvents)
%{
function tsdA = AddEvents(tsd,ISI)
% tsd is the original data; ISI is the CS-US interval from the Phase field
LVto = diff(tsd(:,1))>5; % Logical vector flagging 1st spike in each trial,
% but displaced one event too early because of differencing
LVtoff = [LVto;true]; % logical vector flagging last spike in
% each trial
LVto = [true;LVto]; % adding onset of 1st trial to logical vector and
% displacing flags one event forward to where they belong
TrlOnTms = tsd(LVto,1)-.001; % trial onsets set to 1 ms before 1st spike time
TrlOnD = [TrlOnTms 80*ones(size(TrlOnTms))]; % trial onset data; 80 is the
% event code for trial on
TrlOffTms = tsd(LVtoff,1)+.001; % trial onsets set to 1 ms after last spike
TrlOffD = [TrlOffTms 90*ones(size(TrlOffTms))]; % trial offset data; 90 is
% the event code for trial off

CSonTms = nan(size(TrlOnTms));
USonTms = nan(size(TrlOnTms));
CSoffTms = nan(size(TrlOnTms));
i = 1;
for t = TrlOnTms'
    CSonTms(i) = tsd(find(tsd(:,1)>t & tsd(:,2)==30,1),1)-.001; % subtract
    % 1 ms from 1st pulse time to insure onset precedes 1st pulse
    LV = tsd(:,1)>CSonTms(i) & tsd(:,1)<CSonTms(i)+1; % compassing current CS
    d = tsd(LV,:); 
    r = find(d(:,2)==30,1,'last'); % last pulse of current CS  
    CSoffTms(i) = d(r,1)+.001; % add 1 ms to insure that this event occurs
    % after the final pulse
    USonTms(i) = CSonTms(i)+ISI/1000;
    i=i+1;
end
    
CSonD = [CSonTms 50*ones(size(CSonTms))]; % CS on data; 50 is the event
% code for CS on
CSoffD = [CSoffTms 60*ones(size(CSoffTms))]; % CS off data; 60 is the event
% code for CS off
USonD = [USonTms 70*ones(size(USonTms))]; % US on data; 70 is the event code

tsdA = sortrows([tsd;TrlOnD;TrlOffD;CSonD;CSoffD;USonD]);
%}
% The USoff event was added later--see cell above that follows the
% TSimporteventcodes cell
%% Creating trials

TSdefinetrialtype('PreTrl',[TrlOn CSon])
TStrialstat('SpkTms',@TSparse,'result=time(1);',spike)
TSapplystat('NumSpks','SpkTms',@numel)
TSapplystat('SpkRate',{'NumSpks' 'TrialDuration'},@rdivide)
TSapplystat('ISpkIstats','SpkTms',@SpkStats) % mean, std & cv
%{
function O = SpkStats(D)
ispi = diff(D);
O = [mean(ispi) std(ispi) std(ispi)/mean(ispi)];
%}
TSapplystat('SpikeCountInLast1s',{'EndTime' 'SpkTms'},@SpkCnt)
%{
function Cnt = SpkCnt(T,E,W)
Cnt = sum(T>(E-W)); 
%}
TScombineover('PreSpkRates','SpkRate')
TScombineover('Pre1sSpkCnts','SpikeCountInLast1s')
TSapplystat('FanoFactor','Pre1sSpkCnts',@Fano)
%{
function F = Fano(Cnts)
F = var(Cnts)/mean(Cnts);
%}
TScombineover('FanoFactor_S','FanoFactor')
TScombineover('BtwTrialFanoFactors','FanoFactor_S')
%% Between-trial Fano factor graphic
figure
cdfplot(Experiment.BtwTrialFanoFactors)
hold on
plot([1.15 1.15],ylim,'k--')
title('')
ylabel('Cumulative Fraction of Cells/Subjects')
xlabel('Between Trial Fano Factor (\sigma^2/\mu)')

%% Within-trial Fano factors
TSapplystat('FanoFac','SpkTms',@FanoF)
%{
function O = FanoF(T)
%%
edges = T(1):.1:T(end);
N = histc(T,edges);
lam = numel(T)/(T(end)-T(1)); % estimated firing rate
n = numel(edges); % number of samples
FF = var(N)/mean(N); % estimated Fano Factor
if FF<1
    p = gamcdf(FF,(n-1)/2,2/(n-1));
else
    p = 1-gamcdf(FF,(n-1)/2,2/(n-1));
end    
O = [FF lam n p]; % [Fano Factor lambda n (p if 1)]
%}
TScombineover('FanoFacs','FanoFac')
%% Cumulative Distributions
TSapplystat('','FanoFacs',@TSplotcdfs,'Rows',5,'DataCols',{1},'Xlm',[0 2],...
    'Xlbl','Fano Factor w/i Trial')
for S=1:10
    subplot(5,2,S)
    hold on
    plot([.56 .56],ylim,'k--',[1.57 1.57],ylim,'k--')
end
  
%%
TSdefinetrialtype('CSUS',[CSon USon])
TStrialstat('SpkTms',@TSparse,'result=time(1)-starttime;',spike)
TSapplystat('NumSpks','SpkTms',@numel)
TSapplystat('SpkRate',{'NumSpks' 'TrialDuration'},@rdivide)
TScombineover('CSUSrates','SpkRate')

%% Spike rate diffs
% In 3 subjects, there is 1 more CS "trial" than Pre[CS] "trial".
% Following code finds out where the slippage occurs
S=4;
t=1;
while true
    if Experiment.Subject(S).Session(1).TrialPreTrl.Trial(t).EndTime ~=...
        Experiment.Subject(S).Session(1).TrialTrl.Trial(t).StartTime
        disp(t)
        break
    elseif t==length(Experiment.Subject(S).Session(1).TrialPreTrl.Trial(t).EndTime)
        disp('End of data')
        break
    end
    t = t+1;
end
% For S2, there is a bogus CS presentation starting at Row 139947 and
% ending at Row 140011, with no spikes in the entire trial, only CS pulses,
% and those continue after the US, which they should not do in this
% subject. For S3, there is a trial beginning at Row 243786 that likewise
% begins with a CS on (no pre); it ends at Row 243812. For S4, things go
% screwy at Row 180037, which is shortly before the end and continue screwy
% to the end. All the other cells ("subjects") have equal numbers of PreTrl
% and Trl trials.
%%
%%
Experiment.Subject(2).Session(1).TSDataA(139947:140011,:)=[]; % deleting
% bad data from S2
Experiment.Subject(3).Session(1).TSDataA(243786:243812,:)=[]; % deleting bad
% data from S3
Experiment.Subject(4).Session(1).TSDataA(180037:end,:)=[]; % deleting bad
% data from S4
%% Deleting bad trials from S4
Experiment.Subject(4).Session(1).TrialTrl.Trial(611:end)=[];
Experiment.Subject(4).Session(1).TrialTrl.NumTrials = 610;
Experiment.Subject(4).Session(1).TrialPreTrl.Trial(611:end)=[];
Experiment.Subject(4).Session(1).TrialPreTrl.NumTrials = 610;
Experiment.Subject(4).Session(1).TrialPostTrl.Trial(611:end)=[];
Experiment.Subject(4).Session(1).TrialPostTrl.NumTrials = 610;
%% Deleting bad data from session stat fields for S4
Experiment.Subject(4).Session(1).PreSpkRates(611:end)=[];
Experiment.Subject(4).Session(1).PostSpkRates(611:end)=[];

%%
TSsaveexperiment % after cleaning out bad data and bad results therefrom
%% Spike Rate Diffs
TSapplystat('CSUS_Pre_RateDiffs',{'PreSpkRates' 'CSUSrates' },@minus)
%%
TSsaveexperiment
%% Cumulative Records of rate differences
TSapplystat('','CSUS_Pre_RateDiffs',@TSplotcumrecs,'Rows',5,...
    'Xlbl','Trial','LeftYlbl','Cm \lambda_p_r_e-\lambda_c_s_u_s')
% Adding acquisition trial
for S=1:10
    subplot(5,2,S)
    hold on
    switch S
        case 1
            ylim([-1 2000])
            plot([1 1],ylim,'k--')
        case 2
            ylim([-1 8000]);xlim([0 700])
            plot([406 406],ylim,'k--')
        case 3
            ylim([-3000 15000])
            plot([286 286],ylim,'k--')
        case 4
            ylim([-3000 2000]);xlim([0 700])
            plot([220 220],ylim,'k--')
        case 5
            ylim([-1 7000]);xlim([0 450])
            plot([1 1],ylim,'k--')
        case 6
            xlim([ 0 600]);ylim([-5000 5000])
            plot([345 345],ylim,'k--')
        case 7
            xlim([ 0 650]);ylim([-3000 10000])
            plot([231 231],ylim,'k--')
        case 8
            xlim([ 0 450]);ylim([-1 20000])
            plot([1 1],ylim,'k--')
        case 9
            xlim([ 0 300]);ylim([-1 10000])
            plot([1 1],ylim,'k--')
        case 10
            xlim([ 0 550]);ylim([-2500 500])
            plot([252 252],ylim,'k--')
    end
end
AcqTrlVecs = [1 1;406 406;286 286;220 220;1 1;345 345;231 231;1 1;1 1;252 252];
%% USCSoff trials for cells 8 to 11
TSlimit('Subjects',8:11)
TSdefinetrialtype('USCSoff',[USoff CSoff])
TStrialstat('SpkTms',@TSparse,'result=time(1)-starttime;',spike)
TSapplystat('NumSpks','SpkTms',@numel)
TSapplystat('SpkRate',{'NumSpks' 'TrialDuration'},@rdivide)
TScombineover('USCSoffSpikeRates','SpkRate')
TSapplystat('USCSoff_CSUSrateDiffs',{'USCSoffSpikeRates' 'CSUSrates'},@minus)
%%
TSsaveexperiment
%% The above code in this cell failed to generate a the rate differences for S7
% because the last USCSoff "trial" doesn't exist, so I did it "by hand":
%
Experiment.Subject(7).Session(1).USCSoff_CSUSrateDiffs =...
    Experiment.Subject(7).Session(1).USCSoffSpikeRates -...
    Experiment.Subject(7).Session(1).CSUSrates(1:end-1);
%% Cm Recs rate diffs btw last part (post US) and 1st part of CS
TSapplystat('','USCSoff_CSUSrateDiffs',@TSplotcumrecs,'Xlbl','Trial','LeftYlbl','Cm USCSoffCSUSrateDiff')

%% Raster plots
pre = .2; % duration of pre CS onset interval (in sec)
post = 1; % duration of post US interval
TSlimit('Subjects',1)
TSapplystat('','TSDataA',@Raster,pre,post)
%{
function Raster(tsd,pre,post)
LVcs = tsd(:,2)==50; % flags CS onsets
LVus = tsd(:,2)==70; % flags US onsets
Db = [tsd(LVcs,1)-pre 45*ones(size(tsd(LVcs,1)))]; % to be inserted onset data
De = [tsd(LVus,1)+post 75*ones(size(tsd(LVus,1)))];% to be inserted offset data
tsd = sortrows([tsd;Db;De]);
TSraster(tsd,[45 75],[50 0;70 0;60 0;40 0],['g+';'b+';'rv';'k.'])
xlim([0 1])
title(['Cell' num2str(evalin('caller','sub'))])
%}
%% Raster plots for pub
pre = .3; % duration of pre CS onset interval (in sec)
post = .6; % duration of post US interval
figure

ax1 = subplot(1,2,1);
TSlimit('Subjects',6)
TSapplystat('','TSDataA',@Raster2,pre,post,ax1)

ax2 = subplot(1,2,2);
TSlimit('Subjects',7)
TSapplystat('','TSDataA',@Raster2,pre,post,ax2)


%% Change points
TSapplystat('CSUS_PRE_RateDiffsCPs','CSUS_Pre_RateDiffs',@cp_wrapper,1,2,6)
%%
TSlimit('Subjects',8:11)
TSapplystat('CSUS_Pre_RateDiffCPs','CSUS_PreRateDiffs',@cp_wrapper,1,2,6) Setting 
%% Setting trial number in first row of CSUS_Pre_RateDiffCPs to 2 because
% it is analytic that acquisition cannot occur until after the first trial
for S=1:11
    Experiment.Subject(S).Session(1).CSUS_PRE_RateDiffsCPs(1,1)=2;
end

%% Reversing sign of CP data, so that positive slope indicates acquisition
for S=1:10
    Experiment.Subject(S).Session(1).CSUS_PRE_RateDiffsCPs(2:end,2) = ...
        -Experiment.Subject(S).Session(1).CSUS_PRE_RateDiffsCPs(2:end,2);
end
%% Post US firing
%{
Extracted from Dan-Anders email on July 28, 2017
"Here are the US pulse times in [the] original data files."
436: 0, 2, 4, 6, 8 and 20, 22, 24, 26, 28
457: 0, 2, 4, 6, 8 and 20, 22, 24, 26, 28
474: 0, 2, 4, 6, 8 and 20, 22, 24, 26, 28  
479: 0, 2, 4, 6, 8 and 20, 22, 24, 26, 28
483: 0, 2, 4, 6, 8 and 20, 22, 24, 26, 28
485: 0, 2, 4, 6, 8 and 20, 22, 24, 26, 28
487: 0, 2, 4, 6, 8 and 20, 22, 24, 26, 28
505: 0, 2, 4, 6, 8 and 12, 14, 16, 18, 20
524: 0, 2, 4, 6, 8 and 12, 14, 16, 18, 20
536: 0, 2, 4, 6, 8 and 12, 14, 16, 18, 20
543: 0, 2, 4, 6, 8 and 12, 14, 16, 18, 20
%}
TSdefinetrialtype('US',[USon USoff])
%%
TSlimit('Subjects','all')
TSlimit('Trials','all')
TStrialstat('PostB1spkTms',@InterBurstSpikes,spike)
% latency to 1st spike measured from last pulse in 1st US burst and any
% subsequent inter-spike intervals
%{
function O = InterBurstSpikes(tsd,spikecode)
starttime = tsd(1,1);
LVspk = tsd(:,2)==spikecode; % flags spikes
if sum(LVspk)<1 % no spikes
    O=[];
else
    O = tsd(LVspk,1) - starttime - .008;
    % last pulse in 1st burst of stimulation occurs 8 ms after US onset
    O = O(O>0); % deleting negative spike latencies, which come from spikes
        % recorded during stimulation burst
end
%}
%%
TSapplystat('IspkIs','PostB1spkTms',@IspkInts)
%{
function O = IspkInts(tsd)
if isempty(tsd)
    O=[];
else
    O = [tsd(1,1);diff(tsd(:,1))];
    % latency to 1st spike measured from last pulse in 1st US burst and any
    % subsequent inter-spike intervals
end
%}
%%
maxnum = 5; % maximum # of interspike intervals
TSapplystat('IspIbyPos','IspkIs',@rearrange,maxnum)
% rearranges col vector of interspike intervals into an row vector with a
% fixed length given by the argument maxnum. 1st col are the latencies to
% 1st spike; subsequent columns are the successive interspike intervals
%{
function O = rearrange(ints,cols)
O = nan(1,cols);
if isempty(ints)
    return
elseif length(ints)>cols
    c = cols;
else
    c = length(ints);
end
O(1,1:c) = ints(1:c);
%}
%%
TScombineover('USIspIarray','IspIbyPos') % creates an array at the session
% level in which each row is a trial, the 1st col is the latency to the 1st
% spike (following the 1st stim burst) and subsequent columns give the
% subsequent interspike intervals, with nans in the blank positions. After
% studying the histograms of these 1st latencies and interspike intervals,
% for the first two cells, I concluded that they were mostly determined by
% the stimulation parameters, although that conclusion warrants further
% thought and study. It does not easily explain the 1st spike latencies,
% because these were AFTER the last pulse of the 1st burst, nor, for
% similar reasons does it explain the 1st interspike interval. It does
% explain subsequence interspike intervals and there the spacing at the
% interpulse interval (2ms) is readily apparent. Anyway, on to the post US

%%
TSdefinetrialtype('USoffToTrlOff',[USoff TrlOff])
TSlimit('Subjects',1:7) % the subjects with a 300 ms CS terminating just
% before US onset
TStrialstat('SpkTms',@TSparse,'result=time(1)-starttime;',spike)
% spike times measured from US offset
%%
TSapplystat({'FrstSpkLat' 'SpkTmsAndIspkIs'},'SpkTms',@SpkTmsIspI)
% 3-col array with spike # in the 1st column, spike time in 2nd col and the 
% ensuing inter-spike interval in the third column 
%{
function [O1,O2] = SpkTmsIspI(spktms)
O1 = spktms(1); % latency of 1st spike following US offset
O2 = [(1:size(spktms,1)-1)' spktms(1:end-1) diff(spktms)];
%}
%%
TScombineover('FrstSpkeLat_s','FrstSpkLat')
%%
QV = [0.0500 0.1000 0.2500 0.5000 0.7500 0.9000 0.9500];
TSapplystat('FrstSpkLatQuantiles','FrstSpkeLat_s',@quantile,QV)
%%
% figure
hold on
S=7;
plot(S*ones(7,1),Experiment.Subject(S).Session(1).FrstSpkLatQuantiles,'k*')
%%
TScombineover('FrstSpkLatQntlsCells_S','FrstSpkLatQuantiles')
TScombineover('FrstSpkLatQntlsCells_Ss1_7','FrstSpkLatQntlsCells_S')
%%
edges = 0:.025:2.5;
TSapplystat('Bin25SpkRates','SpkTms',@histc,edges)
%%
edges = 0:.05:3.5;
TSapplystat('Bin50SpkRates','SpkTms',@histc,edges)
edges = 0:.1:3.5;
TSapplystat('Bin100SpkRates','SpkTms',@histc,edges)
edges = 0:.2:3.4;
TSapplystat('Bin200SpkRates','SpkTms',@histc,edges)
%%
% Trial 370 of S2 had only a single interspike interval. For some reason,
% this caused all of the above fields to be row vectors rather than column
% vectors for this one trial. Following code fixes this
TSlimit('Subjects',2);TSlimit('Trials',370)
TSapplystat('Bin25SpkRates','Bin25SpkRates',@transpose)
TSapplystat('Bin50SpkRates','Bin50SpkRates',@transpose)
TSapplystat('Bin100SpkRates','Bin100SpkRates',@transpose)
TSapplystat('Bin200SpkRates','Bin200SpkRates',@transpose)
TSsaveexperiment
%%
TSlimit('Subjects',1:7)
TSlimit('Trials','all')
TScombineover('Bin25SpkRates_s','Bin25SpkRates')

Num25msbins = 141;
TSapplystat('Bin25SpkRates_s','Bin25SpkRates_s',@reshape,[],Num25msbins)
% reshapes the column vector so that each row is a trial and the columns
% are the bins. Thus, there are as many columns (Num25msbins) as there are
% spike counting bins. In general, it is a bad idea to make the output
% field the same as the input field, but in this case it seems harmless
%
TSapplystat('MeansOf25msBins','Bin25SpkRates_s',@mean)
TSapplystat('StdsOf25msBins','Bin25SpkRates_s',@std)
%%
TSapplystat('','MeansOf25msBins',@TSplot,'Xlbl','Successive 25 ms bins',...
    'Ylbl','Mean Spike Count','Xlm',[0 20])
%%
TSlimit('Subjects',1)
TSapplystat('',{'Bin25SpkRates_s'},@MnsWith2stdErr,[0 100])
%{
function MnsWith2stdErr(A,xlm)
mu = mean(A); % bin means
stder = std(A)/sqrt(size(A,1)); % standard errors
b = 1:size(A,2);
figure
errorbar(b,mu,2*stder)
xlabel('Successive 25 ms Bins')
ylabel('Mean # Spikes +/-2 std er')
xlim(xlm)
%}
% This plot suggests a rich structure, as there are several large changes
% in the mean spike count over the span of just 1 or 2 25ms bins. Question
% is how to bring out the structure. Two ideas: 1) Look at the interspike
% interval distributions as a function of the # of spikes post US offset.
% 2) At short successive intervals, look at the distribution of interspike
% intervals in progress. Imagine a line extending up through the raster
% plot at a given post-USoff latency and determining the distribution of
% inter-spike intervals that straddle that line. Then, move the line and
% redo.
%%
TSlimit('Subjects',1:7)
TScombineover('SpkNumTimesIspkI','SpkTmsAndIspkIs','t')
% 4-col array at Session level containing data from the USoffToTrlOff
% trials. 1st col is spike # within a trial; 2nd col is spike time; 3rd col
% is interspike interval; 4th col is Trial number

%%
TSlimit('Subjects',1:7)
HighestSpikeNum = 100;
TSapplystat('USoffTrlOffIspIsBySpkNum','SpkNumTimesIspkI',...
    @IspIsortedBySpkNum,HighestSpikeNum)
%{
function O = IspIsortedBySpkNum(A,HSN)
% O is an array with HSN columns and as many rows as there were trials
% Each column, c, gives the interspike interval commenced by spike#c in
% Trial r
%
% A is a 4-col array: 1st col is spike # within a trial; 2nd col is spike 
% time; 3rd col is interspike interval; 4th col is Trial number

O = nan(A(end),HSN); % Initializing output array; there are as many rows as
% trials and HSN columns, where HSN is the highest spike count to be
% considered
for sp=1:HSN % stepping through successive spikes (columns of O)
    for t = 1:A(end) % stepping through the trials (rows of O)
        LVt = A(:,4)==t; % flags rows in A belonging to Trial t
        nsp = sum(LVt); % # of spikes in Trial t
        if nsp < sp % spike count in Trial t < nsp
            continue % go on to next trial
        else % at least nsp spikes in Trial t
            LVsp = A(:,1)==sp; % flags rows with spike count sp
            try
                O(t,sp) = A(LVt&LVsp,3); % the interspike interval commenced by
                % spike # sp in Trial t
            catch ME
                keyboard
            end
        end
    end
end    
%}
%%

%%
TrlRange = 1:100; % range of trials over which quantiles to be computed
QV = [0.0500 0.1000 0.2500 0.5000 0.7500 0.9000 0.9500]; % quantile p's
TSapplystat('USIspkIQntlsBySpkCntTrls1_100','USoffTrlOffIspIsBySpkNum',...
    @IspkIQs,TrlRange,QV)
%{
function O = IspkIQs(IspkIbySpkCnt,TrlRange,QV)
O = quantile(IspkIbySpkCnt(TrlRange,:),QV);   
%}
%% Multi-panel plot of log quantiles; 1st 100 trials
% tr=TrlRange;
S = 7; % Each setting from 1:7 produces a panel
D = Experiment.Subject(S).Session(1).USIspkIQntlsBySpkCntTrls1_100;
subplot(4,2,S)
semilogy(tr,D(1,:),tr,D(2,:),tr,D(3,:),tr,D(4,:),tr,D(5,:),tr,D(6,:),tr,D(7,:))
set(gca,'FontSize',12,'YTick',[.002 .005 .01 .02 .05 .1],'YTickLabel',...
    {'2' '5' '10' '20' '50' '100'})
switch S
    case 1
        ylim([.004 .12])
    case 2
        ylim([.002 .06])
    case 3
        ylim([.0015 .12])
    case 4
        ylim([.002 .12])
    case 5
        ylim([.001 .12])
    case 6
        ylim([.001 .1])
    case 7
        ylim([.002 .05])
end
if S >5
    xlabel('Spike Count')
end

if mod(S,2)>0
    ylabel('IspkI Qntl (ms), log scale')
end
if S>6
    leg = legend('.05','.1','.25','.5','.75','.9','.95');
    set(leg,'FontSize',14)
end
title(['Cell ' num2str(S) 'USoffToTrlOff, Trls 1:100'],'FontSize',12)

%% CDFs of 1st spk latency and 1st 5 IspkIs; 1st 100 trials
% figure
tr=TrlRange;
S = 7; % Each setting from 1:7 produces a panel
D = Experiment.Subject(S).Session(1).USoffTrlOffIspIsBySpkNum(tr,1:5);
subplot(4,2,S)
cdfplot(Experiment.Subject(S).Session(1).FrstSpkeLat_s)
hold on
cdfplot(D(tr,1))
cdfplot(D(tr,2))
cdfplot(D(tr,3))
cdfplot(D(tr,4))
cdfplot(D(tr,5))
set(gca,'FontSize',12)
xlabel('Latency(s)')
title(['Cell ' num2str(S) '1st 100 trials'],'FontSize',12)
% set(leg,'FontSize',12)
if S<1
    xlim([0 .05])
elseif S==2
    xlim([0 .01])
else
    xlim([0 .02])
end
if S>6
    leg=legend('0-1','1-2','2-3','3-4','4-5','5-6','location','SE');
    set(leg,'FontSize',12)
end

%% Some graphics (Cell 1 only for time being)
figure
cdfplot(A(A(:,1)==1,2)) % cumulative distribution of 1st spike latencies
xlim([0 .025])
title('Cell 1: 1st spike latency','FontSize',18)
xlabel('Latency from US offset (ms)','FontSize',14)

%%
figure
cdfplot(O(:,1))
xlim([0 .025])
xlabel('Interspike Interval (ms)')
hold on
cdfplot(O(:,2))
cdfplot(O(:,3))
cdfplot(O(:,4))
legend('Spk1-2','Spk2-3','Spk3-4','Spk4-5','location','SE')
title('Cell 1','FontSize',18)
%%
figure
cdfplot(O(:,4))
xlim([0 .025])
xlabel('Interspike Interval (ms)')
hold on
cdfplot(O(:,5))
cdfplot(O(:,6))
cdfplot(O(:,7))
legend('Spk4-5','Spk5-6','Spk6-7','Spk7-8','location','SE')
title('Cell 1','FontSize',18)
xlim([0 .1])
%%
figure
cdfplot(O(:,7))
xlim([0 .1])
xlabel('Interspike Interval (ms)')
hold on
cdfplot(O(:,8))
cdfplot(O(:,9))
cdfplot(O(:,10))
cdfplot(O(:,11))
legend('Spk7-8','Spk8-9','Spk9-10','Spk10-11','Spk11-12','location','SE')
title('Cell 1','FontSize',18)

%%
Q = [.05 .1 .25 .5 .75 .9 .95]; % quantiles
QA = quantile(O,Q);
%%
figure
plot(QA(1,:))
hold on
plot(QA(2,:))
plot(QA(3,:))
plot(QA(4,:))
plot(QA(5,:))
plot(QA(6,:))
plot(QA(7,:))
xlabel('Spike Count','FontSize',14)
ylabel('Interspike Interval Quantile (ms)','FontSize',14)
leg=legend('.05','.1','.25','.5','.75','.9','.95','location','N');
set(leg,'FontSize',18)
title('Cell 1','FontSize',18)

%%
figure
semilogy(QA(1,:))
hold on
semilogy(QA(2,:))
semilogy(QA(3,:))
semilogy(QA(4,:))
semilogy(QA(5,:))
semilogy(QA(6,:))
semilogy(QA(7,:))
xlabel('Spike Count','FontSize',14)
ylabel('Interspike Interval Quantile (ms)','FontSize',14)
leg = legend('.05','.1','.25','.5','.75','.9','.95');
title('Cell 1','FontSize',18)
set(leg,'FontSize',18)
% CELL 2 NOTES: Eventually--after 40 or 50 spikes the interspike intervals are log
% normally distributed, with a median of 17 ms, but initially, they come
% much faster (median <5 ms) and they're not log normally distributed; the
% distribution of log(IspkI) is skewed left, with a relatively long tail,
% albeit at the fastest point (3rd spike), 95% of the intervals are less
% than the eventual median
ylim([.002 .12])
set(gca,'YTick',[.002 .005 .01 .02 .05 .1],'YTickLabel',{'.002' '.005' '.01' '.02' '.05' '.1'})

%% Eventual distribution
LV = Experiment.Subject(1).Session(1).SpkNumTimesIspkI(:,1)>50;
% flags spike counts >50, where the eventual distribution obtains
[MUHAT,SIGMAHAT,MUCI,SIGMACI] =...
    normfit(log10(Experiment.Subject(1).Session(1).SpkNumTimesIspkI(LV,3)))
%{
MUHAT = -1.5194    
MUCI =-1.5229
      -1.5158

SIGMAHAT =0.2306
SIGMACI =0.2281
         0.2331   
%}
%
figure
cdfplot(log10(Experiment.Subject(1).Session(1).SpkNumTimesIspkI(LV,3)))
set(gca,'FontSize',12)
hold on
xlabel('log10(IspkIs)')
title('Cell 1 USoff to Trial end, Spike Counts > 50')
xvals = -2.5:.1:0; % for plotting
plot(xvals,normcdf(xvals,MUHAT,SIGMAHAT))
xlim([-2.5 0])
legend('IspIs','Fit(\mu=-1.5194,\sigma=.2306)','location','SE')
%%
figure
edges=.001:.001:.1;
N=histc(Experiment.Subject(1).Session(1).SpkNumTimesIspkI(LV,3),edges);
bar(edges,N,'histc')
set(gca,'FontSize',12)
xlabel('InterSpike Interval (s)')
title('Cell 1: Data from Spike Counts>50 after US offset')
%%
LV3 =Experiment.Subject(1).Session(1).SpkNumTimesIspkI(:,1)==3;
figure
edges=.001:.001:.1;
N=histc(Experiment.Subject(1).Session(1).SpkNumTimesIspkI(LV3,3),edges);
bar(edges,N,'histc')
set(gca,'FontSize',12)
xlabel('InterSpike Interval (s)')
title('Cell 1: Data from Spike Count 3 after US offset')

%% Writing some IspkI files for sending to James
S=7;
str=sprintf('S%dforJames',S);
fid = fopen(str,'w');
D = [Experiment.Subject(S).Session(1).FrstSpkeLat_s(1:100) ...
    Experiment.Subject(S).Session(1).USoffTrlOffIspIsBySpkNum(1:100,1:5)];
fprintf(fid,'%0.4f %0.4f %0.4f %0.4f %0.4f %0.4f\n',D');
fclose(fid);
type(str)

%% Hand-built distribution models
D=S1forJames(:,1);
x=0:.001:.15;
cdfplot(S1forJames(:,1))
xlim([0 .02])
hold on
e=expfit(D);
plot(x,expcdf(x,e))
%
d1=D(D<.0014);
d2=D(D>.0014&D<.0039);
d3=D(D>.0039&D<.01);
d4=D(D>.01&D<.015);
d5=D(D>.015&D<.1);
%
[m1,s1] = normfit(d1);
[m2,s2] = normfit(d2);
[m3,s3] = normfit(d3);
[m4,s4] = normfit(d4);

[m5,s5] = normfit(d5);
%%
plot(x,(length(d1)/100)*normcdf(x,m1,s1)+(length(d2)/100)*normcdf(x,m2,s2)+...
    (length(d3)/100)*normcdf(x,m3,s3)+(length(d4)/100)*normcdf(x,m4,s4)+...
    (length(d5)/100)*normcdf(x,m5,s5));
xlim([0 .075])
set(gca,'FontSize',18)
xlabel('InterSpike Interval (s)')
ylabel('Cumulative Probability')
title('S1 postUS 1st Spike Latency')
legend('Data','Exp. Model','Mixture Model','location','SE')
%%
figure
hist(log10(D),20)

%%
figure
cdfplot(log10(D))
%%
figure
C=-3.2:.1:-.5;
N=hist(log10(D),C);
Do = sort(log10(D));
cdf = .01:.01:1;
[Az,H1,H2]=plotyy(C,N/sum(N),Do,cdf,@stairs);
set(H1,'LineWidth',2)
set(H2,'LineWidth',2)
Xvls =[-3 -2.5 -2 -1.5 -1 -.5];
set(Az(1),'XTick',Xvls,'XTickLabel',{'1' '3' '10' '32' '100' '320'})
xlabel('InterSpike Interval (ms, log scale)')
set(gca,'FontSize',18)
set(Az(2),'FontSize',18)
ylabel('Probability')
ylabel(Az(2),'Cumulative Probability')
title('S1postUS_1stSpkLatDist')

%%
TSapplystat('TSDataB',{'TSData' 'Phase'},@AugTSdata) % creates TSDataA
% (A is for augmented) field at the Session level in which there are CSon,
% CSoff and MidITI events. This enables me to define the appropriate kind
% of trials
%{
function tsdA = AugTSdata(tsd,CSdur)
CSdur = CSdur/1000;
Mtms=tsd(tsd(:,2)==20,1); % Mtm is 0.2 s prior to triggering CS
CSonTms = Mtms+.2;
CSoffTms = CSonTms+CSdur;
ITIs = CSoffTms(2:end)-CSonTms(1:end-1);
HalfITIs=round(1000*ITIs/2)/1000;
LstPostInt=tsd(end,1)-CSoffTms(end);
MidITItms = CSoffTms+[HalfITIs;LstPostInt];
MidITItms=[tsd(1)-.001;MidITItms];
MidITIs = [MidITItms 60*ones(length(MidITItms),1)]; % 60 is code for mid
% ITI
CSons = [CSonTms 30*ones(length(CSonTms),1)]; % 30 is code for CS on
CSoffs=[CSoffTms 50*ones(length(CSoffTms),1)]; % 50 is code for CS off
tsdA = sortrows([tsd;CSons;CSoffs;MidITIs]);
%}

%% Referencing trial times to 0 at CS onset
TSlimit('Subjects','all')
TSapplystat('TSDataB','TSDataA',@BinDanAndersData)
% Creates field at session level containing tsd with times referenced to 0's
% at CS onsets
%{
function BV = BinDanAndersData(tsd)
rTon = find(tsd(:,2)==80); % row #s for trial on
rToff = find(tsd(:,2)==90); % row #s for trial off
rcs = find(tsd(:,2)==50); % row #s for CS on
BV = tsd;
for r = 1:length(rToff)
    BV(rTon(r):rToff(r),1) = BV(rTon(r):rToff(r),1)-BV(rcs(r),1);
end
%}
%% Peri-ISI histograms
figure
bw=.015; % bin width
TSapplystat('PeriCShistogram',{'TSDataB' 'CSUS_PRE_RateDiffsCPs'},@PeriCShist,bw)
% Creates a field at the Session level containing the bins and the
% p_s, the momentary probability of a spike in a 1 ms interval within
% each bin; the bins themselves are bw wide. Also produces the (normalized)
% histograms in a multi-plot figure
%{
function Hst = PeriCShist(D,CPs,bw)
%%
S = evalin('caller','sub');
%
switch S
    case 6
        nb = 5; % number of steps back in the CPs
    case 7
        nb = 3;
    case {(1) (3) (4) (8) (9) (10)}
        nb = 1;
    otherwise
        nb = 2;
end
CPs(1,1)=1;
ron = find(D(:,2)==80); % row #s at trial starts
 % row # where clear pause has emerged
D(1:ron(CPs(end-nb,1)),:)=[]; % deleting data prior to clear pause trials
nt = CPs(end,1) - CPs(end-nb,1); % number of trials from which counts computed
LVs = D(:,2)==40; % flags spikes
EndCSUS = D(find(D(:,2)== 70,1),1)-.005; % trial time when CS terminates
LVdh = D(:,1)>-.3 & D(:,1)<EndCSUS; % flags all the stretches between -0.3s
% and the end of the CS-US interval
edges = -.3:bw:EndCSUS;
N = histc(D(LVdh&LVs,1),edges);
Hst = [edges' N/(nt*bw/.001)]; % p_s in a 1ms bin is the count divided by
% the product of the number of trials over which histogram was computed
% (nt) and the # of 1 ms bins in a histogram bin (bw/.001) 
subplot(5,2,S)
% bar(edges,N,'histc')
bar(edges,N/(nt*bw/.001),'histc')
xlim([-.3 .3])
if S>8
    xlabel('Trial Time (s)')
end
if mod(S,2)>0
    ylabel('p_s')
end
%}
%% Pause onset analysis using binarize spike vectors
TSlimit('Subjects',1:10)
TSapplystat('PutativePsOnsets',{'TSDataB' 'PeriCShistogram'},@pson)
% Creates a 3-col field at the Session level: col 1 is the onset time of a
% putative pause; col 2 is the firing rate prior to that onset minus the
% firing rate after that onset (so positive values indicate a dip in
% firing); col 3 is the weight of the evidence in favor of there being a
% change at that point (NB not necessarily a downward change!)
%{
function OT = pson(tsd,pdf)
bins = -.3:.001:pdf(end-1,1); % bins for binary vector
LVs = tsd(:,2)==40; % flags spikes
ron = find(tsd(:,2)==80); % rows where trials start
roff = find(tsd(:,2)==90); % rows where trials end
OT = nan(length(ron),3);
pdf(end,:)=[]; % deleting last entry, which is always 0
for r = 1:length(ron)
    %%
    Dd=tsd(ron(r):roff(r),1); % times btw trial start & trial end
    LVd=LVs(ron(r):roff(r)); % flags spike times in this same stretch
    Ntst=histc(Dd(LVd),bins); % counts spike times into successive
    %% 1ms wide bins btw -.3 and end of CS
    bv = [bins' Ntst]; % binarized spike vector
    bv(bv(:,2)>1,2)=1; % eliminating the (rare) integers>1
    LVcs = pdf(:,1)>0; % flags portions of pdf within CS
    Mn = min(pdf(LVcs,2)); % smallest p_s within CS
    alphaB = 1; % setting hyperparameters
    betaB = ceil((1-Mn)/Mn); % setting hyperparameters
    alphaA = 1; % setting hyperparameters
    AvCSp = mean(pdf(~LVcs,2)); % setting hyperparameters
    betaA = ceil((1-AvCSp)/AvCSp);% setting hyperparameters
    [CP,Odds] = BernCP(flipud(bv(:,2)),alphaB,betaB,alphaA,betaA,1/length(bv));
    W = log10(Odds);
    cp = bv(length(bins)-CP,1);
    LVps = bv(:,1)>cp;
    lamPs = sum(bv(LVps,2))/(pdf(end,1)-cp);
    lamPre = sum(bv(~LVps,2))/(cp-pdf(1,1));
    OT(r,:) = [cp lamPre-lamPs W];    
end    
%}
%% Graphing pause onset data
AcqTrlVecs = [2 2;406 406;286 286;220 220;2 2;345 345;231 231;2 2;2 2;252 252];
for S=1:10
    OT = Experiment.Subject(S).Session(1).PutativePsOnsets; % pause stats
    tv = 1:length(OT); % trial vector
    Xlm = [0 length(OT)+10];
    figure
    subplot(3,1,1)
    plot(1:length(OT),OT(:,2));
    xlim(Xlm)
    hold on
    plot(xlim,[0 0],'r')
    plot(AcqTrlVecs(S,:),ylim,'k--')
    ylabel('\lambda_p_r_e - \lambda_p_s')
    title(['Cell' num2str(S)])
    
    subplot(3,1,2)
    plot(1:length(OT),OT(:,3));
    xlim(Xlm)
    hold on
    plot(xlim,[0 0],'r')
    plot(AcqTrlVecs(S,:),ylim,'k--')
    ylabel({'Weight of Evidence';'for a Change'})
    title(['Cell' num2str(S)])
    %%
    subplot(3,1,3)
    if S == 1
        crit = .3;
    else
        crit = 1;
    end
    LV = (OT(:,2)>0) & (OT(:,3)>crit); 
    % flags trials with evidence of a pause
    plot(tv(LV),OT(LV,1),'o')
    xlim(Xlm)
    hold on
    plot(xlim,[0 0],'r')
    plot(AcqTrlVecs(S,:),ylim,'k--')
    xlabel('Trial')
    ylabel('Pause Onset')
end
% I adjusted x and y limits by hand after these plots were made and added
% text giving the weight crit
%%
for S=1:10
    figure(S)
    saveas(gcf,['S' num2str(S) 'PauseOnsetStats.fig'])
    saveas(gcf,['S' num2str(S) 'PauseOnsetStats.pdf'])
end
%% Graphing pause onset data
AcqTrlVecs = [2 2;406 406;286 286;220 220;2 2;345 345;231 231;2 2;2 2;252 252];
figure
for S=5:6
    
    OT = Experiment.Subject(S).Session(1).PutativePsOnsets; % pause stats
    tv = 1:length(OT); % trial vector
    Xlm = [0 length(OT)+10];
    if S==5
        subplot(3,2,1)
        plot(1:length(OT),OT(:,2));
        xlim(Xlm)
        ylim([-100 150])
        hold on
        plot(xlim,[0 0],'k--')
        plot(AcqTrlVecs(S,:),ylim,'k--')
        ylabel('\lambda_p_r_e - \lambda_p_s')
        title(['Cell' num2str(S)])

        subplot(3,2,3)
        plot(1:length(OT),OT(:,3));
        xlim(Xlm)
        hold on
        plot(xlim,[0 0],'k--')
        ylim([-.7 5])
        plot(AcqTrlVecs(S,:),ylim,'k--')
        ylabel({'Weight of Evidence';'for a Change'})

        subplot(3,2,5)
        if S == 1
            crit = .3;
        else
            crit = 1;
        end
        LV = (OT(:,2)>0) & (OT(:,3)>crit); 
        % flags trials with evidence of a pause
        plot(tv(LV),OT(LV,1),'.')
        xlim(Xlm)
        ylim([-.15 .3])
        hold on
        plot(xlim,[0 0],'k--')
        plot(AcqTrlVecs(S,:),ylim,'k--')
        xlabel('Trial')
        ylabel('Pause Onset')
    else
        subplot(3,2,2)
        plot(1:length(OT),OT(:,2));
        xlim(Xlm)
        ylim([-100 150])
        hold on
        plot(xlim,[0 0],'k--')
        plot(AcqTrlVecs(S,:),ylim,'k--')
        title(['Cell' num2str(S)])

        subplot(3,2,4)
        plot(1:length(OT),OT(:,3));
        xlim(Xlm)
        ylim([-.7 5])
        hold on
        plot(xlim,[0 0],'k--')
        plot(AcqTrlVecs(S,:),ylim,'k--')

        subplot(3,2,6)
        if S == 1
            crit = .3;
        else
            crit = 1;
        end
        LV = (OT(:,2)>0) & (OT(:,3)>crit); 
        % flags trials with evidence of a pause
        plot(tv(LV),OT(LV,1),'.')
        xlim(Xlm)
        ylim([-.15 .3])
        hold on
        plot(xlim,[0 0],'k--')
        plot(AcqTrlVecs(S,:),ylim,'k--')
        xlabel('Trial')    
    end
end
%%
Experiment.PauseAcquisitionTrials = AcqTrlVecs(:,1);
for S = 1:10
    Experiment.Subject(S).Session(1).PauseAcqTrial = AcqTrlVecs(S,1);
end
%%
TSlimit('Subjects',2:10)
crit = 1;
TSapplystat('PsOnsetLatencies',{'PutativePsOnsets' 'PauseAcqTrial'},@PsOnLats,crit)
%{
function L = PsOnLats(ppo,trl,crit)
LVp = ppo(:,2)>0; % positive (p_pre - p_ps)
LVw = ppo(:,3)>crit;
LVpst = [false(trl,1);true(length(ppo)-trl,1)];
L = ppo(LVpst&LVw&LVp,1);
%}
TSlimit('Subjects',1)
crit = .3;
TSapplystat('PsOnsetLatencies',{'PutativePsOnsets' 'PauseAcqTrial'},@PsOnLats,crit)
%% Cumulative distributions of pause onset latencies
TSlimit('Subjects','all')
TSapplystat('','PsOnsetLatencies',@TSplotcdfs,'Rows',5,'Xlbl','','Xlm',[-.1 .25])