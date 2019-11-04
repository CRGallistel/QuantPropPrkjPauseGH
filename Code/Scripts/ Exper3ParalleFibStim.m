% script m-file FredParalleFib.m
% This code analyzes the data from Fredrik's parallel fiber experiment.
% Much of it is copied from the FredrickExperimentCode.m that analyzed the
% data from his experiment in which the CS was stimulation of the dorsum of
% the paw. The initial portions setting up the Experiment structure are
% copied from the FredExperCode.m file
%{
Quote from Fredrik email: "So, unfortunately it is a bit messy. I did those
experiments over 5 years because they were so difficult to do, so we
improved (or at least we think we did) the protocol along the way. In the
end it did not seem to matter.

ISI 150: 100 Hz, 300 ms CS"
[In the final analysis, I excluded two of these 10 data sets (e1007 & e1032).
In e1007, the eye discerns a pause, but there is a total shut down of firing
for about half a second at CS offset and this defeats the algorithm. In e1032
the eye discerns a very slight decrease in firing during the CS with, on
most trials a clear increase in firing at the time the US is anticipated
and a prounced diminution in firing at CS offset. Alos, the firing looks
oddly even, which inspires my mistrust --CRG]
 
"ISI 200: 100 Hz, two cells 200 ms (isi200_e827 & isi200_e747) CS
and five cells 800 ms CS (the rest)
ISI 300: 50 Hz, 300 ms CS"
%}

cd('/Users/galliste/Documents/Current MSs/PK Pause/FredExper3')
CD = cd;
Subs = sort([885 889 913 9171 9172 988 994 1023 747 762 772 786 ...
    819 820 827 711 7331 7332 736 737 800]); % NB 7331 & 7332 were cells 
%   from the same subject! Ditto for 9171 & 9172
TSinitexperiment('/Users/galliste/Documents/Current MSs/PK Pause/FredExper3/FredParallelFib',...
    101,Subs,'electrodes','Hesslow')
Experiment.Info.LoadFunction='LoadFredrik';
Experiment.Info.InputTimeUnit = 1;
Experiment.Info.OutputTimeUnit=1;
Experiment.Info.FileExtension ='.mat'; % distinguishes data files
Experiment.Info.FilePrefix = ''; % There is no initial character of sequence
% of characters that picks out the data files

TSexperimentbrowser % Calls the browser window that enables browsing around
% the ever growing, mulitlevel Experiment structure
addpath([CD '/TSlib'],[CD '/HelperFunctions'])
addpath([CD '/TSlib/GUI'])
addpath([CD '/TSlib/Support'])
%% Loading data
disp('Choose ''ISI 150 (10 cells)''')
TSloadsessions
% ISI 150 (10 cells), marker is 0.2 sec before CS onset
disp('Choose ''ISI 200 (7 cells)''')
TSloadsessions % Choose ISI 200 (7 cells), marker is 0.2 sec before CS onset
disp('Choose ''ISI 300 (6 cells)''')
TSloadsessions % Choose ISI ISI 300 (6 cells), marker is 0.2 sec before CS onset
TSimporteventcodes('eventcodes.txt') % importing event codes for mark, spike, CSon, CSoff and
%  MidITI events. This command brings up a browse window for you to browse
%  for the file named eventcodes.
%% Creating fields that identify which Subject index #s belong to which CS-US
% training interval
Experiment.CSUS150=[12 13 14 15 16 17 20 21];
Experiment.CSUS200 = [4 5 6 7 9 10 11];
Experiment.CSUS300 = [1 2 3 8 18 19];
%% Putting CS parameters in the appropriate notes fields and noting where
% a "subject" is actually the 2nd cell. (It was easier to modify the file
% names than to write the load function to accommodate the two subjects in
% which two cells were recorded)
for S = [5 6 7 9 10]
    Experiment.Subject(S).SubNotes = '0.8s CS @ 100HZ; 0.2s ISI';
    Experiment.Subject(S).Session(1).SesNotes = '0.8s CS @ 100Hz';
end
for S = [4 11]
    Experiment.Subject(S).SubNotes = '0.2s CS @ 100HZ; 0.2s ISI';
    Experiment.Subject(S).Session(1).SesNotes = '0.2s CS @ 100Hz';
end
%
for S = [12 13 14 15 16 17 20 21]
    Experiment.Subject(S).SubNotes = '0.3s CS @ 100 Hz; 0.15s ISI';
    Experiment.Subject(S).Session(1).SesNotes = '0.3s CS @ 100Hz';    
end
Experiment.Subject(21).SubNotes = '2nd cell from same subject as 9171: 0.3s CS @ 100 Hz; 0.15s ISI';
Experiment.Subject(21).Session(1).SesNotes = '2nd cell from same subject as 9171';
%
for S = [1 2 3 8 18 19]
    Experiment.Subject(S).SubNotes = '0.3s CS @ 50 Hz; 0.3s ISI';
    Experiment.Subject(S).Session(1).SesNotes = '0.3s CS @ 50Hz';    
end
Experiment.Subject(19).SubNotes = '2nd cell from same subject as 7331';
Experiment.Subject(19).Session(1).SesNotes = '2nd cell from same subject as 7331';
%% Begin a day's analysis (comment out this cell's code when running whole script)
TSloadexperiment('/Users/galliste/Documents/Current MSs/PK Pause/FredExper3/FredParallelFib')
TSexperimentbrowser
cd('/Users/galliste/Documents/Current MSs/PK Pause/FredExper3')
addpath('/Users/galliste/Documents/Current MSs/PK Pause/FredExper3/HelperFunctions')
%%
TSapplystat('TSDataA',{'TSData' 'Phase'},@AugTSdata) % creates TSDataA
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
TSsetdata('TSDataA')
%% Raster plots for initial examination of data
% Experiment structure must be in workspace. Three views are shown of each
% cell's raster: 1) the overall view (complete cycles) 2) the view limited
% to 0.2 s before CS onset to 6.5s after CS onset; 3) the view limited to
% 0.2 s before CS onset to 0.8 s after CS onset. You get the successive
% views by hitting the space key (or any key). To exit, activate workspace
% and then hit Cntl C
for S = 1:21 % stepping through the specified subjects
    % 1:Experiment.Subject(S).NumSessions
    for s = 1 % stepping through the cells
        % for a given subject 
        TSlimit('Subjects',S); % data sets to be examined
        TSlimit('Sessions', s)
        CSdur = Experiment.Subject(S).Session(s).Phase; % retrieving
          % duration of the CS
        Xlim1 = [0 6.5]; % first pair of x-axis limits
        Xlim2 = [0 1.1]; % second pair of x-axis limits (closer look)
        TSapplystat('','TSDataA',@DispRaster,CSdur,Xlim1,Xlim2) % execute this cell
        % to see raster plots at a given CS duration. After examining full plot,
        % hit any key to impose first pair of x-axis limits. After examining that
        % view, hit any key to impose 2nd pair. To then go on to next data set,
        % again hit any key.
     %{
      function DispRaster(tsdat,CSdur,Xlm1,Xlm2)
    % displays rasters for a specified CS duration with pause and then Xlm1 and
    % the Xlm2
    if CSdur>1
        CSdur = CSdur/1000;
    end

    TSraster(tsdat,{[20 20] [20 inf]},[20 0;40 0;60 0],['r+';'k.';'r*'])
        % red crosses mark 0.2s prior to CS onset; black dots marks spikes;
        % red asterisks mark middle of intertrial interval
    E = evalin('caller','sub');
    C = evalin('caller','ses');
    if isempty(C)
        title(['Electrode Index ' num2str(E)])
    else
        title(['Electrode Index ' num2str(E), ', Cell ' num2str(C)])
    end
    pause
    hold on
    xlim(Xlm1)
    plot([.2 .2],ylim,'b',[.2+CSdur .2+CSdur],ylim,'b')
    pause
    xlim(Xlm2)
    pause
    % close all
   %}
    end
    % close all
end
TSsaveexperiment
%% Fano Factors
% Spike Counts in 1 s start and end intervals PreCS
TSlimit('Subjects','all')
TSlimit('Sessions','all')
TSdefinetrialtype('MidCSon',[midITI CSon])
TStrialstat('SpkTms',@TSparse,'result=time(1)-starttime;',spike)
TSapplystat('FrstAndLast1sSpkCnt','SpkTms',@SpkCnts) % trial level
%{
function O = SpkCnts(T)
O = [];
if isempty (T) || T(end)-T(1)<2 % only continue if spike train at least 2s long
    return
end
LV1 = T>T(1) & T<=T(1)+1;
LV2 = T>T(end)-1;
O = [sum(LV1)+1 sum(LV2)];
%}
TScombineover('PreSpkCnts','FrstAndLast1sSpkCnt') % session/cell level
TSapplystat('FanoFacsPreCS','PreSpkCnts',@FanoFac) % session/cell level
%{
function O = FanoFac(Cnts)
if isempty(Cnts) || size(Cnts,1)<10
    O=[];
else
    O = [var(Cnts(:,1))/mean(Cnts(:,1)) var(Cnts(:,2))/mean(Cnts(:,2)) size(Cnts,1)];
end % 2 Fano Factors and the n (sample size)
%}
TScombineover('PreCSFanos','FanoFacsPreCS','t') % subject level
TScombineover('PreCSFanoFactors','PreCSFanos','t') % Experiment level

% Spike Counts in 1 s start and end intervals PostCS
TSdefinetrialtype('PostCS',[CSoff midITI])
TStrialstat('SpkTms',@TSparse,'result=time(1)-starttime;',spike)
TSapplystat('FrstAndLast1sSpkCnt','SpkTms',@SpkCnts) % see above cell for helper
TScombineover('PostSpkCnts','FrstAndLast1sSpkCnt') % session level
TSapplystat('FanoFacsPostCS','PostSpkCnts',@FanoFac) % session level; see above cell for helper
TScombineover('PostCSFanos','FanoFacsPostCS','t') % subject level
TScombineover('PostCSFanoFactors','PostCSFanos','t') % Experiment level

% CDFs of logged Fano Factors
figure
subplot(1,2,1)
cdfplot(log10(Experiment.PreCSFanoFactors(:,1)))
hold on
cdfplot(log10(Experiment.PreCSFanoFactors(:,2)))
plot([log10(.47) log10(.47)],ylim,'k--',[log10(1.73) log10(1.73)],ylim,'k--')
xlim([-1 1])
set(gca,'XTick',[-1 -.7 -.3 0 .3 .7 1],'XTickLabel',{'.1' '.2' '.5' '1' '2' '5' '10'})
xlabel('Fano Factor (log scale)')
ylabel('Cumulative Fraction of Cells')
title('Pre')

subplot(1,2,2)
cdfplot(log10(Experiment.PostCSFanoFactors(:,1)))
hold on
cdfplot(log10(Experiment.PostCSFanoFactors(:,2)))
plot([log10(.47) log10(.47)],ylim,'k--',[log10(1.73) log10(1.73)],ylim,'k--')
set(gca,'XTick',[-1 -.7 -.3 0 .3 .7 1],'XTickLabel',{'.1' '.2' '.5' '1' '2' '5' '10'})
xlabel('Fano Factor (log scale)')
title('Post')
saveas(gcf,'Figures/CDFsPre&PostFanoFactors')
disp(char({'This figure gives the distribution of Fano Factors for the';...
    'pre- and post-CS spike counts in 1s windows. The dashed lines';...
    'indicated lower and upper limits expected on the hypothesis that';...
    'the spikes were generated by a stationery Poisson process. The';...
    'fact that 50% factors were greater than the upper limit implies';...
    'that the generative process was non-stationary';''}))

%% Delimiting pauses using BinPsOnOff
TSlimit('Subjects','all')
TSsetdata('TSDataA')
TSdefinetrialtype('MdToMd',[midITI midITI])
BW = .001; % bin width for binary spike vector
% TStrialstat('BinSpkVec',@Bindata2,BW)
% 
TStrialstatMult({'BinSpkVec' 'FrstSpkTm' 'LastSpkTm'},@Bindata2,BW,spike,CSon)
% creates a field in the MdToMd trials whose first column is spike time
% referenced to CS onset at 0 and whose 2nd col is the binary vector
% representation of the spike train. The Errors output is a cell array that
% gives the error messages where the data caused errors. 
%{
function [spikes,Mn,Mx] = Bindata2(tsd,BW)
% A Benoulli vector digitization of the spike train with bins of  width BW
% Syntax     spikes = Bindata(tsd,BW)
LV = tsd(:,2)==40; % flags spikes
LVCSon = tsd(:,2)==30; % flags CS on
spks = tsd(LV,1)-tsd(LVCSon,1); % spike times referenced to 0 at CS on
Mn = min(spks);
Mx = max(spks);
bins = Mn:BW:Mx; % creating the bin edges
BV = histc(spks,bins);
t = bins'+BW/2;
spikes = [t BV];
%}
% 
TScombineover('BinSpkVecS','BinSpkVec','t')
% field at the Session level with all the binary vectors with a 3rd col
% with integers specifying which trial
%% Probability distributions for preCS, CS, and postCS intervals
Edges = [0 .005 .01 .015 .02 .025 .03 .04 inf]; 
TSapplystat({'PreCSispkiPDF' 'CSispkiPDF' 'PostCSispkiPDF' 'Edges'},...
    'TSDataA',@PrbDst,Edges,midITI,CSon,CSoff)
%% 
TSlimit('Subjects',Experiment.CSUS150)
CSdur=.15;
TSapplystat({'PsOn' 'WpsOn' 'PsOff' 'WpsOff' 'lmPre' 'lmDur' 'lmAft' ...
    'MxIspkI' 'LatOfMxIspkI'},'BinSpkVecS',@BinPsOnOff2,CSdur)
% Creates 9 fields at the Session level: PsOn = vector of pause onset
% estimates; WpsOn = vector of the weights of the evidence for those
% onsets; PsOff = vector of estimates of the pause off times; WpsOff =
% vector of the weights of the evidence for those offsets; lmPre = vector
% of the rates of pressing in the interval of duration dur preceding CS
% onset; lmDur = vector of the rates in the middle 40% of the CS; lmAft =
% vector of the rates in the interval of duration dur after the offset of
% the CS; MxIspkI = maximum interspike interval btw PsOn and PsOff;
% LatOfMxIspkI = interval from CSon to onset of maximum interspike interval
% (set to 0 when maximum interspike intervals begins at PsOn and PsOn is
% negative)
%{
function [On,Won,Off,Woff,lb,ld,la,MxIpI,LatMx] = BinPsOnOff2(D,dur)
% Computes the onsets and offsets of the pauses
D(D(:,2)>1,2)=1; % This corrects a single bin with a count of 2 in the
% BinSpkVecS field for S10 in the parallel fiber data. That triggers a
% crash of this function when BernCP is called
r = find(D(:,1)>-.0005&D(:,1)<.0005); % vector of row numbers for CS onset
NmTrls = length(r);
On = nan(NmTrls,1); % initializing On vector
Off = nan(NmTrls,1); % initializing Off vector
Won = zeros(NmTrls,1); % initializing vector for weight of the evidence for
% pause onset
Woff = zeros(NmTrls,1); % inti
lb = nan(NmTrls,1);
ld = nan(NmTrls,1);
la = nan(NmTrls,1);
MxIpI = nan(NmTrls,1);
LatMx = nan(NmTrls,1);
LVs = D(:,2)>0; % flags rows with spikes
LV = D(:,1)>-dur & D(:,1)<3*dur & LVs; % flags rows w spikes btw -dur and
% 3dur
edges = -dur:.2*dur:3*dur;
SpkCnts= histc(D(LV,1),edges);
[Nsdur,loc] = min(SpkCnts(1:end-1)); % lowest spike count within an interval
% equal to 1/5th of CS duration and the bin that contains it 
% 
LVpre = D(:,1)<0 & D(:,1)>-2*dur; % flags pre CS rows over an immediately
% preceding interval equal in duration to the CS
LVpst = D(:,1)>dur & D(:,1)<2*dur; % flags post CS rows over an immediately
% following interval equal in duration to the CS
LVb = D(:,1)>edges(loc); % flags rows greater than lower edge of lowest-count 
%  bin in the histogram
LVe = D(:,1)<edges(loc)+.2*dur; % flags rows less than the upper edge of 
% the lowest-count bin in the histogram
LVmid = LVb&LVe;
Npre = sum(LVpre); % number of pretrial bins 
Nspre = sum(LVpre&LVs); % number of pretrial spikes
if Nspre<1
    Nspre = 1;
end
HypPre = [1 round((Npre-Nspre)/Nspre)]; % hyperparameters of beta prior on
    % pre-CS probability of a spike in a 1 ms bin
if isnan(HypPre(2)) || ~beta(HypPre(1),HypPre(2))>0 
    keyboard
end
Ndur = sum(LVmid); % number of binarize bins in lowest bin of histogram 
if Nsdur<1
    Nsdur = 1;
end
HypDur = [1 round((Ndur-Nsdur)/Nsdur)]; % hyperparameters on the probability
% of a spike in any one ms in the middle 40% of the CS bins
if isnan(HypDur(2)) || ~beta(HypDur(1),HypDur(2))>0
    keyboard
end
Npst = sum(LVpst); % # of post-CS bins
Nspst = sum(LVpst&LVs); % # of spikes in post-CS bins
if Nspst<1
    Nspst = 1;
end
HypPst = [1 round((Npst-Nspst)/Nspst)];
if isnan(HypPst(2)) || ~beta(HypPst(1),HypPst(2))>0
    keyboard
end
nbpers = 1000; % number of bins per second
%
NbinsBack = dur*1000; % number of bins prior to CS onset for pause On
NbinsFwd = dur*700; % number of bins after CS onset for pause On
NbinsToStrt = 300*dur; % number of bins to start of data for pause Off
NbinsPst = 2000*dur; % number of bins to end of data for pause Off
%%
for tt = 1:NmTrls
    %%
    if r(tt)-NbinsBack < 1 % happens in some data sets
        NbinsBack=r(tt)-1;
    end
    Don = flipud(D(r(tt)-NbinsBack:r(tt)+NbinsFwd,2)); % the binary spike
    % vector running backwards in time from 70% of the way through the CS
    % (whence the 700 times a time equal to CS onset minus CS duration in
    % seconds)
    try
    [CPon,OddsD] = BernCP_local(Don,HypDur(1),HypDur(2),HypPre(1),HypPre(2),1/length(Don));
    % the onset cp in # of bins counting backwards
    catch ME
        keyboard
    end
    %
    Doff = D(r(tt)+NbinsToStrt:r(tt)+NbinsPst,2); % data for finding pause
    % offset run from 70% of CS-US interval to 200% of it, that is, from .7
    % of the way through the interval to twice the interval
    
    [CPoff,OddsU] = BernCP_local(Doff,HypDur(1),HypDur(2),HypPst(1),HypPst(2),1/length(Doff));
%     if tt == 12
%         keyboard
%     end
    if sum(Don(1:CPon))>sum(Don(CPon:end))...
            ||sum(Doff(1:CPoff))>sum(Doff(CPoff:end))...
            ||sum(Don)+sum(Doff)<10
        % more spikes after CPon than before (remember: it's backwards here)
        % or during CP than after or too few spikes to make valid decision
        continue % don't record pause onsets & offsets
    end
%     if tt == 15;keyboard;end
    On(tt) = (dur*700-CPon)/1000; % the estimated pause onset time in seconds
    % relative to the CS onset time
    Won(tt) = log10(OddsD); % the weight of the evidence for the onset
    Off(tt) = (dur*300+CPoff)/1000; % the estimated pause offset time in seconds
    % relative to a 0 at CS onset
    Woff(tt) = log10(OddsU); % weight of the evidence for the pause offset
    lb(tt) = sum(D(r(tt)-NbinsBack:r(tt),2))/dur;
    ld(tt) = sum(D(r(tt)+NbinsToStrt:r(tt)+NbinsFwd,2))/dur;
    la(tt) = sum(D(r(tt)+dur*nbpers:r(tt)+2*dur*nbpers,2))/dur;
    % Following code added 04/23/18
    LVtt = D(:,3)==tt; % flags data from current trial
    LVps = D(:,1)>On(tt) & D(:,1)<Off(tt); % flags times where trial time
    % falls within the PsOn and PsOff times for this trial (but across all
    % trials)
    SpkTms = D(LVtt&LVps&LVs,1); % spike times within the pause boundaries on
    % THIS trial
    if isempty(SpkTms)
        MxIpI(tt) = Off(tt)-On(tt);
        LatMx(tt) = 0;
    else
        Tms = [On(tt);SpkTms;Off(tt)];
        IspIs = diff(Tms); % intervals
        [MxIpI(tt),indx] = max(IspIs);
        LatMx(tt) = Tms(indx)-On(tt);
        if LatMx(tt)<0;LatMx(tt)=0;end
    end        
end

function [CP,Odds,LkFun,CL,CPmean] = BernCP_local(Data,AlphaB,BetaB,AlphaA,BetaA,pc)
% Computes maximally likely change point (CP) and odds in favor of a change
% given a binary data string, Data, the parameters of the beta prior
% on the value of p before the change (AlphaB & BetaB), the parameters
% of the beta prior on the p after the change (AlphaA,BetaA), and the
% prior probability of a change after any given trial (pc). This is called
% from Line 515 with AlphaB & BetaB evolving and from Line 597 with them
% set to their initial agnostic values (.5 or 1)
%
%Syntax [CP,Odds,LkFun,CL] = BernCP(Data,AlphaB,BetaB,AlphaA,BetaA,pc)
%
% CL is a cell array giving lower (L) and upper (U) relative likelihood
% limits as follows:
% CL = {[.01L] [.05L] [.1L] [.2L] [.5L] [.5U] [.2U] [.1U] [.05U] [.01U]}
% CPmean is the mean of the marginal likelihood distribution on CP. The CP
% output is this value rounded to the nearest integer

Data = reshape(Data,length(Data),1); % making sure Data is a column vector

Data = Data+0; % converts logic vectors to 0's & 1's

if ~isempty(setdiff(Data,[0 1]))
    disp('Error: Input is not a binary vector')
    return
end

CP = nan; % intializing

Odds = nan;% intializing

LkFun = zeros(length(Data),1);% intializing

CL(1,1:10) = {[]};% intializing

CPmean = nan;

if nargin<6
    display({['Not enough input arguments;'];['6 required: '];...
        ['Vector of binary data'];...
        ['Hyperparameter alpha for prior on p before change'];...
        ['Hyperparameter beta for prior on p before change'];...
        ['Hyperparameter alpha for prior on p after change'];...
        ['Hyperparameter beta for prior on p after change'];...
        ['Probability of a change after any one trial']})
    return
end
    
if ~(beta(AlphaB,BetaB)>0)||~(beta(AlphaA,BetaA)>0)
    
    display('Prior too strong; makes beta(a,b)=0')
    
    return
    
end
%%

p = .01:.01:.99; % row vector of p values

dp = p(2)-p(1); % delta p

L = (1:length(Data))'; % cumulative number of observations (col vector)

Sf = cumsum(Data); % cumulative successes up to and including the nth datum
% (col vector w length = length(Data))

Ff = L - Sf; % cumulative failures up to and including the nth datum
% (col vector w length = the length of the data)

Sr = sum(Data) - Sf; % succeses after nth datum (col vector w length =
% the length of the data))

Lr = ((L(end)-1):-1:0)'; % number of observations after the nth datum
% (length = the length of the data)

Fr = Lr - Sr; % failures after the nth datum (col vector, w length =
% the length of the data)

%% The following relies on the analytic result that
% given a beta prior with parameters, alpha & beta, and observed
% numbers of successes and failures, nS & nF, the posterior marginal
% likelihood of a change hypotheses, given the data and the assumption of
% a beta prior with parameters alpha & beta) is beta(alpha+nS,beta+nF).

LgPstLkf = log(beta(AlphaB + Sf, BetaB + Ff)); % col vector
% giving the log likelihood for the data up to and including datum(n) under
% the prior distribution on p specified by hyperparameters AlphaB and BetaB

LgPstLkr = log(beta(AlphaA + Sr, BetaA + Fr)); % col vector
% giving the log likelihood for the data after datum(n) under the prior
% distribution specified by hyperparameters AlphaA & BetaA


%%
LgLkFun = LgPstLkf + LgPstLkr; % the log likelihood function giving, as a
% function of cp, the log likelihood of a model that assumes a change in p 
% after observation cp.

LgLkFunN = LgLkFun + abs(max(LgLkFun)); % setting peak value of log
% likelihood function to 0, so that when it is transformed into a
% likelihood function prior to integration its values are reasonable. This
% is equivalent to rescaling the likelihood function, which has no effect
% on the relative likelihood or the expectation, which are the two
% quantities of interest

LkFun = exp(LgLkFunN); % the likelihood function
% The last value in this vector is the likelihood of a model
% that assumes a change after the end of the data, in other words, no
% change within the span of the observations. Thus, assuming a uniform
% prior on cp (change equally likely after any trial), the marginal
% likelihood of a change is the mean value of this likelihood function

CPmean = sum((1:length(Data)).*LkFun')/sum(LkFun); % expectation of the
% likelihood function (NB, not its mean value)

RelLkhd = (sum(LkFun(1:end))/(length(Data)))/LkFun(end); % the Bayes
% Factor, the ratio of the posterior likelihoods of the two models (change
% vs no-change). (sum(LkFun(1:end))/(length(Data)) is the mean likelihood,
% which, assuming a uniform prior, is the the marginal likelihood of a
% change; LkFun(end) is the likelihood of no change.
%%

pNC = 1-pc; % Probabilty that p does not change after any given trial

PriorOdds = L(end)*pc/pNC; % prior odds that there is a change
% in a sequence as long as Data = Np_c/(1-p_c)

Odds = PriorOdds*RelLkhd; % Posterior odds equal relative likelihood of the
% two models (change: no-change) times the prior odds in favor of one or
% the other

BF = RelLkhd;
%% Computing CP & relative likelihood limits on its location

CP = round(CPmean);

Mx = max(LkFun(1:end-1)); % maximum of the likelihood function for
% change points before the last datum and location of this maximum
% (the estimate of CP)

CL{1} = find(LkFun>.01*Mx,1)-1;

CL{2} = find(LkFun>.025*Mx,1)-1;

CL{3} = find(LkFun>.1*Mx,1)-1;

CL{4} = find(LkFun>.2*Mx,1)-1;

CL{5} = find(LkFun>.5*Mx,1)-1;

CL{6} = find(LkFun>.5*Mx,1,'last')+1;

CL{7} = find(LkFun>.2*Mx,1,'last')+1;

CL{8} = find(LkFun>.1*Mx,1,'last')+1;

CL{9} = find(LkFun>.025*Mx,1,'last')+1;

CL{10} = find(LkFun>.01*Mx,1,'last')+1;
%}

TSlimit('Subjects',16) % to get proper estimation of pause onsets and
% offsets in this data set, I had to shorten the post-CS interval that is
% hard coded into the pause parameter algorithm to equal the CS duration
% rather than to be twice the CS duration, which is why this calls
% BinPsOnOff3 as a helper rather than BinPsOnOff2
TSapplystat({'PsOn' 'WpsOn' 'PsOff' 'WpsOff' 'lmPre' 'lmDur' 'lmAft' ...
    'MxIspkI' 'LatOfMxIspkI'},'BinSpkVecS',@BinPsOnOff3,CSdur)
%
TSlimit('Subjects',Experiment.CSUS200) % Experiment.CS200electrodes
CSdur=.2;
% TStrialstat('BinSpkVec',@Bindata2,BW)
TSapplystat({'PsOn' 'WpsOn' 'PsOff' 'WpsOff' 'lmPre' 'lmDur' 'lmAft' ...
    'MxIspkI' 'LatOfMxIspkI'},'BinSpkVecS',@BinPsOnOff2,CSdur)
% 
TSlimit('Subjects',Experiment.CSUS300) %
CSdur=.3;
% TStrialstat('BinSpkVec',@Bindata2,BW)
TSapplystat({'PsOn' 'WpsOn' 'PsOff' 'WpsOff' 'lmPre' 'lmDur' 'lmAft' ...
    'MxIspkI' 'LatOfMxIspkI'},'BinSpkVecS',@BinPsOnOff2,CSdur)

TSsaveexperiment
%% Raster plots with Probability Distributions (PDFs)
% Experiment structure must be in workspace. This plotting loop has to be
% modified for each CSUS conditions: Set the CSdur appropriately and set
% the range for SwiPhase appropriately (so that upper end is # of cells in
% that condition
Xlim3 = [0 1.1]; 
CSdur =150; % specify CS duration condition for which one wants plots
for SwiPhase=1:8; % The subject within that condition
    Cond = ['CSUS' num2str(CSdur)];
    for c = 1:Experiment.Subject(Experiment.(Cond)(SwiPhase)).NumSessions
        % stepping through the sessions (cells) for the specified subject in
        % the specified condition
        figure;
        Ax1 = subplot('Position',[.1 .1 .4 .8]);
        Ax2 = subplot('Position',[.6 .675 .3 .225]);
        Ax3 = subplot('Position',[.6 .3825 .3 .225]);
        Ax4 = subplot('Position',[.6 .1 .3 .225]);
        %
        TSlimit('Subjects',Experiment.(Cond)(SwiPhase)); % data set
        %  to be examined 
        Xlim1 = [0 1.5]; % first pair of x-axis limits

        TSapplystat('',{'TSDataA' 'Phase' 'PsOn' 'PsOff'},@DispRaster4,...
                Xlim3,Ax1) % raster plot in left subplot
        %{
        function DispRaster4(tsdat,CSdur,Bs,Es,Xlm1,Ax)
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
        %}
%         title([num2str(Eflds{C}(1:5)) ', Sub ' num2str(S) ', Cell ' num2str(s)])
        TSapplystat('',{'PreCSispkiPDF' 'CSispkiPDF' 'PostCSispkiPDF'},...
            @PlotPrbDists,[Ax2 Ax3 Ax4]) % probability distributions in right
        % subplots
            
        %{
        function PlotPrbDists(preD,csD,postD,Hv)
        axes(Hv(1))
        h=bar(preD,1); % bar plot of prob dist in upper right panel
        ylim([0 .7])
        set(h,'FaceColor',[1 1 1],'LineWidth',2)
        ylabel('Probability','FontSize',12)
        set(gca,'XTick',.5:1:8.5,'XTickLabel',{'0' '5' '10' '15' '20' '25' '30' '40' '>40'})
        text(.2,.6,'PreCS')

        axes(Hv(2))
        h=bar(csD,1);
        ylim([0 .7])
        set(h,'FaceColor',[1 1 1],'LineWidth',2)
        ylabel('Probability','FontSize',12)
        set(gca,'XTick',.5:1:8.5,'XTickLabel',{'0' '5' '10' '15' '20' '25' '30' '40' '>40'})
        text(.2,.6,'CS')

        https://www.zomato.com/new-york-city/cull-pistol-oyster-bar-chelsea/menuaxes(Hv(3))
        h=bar(postD,1);
        ylim([0 .7])
        set(h,'FaceColor',[1 1 1],'LineWidth',2)
        ylabel('Probability','FontSize',12)
        set(gca,'XTick',.5:1:8.5,'XTickLabel',{'0' '5' '10' '15' '20' '25' '30' '40' '>40'})
        xlabel('InterspikeInterval (ms)','FontSize',12)
        text(.2,.6,'PostCS')   
        %}
    end % of plotting a figure for one cell showing raster plot on left and the
    % three probability distributions on the right
end
%% Saving figures; This cell also needs to be modified for each group
for f = 1:8
    S = Experiment.CSUS150(f);
    figure(f)
    saveas(gcf,['Figures/RastersWpauseOnOffs&PDFs/CSUS150rasters&pdfs-fig'...
        num2str(f) 'of 8.pdf'])
end

%% Quartiles of the Pause Stat Distributions for individual cells. 
TSlimit('Subjects','all')
QV = [.25 .5 .75];
TSapplystat('Qs1to3_PsOn','PsOn',@quantile,QV)
TSapplystat('Qs1to3_PsOff','PsOff',@quantile,QV)
TSapplystat('Qs1to3_lmPre','lmPre',@quantile,QV)
TSapplystat('Qs1to3_lmDur','lmDur',@quantile,QV)
TSapplystat('Qs1to3_lmAft','lmAft',@quantile,QV)
TSapplystat('Qs1to3_WpsOn','WpsOn',@quantile,QV)
TSapplystat('Qs1to3_WpsOff','WpsOff',@quantile,QV)
TSapplystat('Qs1to3_MxIspkI','MxIspkI',@quantile,QV)
TSapplystat('Qs1to3_LatOfMxIspkI','LatOfMxIspkI',@quantile,QV)
%% Means and std's of pause onsets and offsets
TSapplystat('PsOnMeanStdCoV','PsOn',@meanstdcov)
%{
function O = meanstdcov(D)
O = [nanmean(D) nanstd(D) nanstd(D)/nanmean(D)];
%}
TSapplystat('PsOffMeanStdCoV','PsOff',@meanstdcov)
%% Pause widths
TSapplystat('PsWidths',{'PsOff' 'PsOn'},@minus)
TSapplystat('Qs1to3_PsWidths','PsWidths',@quantile,QV)
TSapplystat('PsWidthMeanStdCoV','PsWidths',@meanstdcov)

%% Pause parameter correlations
TSapplystat('PsParamCorrelations',{'PsOn' 'PsOff' 'PsWidths' 'MxIspkI'},@PsParamCorr)
% 4x4 array at Session level giving pairwise correlations: A(1,2) is the
% correlation of PsOff w PsOn; A(1,3) is the correlation of PsWidth w PsOn;
% A(1,4) is the correlation of MxIspI w PsOn; A(2,3) is the correlation of
% PsOff with PsWidth; A(2,4) is the correlation of MxIspI w PsOff; A(3,4)
% is the correlation of MxIspI w PsWidth
%{
function Corrs = PsParamCorr(On,Off,Wdth,MxIsI)
LV  = ~isnan(On);
Corrs = corr([On(LV) Off(LV) Wdth(LV) MxIsI(LV)]);
%}
%%
TSapplystat('PsParamsRowVec','PsParamCorrelations',@PProwvec)
% Reduces preceding 4x4 array to a row vector of the 6 correlations
%{
function RV = PProwvec(A)
RV = [A(1,2) A(1,3) A(1,4) A(2,3) A(2,4) A(3,4)];
%}
%% Carrying Pause Stats up to Subject level
TScombineover('PsOn_S','PsOn')
TScombineover('PsOff_S','PsOff')
TScombineover('WpsOn_S','WpsOn')
TScombineover('WpsOff_S','WpsOff')
TScombineover('lmPre_S','lmPre')
TScombineover('lmDur_S','lmDur')
TScombineover('lmAft_S','lmAft')
TScombineover('Qs1to3_PsOn_S','Qs1to3_PsOn')
TScombineover('Qs1to3_PsOff_S','Qs1to3_PsOff')
TScombineover('Qs1to3_lmPre_S','Qs1to3_lmPre')
TScombineover('Qs1to3_lmDur_S','Qs1to3_lmDur')
TScombineover('Qs1to3_lmAft_S','Qs1to3_lmAft')
TScombineover('Qs1to3_WpsOn_S','Qs1to3_WpsOn')
TScombineover('Qs1to3_WpsOn_S','Qs1to3_WpsOn')
TScombineover('MxIspkI_S','MxIspkI')
TScombineover('LatOfMxIspkI_S','LatOfMxIspkI')
TScombineover('Qs1to3_MxIspkI_S','Qs1to3_MxIspkI')
TScombineover('Qs1to3_LatOfMxIspkI_S','Qs1to3_LatOfMxIspkI')
TScombineover('PsOnMeanStdCoV_S','PsOnMeanStdCoV')
TScombineover('PsOffMeanStdCoV_S','PsOffMeanStdCoV')
TScombineover('PsWidthMeanStdCoV_S','PsWidthMeanStdCoV')
TScombineover('Qs1to3_PsWidths_S','Qs1to3_PsWidths')
TScombineover('VecsOfPsParamCorrs','PsParamsRowVec')
%% Carrying Pause Stats up to Group fields at Experiment level
TSlimit('Subjects',Experiment.CSUS150)
TScombineover('G150PsOn','PsOn_S','t')
TScombineover('G150PsOff','PsOff_S','t')
TScombineover('G150WpsOn','WpsOn_S','t')
TScombineover('G150WpsOff','WpsOff_S','t')
TScombineover('G150lmPre','lmPre_S','t')
TScombineover('G150lmDur','lmDur_S','t')
TScombineover('G150lmAft','lmAft_S','t')
TScombineover('G150Qs1to3_PsOn','Qs1to3_PsOn_S','t')
TScombineover('G150Qs1to3_PsOff','Qs1to3_PsOff_S','t')
TScombineover('G150Qs1to3_lmPre','Qs1to3_lmPre_S','t')
TScombineover('G150Qs1to3_lmDur','Qs1to3_lmDur_S','t')
TScombineover('G150Qs1to3_lmAft','Qs1to3_lmAft_S','t')
TScombineover('G150Qs1to3_WpsOn','Qs1to3_WpsOn_S','t')
TScombineover('G150Qs1to3_WpsOn','Qs1to3_WpsOn_S','t')
TScombineover('G150MxIspkIs','MxIspkI_S','t')
TScombineover('G150LatOfMxIspkIs','LatOfMxIspkI_S','t')
TScombineover('G150Qs1to3OfMxIspkIs','Qs1to3_MxIspkI_S','t')
TScombineover('G150Qs1to3OfLatofMxIspkI','Qs1to3_LatOfMxIspkI_S','t')
TScombineover('G150PsOnMeanStdCoV','PsOnMeanStdCoV_S','t')
TScombineover('G150PsOffMeanStdCoV','PsOffMeanStdCoV_S','t')
TScombineover('G150PsWidthMeanStdCoV','PsWidthMeanStdCoV_S','t')
TScombineover('G150Qs1to3_PsWidths','Qs1to3_PsWidths_S')
TScombineover('G150PsParamsCorrelations','VecsOfPsParamCorrs')
%
TSlimit('Subjects',Experiment.CSUS200)
TScombineover('G200PsOn','PsOn_S','t')
TScombineover('G200PsOff','PsOff_S','t')
TScombineover('G200WpsOn','WpsOn_S','t')
TScombineover('G200WpsOff','WpsOff_S','t')
TScombineover('G200lmPre','lmPre_S','t')
TScombineover('G200lmDur','lmDur_S','t')
TScombineover('G200lmAft','lmAft_S','t')
TScombineover('G200Qs1to3_PsOn','Qs1to3_PsOn_S','t')
TScombineover('G200Qs1to3_PsOff','Qs1to3_PsOff_S','t')
TScombineover('G200Qs1to3_lmPre','Qs1to3_lmPre_S','t')
TScombineover('G200Qs1to3_lmDur','Qs1to3_lmDur_S','t')
TScombineover('G200Qs1to3_lmAft','Qs1to3_lmAft_S','t')
TScombineover('G200Qs1to3_WpsOn','Qs1to3_WpsOn_S','t')
TScombineover('G200Qs1to3_WpsOn','Qs1to3_WpsOn_S','t')
TScombineover('G200MxIspkIs','MxIspkI_S','t')
TScombineover('G200LatOfMxIspkIs','LatOfMxIspkI_S','t')
TScombineover('G200Qs1to3OfMxIspkIs','Qs1to3_MxIspkI_S','t')
TScombineover('G200Qs1to3OfLatofMxIspkI','Qs1to3_LatOfMxIspkI_S','t')
TScombineover('G200PsOnMeanStdCoV','PsOnMeanStdCoV_S','t')
TScombineover('G200PsOffMeanStdCoV','PsOffMeanStdCoV_S','t')
TScombineover('G200PsWidthMeanStdCoV','PsWidthMeanStdCoV_S','t')
TScombineover('G200Qs1to3_PsWidths','Qs1to3_PsWidths_S')
TScombineover('G200PsParamsCorrelations','VecsOfPsParamCorrs')
%
TSlimit('Subjects',Experiment.CSUS300)
TScombineover('G300PsOn','PsOn_S','t')
TScombineover('G300PsOff','PsOff_S','t')
TScombineover('G300WpsOn','WpsOn_S','t')
TScombineover('G300WpsOff','WpsOff_S','t')
TScombineover('G300lmPre','lmPre_S','t')
TScombineover('G300lmDur','lmDur_S','t')
TScombineover('G300lmAft','lmAft_S','t')
TScombineover('G300Qs1to3_PsOn','Qs1to3_PsOn_S','t')
TScombineover('G300Qs1to3_PsOff','Qs1to3_PsOff_S','t')
TScombineover('G300Qs1to3_lmPre','Qs1to3_lmPre_S','t')
TScombineover('G300Qs1to3_lmDur','Qs1to3_lmDur_S','t')
TScombineover('G300Qs1to3_lmAft','Qs1to3_lmAft_S','t')
TScombineover('G300Qs1to3_WpsOn','Qs1to3_WpsOn_S','t')
TScombineover('G300Qs1to3_WpsOn','Qs1to3_WpsOn_S','t')
TScombineover('G300MxIspkIs','MxIspkI_S','t')
TScombineover('G300LatOfMxIspkIs','LatOfMxIspkI_S','t')
TScombineover('G300Qs1to3OfMxIspkIs','Qs1to3_MxIspkI_S','t')
TScombineover('G300Qs1to3OfLatofMxIspkI','Qs1to3_LatOfMxIspkI_S','t')
TScombineover('G300PsOnMeanStdCoV','PsOnMeanStdCoV_S','t')
TScombineover('G300PsOffMeanStdCoV','PsOffMeanStdCoV_S','t')
TScombineover('G300PsWidthMeanStdCoV','PsWidthMeanStdCoV_S','t')
TScombineover('G300Qs1to3_PsWidths','Qs1to3_PsWidths_S')
TScombineover('G300PsParamsCorrelations','VecsOfPsParamCorrs')

%% Graphing Pause On & Off and Width their CoVs
figure
subplot(2,3,2)
    D1 = Experiment.G150PsOffMeanStdCoV(:,[1 3]);
    D2 = Experiment.G200PsOffMeanStdCoV(:,[1 3]);
    D3 = Experiment.G300PsOffMeanStdCoV(:,[1 3]);
    plot(repmat(.15,1,length(D1)),D1(:,1),'k*',repmat(.2,1,length(D2)),D2(:,1),'k*',...
        repmat(.3,1,length(D3)),D3(:,1),'k*')
    ylabel('Pause Off Latency (s)')
    xlim([0 .5]); ylim([0 .6])
subplot(2,3,5)
    plot(repmat(.15,1,length(D1)),D1(:,2),'k*',repmat(.2,1,length(D2)),D2(:,2),'k*',...
        repmat(.3,1,length(D3)),D3(:,2),'k*')
    xlim([0 .5]); ylim([0 .6])
    xlabel('CS-US Interval (s)')
%
subplot(2,3,1)
    D1 = Experiment.G150PsOnMeanStdCoV(:,[1 3]);
    D2 = Experiment.G200PsOnMeanStdCoV(:,[1 3]);
    D3 = Experiment.G300PsOnMeanStdCoV(:,[1 3]);
    plot(repmat(.15,1,length(D1)),D1(:,1),'k*',repmat(.2,1,length(D2)),D2(:,1),'k*',...
        repmat(.3,1,length(D3)),D3(:,1),'k*')
    ylabel('Pause On Latency (s)')
    xlim([0 .5]); ylim([0 .6])
subplot(2,3,4)
    plot(repmat(.15,1,length(D1)),D1(:,2),'k*',repmat(.2,1,length(D2)),D2(:,2),'k*',...
        repmat(.3,1,length(D3)),D3(:,2),'k*')
    xlim([0 .5]); ylim([0 .6])
    ylabel('CoV (\sigma/\mu) ')
    xlabel('CS-US Interval (s)')    
subplot(2,3,3)
    D1 = Experiment.G150PsWidthMeanStdCoV(:,[1 3]);
    D2 = Experiment.G200PsWidthMeanStdCoV(:,[1 3]);
    D3 = Experiment.G300PsWidthMeanStdCoV(:,[1 3]);
    plot(repmat(.15,1,length(D1)),D1(:,1),'k*',repmat(.2,1,length(D2)),D2(:,1),'k*',...
        repmat(.3,1,length(D3)),D3(:,1),'k*')
    xlim([0 .5]); ylim([0 .6])
    ylabel('Pause Width (s)')
subplot(2,3,6)   
    plot(repmat(.15,1,length(D1)),D1(:,2),'k*',repmat(.2,1,length(D2)),D2(:,2),'k*',...
        repmat(.3,1,length(D3)),D3(:,2),'k*')
    xlim([0 .5]); ylim([0 .6])
    xlabel('CS-US Interval (s)') 
saveas(gcf,'Figures/PsOnOff&WdthVsCSUSwCoVs.pdf')
%% Graphing Correlations
% A(1,2) is Off vs On; A(1,3) is the correlation of PsWidth w PsOn;
% A(1,4) is the correlation of MxIspI w PsOn; A(2,3) is the correlation of
% PsOff with PsWidth; A(2,4) is the correlation of MxIspI w PsOff; A(3,4)
% is the correlation of MxIspI w PsWidth
figure
subplot(3,1,1)
D1 = Experiment.G150PsParamsCorrelations;
plot(ones(size(D1,1),1),D1(:,1),'k*',2*ones(size(D1,1),1),D1(:,2),'k*',...
    3*ones(size(D1,1),1),D1(:,3),'k*',4*ones(size(D1,1),1),D1(:,4),'k*',...
    5*ones(size(D1,1),1),D1(:,5),'k*',6*ones(size(D1,1),1),D1(:,6),'k*')
% set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{'?v?' 'Wv?' 'Mv?' 'Wv?' 'Mv?' 'MvW'})
set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{char([8595 118 8593]) ...
    char([87 118 8593]) char([77 118 8593]) char([87 118 8595]) ...
    char([77 118 8595]) 'MvW'})
xlim([.5 6.5])
hold on
plot(xlim,[0 0],'k--')
title('CS-US = 0.15')

subplot(3,1,2)
D2 = Experiment.G200PsParamsCorrelations;
plot(ones(size(D2,1),1),D2(:,1),'k*',2*ones(size(D2,1),1),D2(:,2),'k*',...
    3*ones(size(D2,1),1),D2(:,3),'k*',4*ones(size(D2,1),1),D2(:,4),'k*',...
    5*ones(size(D2,1),1),D2(:,5),'k*',6*ones(size(D2,1),1),D2(:,6),'k*')
% set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{'?v?' 'Wv?' 'Mv?' 'Wv?' 'Mv?' 'MvW'})
set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{char([8595 118 8593]) ...
    char([87 118 8593]) char([77 118 8593]) char([87 118 8595]) ...
    char([77 118 8595]) 'MvW'})
xlim([.5 6.5])
hold on
plot(xlim,[0 0],'k--')
title('CS-US = 0.2')
ylabel('Correlation Coefficient')


subplot(3,1,3)
D3 = Experiment.G300PsParamsCorrelations;
plot(ones(size(D3,1),1),D3(:,1),'k*',2*ones(size(D3,1),1),D3(:,2),'k*',...
    3*ones(size(D3,1),1),D3(:,3),'k*',4*ones(size(D3,1),1),D3(:,4),'k*',...
    5*ones(size(D3,1),1),D3(:,5),'k*',6*ones(size(D3,1),1),D3(:,6),'k*')
% set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{'?v?' 'Wv?' 'Mv?' 'Wv?' 'Mv?' 'MvW'})
set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{char([8595 118 8593]) ...
    char([87 118 8593]) char([77 118 8593]) char([87 118 8595]) ...
    char([77 118 8595]) 'MvW'})
xlim([.5 6.5])
hold on
plot(xlim,[0 0],'k--')
title('CS-US = 0.3')
saveas(gcf,'Figures/PauseParameterCorrelations.pdf')
%% CDFs of Quartiles of the Pause Onsets and Offsets
figure 
subplot(3,2,1)
    cdfplot(Experiment.G150Qs1to3_PsOn(:,1))
    hold on
    cdfplot(Experiment.G150Qs1to3_PsOn(:,2))
    cdfplot(Experiment.G150Qs1to3_PsOn(:,3))
    xlim([-.1 .45])
    title('CS->US = 0.15')
    xlabel('')
    ylabel('')
    legend('Q1','Mdn','Q3','location','SE')
subplot(3,2,2)
    cdfplot(Experiment.G150Qs1to3_PsOff(:,1))
    hold on
    cdfplot(Experiment.G150Qs1to3_PsOff(:,2))
    cdfplot(Experiment.G150Qs1to3_PsOff(:,3))
    xlim([-.1 .45])
    xlabel('')
    title('')
    ylabel('')
subplot(3,2,3)
    cdfplot(Experiment.G200Qs1to3_PsOn(:,1))
    hold on
    cdfplot(Experiment.G200Qs1to3_PsOn(:,2))
    cdfplot(Experiment.G200Qs1to3_PsOn(:,3))
    xlim([-.1 .45])
    title('CS->US = 0.20')
    xlabel('')
    ylabel('Cumulative Fraction of Cells')
subplot(3,2,4)
    cdfplot(Experiment.G200Qs1to3_PsOff(:,1))
    hold on
    cdfplot(Experiment.G200Qs1to3_PsOff(:,2))
    cdfplot(Experiment.G200Qs1to3_PsOff(:,3))
    xlim([-.1 .45])
    xlabel('')
    title('')
    ylabel('')
subplot(3,2,5)
    cdfplot(Experiment.G300Qs1to3_PsOn(:,1))
    hold on
    cdfplot(Experiment.G300Qs1to3_PsOn(:,2))
    cdfplot(Experiment.G300Qs1to3_PsOn(:,3))
    xlim([-.1 .45])
    title('CS->US = 0.30')
    xlabel('Pause Onset Latency (s)')
    ylabel('')
subplot(3,2,6)
    cdfplot(Experiment.G300Qs1to3_PsOff(:,1))
    hold on
    cdfplot(Experiment.G300Qs1to3_PsOff(:,2))
    cdfplot(Experiment.G300Qs1to3_PsOff(:,3))
    xlabel('Pause Offset Latency (s)')
    xlim([-.1 .45]) 
    ylabel('')
    title('')
saveas(gcf,'Figures/CDFsPsOn&OffQuartiles.pdf')
%% CDFs of Pause Width Quartiles by Group
figure 
subplot(3,1,1)
    cdfplot(Experiment.G150Qs1to3_PsWidths(:,1))
    hold on
    cdfplot(Experiment.G150Qs1to3_PsWidths(:,2))
    cdfplot(Experiment.G150Qs1to3_PsWidths(:,3))
    xlim([-.1 .5])
    title('CS->US = 0.15')
    xlabel('')
    ylabel('')
    legend('Q1','Median','Q3','location','SE')

subplot(3,1,2)
    cdfplot(Experiment.G200Qs1to3_PsWidths(:,1))
    hold on
    cdfplot(Experiment.G200Qs1to3_PsWidths(:,2))
    cdfplot(Experiment.G200Qs1to3_PsWidths(:,3))
    xlim([-.1 .5])
    title('CS->US = 0.20')
    xlabel('')
    ylabel('Cumulative Fraction of Cells')
subplot(3,1,3)
    cdfplot(Experiment.G300Qs1to3_PsWidths(:,1))
    hold on
    cdfplot(Experiment.G300Qs1to3_PsWidths(:,2))
    cdfplot(Experiment.G300Qs1to3_PsWidths(:,3))
    xlim([-.1 .5])
    title('CS->US = 0.30')
    xlabel('Pause Width (s)')
    ylabel('')


% saveas(gcf,'Figures/CDFsPsWidthsByCSUSgroup.pdf')    
%% Graphing longest interspike-interval & latency to it
figure
subplot(3,2,1)
    cdfplot(Experiment.G150Qs1to3OfMxIspkIs(:,1))
    hold on
    cdfplot(Experiment.G150Qs1to3OfMxIspkIs(:,2))
    cdfplot(Experiment.G150Qs1to3OfMxIspkIs(:,3))
    xlim([0 .6])
    xlabel('')
    ylabel('')
    title('CSUS150')
    legend('Q1','Median','Q3','location','SE')
subplot(3,2,2)
    cdfplot(Experiment.G150Qs1to3OfLatofMxIspkI(:,1))
    hold on
    cdfplot(Experiment.G150Qs1to3OfLatofMxIspkI(:,2))
    cdfplot(Experiment.G150Qs1to3OfLatofMxIspkI(:,3))
    xlabel('')
    ylabel('')
    title('')
subplot(3,2,3)
    cdfplot(Experiment.G200Qs1to3OfMxIspkIs(:,1))
    hold on
    cdfplot(Experiment.G200Qs1to3OfMxIspkIs(:,2))
    cdfplot(Experiment.G200Qs1to3OfMxIspkIs(:,3))
    xlim([0 .6])
    ylabel('Cumulative Fraction of Cells')
    xlabel('')
    title('CSUS200')
subplot(3,2,4)
    cdfplot(Experiment.G200Qs1to3OfLatofMxIspkI(:,1))
    hold on
    cdfplot(Experiment.G200Qs1to3OfLatofMxIspkI(:,2))
    cdfplot(Experiment.G200Qs1to3OfLatofMxIspkI(:,3))
    xlabel('')
    ylabel('')
    title('')
subplot(3,2,5)
    cdfplot(Experiment.G300Qs1to3OfMxIspkIs(:,1))
    hold on
    cdfplot(Experiment.G300Qs1to3OfMxIspkIs(:,2))
    cdfplot(Experiment.G300Qs1to3OfMxIspkIs(:,3))
    xlim([0 .6])
    xlabel('Longest Interspike Interval (s)')
    ylabel('')
    title('CSUS300')
subplot(3,2,6)
    cdfplot(Experiment.G300Qs1to3OfLatofMxIspkI(:,1))
    hold on
    cdfplot(Experiment.G300Qs1to3OfLatofMxIspkI(:,2))
    cdfplot(Experiment.G300Qs1to3OfLatofMxIspkI(:,3))
    xlabel('Latency To Onset (s)')
    ylabel('')
    title('')
    
TSsaveexperiment % conclusion of analysis
%%
mainProg = readtextfile('FredParalleFibFinalAnalysis.m');
Experiment.Script=mainProg;