% script m-file Exper2PawStim.m
% "Exper2PawStim.m" is the file name for this script. This script
% contains all of the code for analyzing the data and generating the
% figures
CD = cd;
if ~exist('TSlib','dir')
    disp(char({'';'The TSlib folder is not on Matlab'' search path.'...
        'Move it and the HelperFunctions folder into the current';...
        'directory and run the script again';''}))
    return
elseif ~exist('HelperFunctions','dir')
    disp(char({'';'The HelperFunctions folder is not on Matlab''s';...
        'search path. Move it into the current directory and';...
        'run the script again';''}))
    return
end
addpath([CD '/TSlib'],[CD '/HelperFunctions'])
addpath([CD '/TSlib/GUI'])
addpath([CD '/TSlib/Support'])
if ~all([exist('ISI 150','file') exist('ISI 200','file') exist('ISI 300','file') ...
        exist('ISI 400','file') exist('ISI 450','file')])
    disp(char({'';'Some or all of the 5 raw data folders are not in the';...
        'current directory. Move them into it and rerun the script';''}))
    return
end
warning('off','all')
%% Cell 1: Creating the experiment structure
% Create the Experiment structure into which everything will go. A structure
% is a kind of "variable" in Matlab. 'Variable' is in scare quotes because
% this "variable" is really a hierarchical data structure into which every
% kind of data and metadata is put. In the Gallistel-King TSlib system,
% this structure is what "keeps it all together". When the experiment has
% been accepted for publication and one wants to archive it, even this
% script can itself be stored in the structure. In doing that, one would
% convert the script to a function and embed all of the helper functions in
% it. This "variable" is always saved as its own file. Thus, to upload to a
% public data base everything about the experiment--the raw data, the
% process-control code, the data-analysis and figure-generation
% table-generation code, and all the metadata, one would upload a single
% file.

% ID numbers for the subjects in the experiment, that is, the hopefully
% unique identification numbers (in effect social security numbers) that
% they were assigned when they entered the lab or perhaps only when they
% entered the experiment
Subs = [ 929.00 936.00 951.00 959.00 961.00 965.00 970.00 971.00 972.00 ...
    975.00 976.00 980.00 984.00 1018.00 1037.00 1041.00 1048.00 1049.00 ...
 1052.00 1058.00 1060.00 1062.00 1069.00 1071.00 1073.00 1075.00 1077.00 ...
 1082.00 1084.00 1087.00 1090.00 1105.00 1109.00 1110.00 1111.00 1113.00 ...
  1115.00 1116.00 1118.00 1119.00 1120.00 1121.00 1123.00 1125.00 1126.00 ...
  1131.00 1135.00 1138.00 1140.00 1149.00 1150.00 1157 1159 1161];

%
TSinitexperiment([CD '/PurkinjeCellPauseExper2'],100,Subs,'electrodes','Hesslow')
Experiment.Info.LoadFunction='LoadFredrik';
% 'LoadFredrik' is the name of the first helper function. A helper function
% helps a principal function do its job.The principal function in this
% case--the function that calls on this helper--is the TSloadsessions
% function, the function that loads session by session raw data into the
% TSData field of a Session field within the hierarchical Experiment
% structure. Whenever the user calls TSloadsessions, it goes to this field
% in the Experiment structure in order to get the name of the helper
% function that understands the structure of the data to be
% loaded. Data structures are usually lab-specific, which makes it
% impossible to create a session-loading function that loads data from no-
% matter-whose lab. For each lab that uses the TSsystem, a load function
% must be written, by someone who understands the structure of the raw data
% files in that lab. Most helper functions are simple; often only 1 to 10
% lines of code, but they can be arbitrarily long and complex. Load
% functions tend to be long, complex and hard to understand by someone
% unfamiliar with the structure of the raw data files in a given lab. We
% provide separately several heavily annotated load functions to give the
% new user examples of what needs to go into a load function.
Experiment.Info.InputTimeUnit = 1;%  The number of seconds in the time unit
%  in the raw data files that are to be loaded (e.g., .01 in data where the
% event clock ticks every 10 ms and the recorded times are the number of
% ticks up to the given event time).

Experiment.Info.OutputTimeUnit=1; % The number of seconds in the time unit 
% that one wants to work with in analyzing the data. If one wanted to work
% with times in minutes, this number would be 60.

Experiment.Info.FileExtension ='.mat'; % The extension by which data files
% can be uniquely recognized. Even folders that apparently contain only
% data files usually have system-generated hidden files. Trying to load one
% of these files will crash the load function.

Experiment.Info.FilePrefix = ''; % There is no initial character of sequence
% of characters that picks out the data files

TSexperimentbrowser % Calls the browser window that enables browsing around
% the ever growing, mulitlevel Experiment structure

if exist('ISI 150','dir') && exist('ISI 200','dir') && exist('ISI 300','dir')...
        && exist('ISI 400','dir') && exist('ISI 450','dir')
    TSloadsessions('ISI 150') % Choose ISI 150 (n=22), marker is 0.2 sec before CS onset
    TSloadsessions('ISI 200') % Choose ISI 200 (n=22), marker is 0.2 sec before CS onset
    TSloadsessions('ISI 300') % Choose ISI 300 (n=75), marker is 0.2 sec before CS onset
    TSloadsessions('ISI 400') % Choose ISI 400 (n=5), marker is 0.2 sec before CS onset
    TSloadsessions('ISI 450') % Choose ISI 450 (n=3), marker is 0.2 sec before CS onset
    % Once sessions have been loaded into the Experiment, the raw data exists
    % in two forms: 1) the file generated by the process-control software at
    % the time the experiment was run; 2) the contents of the TSData field for
    % that Session and that Subject within the hierarchical Experiment
    % structure. Subject and Session are both indexed field headings; thus,
    % there is Subject(1), Subject(2),...,Subject(n); and under Subject(n), 
    % there is Session(1), Session(2), ..., Session(m). The raw data for
    % Session (m) of Subject(n) live in the field TSData under Subject(n) and
    % Session(m) of the hierarchy. All of the data analysis operations apply to
    % the raw data in this TSData field (or in augmented copies of it, which
    % may be created in the course of the analysis). The original raw data need
    % never be accessed again, although copies of them should, of course, be
    % kept in case doubt should ever arise about whether the data in the TSData
    % field are in fact a faithful copy of the original raw data.
else
    disp(char({'';'The TSloadsessions commands in Cell 1assume that the raw';...
    'data are in directories named ISI 150, ISI 200, ISI 300, ISI 400 & ISI 450';...
    'and that those directories are subdirectories of the current directory.';...
    'If the raw data have been put somewhere else, then call these 5 TSsession';...
    'commands without an argument, in which case the command opens a browser';...
    'window for you to find the folder that contains the data to be loaded.';...
    'The commands begin at Line 74 of the script.';''}))
    return
end
% IN THIS ELECTROPYSIOLOGICAL DATA SET, 'SESSION' AND 'SUBJECT' CORRESPOND 
% TO CELL AND SUBJECT, respectively. For most subjects, there was only
% one cell; but for a few, there were as many as three cells.
if exist('eventcodes.txt','file')
    TSimporteventcodes('eventcodes.txt') % importing event codes for mark, spike, CSon, CSoff and
    %  MidITI events. This command brings up a browse window for you to browse
    %  for the text file named eventcodes, which is simple, experiment-specific
    %  text file created by the user. In this electrophyisological experiment,
    %  there are only three events plus two pseudoevents. The events are 'CSon',
    %  'CSoff' and 'spike'. The pseudo-events are 'marker' and 'midITI'. The
    %  marker event was generated at the time the experiment was run; it marks 
    %  point 200ms prior to CS onset. The midITI event marks the point that is
    %  the half way through each intertrial interval. This pseudoevent is inserted
    %  into the TSData field by code found at Line 311. The text file
    %  containing the event codes for this experiment looks like this:
    %{
    marker = 20
    spike = 40
    CSon = 30
    CSoff = 50
    midITI = 60
    %}
    % The TSimporteventcodes function reads this into a field named EvenCodes 
    % atthe Experimentlevel of the hierarchy, the contents of which look like
    % this:
    %{
    marker: 20
    spike: 40
    CSon: 30
    CSoff: 50
    midITI 60
    %}
    % This field serves as a dictionary that translates between the textual
    % event codes understood by humans (e.g., "spike") and the numerical event
    % codes understood by Matlab. Code in the TSsystem refers to events by
    % their readily understood names. When that code is executed, however,
    % Matlab reads the EventCodes field to get the corresponding numbers, the
    % numbers by which those events are identified in the raw data. It is
    % impossible to overstate how much more intelligible the data-analysis code
    % becomes when the events operated on are identified by intelligible
    % names rather than by numbers
else
    disp(char({'';'The TSimportevent codes expected to find a file named eventcodes.txt';...
        'in the current directory. If you call the TSimporteventcodes function w/o';...
        'an argument, it will open a file browser window, allowing you to browse';...
        'for that file. If you cannot find it, see script comment near end of Cell 1';...
        'for what the simple file should contain.';''}))
    return
end
TSsaveexperiment % saves experiment in a file with the name supplied in the
% TSinitexperiment command near the top of this cell
%% Cell 2: Daily START
% cd '/Users/galliste/Google Drive/ForDavid' % Commented out because it is
% specific to my computer. In building a script (analysis), I always have
% this cell to speed up getting started each time I return to the work

% TSloadexperiment('FredExpMarchExmpl2') % commented out so that when whole
% script is run, it does not reload the just-saved structure and does not
% query whether to overwrite Experiment structure already in workspace

TSexperimentbrowser % opens the browser that allows one to roam around the
% experiment structure.

%% Cell 3a: Creating fields that are used to group subjects in accord with
% the experimental condition, that is to say in this case, in accord with
% the CS-US interval used in training these decerebrate ferret subjects

TScombineover('CSdur','Phase'); % raising CS duration info to Subject level
% The load function entered this info into the Phase field at the Session
% level when the raw data were loaded. 'Phase' here is synonymous with
% 'condition'; it specifies the duration of the CS-US interval used in
% training a subject


TScombineover('CSdurByElectrode','CSdur','t') % Creates a field named
% 'CSdurByElectrode' at the Experiment level, which field contains a 2-col
% array, with CS durations in 1st column and electrode in 2nd col. 'Subject'
% and 'electrode' have the same referent here, because there was only one
% penetration per (decerebrate) subject. To create this field,
% TScombineover reads the fields named 'CSdur' at the Subject level

%
TSapplystat({'CS150electrodes' 'CS200electrodes' 'CS300electrodes' 'CS400electrodes' 'CS450electrodes'},...
    'CSdurByElectrode',@CSdurElectrodes)
% This command creates five fields at the Experiment level. The names of
% these five fields are enclosed in single quotes as the first argument of
% the function (before the first comma). The second argument of the
% TSapplystat command ('CSdurByElectrode') is the name of the field that is
% read in order to obtain the necessary information for creating this new
% fields. CSdurElectrodes is another helper function; here is the code it
% contains:
%{
function [r150,r200,r300,r400,r450] = CSdurElectrodes(D)
% creates logical vectors flagging electrodes in a given CS duration group
r150 = unique(D(D(:,1)==150,2))'; % col 2 for all rows with 150 in 1st col
r200 = unique(D(D(:,1)==200,2))';
r300 = unique(D(D(:,1)==300,2))';
r400 = unique(D(D(:,1)==400,2))';
r450 = unique(D(D(:,1)==450,2))';
%}
% The contents of these fields will be referred to whenever an analysis is
% to be limited to the data coming from a particular CS-US training
% interval, because these fields contain the index numbers of the subjects
% trained with a given CS-US interval

%% Cell 3b: Fields specifying the US-US intervals and whether or not the spike
% recorder ran continuously, that is, during the intertrial intervals. This
% info was not contained in the metadata. It was discovered by examining
% raster plots.

Experiment.USUSCycleDur6s=[1:3 17 18 22 24:26 28]; % subjects trained
% with 6s US-US intervals. 

Experiment.USUSCycleDur16s = [4:16 19:21 23 27 29:55]; % subjects trained
% with 16s US-US intervals

Experiment.DataThroughoutCycle = [1 1;9 1;14 1;16 1;17 1;18 1;19 1;21 1;...
    36 1;36 2;36 3;54 1];
    
%% Cell 4: Raster plots for initial examination of data
% Three views are shown of each cell's raster: 1) the overall view
% (complete cycles) 2) the view limited to 0.2 s before CS onset to 2s
% after CS onset; 3) the view limited to 0.2 s before CS onset to 1.2s
% after CS onset. You get the successive views by hitting the space key (or
% any key)
disp(char({'';'The for loops in Cell 4 of the script show 3 successive views';...
    'of the spike raster for successive subjects and successive cells';...
    'within subjects;. The first view is the whole raster; the 2nd zooms';...
    'in for a close look by limiting the x axis extent to [0 6.5s]';...
    'The 3rd view shows only 0 to 1.2s and it shows the onset and';...
    'and offset of the CS (vertical blue lines). You progress from';...
    'view to view by hitting any key (e.g., the space key).';...
    '     After the data for each subject,it asks you whether';...
    'you want to continue to the next subject. Hit ''y'' or';...
    'return to continue; hit ''n'' to terminate the sequence';''}))
pause(10)
N = 1;
while N<Experiment.NumSubjects
    for S = 1:Experiment.NumSubjects % stepping through the specified
        % subjects; 
        for s = 1:Experiment.Subject(S).NumSessions % stepping through the cells
            % for a given subject
            TSlimit('Subjects',S); % data sets to be examined
            TSlimit('Sessions', s)
            CSdur = Experiment.Subject(S).Session(s).Phase; % retrieving
              % duration of the CS
            Xlim1 = [0 6.5]; % first pair of x-axis limits
            Xlim2 = [0 1.2]; % second pair of x-axis limits (closer look)
            TSapplystat('','TSData',@DispRaster,CSdur,Xlim1,Xlim2) % execute this cell
            % to see raster plots at a given CS duration. After examining full plot,
            % hit any key to impose first pair of x-axis limits. After examining that
            % view, hit any key to impose 2nd pair. To then go on to next data set,
            % again hit any key. DispRaster is a helper function and here is
            % its code
            %{
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
            %}
            term=input('Continue to next subject? (y/n) ','s');
            if strcmp('n',term)
                N = 1000;
                TSlimit('Sessions','all')
                break
            end
        end
        if strcmp('n',term)
            TSlimit('Subjects','all')
            break
        end
    end
end

%% Cell 5: Maximum and minimum interspike intervals
%  These enable recognition of those cells where recording was continuous
% vs those where it was turned off during the intertrial interval
TSlimit('Sessions','all') 
TSsessionstat({'MxISpI' 'MinISpI'},@MxMin)
% Creates two fields at the Session/cell level, named 'MxISpI' and 'MinISpI',
% containing, respectively, the maximum and minimum interspike interval.
% MxMin is a helper function, and here is its code:
%{
function [Mx,Min] = MxMin(tsd)
ispkis = diff(tsd(:,1));
Mx = max(ispkis);
Min = min(ispkis);
%}
TScombineover('MxISpkI_El','MxISpI') % raises maximum interspike interval
% from the field at the session level ('MxISpI') to a field at the Subject
% level ('MxISpkI_El')
TScombineover('MxISpkI_Ex','MxISpkI_El','t') % ditto for min inter-spike interval

%% Cell 6: Inserting events into the raw data
% The original raw data recorded only one event (spike) and a pseudoevent
% (marker) marking a point 0.2 s prior to CS onset. The following code uses
% the location of the marker pseudoevent plus metadata knowledge of the CS-US
% interval for a given Subject to insert two real events into the raw data,
% namely CSon and CSoff, plus another pseudoevent, midITI, marking the
% point in the middle of the intertrial interval. It puts the augmented
% "raw" data in a new field TSDataA. It COULD put the augmented raw data
% back into the TSData field, but best practice is to leave that field
% unaltered by any processing
TSapplystat('TSDataA',{'TSData' 'Phase'},@AugTSdata) % creates TSDataA
% (A is for augmented) field at the Session level in which there are CSon,
% CSoff and MidITI events. This enables me to define the appropriate kind
% of trials
%{
function tsdA = AugTSdata(tsd,CSdur)
CSdur = CSdur/1000;
Mtms=tsd(tsd(:,2)==20,1);
CSonTms = Mtms+.2;
CSoffTms = CSonTms+CSdur;
ITIs = CSoffTms(2:end)-CSonTms(1:end-1);
HalfITIs=round(1000*ITIs/2)/1000;
LstPostInt=tsd(end,1)-CSoffTms(end);
MidITItms = CSoffTms+[HalfITIs;LstPostInt];
MidITItms=[tsd(1)-.001;MidITItms];
MidITIs = [MidITItms 60*ones(length(MidITItms),1)];
CSons = [CSonTms 30*ones(length(CSonTms),1)];
CSoffs=[CSoffTms 50*ones(length(CSoffTms),1)];
tsdA = sortrows([tsd;CSons;CSoffs;MidITIs]);
%}

TSsetdata('TSDataA') % makes the newly created augmented raw field the
% active data field
%% Cell 7: Defining three trial types: MidCSon, CS, and PostCS. These trial
% types partition each trial into a pre-CS portion, a CS portion (during
% the CS, which was electrical stimulation of the dorsum of the paw), and a
% post-CS portion, allowing us to compute statistics separately for each
% portion on each trial
TSdefinetrialtype('MidCSon',[midITI CSon]) % defines as a "trial type" the
% interval from the pseudoevent that marks the middle of the intertrial
% interval to the CS onset. When one defines a trial type, it becomes a
% field at the session level. This field will contain indexed trials of
% that type as soon as one makes a call to the TStrialstat function.

TStrialstat('SpkTms',@TSparse,'result=time(1)-starttime;',spike)
% This command creates indexed Trials [Trial(1), Trial(2),...,Trial(n)
% underneath the above defined trial type. n is the number of occurrences
% of the interval specified by the TSdefinetrialtype command, in this case,
% the interval from the midITI marker to CS onset. Each of the indexed
% trials has a number of fields that are automatically created by this
% command the first time it is called for a given trial type, to wit:
% 'StartTime', 'EndTime', 'TrialDuration','sloc'(line number in raw data
% where trial starts),'eloc' (line number where it ends); and it then
% creates a field with the name specified by the user in the call to this
% function, in this case 'SpkTmns'. And it gives a bit of code that
% computes the statistic that goes in this field, which, in this case, is
% the times within each trial--the times referenced to a 0 at the intertrial
% interval marker for that trial--at which spikes occur
%
TSapplystat('NumSpks','SpkTms',@numel) % Creates a field at the Trial leve
% containing the # of spikes

TSapplystat('Rate',{'SpkTms' 'TrialDuration'},@spkspersec)% trial level
% Creates a field at the trial level containing the rate of firing 
%{
function r=spkspersec(spktms,td)
N = numel(spktms);
D=td-spktms(1); % rate computed only once the 1st spike in the trial occurs
% because that was when the recorder was turned on prior to CS delivery
if N==0
    r=0;
else
    r=N/D;
end
%}
TSapplystat('IspkIs','SpkTms',@diff) % Creates a field at the trial level
% containing the interspike intervals
TSapplystat('MinIspkI','IspkIs',@min) % field containing the minimum ispki
TSapplystat('MaxIspkI','IspkIs',@max) % field containing the max ispki
TSapplystat('IspkIquantiles','IspkIs',@quantile,...
    [.125 .25 .375 .5 .625 .75 .875]) % field containing quantiles of the
% ispki distribution
TScombineover('PreCSrates','Rate') % Raising the statistics on pre-CS firing
% rates to a field at the Cell (Session) level
TScombineover('PreCSispis','IspkIs') % ditto for the interspike intervals
TScombineover('PreCSminIspkIs','MinIspkI') % ditto for min ispki
TScombineover('PreCSmaxIspkIs','MaxIspkI') % ditto for max ispki
TScombineover('PreCSispkiQuants','IspkIquantiles') % ditto for quantiles

%
TSdefinetrialtype('CS',[CSon CSoff]) % defines a 2nd trial type. Defining
% a new trial type makes it the active trial type.  TStrialstat commands
% always operate only on whatever is the active trial type

TStrialstat('SpkTms',@TSparse,'result=time(1)-starttime;',spike) % see above
% for 1st trial type

% Same as for 1st trial type
TSapplystat('NumSpks','SpkTms',@numel) % trial level
TSapplystat('Rate',{'NumSpks' 'TrialDuration'},@rdivide)% trial level
TSapplystat('IspkIs','SpkTms',@diff) % trial level
TSapplystat('MinIspkI','IspkIs',@min) % trial level
TSapplystat('MaxIspkI','IspkIs',@max) % trial level
TSapplystat('IspkIquantiles','IspkIs',@quantile,...
    [.125 .25 .375 .5 .625 .75 .875]) % trial level

% Same as for 1st trial type: raising results to the Cell level, with
% separate fields at that level depending on which trial type the
% statistics come from
TScombineover('CSrates','Rate') % cell level firing rates during CS
TScombineover('CSispis','IspkIs') % ditto for ispki's
TScombineover('CSminIspkIs','MinIspkI') % ditto for min ispki
TScombineover('CSmaxIspkIs','MaxIspkI') % ditto for max ispki
TScombineover('CSispkiQuants','IspkIquantiles') % ditto for quantiles


% Same statistics for a third trial type
TSdefinetrialtype('PostCS',[CSoff midITI]) % defining a 3rd trial type
TStrialstat('SpkTms',@TSparse,'result=time(1)-starttime;',spike) % see above
TSapplystat('NumSpks','SpkTms',@numel) % trial level
TSapplystat('Rate','SpkTms',@postspkspersec)% trial level
%{
function r=postspkspersec(spktms)
N = numel(spktms);
D=spktms(end);
if N==0
    r=0;
else
    r=N/D;
end
%}

TSapplystat('IspkIs','SpkTms',@diff) % trial level
TSapplystat('MinIspkI','IspkIs',@min) % trial level
TSapplystat('MaxIspkI','IspkIs',@max) % trial level
TSapplystat('IspkIquantiles','IspkIs',@quantile,...
    [.125 .25 .375 .5 .625 .75 .875])
TScombineover('PostCSrates','Rate') % cell level
TScombineover('PostCSispis','IspkIs') % cell level
TScombineover('PostCSminIspkIs','MinIspkI') % cell level
TScombineover('PostCSmaxIspkIs','MaxIspkI') % cell level
TScombineover('PostCSispkiQuants','IspkIquantiles') % cell level


%% Cell 8a: Some Basic Graphs (these were not used in publication; they were just 
% to get "the lay of the land" and spot weird/suspicious data. For each
% cell, they plot the cumulative distributions of the firing rates during
% the pre-CS intervals (black curves), during the CS (red) and during the
% post-CS interval (green). Thus one can see at a glance whether the firing
% rate during the CS was less than during the pre- and post-CS intervals
% (red curve will lie to left of black & green) and whether the firing rate
% differed in the post-CS interval from what it was in the pre-CS interval.
% This cell both creates and saves the graphs, with separate figures (or
% sequences of figures) for each training condition (each CS duration). The
% figures are saved in a directory it creates named
% 'CDFsPreDurPstFiringRate' inside the Figures subdirectory
mkdir('Figures')

disp(char({'';'Creating figures with multiple panels. Each panel shows';...
    'the cumulative distribution functions for the pre- during- & post-CS';...
    'firing rates for one cell (black, red & green cruves, respectively)';...
    'The figures are saved in a directory named CDFsPreDurPstFiringRate';''}))
pause(5)
mkdir('Figures/CDFsPreDurPstFiringRate')
TSlimit('Subjects',Experiment.CS150electrodes)
TSlimit('Sessions','all')
TSapplystat('',{'PreCSrates' 'CSrates' 'PostCSrates'},@TSplotcdfs,'Rows',6,...
    'Xlbl','Spikes/s','Handle','H15') % cumulative distributions of the pre-CS rates in
% subjects trained with 0.15s CS-US interval. Multipanel graph, with each
% panel for a different subject within that condition
L = length(H15); % # of figures to be saved
for f = 1:L % stepping through the just created figure(s)
    set(H15(f),'Name',sprintf('CSdur=.15s; fig %d of %d',f,L))
    subplot(6,2,1)
    legend('Pre','Dur','Pst','location','SE')
    saveas(gcf,sprintf('Figures/CDFsPreDurPstFiringRate/CS150fig%dof%d',f,L))
end
pause(5)
close all
%% Cell 8b
TSlimit('Subjects',Experiment.CS200electrodes)
TSlimit('Sessions','all')
TSapplystat('',{'PreCSrates' 'CSrates' 'PostCSrates'},@TSplotcdfs,'Rows',6,...
    'Xlbl','Spikes/s','Handle','H20') % cumulative distributions of the pre-CS rates in
% subjects trained with 0.2s CS-US interval. Multipanel graph, with each
% panel for a different subject within that condition
L = length(H20); % # of figures to be saved
for f = 1:L % stepping through the just created figure(s)
    figure(H20(f))
    set(gcf,'Name',sprintf('CSdur=.2s; fig %d of %d',f,L))
    subplot(6,2,1)
    legend('Pre','Dur','Pst','location','SE')
    saveas(gcf,sprintf('Figures/CDFsPreDurPstFiringRate/CS200fig%dof%d',f,L))
end
pause(15)
close all
%% Cell 8c
TSlimit('Subjects',Experiment.CS300electrodes)
TSlimit('Sessions','all')
TSapplystat('',{'PreCSrates' 'CSrates' 'PostCSrates'},@TSplotcdfs,'Rows',6,...
    'Xlbl','Spikes/s','Handle','H30') % cumulative distributions of the pre-CS rates in
% subjects trained with 0.3s CS-US interval. Multipanel graph, with each
% panel for a different subject within that condition
L = length(H30); % # of figures to be saved
for f = 1:L % stepping through the just created figure(s)
    figure(H30(f))
    set(gcf,'Name',sprintf('CSdur=.3s; fig %d of %d',f,L))
    subplot(6,2,1)
    legend('Pre','Dur','Pst','location','SE')
    saveas(gcf,sprintf('Figures/CDFsPreDurPstFiringRate/CS300fig%dof%d',f,L))
end
pause(60)
close all
%% Cell 8d
TSlimit('Subjects',Experiment.CS400electrodes)
TSlimit('Sessions','all')
TSapplystat('',{'PreCSrates' 'CSrates' 'PostCSrates'},@TSplotcdfs,'Rows',6,...
    'Xlbl','Spikes/s','Handle','H40') % cumulative distributions of the pre-CS rates in
% subjects trained with 0.4s CS-US interval. Multipanel graph, with each
% panel for a different subject within that condition
L = length(H40); % # of figures to be saved
for f = 1:L % stepping through the just created figure(s)
    figure(H40(f))
    set(gcf,'Name',sprintf('CSdur=.4s; fig %d of %d',f,L))
    subplot(6,2,1)
    legend('Pre','Dur','Pst','location','SE')
    saveas(gcf,sprintf('Figures/CDFsPreDurPstFiringRate/CS400fig%dof%d',f,L))
end
pause(5)
close all
%% Cell 8e
TSlimit('Subjects',Experiment.CS450electrodes)
TSlimit('Sessions','all')
TSapplystat('',{'PreCSrates' 'CSrates' 'PostCSrates'},@TSplotcdfs,'Rows',6,...
    'Xlbl','Spikes/s','Handle','H45') % cumulative distributions of the pre-CS rates in
% subjects trained with 0.45s CS-US interval. Multipanel graph, with each
% panel for a different subject within that condition
L = length(H45); % # of figures to be saved
for f = 1:L % stepping through the just created figure(s)
    figure(H45(f))
    set(gcf,'Name',sprintf('CSdur=.45s; fig %d of %d',f,L))
    subplot(6,2,1)
    legend('Pre','Dur','Pst','location','SE')
    saveas(gcf,sprintf('Figures/CDFsPreDurPstFiringRate/CS450fig%dof%d',f,L))
end
pause(5)
close all

%% Cell 9: Differences between pre-CS and CS firing rates 
TSlimit('Subjects','all')
TSlimit('Sessions','all')
TSapplystat('CSratesMinusPreCSrates',{'CSrates' 'PreCSrates'},@minus)
TSapplystat('CSratesMinusPostCSrates',{'CSrates' 'PostCSrates'},@minus)

%% Cell 10: Probability distributions for preCS, CS, and postCS interspike
% intervals
TSlimit('Subjects','all')
Edges = [0 .005 .01 .015 .02 .025 .03 .04 inf]; 
TSapplystat({'PreCSispkiPDF' 'CSispkiPDF' 'PostCSispkiPDF' 'Edges'},...
    'TSDataA',@PrbDst,Edges,midITI,CSon,CSoff)
%{
function [Pre,CS,Post,E]=PrbDst(tsd,Edges,midITI,CSon,CSoff)
% computes three empirical discrete probability distribution functions, one 
% from the preCS data, one from the CS data, and one from the postCS data, 
% using the vector Edges to define the intervals
%%
mid = find(tsd(:,2)==midITI); % rows where trials notionally begin & end
CSb=find(tsd(:,2)==CSon); % CSon rows
CSe=find(tsd(:,2)==CSoff); % CSoff rows
PreIspI = [];
CSispi = [];
PostIspI = [];
%
for t=1:length(CSb) % stepping through the 20 trials
    ispi1 = diff(tsd(mid(t)+1:CSb(t)-1,1)); % pre ispki's
    ispi2 = diff(tsd(CSe(t):mid(t+1)-1,1)); % post ispki's
    ispi3 = diff(tsd(CSb(t):CSe(t),1));
    PreIspI(end+1:end+length(ispi1),1)=ispi1;
    CSispi(end+1:end+length(ispi3),1) = ispi3;
    PostIspI(end+1:end+length(ispi2),1)=ispi2;
end
Pre = histc(PreIspI,Edges);
CS = histc(CSispi,Edges);
Post = histc(PostIspI,Edges);
Pre = Pre/sum(Pre); % normalized
CS = CS/sum(CS); % normalized
Post = Post/sum(Post); % normalized
E=Edges;   
%}    
%% Cell 11: Histograms of interspike intervals for the 12 cells w
% continuous recording
CSs = Experiment.DataThroughoutCycle; % 2-col array, with subject in 1st
% col and cell (session) in 2nd
E=logspace(-3,-.4); % bin edges
plt=1;
sbplt=1;
for c = 1:length(CSs)
    D = diff(Experiment.Subject(CSs(c,1)).Session(CSs(c,2)).TSDataA(:,1));
    % interspike intervals (very slightly contaminated by non-spike events)
    if mod(plt,8)==1
        figure
        sbplt = 1;
    end
    subplot(4,2,sbplt)
    hist(D,E)
    xlim([0 .2])
    Ylm = ylim;
    text(.04,.8*(Ylm(2)),['E' num2str(CSs(c,1)) ', cell' num2str(CSs(c,2))])
    if sbplt>6
        xlabel('ISpkI (s, loged bin widths)','FontSize',14)
    end
    if mod(sbplt,2)==1
        ylabel('Spike Counts','FontSize',14)
    end
    plt=plt+1;
    sbplt=sbplt+1;    
end
disp(char({'';'There were 12 cells for which the recorder was not turned off';...
    'during the intertrial interval. These provide the best data in re';...
    'the distribution of those intervals. These two figures show the';...
    'histograms of the interspike intervals for each of those 12 cells';...
    'They are saved in a folder named Interpsike Int Histograms';''}))
pause(10)
mkdir('Figures/Interpsike Int Histograms')
for f=1:2
    saveas(gcf,sprintf('Figures/Interpsike Int Histograms/IspkIhistogramsFig%dof2',f))
    close
end
%% Cell 12a: Computing & Graphing Fano Factors
% Spike Counts in 1 s start and end intervals PreCS
TSlimit('Subjects','all')
TSlimit('Sessions','all')
TSsettrialtype('MidCSon')
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
TSsettrialtype('PostCS')
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
    'fact that most factors were greater than the upper limit implies';...
    'that the generative process was non-stationary';''}))

%% Cell 12b: Checking exponential fits for the few subjects whose Fano Factors
% are within the Poisson range
% Following Subjects & Sessions have Fano factors within the
% Poisson range
%{
 Pre
S  s
1  1
8  2
21 1
 Post
1  1
3  1
8  2
21 1
33 1
34 6
%}
% Examining fits of exponetial for those few subjects w Fano Factors in
% Poisson range
figure
subplot(3,1,1)
D = [Experiment.Subject(1).Session(1).PreCSispis;...
    Experiment.Subject(1).Session(1).PostCSispis]-.004;
D(D<0)=[];
S1mu = expfit(D);
cdfplot(D)
x = linspace(0,.07);
hold on;plot(x,expcdf(x,S1mu),'k--')
xlim([0 .1])
title('S1,Cell 1')
ylabel('')
xlabel('')

subplot(3,1,2)
D = [Experiment.Subject(8).Session(2).PreCSispis;...
    Experiment.Subject(8).Session(2).PostCSispis]-.005;
D(D<0)=[];
S8mu = expfit(D);
cdfplot(D)
x = linspace(0,.07);
hold on;plot(x,expcdf(x,S8mu),'k--')
xlim([0 .1])
title('S8,Cell 2')
ylabel('Cumulative Fraction of Inter-Spike Intervals')
xlabel('')

subplot(3,1,3)
D = [Experiment.Subject(21).Session(1).PreCSispis;...
    Experiment.Subject(21).Session(1).PostCSispis]-.007;
D(D<0)=[];
S21mu = expfit(D);
cdfplot(D)
x = linspace(0,.07);
hold on;plot(x,expcdf(x,S21mu),'k--')
xlim([0 .1])
title('21,Cell 1')
xlabel('Inter-Spike Interval (s)')
ylabel('')
saveas(gcf,'Figures/ExpoFits')
% Fits are not good. Empirical cdf rises too sharply and has too long a tail 
disp(char({'';'For those few subjects with plausible Poisson Fano Factors,';...
    'this figure compares the empirical distributions of the interspike';...
    'intervals (solid curves) to the best-fitting exponentials (dashed';...
    'curves). The fits are poor; the empirical distributions rise more';...
    'abruptly than the best-fit exponentials and have longer tails.';''}))    
pause(10)
close all

%% Cell 13: Binarizing the data in preparation for using binary CP code to get pause stats
TSlimit('Subjects','all');
TSlimit('Sessions','all');
TSdefinetrialtype('MdToMd',[midITI midITI])
BW = .001; % bin width for binary spike vector
% 
TStrialstatMult({'BinSpkVec' 'FrstSpkTm' 'LastSpkTm'},@Bindata2,BW,spike,CSon)
% creates a field named 'BinSpkVec' in the MdToMd trials whose first column
% is spike time referenced to CS onset at 0 and whose 2nd col is the 
% binary vector representation of the spike train. Also fields giving the
% first and last spike times. Below is the code for the helper function
%{
function [spikes,Mn,Mx] = Bindata2(tsd,BW,spike,CSon)
% A Benoulli vector digitization of the spike train with bins of  width BW
% Syntax     spikes = Bindata2(tsd,BW,spike,CSon)
LV = tsd(:,2)==spike; % flags spikes
LVCSon = tsd(:,2)==CSon; % flags CS on
spks = tsd(LV,1)-tsd(LVCSon,1); % spike times referenced to 0 at CS on
Mn = min(spks);
Mx = max(spks);
bins = Mn:BW:Mx; % creating the bin edges
BV = histc(spks,bins);
t = bins'+BW/2;
spikes = [t BV];
%}

%  Carrying binary spike vectors up to the Cell (Session level) level
TScombineover('BinSpkVecS','BinSpkVec','t')
% field at the Session level with all the binary vectors with a 3rd col
% with integers specifying which trial

%% Cell 14a: Computing pause stats. This is THE KEY CELL
TSlimit('Subjects',Experiment.CS200electrodes);
TSlimit('Sessions','all')
CSdur = .2; % CS duration
TSapplystat({'PsOn' 'WpsOn' 'PsOff' 'WpsOff' 'lmPre' 'lmDur' 'lmAft' ...
    'LngstISpI' 'LatToLngstISpI'},'BinSpkVecS',@BinPsOnOff2,CSdur)
% Takes the data in the BinSpkVecS  at Session level along with the CSdur and
% creates 9 fields at the Session level: PsOn = vector of pause onset
% estimates; WpsOn = vector of the weights of the evidence for those
% onsets; PsOff = vector of estimates of the pause off times; WpsOff =
% vector of the weights of the evidence for those offsets; lmPre = vector
% of the rates of pressing in the interval of duration dur preceding CS
% onset; lmDur = vector of the rates in the middle 40% of the CS; lmAft =
% vector of the rates in the interval of duration dur after the offset of
% the CS; LngstISpI = longest interspike interval within the pause;
% LatToLngstISpI = interval from CSon to PsOn. Because this helper function
% is a critical part of the analysis that would not be easily reconstructed,
% its code is reproduced here. At one point, it calls yet another custom
% function, BernCP_local. That custom function is, however, embedded in
% this function, so all the key code is here and cannot become separated
% from this script
%{
function [On,Won,Off,Woff,lb,ld,la,MxIpI,LatMx] = BinPsOnOff2(D,dur)
% Computes the onsets and offsets of the pauses
D(D(:,2)>1,2)=1; % This corrects a single bin with a count of 2 in the
% BinSpkVecS field for S10 in the parallel fiber data. That triggers a
% crash of this function when BernCP is called
r = find(D(:,1)>-.0005 & D(:,1)<.0005); % vector of row numbers for CS onset
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
NbinsBack = floor(dur*1000); % number of bins prior to CS onset for pause On
NbinsFwd = ceil(edges(loc+1)*1000); % positions start of backward search for CS
% onset at forward edge of minimum bin in peristimulus histogram
NbinsToStrt = floor(edges(loc)*1000); % positions start of forward search for pause
% offset at backward edge of minimum bin in peristimulus histogram
NbinsPst = ceil(3000*dur); % number of bins to end of data for pause Off
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
    if r(tt)+NbinsPst<length(D) 
        Doff = D(r(tt)+NbinsToStrt:r(tt)+NbinsPst,2); % data for finding 
        % pause offset run from backward edge of minimum bin to three
        % times the duration of the CS
    else % don't go past end of data
        Doff = D(r(tt)+NbinsToStrt:end,2);
    end
    
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
    On(tt) = edges(loc+1)-CPon/1000; % the estimated pause onset time in
    % seconds relative to the CS onset time (remember CP is backwards from
    % forward edge of minimum bin
    Won(tt) = log10(OddsD); % the weight of the evidence for the onset
    Off(tt) = edges(loc)+CPoff/1000; % the estimated pause offset time in
    % seconds relative to a 0 at CS onset
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
    display({'Not enough input arguments;';'6 required: ';...
        'Vector of binary data';...
        'Hyperparameter alpha for prior on p before change';...
        'Hyperparameter beta for prior on p before change';...
        'Hyperparameter alpha for prior on p after change';...
        'Hyperparameter beta for prior on p after change';...
        'Probability of a change after any one trial'})
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
%% Cell 14b Repeating above computation of pause parameters for subjects
% in other CS-US duration conditions. Notice how compact this code is. 
% That is the huge advantage of modularization
TSlimit('Subjects',Experiment.CS150electrodes) % Experiment.CS150electrodes
CSdur=.15;
TSapplystat({'PsOn' 'WpsOn' 'PsOff' 'WpsOff' 'lmPre' 'lmDur' 'lmAft' ...
    'LngstISpI' 'LatToLngstISpI'},'BinSpkVecS',@BinPsOnOff2,CSdur)

TSlimit('Subjects',Experiment.CS300electrodes) % 
TSlimit('Sessions','all')% 
CSdur=.3;
TSapplystat({'PsOn' 'WpsOn' 'PsOff' 'WpsOff' 'lmPre' 'lmDur' 'lmAft' ...
    'LngstISpI' 'LatToLngstISpI'},'BinSpkVecS',@BinPsOnOff2,CSdur)

TSlimit('Subjects',Experiment.CS400electrodes) %
CSdur=.4;
TSsettrialtype('MdToMd')
TSapplystat({'PsOn' 'WpsOn' 'PsOff' 'WpsOff' 'lmPre' 'lmDur' 'lmAft' ...
    'LngstISpI' 'LatToLngstISpI'},'BinSpkVecS',@BinPsOnOff2,CSdur)

TSlimit('Subjects',Experiment.CS450electrodes)
CSdur=.45;
TSapplystat({'PsOn' 'WpsOn' 'PsOff' 'WpsOff' 'lmPre' 'lmDur' 'lmAft' ...
    'LngstISpI' 'LatToLngstISpI'},'BinSpkVecS',@BinPsOnOff2,CSdur)
TSsaveexperiment
%% Cell 15: Raster plots with Probability Distributions (PDFs)
% These complex plots enable us to examine the extent to which the above
% algorithm has correctly identified pause onsets and offsets and, at the
% same time, to see the resulting statistics regarding the distribution of
% interspike intervals pre-CS, during the CS, and post-CS. All 100+ of
% these plots will go into the supplementary material for the published
% paper, so that skeptical readers can satisfy themselves that our
% algorithm reliably identified pause onsets and offsets. We do not run
% this as a for loop, because doing so generates so many plots that it
% brings Matlab to its knees. Notice, however, that TSapplystat is
% called twice, each time with a different helper function that does the
% actual plotting. One helper function makes the raster plot; the other
% makes the probability distribution plots to the right of the raster plot
mkdir('Figures/RastersWthPsOn&OFF&IspkIpdfs')
flds=fieldnames(Experiment);
Eflds=flds(15:19);
Xlim3 = [0 1.1]; % the x-axis limits for the raster plots
for C = 1:length(Eflds) % stepping through the CS-US interval conditions
    for S = Experiment.(Eflds{C})
        TSlimit('Subjects',S)
        for s = 1:Experiment.Subject(S).NumSessions % stepping through the
            % cells in a given subject
            figure;
            set(gcf,'Name',[num2str(Eflds{C}(1:5)) ', Sub ' num2str(S) ', Cell ' num2str(s)])
            Ax1 = subplot('Position',[.1 .1 .4 .8]);
            Ax2 = subplot('Position',[.6 .675 .3 .225]);
            Ax3 = subplot('Position',[.6 .3825 .3 .225]);
            Ax4 = subplot('Position',[.6 .1 .3 .225]);
            
            TSlimit('Sessions',s)
            TSapplystat('',{'TSDataA' 'Phase' 'PsOn' 'PsOff'},@DispRaster3,...
                Xlim3,Ax1) % raster plot in left subplot
            %{
            function DispRaster3(tsdat,CSdur,Bs,Es,Xlm1,Ax)
            % displays rasters for a specified CS duration with Xlm1 and CS delimited
            % and pause beginnings and endings marked. Bs is the 'PsOn' field
            % and Es the 'PsOff' field at the Session level
            CSdur = CSdur/1000;
            % figure
            TSraster(tsdat,{[20 20] [20 inf]},[20 0;40 0],Ax,['r+';'k.'])
            %     E = evalin('caller','sub');
            %     C = evalin('caller','ses');
            %     title(['Elect Indx ' num2str(E) '; Cell Indx ' num2str(C)])
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
            title([num2str(Eflds{C}(1:5)) ', Sub ' num2str(S) ', Cell ' num2str(s)])
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

            axes(Hv(3))
            h=bar(postD,1);
            ylim([0 .7])
            set(h,'FaceColor',[1 1 1],'LineWidth',2)
            ylabel('Probability','FontSize',12)
            set(gca,'XTick',.5:1:8.5,'XTickLabel',{'0' '5' '10' '15' '20' '25' '30' '40' '>40'})
            xlabel('InterspikeInterval (ms)','FontSize',12)
            text(.2,.6,'PostCS') 
            %} 
            saveas(gcf,['Figures/RastersWthPsOn&OFF&IspkIpdfs/' get(gcf,'Name')]) 
        end % of plotting a figure for one cell showing raster plot on left and the
        % three probability distributions on the right
          pause(5) % to examine figures for a given subject
          close all
    end % one subject   
end % one CS-US condition
    
%% Cell 16: Computing Pause Widths
TSlimit('Subjects','all')
TSlimit('Sessions','all')
TSapplystat('PsWidth',{'PsOff' 'PsOn'},@minus)
 
%% Cell 17: Means, std's and coefficients of variation of pause onsets,
% offsets and widths
TSapplystat('PsOnMeanStdCoV','PsOn',@meanstdcov)
%{
function O = meanstdcov(D)
O = [nanmean(D) nanstd(D) nanstd(D)/nanmean(D)];
%}
TSapplystat('PsOffMeanStdCoV','PsOff',@meanstdcov)
TSapplystat('PsWidthMeanStdCoV','PsWidth',@meanstdcov)

%% Cell 18: Pause parameter correlations
TSapplystat('PsParamCorrelations',{'PsOn' 'PsOff' 'PsWidth' 'LngstISpI'},@PsParamCorr)
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
%
TSapplystat('PsParamsRowVec','PsParamCorrelations',@PProwvec)
% Reduces preceding 4x4 array to a row vector of the 6 correlations
%{
function RV = PProwvec(A)
RV = [A(1,2) A(1,3) A(1,4) A(2,3) A(2,4) A(3,4)];
%}
%% Cell 19: Quartiles of the Pause Stat Distributions for individual cells
TSlimit('Subjects','all')
TSlimit('Sessions','all')
QV = [.25 .5 .75]; % vector of desired quantiles
TSapplystat('Qs1to3_PsOn','PsOn',@quantile,QV) % quantile is the helper
% function (native to Matlab); the vector of desired quantiles is passed to
% it along with the data from the PsOn field
TSapplystat('Qs1to3_PsOff','PsOff',@quantile,QV)
TSapplystat('Qs1to3_lmPre','lmPre',@quantile,QV)
TSapplystat('Qs1to3_lmDur','lmDur',@quantile,QV)
TSapplystat('Qs1to3_lmAft','lmAft',@quantile,QV)
TSapplystat('Qs1to3_WpsOn','WpsOn',@quantile,QV)
TSapplystat('Qs1to3_WpsOff','WpsOff',@quantile,QV)
TSapplystat('Qs1to3_PsWidth','PsWidth',@quantile,QV)
TSapplystat('Qs1to3_LngstIspI','LngstISpI',@quantile,QV)
TSapplystat('Qs1to3_LatToLngstIspI','LatToLngstISpI',@quantile,QV)

%% Cell 20: Carrying pause parameters up to Subject level
TScombineover('PsOn_S','PsOn','t')
TScombineover('PsOff_S','PsOff','t')
TScombineover('WpsOn_S','WpsOn','t')
TScombineover('WpsOff_S','WpsOff','t')
TScombineover('lmPre_S','lmPre','t')
TScombineover('lmDur_S','lmDur','t')
TScombineover('lmAft_S','lmAft','t')
TScombineover('Qs1to3_PsOn_S','Qs1to3_PsOn','t')
TScombineover('Qs1to3_PsOff_S','Qs1to3_PsOff','t')
TScombineover('Qs1to3_lmPre_S','Qs1to3_lmPre','t')
TScombineover('Qs1to3_lmDur_S','Qs1to3_lmDur','t')
TScombineover('Qs1to3_lmAft_S','Qs1to3_lmAft','t')
TScombineover('Qs1to3_WpsOn_S','Qs1to3_WpsOn','t')
TScombineover('Qs1to3_WpsOff_S','Qs1to3_WpsOff','t')
TScombineover('Qs1to3_PsWidth_S','Qs1to3_PsWidth')
TScombineover('Qs1to3_LngstIspI_S','Qs1to3_LngstIspI')
TScombineover('Qs1to3_LatToLngstIspI_S','Qs1to3_LatToLngstIspI')
TScombineover('PsOnMeanStdCoV_S','PsOnMeanStdCoV')
TScombineover('PsOffMeanStdCoV_S','PsOffMeanStdCoV')
TScombineover('VecsOfPsParamCorrs','PsParamsRowVec')
TScombineover('PsWidthMeanStdCoV_S','PsWidthMeanStdCoV')

%% Cell 21: Carrying pause parameters up to Group fields at Experiment level
TSlimit('Subjects',Experiment.CS150electrodes)
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
TScombineover('G150Qs1to3_WpsOff','Qs1to3_WpsOff_S','t')
TScombineover('G150Qs1to3_PsWidth','Qs1to3_PsWidth_S','t')
TScombineover('G150Qs1to3_LngstIspI','Qs1to3_LngstIspI_S','t')
TScombineover('G150Qs1to3_LaTtoLngstIspI','Qs1to3_LatToLngstIspI_S','t')
TScombineover('G150PsOnMeanStdCoV','PsOnMeanStdCoV_S','t')
TScombineover('G150PsOffMeanStdCoV','PsOffMeanStdCoV_S','t')
TScombineover('G150PsParamsCorrelations','VecsOfPsParamCorrs')
TScombineover('G150PsWidthMeanStdCoV','PsWidthMeanStdCoV_S')
% 
TSlimit('Subjects',Experiment.CS200electrodes)
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
TScombineover('G200Qs1to3_WpsOff','Qs1to3_WpsOff_S','t')
TScombineover('G200Qs1to3_PsWidth','Qs1to3_PsWidth_S','t')
TScombineover('G200Qs1to3_LngstIspI','Qs1to3_LngstIspI_S','t')
TScombineover('G200Qs1to3_LaTtoLngstIspI','Qs1to3_LatToLngstIspI_S','t')
TScombineover('G200PsOnMeanStdCoV','PsOnMeanStdCoV_S','t')
TScombineover('G200PsOffMeanStdCoV','PsOffMeanStdCoV_S','t')
TScombineover('G200PsParamsCorrelations','VecsOfPsParamCorrs')
TScombineover('G200PsWidthMeanStdCoV','PsWidthMeanStdCoV_S')
% 
TSlimit('Subjects',Experiment.CS300electrodes)
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
TScombineover('G300Qs1to3_WpsOff','Qs1to3_WpsOff_S','t')
TScombineover('G300Qs1to3_PsWidth','Qs1to3_PsWidth_S','t')
TScombineover('G300Qs1to3_LngstIspI','Qs1to3_LngstIspI_S','t')
TScombineover('G300Qs1to3_LaTtoLngstIspI','Qs1to3_LatToLngstIspI_S','t')
TScombineover('G300PsOnMeanStdCoV','PsOnMeanStdCoV_S','t')
TScombineover('G300PsOffMeanStdCoV','PsOffMeanStdCoV_S','t')
TScombineover('G300PsParamsCorrelations','VecsOfPsParamCorrs')
TScombineover('G300PsWidthMeanStdCoV','PsWidthMeanStdCoV_S')
% 
TSlimit('Subjects',Experiment.CS400electrodes)
TScombineover('G400PsOn','PsOn_S','t')
TScombineover('G400PsOff','PsOff_S','t')
TScombineover('G400WpsOn','WpsOn_S','t')
TScombineover('G400WpsOff','WpsOff_S','t')
TScombineover('G400lmPre','lmPre_S','t')
TScombineover('G400lmDur','lmDur_S','t')
TScombineover('G400lmAft','lmAft_S','t')
TScombineover('G400Qs1to3_PsOn','Qs1to3_PsOn_S','t')
TScombineover('G400Qs1to3_PsOff','Qs1to3_PsOff_S','t')
TScombineover('G400Qs1to3_lmPre','Qs1to3_lmPre_S','t')
TScombineover('G400Qs1to3_lmDur','Qs1to3_lmDur_S','t')
TScombineover('G400Qs1to3_lmAft','Qs1to3_lmAft_S','t')
TScombineover('G400Qs1to3_WpsOn','Qs1to3_WpsOn_S','t')
TScombineover('G400Qs1to3_WpsOff','Qs1to3_WpsOff_S','t')
TScombineover('G400Qs1to3_PsWidth','Qs1to3_PsWidth_S','t')
TScombineover('G400Qs1to3_LngstIspI','Qs1to3_LngstIspI_S','t')
TScombineover('G400Qs1to3_LaTtoLngstIspI','Qs1to3_LatToLngstIspI_S','t')
TScombineover('G400PsOnMeanStdCoV','PsOnMeanStdCoV_S','t')
TScombineover('G400PsOffMeanStdCoV','PsOffMeanStdCoV_S','t')
TScombineover('G400PsParamsCorrelations','VecsOfPsParamCorrs')
TScombineover('G400PsWidthMeanStdCoV','PsWidthMeanStdCoV_S')
%
TSlimit('Subjects',Experiment.CS450electrodes)
TScombineover('G450PsOn','PsOn_S','t')
TScombineover('G450PsOff','PsOff_S','t')
TScombineover('G450WpsOn','WpsOn_S','t')
TScombineover('G450WpsOff','WpsOff_S','t')
TScombineover('G450lmPre','lmPre_S','t')
TScombineover('G450lmDur','lmDur_S','t')
TScombineover('G450lmAft','lmAft_S','t')
TScombineover('G450Qs1to3_PsOn','Qs1to3_PsOn_S','t')
TScombineover('G450Qs1to3_PsOff','Qs1to3_PsOff_S','t')
TScombineover('G450Qs1to3_lmPre','Qs1to3_lmPre_S','t')
TScombineover('G450Qs1to3_lmDur','Qs1to3_lmDur_S','t')
TScombineover('G450Qs1to3_lmAft','Qs1to3_lmAft_S','t')
TScombineover('G450Qs1to3_WpsOn','Qs1to3_WpsOn_S','t')
TScombineover('G450Qs1to3_WpsOff','Qs1to3_WpsOff_S','t')
TScombineover('G450Qs1to3_PsWidth','Qs1to3_PsWidth_S','t')
TScombineover('G450Qs1to3_LngstIspI','Qs1to3_LngstIspI_S','t')
TScombineover('G450Qs1to3_LaTtoLngstIspI','Qs1to3_LatToLngstIspI_S','t')
TScombineover('G450PsOnMeanStdCoV','PsOnMeanStdCoV_S','t')
TScombineover('G450PsOffMeanStdCoV','PsOffMeanStdCoV_S','t')
TScombineover('G450PsParamsCorrelations','VecsOfPsParamCorrs')
TScombineover('G450PsWidthMeanStdCoV','PsWidthMeanStdCoV_S')

%% Cell 22: Plotting Pause Onset & Offset Latencies and Pause Width and
% their CoVs
figure
subplot(2,3,2)
    D1 = Experiment.G150PsOffMeanStdCoV(:,[1 3]);
    D2 = Experiment.G200PsOffMeanStdCoV(:,[1 3]);
    D3 = Experiment.G300PsOffMeanStdCoV(:,[1 3]);
    D4 = Experiment.G400PsOffMeanStdCoV(:,[1 3]);
    D5 = Experiment.G450PsOffMeanStdCoV(:,[1 3]);
    plot(repmat(.15,1,length(D1)),D1(:,1),'k*',repmat(.2,1,length(D2)),D2(:,1),'k*',...
        repmat(.3,1,length(D3)),D3(:,1),'k*',repmat(.4,1,length(D4)),D4(:,1),'k*',...
        repmat(.45,1,length(D5)),D5(:,1),'k*')
    ylabel('Pause Off Latency (s)')
    xlim([0 .5]); ylim([0 .6])
subplot(2,3,5)
    plot(repmat(.15,1,length(D1)),D1(:,2),'k*',repmat(.2,1,length(D2)),D2(:,2),'k*',...
        repmat(.3,1,length(D3)),D3(:,2),'k*',repmat(.4,1,length(D4)),D4(:,2),'k*',...
        repmat(.45,1,length(D5)),D5(:,2),'k*')
    xlim([0 .5]); ylim([0 .6])
    xlabel('CS-US Interval (s)')
%
subplot(2,3,1)
    D1 = Experiment.G150PsOnMeanStdCoV(:,[1 3]);
    D2 = Experiment.G200PsOnMeanStdCoV(:,[1 3]);
    D3 = Experiment.G300PsOnMeanStdCoV(:,[1 3]);
    D4 = Experiment.G400PsOnMeanStdCoV(:,[1 3]);
    D5 = Experiment.G450PsOnMeanStdCoV(:,[1 3]);
    plot(repmat(.15,1,length(D1)),D1(:,1),'k*',repmat(.2,1,length(D2)),D2(:,1),'k*',...
        repmat(.3,1,length(D3)),D3(:,1),'k*',repmat(.4,1,length(D4)),D4(:,1),'k*',...
        repmat(.45,1,length(D5)),D5(:,1),'k*')
    ylabel('Pause On Latency (s)')
    xlim([0 .5]); ylim([0 .6])
subplot(2,3,4)
    plot(repmat(.15,1,length(D1)),D1(:,2),'k*',repmat(.2,1,length(D2)),D2(:,2),'k*',...
        repmat(.3,1,length(D3)),D3(:,2),'k*',repmat(.4,1,length(D4)),D4(:,2),'k*',...
        repmat(.45,1,length(D5)),D5(:,2),'k*')
    xlim([0 .5]); ylim([0 .6])
    ylabel('CoV (\sigma/\mu) ')
    xlabel('CS-US Interval (s)')    
subplot(2,3,3)
    D1 = Experiment.G150PsWidthMeanStdCoV(:,[1 3]);
    D2 = Experiment.G200PsWidthMeanStdCoV(:,[1 3]);
    D3 = Experiment.G300PsWidthMeanStdCoV(:,[1 3]);
    D4 = Experiment.G400PsWidthMeanStdCoV(:,[1 3]);
    D5 = Experiment.G450PsWidthMeanStdCoV(:,[1 3]);
    plot(repmat(.15,1,length(D1)),D1(:,1),'k*',repmat(.2,1,length(D2)),D2(:,1),'k*',...
        repmat(.3,1,length(D3)),D3(:,1),'k*',repmat(.4,1,length(D4)),D4(:,1),'k*',...
        repmat(.45,1,length(D5)),D5(:,1),'k*')
    xlim([0 .5]); ylim([0 .6])
    ylabel('Pause Width (s)')
subplot(2,3,6)   
    plot(repmat(.15,1,length(D1)),D1(:,2),'k*',repmat(.2,1,length(D2)),D2(:,2),'k*',...
        repmat(.3,1,length(D3)),D3(:,2),'k*',repmat(.4,1,length(D4)),D4(:,2),'k*',...
        repmat(.45,1,length(D5)),D5(:,2),'k*')
    xlim([0 .5]); ylim([0 .6])
    xlabel('CS-US Interval (s)')

saveas(gcf,'Figures/PsOnOffWdth&CoVs')
    % The on latency, the off latency, and the pause width all increase
    % with the CS-US interval during training, while the CoVs remain
    % constant
pause(20)
close
%% Cell 23: Plotting Correlations
% A(1,2) is Off vs On; A(1,3) is the correlation of PsWidth w PsOn;
% A(1,4) is the correlation of MxIspI w PsOn; A(2,3) is the correlation of
% PsOff with PsWidth; A(2,4) is the correlation of MxIspI w PsOff; A(3,4)
% is the correlation of MxIspI w PsWidth
figure
subplot(5,1,1)
D1 = Experiment.G150PsParamsCorrelations;
plot(ones(size(D1,1),1),D1(:,1),'k*',2*ones(size(D1,1),1),D1(:,2),'k*',...
    3*ones(size(D1,1),1),D1(:,3),'k*',4*ones(size(D1,1),1),D1(:,4),'k*',...
    5*ones(size(D1,1),1),D1(:,5),'k*',6*ones(size(D1,1),1),D1(:,6),'k*')
set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{char([8595 118 8593]) ...
    char([87 118 8593]) char([77 118 8593]) char([87 118 8595]) ...
    char([77 118 8595]) 'MvW'})
% When this is run, the question marks should be replaced with up and down
% arrows; Matlab's editor converts them to question marks when the code is
% saved. The sequence of conversions is: down, up, up, up, down, down. The
% corresponding sequence of 6 correlations is: 1) offset vs onset; 2) Width
% vs Onset; 3) Max ispki vs Onset; 4) Width vs Offset; 5) Max inspki vs
% Offset; 6) Max ispki vs Width. I have replaced the question marks with
% the numerical codes for these arrows. Unfortunately, that code is font
% and/or platform specific
xlim([.5 6.5])
hold on
plot(xlim,[0 0],'k--')
title('CS-US = 0.15')
%
subplot(5,1,2)
D2 = Experiment.G200PsParamsCorrelations;
plot(ones(size(D2,1),1),D2(:,1),'k*',2*ones(size(D2,1),1),D2(:,2),'k*',...
    3*ones(size(D2,1),1),D2(:,3),'k*',4*ones(size(D2,1),1),D2(:,4),'k*',...
    5*ones(size(D2,1),1),D2(:,5),'k*',6*ones(size(D2,1),1),D2(:,6),'k*')
set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{char([8595 118 8593]) ...
    char([87 118 8593]) char([77 118 8593]) char([87 118 8595]) ...
    char([77 118 8595]) 'MvW'})
% in re '?', see first subplot for explanation
xlim([.5 6.5])
hold on
plot(xlim,[0 0],'k--')
title('CS-US = 0.2')

subplot(5,1,3)
D3 = Experiment.G300PsParamsCorrelations;
plot(ones(size(D3,1),1),D3(:,1),'k*',2*ones(size(D3,1),1),D3(:,2),'k*',...
    3*ones(size(D3,1),1),D3(:,3),'k*',4*ones(size(D3,1),1),D3(:,4),'k*',...
    5*ones(size(D3,1),1),D3(:,5),'k*',6*ones(size(D3,1),1),D3(:,6),'k*')
set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{char([8595 118 8593]) ...
    char([87 118 8593]) char([77 118 8593]) char([87 118 8595]) ...
    char([77 118 8595]) 'MvW'})
% in re '?', see first subplot for explanation
xlim([.5 6.5])
hold on
plot(xlim,[0 0],'k--')
title('CS-US = 0.3')
ylabel('Correlation Coefficient')

subplot(5,1,4)
D4 = Experiment.G400PsParamsCorrelations;
plot(ones(size(D4,1),1),D4(:,1),'k*',2*ones(size(D4,1),1),D4(:,2),'k*',...
    3*ones(size(D4,1),1),D4(:,3),'k*',4*ones(size(D4,1),1),D4(:,4),'k*',...
    5*ones(size(D4,1),1),D4(:,5),'k*',6*ones(size(D4,1),1),D4(:,6),'k*')
set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{char([8595 118 8593]) ...
    char([87 118 8593]) char([77 118 8593]) char([87 118 8595]) ...
    char([77 118 8595]) 'MvW'})
% in re '?', see 1st subplot for explanation
xlim([.5 6.5])
hold on
plot(xlim,[0 0],'k--')
title('CS-US = 0.4')

subplot(5,1,5)
D5 = Experiment.G450PsParamsCorrelations;
plot(ones(size(D5,1),1),D5(:,1),'k*',2*ones(size(D5,1),1),D5(:,2),'k*',...
    3*ones(size(D5,1),1),D5(:,3),'k*',4*ones(size(D5,1),1),D5(:,4),'k*',...
    5*ones(size(D5,1),1),D5(:,5),'k*',6*ones(size(D5,1),1),D5(:,6),'k*')
set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{char([8595 118 8593]) ...
    char([87 118 8593]) char([77 118 8593]) char([87 118 8595]) ...
    char([77 118 8595]) 'MvW'})
% in re '?', see 1st subplot for explanation
xlim([.5 6.5])
hold on
plot(xlim,[0 0],'k--')
title('CS-US = 0.45')
saveas(gcf,'Figures/PauseParameterCorrelations')
disp(char({'';'The x axis labels of the plots in this figure use up and down';...
    'arrows to indicate pause onset & offset, respectively. Unfortunately';...
    'the codes for these are platform & even font specific, so they may not';...
    'reproduce correctly on non-Mackintosh platforms and perhaps even on';...
    'other Mackintoshs. The middle character in each of the 6 lables is';...
    '''v'' (for versus) The first label on each axis should have a down';...
    'arrow (for pause offset) before the ''v'' with an up arrow for';...
    '(pause onset) after the ''v''. The 2nd label is for Pause Width (''W'')';...
    ' versus pause onset; the 3rd for Maximum IspkI (''M'') versus pause';...
    'onset; the 4th for Width versus offset; the 5th for ''M'' vs offset;';...
    'and the last is for Maximum IspkI vs Width.';''}))
pause(20)
close

%% Cell 24: CDFs of Group quartile pause-onset & pause-offset stats
figure % Pause Onsets & Offsets
subplot(5,2,1)
H=cdfplot(Experiment.G150Qs1to3_PsOn(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G150Qs1to3_PsOn(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G150Qs1to3_PsOn(:,3));
set(H,'LineStyle','--','Color','k')
xlim([-.01 .3])
xlabel('');ylabel('')
title('')
legend('median','Q1 & Q3','location','NE')

subplot(5,2,3)
H=cdfplot(Experiment.G200Qs1to3_PsOn(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G200Qs1to3_PsOn(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G200Qs1to3_PsOn(:,3));
xlim([-.01 .3])

set(H,'LineStyle','--','Color','k')
xlabel('');ylabel('')
title('')

subplot(5,2,5)
H=cdfplot(Experiment.G300Qs1to3_PsOn(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G300Qs1to3_PsOn(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G300Qs1to3_PsOn(:,3));
xlim([-.01 .3])

set(H,'LineStyle','--','Color','k')
xlabel('');ylabel('Cumulative Fraction of Cells')
title('')

subplot(5,2,7)
H=cdfplot(Experiment.G400Qs1to3_PsOn(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G400Qs1to3_PsOn(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G400Qs1to3_PsOn(:,3));
xlim([-.01 .3])
set(H,'LineStyle','--','Color','k')
xlabel('');ylabel('')
title('')

subplot(5,2,9)
H=cdfplot(Experiment.G450Qs1to3_PsOn(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G450Qs1to3_PsOn(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G450Qs1to3_PsOn(:,3));
set(H,'LineStyle','--','Color','k')
xlim([-.01 .3])
xlabel('Pause Onset Latency (s)');ylabel('')
title('')

subplot(5,2,2)
H=cdfplot(Experiment.G150Qs1to3_PsOff(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G150Qs1to3_PsOff(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G150Qs1to3_PsOff(:,3));
set(H,'LineStyle','--','Color','k')
xlim([-.1 .6])
xlabel('');ylabel('')
title('')

subplot(5,2,4)
H=cdfplot(Experiment.G200Qs1to3_PsOff(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G200Qs1to3_PsOff(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G200Qs1to3_PsOff(:,3));
xlim([-.1 .6])
set(H,'LineStyle','--','Color','k')
xlabel('');ylabel('')
title('')

subplot(5,2,6)
H=cdfplot(Experiment.G300Qs1to3_PsOff(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G300Qs1to3_PsOff(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G300Qs1to3_PsOff(:,3));
xlim([-.1 .6])
set(H,'LineStyle','--','Color','k')
xlabel('');ylabel('Cumulative Fraction of Cells')
title('')

subplot(5,2,8)
H=cdfplot(Experiment.G400Qs1to3_PsOff(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G400Qs1to3_PsOff(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G400Qs1to3_PsOff(:,3));
xlim([-.1 .6])
set(H,'LineStyle','--','Color','k')
xlabel('');ylabel('')
title('')

subplot(5,2,10)
H=cdfplot(Experiment.G450Qs1to3_PsOff(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G450Qs1to3_PsOff(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G450Qs1to3_PsOff(:,3));
set(H,'LineStyle','--','Color','k')
xlim([-.1 .6])
xlabel('Pause Offset Latency (s)');ylabel('')
title('')
saveas(gcf,'Figures/CDFsPsOn&OffQuartilesByGroup')
disp(char({'';'The panels in this figure plot the cumulative distributions of';...
    'the quartiles of the distributions of pause onset and offset latencies,';...
    'with each panel for a different CS-US training condition. As you look';...
    'down the panels, you see the systematic rightward shifts in these distributions';''}))
pause(20)
close

%% Cell 25: CDFs of Group quartiles of pause width and maximum within-pause
% IspkI
figure % 
subplot(5,2,1)
H=cdfplot(Experiment.G150Qs1to3_PsWidth(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G150Qs1to3_PsWidth(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G150Qs1to3_PsWidth(:,3));
set(H,'LineStyle','--','Color','k')
xlim([0 .6])
xlabel('');ylabel('')
title('')
legend('median','Q1 & Q3','location','NE')

subplot(5,2,3)
H=cdfplot(Experiment.G200Qs1to3_PsWidth(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G200Qs1to3_PsWidth(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G200Qs1to3_PsWidth(:,3));
xlim([0 .6])
set(H,'LineStyle','--','Color','k')
xlabel('');ylabel('')
title('')

subplot(5,2,5)
H=cdfplot(Experiment.G300Qs1to3_PsWidth(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G300Qs1to3_PsWidth(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G300Qs1to3_PsWidth(:,3));
xlim([0 .6])
set(H,'LineStyle','--','Color','k')
xlabel('');ylabel('Cumulative Fraction of Cells')
title('')

subplot(5,2,7)
H=cdfplot(Experiment.G400Qs1to3_PsWidth(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G400Qs1to3_PsWidth(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G400Qs1to3_PsWidth(:,3));
xlim([0 .6])
set(H,'LineStyle','--','Color','k')
xlabel('');ylabel('')
title('')

subplot(5,2,9)
H=cdfplot(Experiment.G450Qs1to3_PsWidth(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G450Qs1to3_PsWidth(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G450Qs1to3_PsWidth(:,3));
set(H,'LineStyle','--','Color','k')
xlim([0 .6])
xlabel('Pause Width (s)');ylabel('')
title('')

subplot(5,2,2)
H=cdfplot(Experiment.G150Qs1to3_LngstIspI(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G150Qs1to3_LngstIspI(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G150Qs1to3_LngstIspI(:,3));
set(H,'LineStyle','--','Color','k')
xlim([-.1 .6])
xlabel('');ylabel('')
title('')

subplot(5,2,4)
H=cdfplot(Experiment.G200Qs1to3_LngstIspI(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G200Qs1to3_LngstIspI(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G200Qs1to3_LngstIspI(:,3));
xlim([-.1 .6])
set(H,'LineStyle','--','Color','k')
xlabel('');ylabel('')
title('')

subplot(5,2,6)
H=cdfplot(Experiment.G300Qs1to3_LngstIspI(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G300Qs1to3_LngstIspI(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G300Qs1to3_LngstIspI(:,3));
xlim([-.1 .6])
set(H,'LineStyle','--','Color','k')
xlabel('');ylabel('Cumulative Fraction of Cells')
title('')

subplot(5,2,8)
H=cdfplot(Experiment.G400Qs1to3_LngstIspI(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G400Qs1to3_LngstIspI(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G400Qs1to3_LngstIspI(:,3));
xlim([-.1 .6])
set(H,'LineStyle','--','Color','k')
xlabel('');ylabel('')
title('')

subplot(5,2,10)
H=cdfplot(Experiment.G450Qs1to3_LngstIspI(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G450Qs1to3_LngstIspI(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G450Qs1to3_LngstIspI(:,3));
set(H,'LineStyle','--','Color','k')
xlim([-.1 .6])
xlabel('Longest IspI (s)');ylabel('')
title('')
saveas(gcf,'Figures/CDFsQrtlsPsWdth&LngstIspkI')
disp(char({'';'The panels in this figure plot the cumulative distributions of';...
    'the quartiles of the distributions of pause widths and maximum within-pause IspkI,';...
    'with each panel for a different CS-US training condition. As you look';...
    'down the panels, you see the systematic rightward shifts in these distributions';''}))
pause(20)
close

%% Cell 26: CDFs of latency to onset of longest IspkI
figure
subplot(5,1,1)
H=cdfplot(Experiment.G150Qs1to3_LaTtoLngstIspI(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G150Qs1to3_LaTtoLngstIspI(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G150Qs1to3_LaTtoLngstIspI(:,3));
set(H,'LineStyle','--','Color','k')
xlim([0 .2])
xlabel('');ylabel('')
title('')
legend('median','Q1 & Q3','location','NE')

subplot(5,1,2)
H=cdfplot(Experiment.G200Qs1to3_LaTtoLngstIspI(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G200Qs1to3_LaTtoLngstIspI(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G200Qs1to3_LaTtoLngstIspI(:,3));
set(H,'LineStyle','--','Color','k')
xlim([0 .2])
xlabel('');ylabel('')
title('')


subplot(5,1,3)
H=cdfplot(Experiment.G300Qs1to3_LaTtoLngstIspI(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G300Qs1to3_LaTtoLngstIspI(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G300Qs1to3_LaTtoLngstIspI(:,3));
set(H,'LineStyle','--','Color','k')
xlim([0 .2])
xlabel('');ylabel('')
title('')


subplot(5,1,4)
H=cdfplot(Experiment.G400Qs1to3_LaTtoLngstIspI(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G400Qs1to3_LaTtoLngstIspI(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G400Qs1to3_LaTtoLngstIspI(:,3));
set(H,'LineStyle','--','Color','k')
xlim([0 .2])
xlabel('');ylabel('')
title('')


subplot(5,1,5)
H=cdfplot(Experiment.G450Qs1to3_LaTtoLngstIspI(:,2));
set(H,'LineWidth',2)
hold on
H=cdfplot(Experiment.G450Qs1to3_LaTtoLngstIspI(:,1));
set(H,'LineStyle','--','Color','k')
H=cdfplot(Experiment.G450Qs1to3_LaTtoLngstIspI(:,3));
set(H,'LineStyle','--','Color','k')
xlim([0 .2])
xlabel('Latency to Longest IspI (s)');ylabel('')
title('')
saveas(gcf,'Figures/CDFsQrtlsLatLngstIspkI')
pause(20)
close

%% Cell 27: Computing & plotting interspike intervals looking backward 
% from pause onsets to CS onsets
TSlimit('Subjects',[Experiment.CS300electrodes Experiment.CS400electrodes ...
    Experiment.CS450electrodes])
TSlimit('Sessions','all')
TSapplystat('BkwdIspkInts',{'BinSpkVecS','PsOn'},@BkwdISpIs)
% % creates field at Session level, giving interspike intervals looking
% backward from PsOn to CS onset: 2 cols: 1st is backward count of
% interspike intervals (-1, -2, etc); 2nd is the interspike intervals
%{
function O=BkwdISpIs(bsv,PsOns)
O = double.empty(0,2); % initializing output array
rws = find(bsv(:,1)>-.0005 & bsv(:,1)<.0005); % CS onset rows
rnv = (1:length(bsv))';
LVs = bsv(:,2)>0; % flags spike rows
i=1;
for br = rws' % stepping through the trials
    if PsOns(i)>0
        er = br+round(PsOns(i)/.001); % row # at pause onset
        LVt = rnv>br&rnv<er; % flags rows btw CS onset and pause onset
        spktms = bsv(LVt&LVs,1);
        n = length(spktms);
        if n>1 % more than 1 spike btw CS onset & pause onset
            O = [O;[(-1:-1:-n+1)' flipud(diff(spktms))]];
        else
            i = i+1;
            continue
        end 
    else
        i=i+1;
        continue
    end
    i=i+1;
end    
%}
%%
% Plotting the backward looking interspike intervals
mkdir('Figures/RetrospectiveIspIs')
TSlimit('Subjects',Experiment.CS300electrodes)
TSapplystat('','BkwdIspkInts',@TSplot,'Scat','.','Xlbl','IspkIs bck from PsOn',...
    'Ylbl','IspkI (s)')
for f=1:9
    figure(f)
    saveas(gcf,['Figures/RetrospectiveIspIs/CSdur300f' num2str(f) 'of9'])
    close
end
%    
TSlimit('Subjects',Experiment.CS400electrodes)
TSapplystat('','BkwdIspkInts',@TSplot,'Scat','.','Xlbl','IspkIs bck from PsOn',...
    'Ylbl','IspkI (s)')
saveas(gcf,['Figures/RetrospectiveIspIs/CSdur400'])

%    
TSlimit('Subjects',Experiment.CS450electrodes)
TSapplystat('','BkwdIspkInts',@TSplot,'Scat','.','Xlbl','IspkIs bck from PsOn',...
    'Ylbl','IspkI (s)')
saveas(gcf,['Figures/RetrospectiveIspIs/CSdur450'])%%
close all

%% Cell 28: Linear regressions on backward interspike intervals
TSlimit('Subjects',[Experiment.CS300electrodes Experiment.CS400electrodes ...
    Experiment.CS450electrodes])
TSapplystat('LinRegBkInts','BkwdIspkInts',@LRbkInts)
%{
function O = LRbkInts(D)
O=double.empty(0,3); %initializing output
if length(D)<5 || min(D(:,1))>-3
    return
else
    [FO G] = fit(D(:,1),D(:,2),'poly1');
    CI = confint(FO);
    O = [FO.p1 CI(1,1) G.rsquare]; % slope of regression, lower limit on
    % slope confidence interval; amt variance accounted for
end
%}
TScombineover('LinRegBkInts_S','LinRegBkInts','t') % raising reg rslts to Subject level
TScombineover('LinRegBkInts_Exp','LinRegBkInts_S','t') % raising to Experiment level
%%
figure
subplot(3,1,1)
    cdfplot(Experiment.LinRegBkInts_Exp(:,1))
    xlabel('Slope');ylabel('Cum Frac of Regs')
    title('')
subplot(3,1,2)
    cdfplot(Experiment.LinRegBkInts_Exp(:,2))
    xlabel('Lwr Limit on Slope');ylabel('Cum Frac of Regs')
    title('')
subplot(3,1,3)
    cdfplot(Experiment.LinRegBkInts_Exp(:,3))
    xlabel('Variance Explained');ylabel('Cum Frac of Regs')
    title('')
saveas(gcf,'Figures/CDFsOfBkwdRegParams')
disp(char({'';'This figure gives the cumulative distributions of the results';...
    'of the linear regressions of IspkI on # of IspkIs counting backward';...
    'from pause onset to the CS onset. The slopes tend to be positive';...
    '(top panel), but the lower confidence limit on the slope is negative';...
    'in 80% of the regressions (middle panel), and less than 5% variance';...
    'is  accounted for in 70% of the cases (bottom panel)';''}))
pause(10)
close
TSsaveexperiment
TSloadscript('PurkinjePauseExperiment') % loads the script into Experiment.Script
% field
TSsaveexperiment
TSloadscript('PurkinjePauseExperiment') % loads script into Experiment.Script field
TSsaveexperiment
%% Max trial-trial intervals
for S = 1:Experiment.NumSubjects
    for s=1:Experiment.Subject(S).NumSessions
        mxdr = 0;
        for t = 1:Experiment.Subject(S).Session(s).TrialMdToMd.NumTrials
            if Experiment.Subject(S).Session(s).TrialMdToMd.Trial(t).TrialDuration>mxdr
                mxdr=Experiment.Subject(S).Session(s).TrialMdToMd.Trial(t).TrialDuration;
            end
        end
        format bank
        fprintf('\n%d,%d,%d: %2.1f',S,s,t,mxdr)
    end
end
            