% Load Data
clear all
clc

% datafile = '/home/melisagumus/Documents/MATLAB/CA1_SimTracker/pyr/pyr_36884_1000/mytrace_36884_syns.dat';
%
% ls('/home/melisagumus/Documents/MATLAB/CA1_SimTracker/pyr/pyr*1000');

f = fullfile('/home','melisagumus','Documents', ...
    'MATLAB','CA1_SimTracker','pyr',{...
    'pyr_29097_1000';...
    'pyr_36884_1000';...
    'pyr_52458_1000';...
    'pyr_68032_1000';...
    'pyr_83606_1000';...
    'pyr_99180_1000';...
    'pyr_106967_1000';...
    'pyr_114754_1000';...
    'pyr_153689_1000';...
    'pyr_177050_1000';...
    'pyr_200411_1000';...
    'pyr_254920_1000';...
    'pyr_286068_1000';...
    'pyr_301642_1000';...
    'pyr_325003_1000'...
    },{...
    'mytrace_29097_syns.dat';...
    'mytrace_36884_syns.dat';...
    'mytrace_52458_syns.dat';...
    'mytrace_68032_syns.dat';...
    'mytrace_83606_syns.dat';...
    'mytrace_99180_syns.dat';...
    'mytrace_106967_syns.dat';...
    'mytrace_114754_syns.dat';...
    'mytrace_153689_syns.dat';...
    'mytrace_177050_syns.dat';...
    'mytrace_200411_syns.dat';...
    'mytrace_254920_syns.dat';...
    'mytrace_286068_syns.dat';...
    'mytrace_301642_syns.dat';...
    'mytrace_325003_syns.dat'...
    });

%% Alternative - does not work
clear all
clc

start_path = fullfile(matlabroot, '/home','melisagumus','Documents', ...
    'MATLAB','CA1_SimTracker','pyr');
topLevelFolder = uigetdir(start_path);
if topLevelFolder == 0
    return;
end

% To get list of all subfolders
allSubFolders = genpath(topLevelFolder);
remain = allSubFolders;
listOfFolderNames = {};
while true
    [singleSubFolder, remain] = strtok(remain, ';');
    if isempty(singleSubFolder)
        break;
    end
    listOfFolderNames = [listOfFolderNames singleSubFolder];
end
numberOfFolders = length(listOfFolderNames);

% process text files in the folders
for k = 1:numberOfFolders
    thisFolder = listOfFolderNames{k};
    fprintf('Processing folder %s\n', thisFolder);

    filePattern = sprintf('mytrace %s/*.dat', thisFolder);
    baseFileNames = dir(filePattern);
    numberOfFiles = length(baseFileNames);

    if numberOfFiles >= 1
        for f = 1:numberOfFiles
            fullFileName = fullfile(thisFolder, baseFileNames{f}.name);
            fprint('Processing text files %s\n', fullFileName);
        end
    else
        fprintf('Folder %s has no text files in it.\n', thisFolder);
    end
end


%%
%cd C:\Users\Melisa\Documents\MATLAB\CA1_SimTracker\pyr\
%files = dir('pyr*1000\mytrace*syns.dat');
%pyr_names = f.names;

alldata = [];
for m = 1:1:15
    temp_data = readtable(f{m},'Delimiter','\t');
    temp_data = table2array(temp_data);
    alldata = [alldata temp_data];
end

data = mat2cell(alldata, 40000, ...
    [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]);

%% Creates a big table consists of inputs from AAC BiC PYR and BC in order

M = [];
current_AAC = [];
current_BC = [];
current_PYR = [];
current_BiC = [];
for m = 1:15  % number of cells
    for k = 2:12  % number of input
        if k ==2
            temp_current_AAC = data{m}(:,k);
            [pks, locs] = findpeaks(data{m}(:,k),'MinPeakDistance',4000); % peak detection with interval based threshold
            temp_AAC = data{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_AAC(element,:) = 0;
            end
            AAC = temp_AAC;
        elseif k == 3
            temp_current_BiC = data{m}(:,k);
            [pks, locs] = findpeaks(data{m}(:,k),'MinPeakDistance',4000); % peak detection
            temp_BiC = data{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_BiC(element,:) = 0;
            end
            BiC = temp_BiC;
        elseif k == 8
            temp_current_PYR = data{m}(:,k);
            [pks, locs] = findpeaks(-data{m}(:,k),'MinPeakDistance',4000); % peak detection
            temp_PYR = data{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_PYR(element,:) = 0;
            end
            PYR = temp_PYR;
        elseif k == 9
            temp_current_BC = data{m}(:,k);
            [pks, locs] = findpeaks(data{m}(:,k),'MinPeakDistance',4000); % peak detection
            %findpeaks(data{m}(:,k),M
            temp_BC = data{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_BC(element,:) = 0;
            end
            BC = temp_BC;
        end
    end
    current_AAC = [current_AAC temp_current_AAC];
    current_PYR = [current_PYR temp_current_PYR];
    current_BiC = [current_BiC temp_current_BiC];
    current_BC = [current_BC temp_current_BC];
    M = [M AAC BiC PYR BC];
end

allcells = mat2cell(M, 40000, ...
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]);

%% Graph EPSCs on PYR

figure
for i = 1:1:10
    temp = allcells{i}(:,3);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(5,2,i)
    histfit(temp)
    hold on
    title (['Pyr' num2str(i)]) % name each graph with the corresponding pyr #
end
suptitle('EPSC Distribution on PYR Cells')

%% Graph IPSCs from AAC onto PYR

figure
for i = 1:1:10
    temp = allcells{i}(:,1);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(5,2,i)
    histfit(temp,100)
    xlim([-2 2])
    hold on
    title (['Pyr' num2str(i)]) % name each graph with the corresponding pyr #
end
suptitle('IPSC Distribution from AAC onto PYR')

%% Graph IPSCs from BiC onto PYR

figure
for i = 1:1:10
    temp = allcells{i}(:,2);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(5,2,i)
    histfit(temp,50)
    xlim([-2 2])
    hold on
    title (['Pyr' num2str(i)]) % name each graph with the corresponding pyr #
end
suptitle('IPSC Distribution from BiC onto PYR')

%% Graph IPSCs from BC onto PYR

figure
for i = 1:1:10
    temp = allcells{i}(:,2);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(5,2,i)
    histfit(temp,50)
    xlim([-2 2])
    hold on
    title (['Pyr' num2str(i)]) % name each graph with the corresponding pyr #
end
suptitle('IPSC Distribution from BC onto PYR')


%% EPSCs from PYR

EPSC = [];
epsc = [];
for i = 1:1:15 % number of PYR cells
    pks_epsc = allcells{i}(:,3);
    pks_epsc(pks_epsc==0)=[];
    epsc_mean = mean(pks_epsc);
    epsc_std = std(pks_epsc);
    epsc = [epsc_mean;epsc_std];
    EPSC = [EPSC epsc];
end

EPSC = array2table(EPSC);
EPSC.Properties.VariableNames = {'pyr1'...
    'pyr2' 'pyr3' 'pyr4' 'pyr5' 'pyr6'...
    'pyr7' 'pyr8' 'pyr9' 'pyr10' 'pyr11'...
    'pyr12' 'pyr13' 'pyr14' 'pyr15'};

subplot(2,2,1)
EPSC_mean = table2array(EPSC(1,:));
EPSC_std = table2array(EPSC(2,:));
x = linspace(0,15,length(EPSC_mean));
scatter(x,EPSC_mean,'black','filled');
xlabel('Individual PYR Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak EPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,EPSC_mean,EPSC_std,'b','LineStyle','none')
title('Mean Peak EPSC onto PYR cells','FontSize',15,'FontWeight','bold')

%% IPSCs only from AAC

IPSC_AAC = [];
ipsc_AAC = [];
for i = 1:1:15 % number of PYR cells
    pks_ipsc_AAC = allcells{i}(:,1);
    pks_ipsc_AAC(pks_ipsc_AAC==0)=[];
    ipsc_AAC_mean = mean(pks_ipsc_AAC);
    ipsc_AAC_std = std(pks_ipsc_AAC);
    ipsc_AAC = [ipsc_AAC_mean;ipsc_AAC_std];
    IPSC_AAC = [IPSC_AAC ipsc_AAC];
end

IPSC_AAC = array2table(IPSC_AAC);
IPSC_AAC.Properties.VariableNames = {'pyr1'...
    'pyr2' 'pyr3' 'pyr4' 'pyr5' 'pyr6'...
    'pyr7' 'pyr8' 'pyr9' 'pyr10' 'pyr11'...
    'pyr12' 'pyr13' 'pyr14' 'pyr15'};

subplot(2,2,2)
IPSC_AAC_mean = table2array(IPSC_AAC(1,:));
IPSC_AAC_std = table2array(IPSC_AAC(2,:));
x = linspace(0,15,length(IPSC_AAC_mean));
scatter(x,IPSC_AAC_mean,'black','filled');
xlabel('Individual PYR Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC from AAC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_AAC_mean,IPSC_AAC_std,'b','LineStyle','none')
title('Mean Peak IPSC from AAC onto PYR cells','FontSize',15,'FontWeight','bold')

%% IPSCs only from BiC

IPSC_BiC = [];
ipsc_BiC = [];
for i = 1:1:15 % number of PYR cells
    pks_ipsc_BiC = allcells{i}(:,2);
    pks_ipsc_BiC(pks_ipsc_BiC==0)=[];
    ipsc_BiC_mean = mean(pks_ipsc_BiC);
    ipsc_BiC_std = std(pks_ipsc_BiC);
    ipsc_BiC = [ipsc_BiC_mean;ipsc_BiC_std];
    IPSC_BiC = [IPSC_BiC ipsc_BiC];
end

IPSC_BiC = array2table(IPSC_BiC);
IPSC_BiC.Properties.VariableNames = {'pyr1'...
    'pyr2' 'pyr3' 'pyr4' 'pyr5' 'pyr6'...
    'pyr7' 'pyr8' 'pyr9' 'pyr10' 'pyr11'...
    'pyr12' 'pyr13' 'pyr14' 'pyr15'};

subplot(2,2,3)
IPSC_BiC_mean = table2array(IPSC_BiC(1,:));
IPSC_BiC_std = table2array(IPSC_BiC(2,:));
x = linspace(0,15,length(IPSC_BiC_mean));
scatter(x,IPSC_BiC_mean,'black','filled');
xlabel('Individual PYR Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC from BiC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BiC_mean,IPSC_BiC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BiC onto PYR cells','FontSize',15,'FontWeight','bold')

%% IPSCs only from BC

IPSC_BC = [];
ipsc_BC = [];
for i = 1:1:15 % number of PYR cells
    pks_ipsc_BC = allcells{i}(:,4);
    pks_ipsc_BC(pks_ipsc_BC==0)=[];
    ipsc_BC_mean = mean(pks_ipsc_BC);
    ipsc_BC_std = std(pks_ipsc_BC);
    ipsc_BC = [ipsc_BC_mean;ipsc_BC_std];
    IPSC_BC = [IPSC_BC ipsc_BC];
end

IPSC_BC = array2table(IPSC_BC);
IPSC_BC.Properties.VariableNames = {'pyr1'...
    'pyr2' 'pyr3' 'pyr4' 'pyr5' 'pyr6'...
    'pyr7' 'pyr8' 'pyr9' 'pyr10' 'pyr11'...
    'pyr12' 'pyr13' 'pyr14' 'pyr15'};

subplot(2,2,4)
IPSC_BC_mean = table2array(IPSC_BC(1,:));
IPSC_BC_std = table2array(IPSC_BC(2,:));
x = linspace(0,15,length(IPSC_BC_mean));
scatter(x,IPSC_BC_mean,'black','filled');
xlabel('Individual PYR Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC from BC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BC_mean,IPSC_BC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BC onto PYR cells','FontSize',15,'FontWeight','bold')


%% all ipsc current onto PYR gathered

% Sum all ipsc currents
all_ipsc = [];
for i = 1:1:15
    tot_cur_ipsc = current_AAC(:,i) + current_BiC(:,i) + current_BC(:,i);
    all_ipsc = [all_ipsc tot_cur_ipsc];
end

% Find the peaks of the summed ipsc currents
peaks_all_PV = [];
for k = 1:1:15
    [pks, locs] = findpeaks(all_ipsc(:,k),'MinPeakDistance',4000); % peak detection
    temp_cur = all_ipsc(:,k);
    allrows = (1:40000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    peaks_all = temp_cur;
    peaks_all_PV = [peaks_all_PV peaks_all];
end

%% All ipsc currents together onto PYR - graph and table

IPSC_all = [];
ipsc_all = [];
for i = 1:1:15 % number of PYR cells
    pks_ipsc_all = peaks_all_PV(:,i);
    pks_ipsc_all(pks_ipsc_all == 0) = [];
    ipsc_all_mean = mean(pks_ipsc_all);
    ipsc_all_std = std(pks_ipsc_all);
    ipsc_all = [ipsc_all_mean;ipsc_all_std];
    IPSC_all = [IPSC_all ipsc_all];
end

IPSC_all = array2table(IPSC_all);
IPSC_all.Properties.VariableNames = {'pyr1'...
    'pyr2' 'pyr3' 'pyr4' 'pyr5' 'pyr6'...
    'pyr7' 'pyr8' 'pyr9' 'pyr10' 'pyr11'...
    'pyr12' 'pyr13' 'pyr14' 'pyr15'};


IPSC_all_mean = table2array(IPSC_all(1,:));
IPSC_all_std = table2array(IPSC_all(2,:));
x = linspace(0,15,length(IPSC_all_mean));
figure
scatter(x,IPSC_all_mean,'black','filled');
xlabel('Individual PYR Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_all_mean,IPSC_all_std,'b','LineStyle','none')
title('Mean Peak IPSC from BC, BiC and AAC onto PYR cells','FontSize',15,'FontWeight','bold')

%% AAC and BC ipsc current onto PYR gathered

% Sum all ipsc currents
AAC_BC_ipsc = [];
for i = 1:1:15
    AAC_BC_cur_ipsc = current_AAC(:,i) + current_BC(:,i);
    AAC_BC_ipsc = [AAC_BC_ipsc AAC_BC_cur_ipsc];
end

% Find the peaks of the summed ipsc currents
peaks_AAC_BC = [];
for k = 1:1:15
    [pks, locs] = findpeaks(all_ipsc(:,k),'MinPeakDistance',4000); % peak detection
    temp_cur = AAC_BC_ipsc(:,k);
    allrows = (1:40000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    pk_AAC_BC = temp_cur;
    peaks_AAC_BC = [peaks_AAC_BC pk_AAC_BC];
end

%% AAC and BC ipsc currents together onto PYR - graph and table

IPSC_AAC_BC = [];
ipsc_AAC_BC = [];
for i = 1:1:15 % number of PYR cells
    pks_ipsc_AAC_BC = peaks_AAC_BC(:,i);
    pks_ipsc_AAC_BC(pks_ipsc_AAC_BC == 0) = [];
    ipsc_AAC_BC_mean = mean(pks_ipsc_AAC_BC);
    ipsc_AAC_BC_std = std(pks_ipsc_AAC_BC);
    ipsc_AAC_BC = [ipsc_AAC_BC_mean;ipsc_AAC_BC_std];
    IPSC_AAC_BC = [IPSC_AAC_BC ipsc_AAC_BC];
end

IPSC_AAC_BC = array2table(IPSC_AAC_BC);
IPSC_AAC_BC.Properties.VariableNames = {'pyr1'...
    'pyr2' 'pyr3' 'pyr4' 'pyr5' 'pyr6'...
    'pyr7' 'pyr8' 'pyr9' 'pyr10' 'pyr11'...
    'pyr12' 'pyr13' 'pyr14' 'pyr15'};


IPSC_AAC_BC_mean = table2array(IPSC_AAC_BC(1,:));
IPSC_AAC_BC_std = table2array(IPSC_AAC_BC(2,:));
x = linspace(0,15,length(IPSC_AAC_BC_mean));
figure
scatter(x,IPSC_AAC_BC_mean,'black','filled');
xlabel('Individual PYR Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_all_mean,IPSC_all_std,'b','LineStyle','none')
title('Mean Peak IPSC from BC and AAC onto PYR cells','FontSize',15,'FontWeight','bold')


%% E/I Ratios on PYR Cells

IPSC_AAC = table2array(IPSC_AAC);
IPSC_BiC = table2array(IPSC_BiC);
IPSC_BC = table2array(IPSC_BC);
IPSC_all= table2array(IPSC_all);
IPSC_AAC_BC = table2array(IPSC_AAC_BC);
EPSC = table2array(EPSC);

Ratios_PYR = [];
E_I_AAC = abs(EPSC(1,:)./IPSC_AAC(1,:))';
E_I_BC = abs(EPSC(1,:)./IPSC_BC(1,:))';
E_I_BiC = abs(EPSC(1,:)./IPSC_BiC(1,:))';
E_I_AAC_BC = abs(EPSC(1,:)./IPSC_AAC_BC(1,:))';
E_I_all = abs(EPSC(1,:)./IPSC_all(1,:))';

pyr = 1:15;
Ratios_PYR = [pyr' E_I_AAC E_I_BC E_I_BiC E_I_AAC_BC E_I_all];
Ratios_PYR = array2table(Ratios_PYR);

Ratios_PYR.Properties.VariableNames = {'pyr_no' 'Ratio_AAC_PYR'...
    'Ratio_BC_PYR' 'Ratio_BiC_PYR' 'Ratio_AAC_BC_on_PYR' 'Ratio_All_PYR'};

%% Display the table as a figure

uitable('Data',Ratios_PYR{:,:},'ColumnName',Ratios_PYR.Properties.VariableNames,...
    'Units', 'Normalized', 'Position',[0, 0, 1, 1]);


%% Raster Plot for one neuron - NOT FIXED

rasterfile = fullfile('/home','melisagumus','Documents', ...
                        'MATLAB','CA1_SimTracker','pyr',{...
                        'pyr_36884_1000';...
                        'pyr_68032_1000';...
                        'pyr_83606_1000';...
                        'pyr_99180_1000';...
                        'pyr_106967_1000';...
                        'pyr_114754_1000';...
                        'pyr_200411_1000';...
                        'pyr_254920_1000';...
                        'pyr_286068_1000';...
                        'pyr_301642_1000'...
                        },{...
                        'spikeraster.dat';...
                        'spikeraster.dat';...
                        'spikeraster.dat';...
                        'spikeraster.dat';...
                        'spikeraster.dat';...
                        'spikeraster.dat';...
                        'spikeraster.dat';...
                        'spikeraster.dat';...
                        'spikeraster.dat';...
                        'spikeraster.dat'...
                        });

%rasterfile = ('/home/melisagumus/Documents/MATLAB/CA1_SimTracker/pyr/pyr_36884_1000/spikeraster.dat');

raster = readtable(rasterfile{1},'Delimiter','\t');
num_spikes = length(raster.Var2);
var2 = (raster.Var2)';

plot([var2; var2], [ones(1,num_spikes); ones(1, num_spikes) + 1], 'b')
%%
num_spikes = raster.Var1';
plot([var2; var2], [num_spikes; num_spikes], 'b')


%%
t = raster(:,1);
nspikes = numel(t);

for ii = 1:nspikes
    line([t(ii) t(ii)]);%[8474 8475], 'Color', 'k');
end

%%


%%
pd = fitdist(AAC, 'Normal');
mean_cell = mean(pd);
x = 0:10:10000;
y = pdf(pd,x);

histogram(AAC,'Normalization','pdf')
ylim([0 4])
xlim([-5 5])
xlabel('EPSC')

%%
