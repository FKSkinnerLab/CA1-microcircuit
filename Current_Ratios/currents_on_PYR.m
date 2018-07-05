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
test = [];
peak_AAC = [];
peak_BC = [];
peak_PYR = [];
peak_BiC = [];
for m = 1:15  % number of cells 
    for k = 2:12  % number of input 
        if k ==2
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
for i = 1:1:10 % number of PYR cells
    ipsc_bc = mean(allcells{i}(:,4));
    IPSC_BC = [IPSC_BC ipsc_bc];
end 

IPSC_BC = array2table(IPSC_BC);
IPSC_BC.Properties.VariableNames = {'pyr1'...
    'pyr2' 'pyr3' 'pyr4' 'pyr5' 'pyr6' 'pyr7' 'pyr8' 'pyr9' 'pyr10'};

%% all current onto PYR gathered



%% Raster Plot for one neuron

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
