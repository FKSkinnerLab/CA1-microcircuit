% Load Data 
clear all
clc 
 
f = fullfile('/home','melisagumus','Documents', ...
    'MATLAB','CA1_SimTracker','BiC',{...
    'BiC_1470_1000';...
    'BiC_1580_1000';...
    'BiC_1635_1000';...
    'BiC_1855_1000';...
    'BiC_2020_1000';...
    'BiC_2130_1000';...
    'BiC_2350_1000';...
    'BiC_2460_1000';...
    'BiC_2570_1000';...
    'BiC_2900_1000';...
    'BiC_3010_1000';...
    'BiC_3175_1000';...
    'BiC_3340_1000';...
    'BiC_3505_1000';...
    'BiC_3615_1000'...
    },{...
    'mytrace_1470_syns.dat';...
    'mytrace_1580_syns.dat';...
    'mytrace_1635_syns.dat';...
    'mytrace_1855_syns.dat';...
    'mytrace_2020_syns.dat';...
    'mytrace_2130_syns.dat';...
    'mytrace_2350_syns.dat';...
    'mytrace_2460_syns.dat';...
    'mytrace_2570_syns.dat';...
    'mytrace_2900_syns.dat';...
    'mytrace_3010_syns.dat';...
    'mytrace_3175_syns.dat';...
    'mytrace_3340_syns.dat';...
    'mytrace_3505_syns.dat';...
    'mytrace_3615_syns.dat'...
    });

%% Alternative - does not work
clear all 
clc

start_path = fullfile(matlabroot, '/home','melisagumus','Documents', ...
    'MATLAB','CA1_SimTracker','BiC');
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
%cd C:\Users\Melisa\Documents\MATLAB\CA1_SimTracker\AAC\
%files = dir('AAC*1000\mytrace*syns.dat');
%pyr_names = f.names;

alldata = [];
for m = 1:1:15
    temp_data = readtable(f{m},'Delimiter','\t');
    temp_data = table2array(temp_data);
    alldata = [alldata temp_data];
end

data = mat2cell(alldata, 40000, ...
    [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]);

%% Creates a big table consists of inputs from BiC PYR and BC in order

M = [];
current_BC = [];
current_PYR = [];
current_BiC = [];
current_cck = [];
current_ivy = [];
current_ngf = [];
current_olm = [];
current_sca = [];
for m = 1:15  % number of cells 
    for k = 2:12  % number of input 
        if k == 3
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
        elseif k == 4
            temp_current_cck = data{m}(:,k);
        elseif k == 5
            temp_current_ivy = data{m}(:,k);
        elseif k == 6
            temp_current_ngf = data{m}(:,k);
        elseif k == 7
            temp_current_olm = data{m}(:,k);
        elseif k == 10
            temp_current_sca = data{m}(:,k);
        end 
    end
    current_PYR = [current_PYR temp_current_PYR];
    current_BiC = [current_BiC temp_current_BiC];
    current_BC = [current_BC temp_current_BC];
    current_cck = [current_cck temp_current_cck];
    current_ivy = [current_ivy temp_current_ivy];
    current_ngf = [current_ngf temp_current_ngf];
    current_olm = [current_olm temp_current_olm];
    current_sca = [current_sca temp_current_sca];
    M = [M BiC PYR BC];
end 

allcells = mat2cell(M, 40000, ...
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]); 

%% Graph EPSCs on BiC

figure 
for i = 1:1:15
    temp = allcells{i}(:,2);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(5,3,i)
    histfit(temp)
    hold on
    title (['BiC' num2str(i)]) % name each graph with the corresponding pyr #
end 
suptitle('EPSC Distribution on BiC')

%% Graph IPSCs from BiC onto BiC

figure 
for i = 1:1:15
    temp = allcells{i}(:,1);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(5,3,i)  
    histfit(temp,50)
    %xlim([-2 2])
    hold on 
    title (['BiC' num2str(i)]) % name each graph with the corresponding pyr #
end 
suptitle('IPSC Distribution from BiC onto BiC')

%% Graph IPSCs from BC onto BiC

figure 
for i = 1:1:15
    temp = allcells{i}(:,3);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(5,3,i)  
    histfit(temp,50)
    %xlim([-2 2])
    hold on 
    title (['AAC' num2str(i)]) % name each graph with the corresponding pyr #
end 
suptitle('IPSC Distribution from BC onto BiC')


%% EPSCs from PYR

EPSC = [];
epsc = [];
for i = 1:1:15 % number of PYR cells
    pks_epsc = allcells{i}(:,2);
    pks_epsc(pks_epsc==0)=[];
    epsc_mean = mean(pks_epsc);
    epsc_std = std(pks_epsc);
    epsc = [epsc_mean;epsc_std];
    EPSC = [EPSC epsc];
end 

EPSC = array2table(EPSC);
EPSC.Properties.VariableNames = {'BiC1'...
    'BiC2' 'BiC3' 'BiC4' 'BiC5' 'BiC6'...
    'BiC7' 'BiC8' 'BiC9' 'BiC10' 'BiC11'...
    'BiC12' 'BiC13' 'BiC14' 'BiC15'};

subplot(3,1,1)
EPSC_mean = table2array(EPSC(1,:));
EPSC_std = table2array(EPSC(2,:));
x = linspace(0,15,length(EPSC_mean));
scatter(x,EPSC_mean,'black','filled');
xlabel('Individual BiCs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak EPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,EPSC_mean,EPSC_std,'b','LineStyle','none')
title('Mean Peak EPSC onto BiC','FontSize',15,'FontWeight','bold')

%%
uitable('Data',EPSC{:,:},'ColumnName',EPSC.Properties.VariableNames,...
    'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

%% IPSCs only from BiC

IPSC_BiC = [];
ipsc_BiC = [];
for i = 1:1:15 % number of PYR cells
    pks_ipsc_BiC = allcells{i}(:,1);
    pks_ipsc_BiC(pks_ipsc_BiC==0)=[];
    ipsc_BiC_mean = mean(pks_ipsc_BiC);
    ipsc_BiC_std = std(pks_ipsc_BiC);
    ipsc_BiC = [ipsc_BiC_mean;ipsc_BiC_std];
    IPSC_BiC = [IPSC_BiC ipsc_BiC];
end 

IPSC_BiC = array2table(IPSC_BiC);
IPSC_BiC.Properties.VariableNames = {'BiC1'...
    'BiC2' 'BiC3' 'BiC4' 'BiC5' 'BiC6'...
    'BiC7' 'BiC8' 'BiC9' 'BiC10' 'BiC11'...
    'BiC12' 'BiC13' 'BiC14' 'BiC15'};

subplot(3,1,2)
IPSC_BiC_mean = table2array(IPSC_BiC(1,:));
IPSC_BiC_std = table2array(IPSC_BiC(2,:));
x = linspace(0,15,length(IPSC_BiC_mean));
scatter(x,IPSC_BiC_mean,'black','filled');
xlabel('Individual BiCs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC from BiC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BiC_mean,IPSC_BiC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BiC onto BiC','FontSize',15,'FontWeight','bold')

%%
uitable('Data',IPSC_BiC{:,:},'ColumnName',IPSC_BiC.Properties.VariableNames,...
    'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

%% IPSCs only from BC 

IPSC_BC = [];
ipsc_BC = [];
for i = 1:1:15 % number of PYR cells
    pks_ipsc_BC = allcells{i}(:,3);
    pks_ipsc_BC(pks_ipsc_BC==0)=[];
    ipsc_BC_mean = mean(pks_ipsc_BC);
    ipsc_BC_std = std(pks_ipsc_BC);
    ipsc_BC = [ipsc_BC_mean;ipsc_BC_std];
    IPSC_BC = [IPSC_BC ipsc_BC];
end 

IPSC_BC = array2table(IPSC_BC);
IPSC_BC.Properties.VariableNames = {'BiC1'...
    'BiC2' 'BiC3' 'BiC4' 'BiC5' 'BiC6'...
    'BiC7' 'BiC8' 'BiC9' 'BiC10' 'BiC11'...
    'BiC12' 'BiC13' 'BiC14' 'BiC15'};

subplot(3,1,3)
IPSC_BC_mean = table2array(IPSC_BC(1,:));
IPSC_BC_std = table2array(IPSC_BC(2,:));
x = linspace(0,15,length(IPSC_BC_mean));
scatter(x,IPSC_BC_mean,'black','filled');
xlabel('Individual BiCs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC from BiC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BC_mean,IPSC_BC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BC onto BiC','FontSize',15,'FontWeight','bold')
%%
uitable('Data',IPSC_BC{:,:},'ColumnName',IPSC_BC.Properties.VariableNames,...
    'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

%% BC and BiC ipsc current onto BiC gathered

% Sum all ipsc currents
all_ipsc = [];
all_ipsc_together = [];

for i = 1:1:15
    tot_cur_ipsc =  current_BiC(:,i) + current_BC(:,i);
    tot_ipsc_together = current_BiC(:,i)...
        +current_BC(:,i)...
        +current_cck(:,i)...
        +current_ivy(:,i)...
        +current_ngf(:,i)...
        +current_olm(:,i)...
        +current_sca(:,i);
    
    all_ipsc = [all_ipsc tot_cur_ipsc];
    all_ipsc_together = [all_ipsc_together tot_ipsc_together];

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

peaks_all_PV_together = [];
for k = 1:1:15
    [pks, locs] = findpeaks(all_ipsc_together(:,k),'MinPeakDistance',4000); % peak detection
    temp_cur_together = all_ipsc_together(:,k);
    allrows = (1:40000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    peaks_all_together = temp_cur_together;
    peaks_all_PV_together = [peaks_all_PV_together peaks_all_together];
end

%% BC and BiC ipsc currents together onto BiC - graph and table

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
IPSC_all.Properties.VariableNames = {'BiC1'...
    'BiC2' 'BiC3' 'BiC4' 'BiC5' 'BiC6'...
    'BiC7' 'BiC8' 'BiC9' 'BiC10' 'BiC11'...
    'BiC12' 'BiC13' 'BiC14' 'BiC15'};
 
IPSC_all_mean = table2array(IPSC_all(1,:));
IPSC_all_std = table2array(IPSC_all(2,:));
x = linspace(0,15,length(IPSC_all_mean));
figure
scatter(x,IPSC_all_mean,'black','filled');
xlabel('Individual BiCs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_all_mean,IPSC_all_std,'b','LineStyle','none')
title('Mean Peak IPSC from BC and BiC onto BiC','FontSize',15,'FontWeight','bold')

%%
uitable('Data',IPSC_all{:,:},'ColumnName',IPSC_all.Properties.VariableNames,...
    'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

%% All ipsc currents together onto BC - graph and table

IPSC_all_together = [];
ipsc_all_together = [];
for i = 1:1:15 % number of PYR cells
    pks_ipsc_all_together = peaks_all_PV_together(:,i);
    pks_ipsc_all_together(pks_ipsc_all_together == 0) = [];
    ipsc_all_mean_together = mean(pks_ipsc_all_together);
    ipsc_all_std_together = std(pks_ipsc_all_together);
    ipsc_all_together = [ipsc_all_mean_together;ipsc_all_std_together];
    IPSC_all_together = [IPSC_all_together ipsc_all_together];
end 

IPSC_all_together = array2table(IPSC_all_together);
IPSC_all_together.Properties.VariableNames = {'BiC1'...
    'BiC2' 'BiC3' 'BiC4' 'BiC5' 'BiC6'...
    'BiC7' 'BiC8' 'BiC9' 'BiC10' 'BiC11'...
    'BiC12' 'BiC13' 'BiC14' 'BiC15'};
 
IPSC_all_mean_together = table2array(IPSC_all_together(1,:));
IPSC_all_std_together = table2array(IPSC_all_together(2,:));
x = linspace(1,15,length(IPSC_all_mean_together));
figure
scatter(x,IPSC_all_mean_together,'black','filled');
axis([x], [0 16])
xlabel('Individual BiC','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_all_mean_together,IPSC_all_std_together,'b','LineStyle','none')
title('Mean Peak IPSC from all inhibitory cells onto BiC','FontSize',15,'FontWeight','bold')

%%
uitable('Data',IPSC_all_together{:,:},'ColumnName',IPSC_all_together.Properties.VariableNames,...
    'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

%% E/I Ratios on BiC Cells

IPSC_BiC = table2array(IPSC_BiC);
IPSC_BC = table2array(IPSC_BC);
IPSC_all= table2array(IPSC_all);  % all refers to BC and BiC together
IPSC_all_together= table2array(IPSC_all_together);

EPSC = table2array(EPSC);

%%
Ratios_BiC = [];
E_I_BC = abs(EPSC(1,:)./IPSC_BC(1,:))';
E_I_BiC = abs(EPSC(1,:)./IPSC_BiC(1,:))';
E_I_all = abs(EPSC(1,:)./IPSC_all(1,:))'; 
E_I_all_together = abs(EPSC(1,:)./IPSC_all_together(1,:))';


%%
aac = 1:15;
Ratios_BiC = [Ratios_BiC aac' E_I_BC E_I_BiC E_I_all E_I_all_together];
Ratios_BiC = array2table(Ratios_BiC);
 
Ratios_BiC.Properties.VariableNames = {'BiC_no' 'Ratio_BC_on_BiC'...
    'Ratio_BiC_on_BiC' 'Ratio_BC_BiC_on_BiC' 'All_ipsc_onto_BiC'};

%% Display the table as a figure

uitable('Data',Ratios_BiC{:,:},'ColumnName',Ratios_BiC.Properties.VariableNames,...
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
