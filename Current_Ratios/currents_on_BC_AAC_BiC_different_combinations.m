%% Calculate Excitatory/Inhibitory Ratios onto BCs, BiCs and AACs
%  Melisa Gumus
%  May 2018 

%% Load Data From Netclamp Results
clear all
close all
clc
g = fullfile('/home','melisagumus','Documents', ...
    'MATLAB','CA1_SimTracker','pvbasket',{...
    'pvbasket_332810_1000';...
    'pvbasket_333500_1000';...
    'pvbasket_333776_1000';...
    'pvbasket_334466_1000';...
    'pvbasket_335018_1000';...
    'pvbasket_335432_1000';...
    'pvbasket_335846_1000';...
    'pvbasket_336260_1000';...
    'pvbasket_332948_1000';...
    'pvbasket_333086_1000';...
    'pvbasket_333224_1000';...
    'pvbasket_333638_1000';...
    'pvbasket_333914_1000';...
    'pvbasket_338192_1000';...
    'pvbasket_338054_1000'...
    },{...
    'mytrace_332810_syns.dat';...
    'mytrace_333500_syns.dat';...
    'mytrace_333776_syns.dat';...
    'mytrace_334466_syns.dat';...
    'mytrace_335018_syns.dat';...
    'mytrace_335432_syns.dat';...
    'mytrace_335846_syns.dat';...
    'mytrace_336260_syns.dat';...
    'mytrace_332948_syns.dat';...
    'mytrace_333086_syns.dat';...
    'mytrace_333224_syns.dat';...
    'mytrace_333638_syns.dat';...
    'mytrace_333914_syns.dat';...
    'mytrace_338192_syns.dat';...
    'mytrace_338054_syns.dat'...
    });

f = fullfile('/home','melisagumus','Documents', ...
    'MATLAB','CA1_SimTracker','AAC',{...
    'AAC_0_1000';...
    'AAC_36_1000';...
    'AAC_180_1000';...
    'AAC_288_1000';...
    'AAC_360_1000';...
    'AAC_468_1000';...
    'AAC_576_1000';...
    'AAC_720_1000';...
    'AAC_828_1000';...
    'AAC_900_1000';...
    'AAC_1008_1000';...
    'AAC_1152_1000';...
    'AAC_1224_1000';...
    'AAC_1332_1000';...
    'AAC_1404_1000'...
    },{...
    'mytrace_0_syns.dat';...
    'mytrace_36_syns.dat';...
    'mytrace_180_syns.dat';...
    'mytrace_288_syns.dat';...
    'mytrace_360_syns.dat';...
    'mytrace_468_syns.dat';...
    'mytrace_576_syns.dat';...
    'mytrace_720_syns.dat';...
    'mytrace_828_syns.dat';...
    'mytrace_900_syns.dat';...
    'mytrace_1008_syns.dat';...
    'mytrace_1152_syns.dat';...
    'mytrace_1224_syns.dat';...
    'mytrace_1332_syns.dat';...
    'mytrace_1404_syns.dat'...
    });

h = fullfile('/home','melisagumus','Documents', ...
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

%% Load Data From Netclamp Results
clear all
close all
clc

g = fullfile('/Users','macklabadmin','Documents','other','pvbasket',{...
    'pvbasket_332810_1000';...
    'pvbasket_333500_1000';...
    'pvbasket_333776_1000';...
    'pvbasket_334466_1000';...
    'pvbasket_335018_1000';...
    'pvbasket_335432_1000';...
    'pvbasket_335846_1000';...
    'pvbasket_336260_1000';...
    'pvbasket_332948_1000';...
    'pvbasket_333086_1000';...
    'pvbasket_333224_1000';...
    'pvbasket_333638_1000';...
    'pvbasket_333914_1000';...
    'pvbasket_338192_1000';...
    'pvbasket_338054_1000'...
    },{...
    'mytrace_332810_syns.dat';...
    'mytrace_333500_syns.dat';...
    'mytrace_333776_syns.dat';...
    'mytrace_334466_syns.dat';...
    'mytrace_335018_syns.dat';...
    'mytrace_335432_syns.dat';...
    'mytrace_335846_syns.dat';...
    'mytrace_336260_syns.dat';...
    'mytrace_332948_syns.dat';...
    'mytrace_333086_syns.dat';...
    'mytrace_333224_syns.dat';...
    'mytrace_333638_syns.dat';...
    'mytrace_333914_syns.dat';...
    'mytrace_338192_syns.dat';...
    'mytrace_338054_syns.dat'...
    });

f = fullfile('/Users','macklabadmin','Documents','other','AAC',{...
    'AAC_0_1000';...
    'AAC_36_1000';...
    'AAC_180_1000';...
    'AAC_288_1000';...
    'AAC_360_1000';...
    'AAC_468_1000';...
    'AAC_576_1000';...
    'AAC_720_1000';...
    'AAC_828_1000';...
    'AAC_900_1000';...
    'AAC_1008_1000';...
    'AAC_1152_1000';...
    'AAC_1224_1000';...
    'AAC_1332_1000';...
    'AAC_1404_1000'...
    },{...
    'mytrace_0_syns.dat';...
    'mytrace_36_syns.dat';...
    'mytrace_180_syns.dat';...
    'mytrace_288_syns.dat';...
    'mytrace_360_syns.dat';...
    'mytrace_468_syns.dat';...
    'mytrace_576_syns.dat';...
    'mytrace_720_syns.dat';...
    'mytrace_828_syns.dat';...
    'mytrace_900_syns.dat';...
    'mytrace_1008_syns.dat';...
    'mytrace_1152_syns.dat';...
    'mytrace_1224_syns.dat';...
    'mytrace_1332_syns.dat';...
    'mytrace_1404_syns.dat'...
    });

h = fullfile('/Users','macklabadmin','Documents','other','BiC',{...
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


%% Write Data on Matrix
alldataBC = [];
for m = 1:1:15
    temp_data = readtable(g{m},'Delimiter','\t');
    temp_data = table2array(temp_data);
    alldataBC = [alldataBC temp_data];
end

dataBC = mat2cell(alldataBC, 40000, ...
    [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]);

alldataAAC = [];
for m = 1:1:15
    temp_data = readtable(f{m},'Delimiter','\t');
    temp_data = table2array(temp_data);
    alldataAAC = [alldataAAC temp_data];
end

dataAAC = mat2cell(alldataAAC, 40000, ...
    [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]);

alldataBiC = [];
for m = 1:1:15
    temp_data = readtable(h{m},'Delimiter','\t');
    temp_data = table2array(temp_data);
    alldataBiC = [alldataBiC temp_data];
end

dataBiC = mat2cell(alldataBiC, 40000, ...
    [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]);

%% Creates a big table consists of inputs from BiC, PYR, and BC... onto BC

M_on_BC = [];
current_BC_on_BC = [];
current_PYR_on_BC = [];
current_BiC_on_BC = [];
current_cck = [];
current_ivy = [];
current_ngf = [];
current_olm = [];
current_sca = [];

for m = 1:15  % number of cells 
    for k = 2:12  % number of input 
        if k == 3
            temp_current_BiC = dataBC{m}(:,k);
            [pks, locs] = findpeaks(dataBC{m}(:,k),'MinPeakDistance',3000); % peak detection
            temp_BiC = dataBC{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_BiC(element,:) = 0;
            end 
            BiC_on_BC = temp_BiC;
        elseif k == 8
            temp_current_PYR = dataBC{m}(:,k);
            [pks, locs] = findpeaks(-dataBC{m}(:,k),'MinPeakDistance',3000); % peak detection
            temp_PYR = dataBC{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_PYR(element,:) = 0;
            end 
            PYR_on_BC = temp_PYR;
        elseif k == 9
            temp_current_BC = dataBC{m}(:,k);
            [pks, locs] = findpeaks(dataBC{m}(:,k),'MinPeakDistance',3000); % peak detection
            %findpeaks(data{m}(:,k),M
            temp_BC = dataBC{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_BC(element,:) = 0;
            end 
            BC_on_BC = temp_BC;
        elseif k == 4
            temp_current_cck = dataBC{m}(:,k);
        elseif k == 5
            temp_current_ivy = dataBC{m}(:,k);
        elseif k == 6
            temp_current_ngf = dataBC{m}(:,k);
        elseif k == 7
            temp_current_olm = dataBC{m}(:,k);
        elseif k == 10
            temp_current_sca = dataBC{m}(:,k);
        
        end 
        
    end
    current_PYR_on_BC = [current_PYR_on_BC temp_current_PYR];
    current_BiC_on_BC = [current_BiC_on_BC temp_current_BiC];
    current_BC_on_BC = [current_BC_on_BC temp_current_BC];
    current_cck_on_BC = [current_cck temp_current_cck];
    current_ivy_on_BC = [current_ivy temp_current_ivy];
    current_ngf_on_BC = [current_ngf temp_current_ngf];
    current_olm_on_BC = [current_olm temp_current_olm];
    current_sca_on_BC = [current_sca temp_current_sca];
    
    M_on_BC = [M_on_BC BiC_on_BC PYR_on_BC BC_on_BC];
    
end 

allcellsBC = mat2cell(M_on_BC, 40000, ...
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]); 

%% Creates a big table consists of inputs from BiC, PYR, and BC... onto AAC

M = [];
current_BC_on_AAC = [];
current_PYR_on_AAC = [];
current_BiC_on_AAC = [];
for m = 1:15  % number of cells 
    for k = 2:12  % number of input 
        if k == 3
            temp_current_BiC = dataAAC{m}(:,k);
            [pks, locs] = findpeaks(dataAAC{m}(:,k),'MinPeakDistance',3000); % peak detection
            temp_BiC = dataAAC{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_BiC(element,:) = 0;
            end 
            BiC_on_AAC = temp_BiC;
        elseif k == 8
            temp_current_PYR = dataAAC{m}(:,k);
            [pks, locs] = findpeaks(-dataAAC{m}(:,k),'MinPeakDistance',3000); % peak detection
            temp_PYR = dataAAC{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_PYR(element,:) = 0;
            end 
            PYR_on_AAC = temp_PYR;
        elseif k == 9
            temp_current_BC = dataAAC{m}(:,k);
            [pks, locs] = findpeaks(dataAAC{m}(:,k),'MinPeakDistance',3000); % peak detection
            %findpeaks(data{m}(:,k),M
            temp_BC = dataAAC{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_BC(element,:) = 0;
            end 
            BC_on_AAC = temp_BC;
            elseif k == 4
                temp_current_cck = dataAAC{m}(:,k);
            elseif k == 5
                temp_current_ivy = dataAAC{m}(:,k);
            elseif k == 6
                temp_current_ngf = dataAAC{m}(:,k);
            elseif k == 7
                temp_current_olm = dataAAC{m}(:,k);
            elseif k == 10
                temp_current_sca = dataAAC{m}(:,k); 
        end 
    end
    current_PYR_on_AAC = [current_PYR_on_AAC temp_current_PYR];
    current_BiC_on_AAC = [current_BiC_on_AAC temp_current_BiC];
    current_BC_on_AAC = [current_BC_on_AAC temp_current_BC];
    M = [M BiC_on_AAC PYR_on_AAC BC_on_AAC];
    current_cck_on_AAC = [current_cck temp_current_cck];
    current_ivy_on_AAC = [current_ivy temp_current_ivy];
    current_ngf_on_AAC = [current_ngf temp_current_ngf];
    current_olm_on_AAC = [current_olm temp_current_olm];
    current_sca_on_AAC = [current_sca temp_current_sca];
end 

allcellsAAC = mat2cell(M, 40000, ...
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]); 

%% Creates a big table consists of inputs from BiC, PYR, and BC... onto BiC

M_BiC = [];
current_BC_on_BiC = [];
current_PYR_on_BiC = [];
current_BiC_on_BiC = [];
for m = 1:15  % number of cells 
    for k = 2:12  % number of input 
        if k == 3
            temp_current_BiC = dataBiC{m}(:,k);
            [pks, locs] = findpeaks(dataBiC{m}(:,k),'MinPeakDistance',3000); % peak detection
            temp_BiC = dataBiC{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_BiC(element,:) = 0;
            end 
            BiC_on_BiC = temp_BiC;
        elseif k == 8
            temp_current_PYR = dataBiC{m}(:,k);
            [pks, locs] = findpeaks(-dataBiC{m}(:,k),'MinPeakDistance',3000); % peak detection
            temp_PYR = dataBiC{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_PYR(element,:) = 0;
            end 
            PYR_on_BiC = temp_PYR;
        elseif k == 9
            temp_current_BC = dataBiC{m}(:,k);
            [pks, locs] = findpeaks(dataBiC{m}(:,k),'MinPeakDistance',3000); % peak detection
            %findpeaks(data{m}(:,k),M
            temp_BC = dataBiC{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_BC(element,:) = 0;
            end
            BC_on_BiC = temp_BC;
            elseif k == 4
                temp_current_cck = dataAAC{m}(:,k);
            elseif k == 5
                temp_current_ivy = dataAAC{m}(:,k);
            elseif k == 6
                temp_current_ngf = dataAAC{m}(:,k);
            elseif k == 7
                temp_current_olm = dataAAC{m}(:,k);
            elseif k == 10
                temp_current_sca = dataAAC{m}(:,k); 
        end 
    end
    current_PYR_on_BiC = [current_PYR_on_BiC temp_current_PYR];
    current_BiC_on_BiC = [current_BiC_on_BiC temp_current_BiC];
    current_BC_on_BiC = [current_BC_on_BiC temp_current_BC];
    M_BiC = [M_BiC BiC_on_BiC PYR_on_BiC BC_on_BiC];
    current_cck_on_BiC = [current_cck temp_current_cck];
    current_ivy_on_BiC = [current_ivy temp_current_ivy];
    current_ngf_on_BiC = [current_ngf temp_current_ngf];
    current_olm_on_BiC = [current_olm temp_current_olm];
    current_sca_on_BiC = [current_sca temp_current_sca];
end 

allcellsBiC = mat2cell(M_BiC, 40000, ...
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]); 

%% Find Mean EPSCs and SD from PYR onto BC

EPSC_BC = [];
epsc_BC = [];
for i = 1:1:15 
    pks_epsc = allcellsBC{i}(:,2);
    pks_epsc(pks_epsc==0)=[];
    epsc_mean = mean(pks_epsc);
    epsc_std = std(pks_epsc);
    epsc_BC = [epsc_mean;epsc_std];
    EPSC_BC = [EPSC_BC epsc_BC];
end 

EPSC_BC_table = EPSC_BC;
num = (1:15)';
EPSC_BC_table = array2table(EPSC_BC_table');
EPSC_BC_table.num = num;
EPSC_BC_table = [EPSC_BC_table(:,end) EPSC_BC_table(:,1) EPSC_BC_table(:,2)];

EPSC_BC_table.Properties.VariableNames = {'BC_Number', 'Mean_Peak', 'Standard_Deviation'};

subplot(3,1,1)
EPSC_BC_mean = EPSC_BC(1,:);
EPSC_BC_std = EPSC_BC(2,:);
x = linspace(0,14,length(EPSC_BC_mean));
scatter(x,EPSC_BC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dEPSC_BC_mean=0.1;
text(x+dx, EPSC_BC_mean+dEPSC_BC_mean, c);
xlabel('Individual BCs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak EPSCs','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,EPSC_BC_mean,EPSC_BC_std,'b','LineStyle','none')
title('Mean Peak EPSCs onto BCs','FontSize',15,'FontWeight','bold')


%% Find Mean EPSCs and SD onto AAC

EPSC_on_AAC = [];
epsc_on_AAC = [];
for i = 1:1:15
    pks_epsc = allcellsAAC{i}(:,2);
    pks_epsc(pks_epsc==0)=[];
    epsc_on_AAC_mean = mean(pks_epsc);
    epsc_on_AAC_std = std(pks_epsc);
    epsc_on_AAC = [epsc_on_AAC_mean;epsc_on_AAC_std];
    EPSC_on_AAC = [EPSC_on_AAC epsc_on_AAC];
end 

EPSC_on_AAC_table = EPSC_on_AAC;
num = (1:15)';
EPSC_on_AAC_table = array2table(EPSC_on_AAC_table');
EPSC_on_AAC_table.num = num;
EPSC_on_AAC_table = [EPSC_on_AAC_table(:,end) EPSC_on_AAC_table(:,1) EPSC_on_AAC_table(:,2)];

EPSC_on_AAC_table.Properties.VariableNames = {'AAC_Number', 'Mean_Peak', 'Standard_Deviation'};

subplot(3,1,2)
EPSC_on_AAC_mean = EPSC_on_AAC(1,:);
EPSC_on_AAC_std = EPSC_on_AAC(2,:);
x = linspace(0,14,length(EPSC_on_AAC_mean));
scatter(x,EPSC_on_AAC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dEPSC_on_AAC_mean=0.1;
text(x+dx, EPSC_on_AAC_mean+dEPSC_on_AAC_mean, c);
xlabel('Individual AACs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak EPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,EPSC_on_AAC_mean,EPSC_on_AAC_std,'b','LineStyle','none')
title('Mean Peak EPSC onto AACs','FontSize',15,'FontWeight','bold')

%% Find Mean EPSCs and SD onto BiC

EPSC_on_BiC = [];
epsc_on_BiC = [];
for i = 1:1:15
    pks_epsc = allcellsBiC{i}(:,2);
    pks_epsc(pks_epsc==0)=[];
    epsc_on_BiC_mean = mean(pks_epsc);
    epsc_on_BiC_std = std(pks_epsc);
    epsc_on_BiC = [epsc_on_BiC_mean;epsc_on_BiC_std];
    EPSC_on_BiC = [EPSC_on_BiC epsc_on_BiC];
end 

EPSC_on_BiC_table = EPSC_on_BiC;
num = (1:15)';
EPSC_on_BiC_table = array2table(EPSC_on_BiC_table');
EPSC_on_BiC_table.num = num;
EPSC_on_BiC_table = [EPSC_on_BiC_table(:,end) EPSC_on_BiC_table(:,1) EPSC_on_BiC_table(:,2)];

EPSC_on_BiC_table.Properties.VariableNames = {'BiC_Number', 'Mean_Peak', 'Standard_Deviation'};

subplot(3,1,3)
EPSC_on_BiC_mean = EPSC_on_BiC(1,:);
EPSC_on_BiC_std = EPSC_on_BiC(2,:);
x = linspace(0,14,length(EPSC_on_BiC_mean));
scatter(x,EPSC_on_BiC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dEPSC_on_BiC_mean=0.1;
text(x+dx, EPSC_on_BiC_mean+dEPSC_on_BiC_mean, c);
xlabel('Individual BiCs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak EPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,EPSC_on_BiC_mean,EPSC_on_BiC_std,'b','LineStyle','none')
title('Mean Peak EPSC onto BiCs','FontSize',15,'FontWeight','bold')


%% IPSCs Only from BiC onto BC

IPSC_BiC_on_BC = [];
ipsc_BiC_on_BC = [];
for i = 1:1:15
    pks_ipsc_BiC = allcellsBC{i}(:,1);
    pks_ipsc_BiC(pks_ipsc_BiC==0)=[];
    ipsc_BiC_mean = mean(pks_ipsc_BiC);
    ipsc_BiC_std = std(pks_ipsc_BiC);
    ipsc_BiC_on_BC = [ipsc_BiC_mean;ipsc_BiC_std];
    IPSC_BiC_on_BC = [IPSC_BiC_on_BC ipsc_BiC_on_BC];
end 

IPSC_BiC_on_BC_table = IPSC_BiC_on_BC;
num = (1:15)';
IPSC_BiC_on_BC_table = array2table(IPSC_BiC_on_BC_table');
IPSC_BiC_on_BC_table.num = num;
IPSC_BiC_on_BC_table = [IPSC_BiC_on_BC_table(:,end) IPSC_BiC_on_BC_table(:,1) IPSC_BiC_on_BC_table(:,2)];

IPSC_BiC_on_BC_table.Properties.VariableNames = {'BC_Number', 'Mean_Peak', 'Standard_Deviation'};


subplot(3,1,1)
IPSC_BiC_on_BC_mean = IPSC_BiC_on_BC(1,:);
IPSC_BiC_on_BC_std = IPSC_BiC_on_BC(2,:);
x = linspace(0,14,length(IPSC_BiC_on_BC_mean));
scatter(x,IPSC_BiC_on_BC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_BiC_on_BC_mean=0.1;
text(x+dx, IPSC_BiC_on_BC_mean+dIPSC_BiC_on_BC_mean, c);
xlabel('Individual BCs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BiC_on_BC_mean,IPSC_BiC_on_BC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BiCs onto BCs','FontSize',15,'FontWeight','bold')

%% IPSCs only from BiC on AAC

IPSC_BiC_on_AAC = [];
ipsc_BiC_on_AAC = [];
for i = 1:1:15
    pks_ipsc_BiC = allcellsAAC{i}(:,1);
    pks_ipsc_BiC(pks_ipsc_BiC==0)=[];
    ipsc_BiC_on_AAC_mean = mean(pks_ipsc_BiC);
    ipsc_BiC_on_AAC_std = std(pks_ipsc_BiC);
    ipsc_BiC_on_AAC = [ipsc_BiC_on_AAC_mean;ipsc_BiC_on_AAC_std];
    IPSC_BiC_on_AAC = [IPSC_BiC_on_AAC ipsc_BiC_on_AAC];
end 

IPSC_BiC_on_AAC_table = IPSC_BiC_on_AAC;
num = (1:15)';
IPSC_BiC_on_AAC_table = array2table(IPSC_BiC_on_AAC_table');
IPSC_BiC_on_AAC_table.num = num;
IPSC_BiC_on_AAC_table = [IPSC_BiC_on_AAC_table(:,end) IPSC_BiC_on_AAC_table(:,1) IPSC_BiC_on_AAC_table(:,2)];

IPSC_BiC_on_AAC_table.Properties.VariableNames = {'AAC_Number', 'Mean_Peak', 'Standard_Deviation'};

subplot(3,1,2)
IPSC_BiC_on_AAC_mean = IPSC_BiC_on_AAC(1,:);
IPSC_BiC_on_AAC_std = IPSC_BiC_on_AAC(2,:);
x = linspace(0,14,length(IPSC_BiC_on_AAC_mean));
scatter(x,IPSC_BiC_on_AAC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_BiC_on_AAC_mean=0.1;
text(x+dx, IPSC_BiC_on_AAC_mean+dIPSC_BiC_on_AAC_mean, c);
xlabel('Individual AACs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC from BiC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BiC_on_AAC_mean,IPSC_BiC_on_AAC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BiCs onto AACs','FontSize',15,'FontWeight','bold')

%% IPSCs only from BiC on BiC

IPSC_BiC_on_BiC = [];
ipsc_BiC_on_BiC = [];
for i = 1:1:15
    pks_ipsc_BiC = allcellsBiC{i}(:,1);
    pks_ipsc_BiC(pks_ipsc_BiC==0)=[];
    ipsc_BiC_on_BiC_mean = mean(pks_ipsc_BiC);
    ipsc_BiC_on_BiC_std = std(pks_ipsc_BiC);
    ipsc_BiC_on_BiC = [ipsc_BiC_on_BiC_mean;ipsc_BiC_on_BiC_std];
    IPSC_BiC_on_BiC = [IPSC_BiC_on_BiC ipsc_BiC_on_BiC];
end 

IPSC_BiC_on_BiC_table = IPSC_BiC_on_BiC;
num = (1:15)';
IPSC_BiC_on_BiC_table = array2table(IPSC_BiC_on_BiC_table');
IPSC_BiC_on_BiC_table.num = num;
IPSC_BiC_on_BiC_table = [IPSC_BiC_on_BiC_table(:,end) IPSC_BiC_on_BiC_table(:,1) IPSC_BiC_on_BiC_table(:,2)];

IPSC_BiC_on_BiC_table.Properties.VariableNames = {'BiC_Number', 'Mean_Peak', 'Standard_Deviation'};

subplot(3,1,3)
IPSC_BiC_on_BiC_mean = IPSC_BiC_on_BiC(1,:);
IPSC_BiC_on_BiC_std = IPSC_BiC_on_BiC(2,:);
x = linspace(0,14,length(IPSC_BiC_on_BiC_mean));
scatter(x,IPSC_BiC_on_BiC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_BiC_on_BiC_mean=0.1;
text(x+dx, IPSC_BiC_on_BiC_mean+dIPSC_BiC_on_BiC_mean, c);
xlabel('Individual BiCs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC from BiCs','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BiC_on_BiC_mean,IPSC_BiC_on_BiC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BiCs onto BiCs','FontSize',15,'FontWeight','bold')

%% IPSCs only from BC on BC

IPSC_BC_on_BC = [];
ipsc_BC_on_BC = [];
for i = 1:1:15
    pks_ipsc_BC = allcellsBC{i}(:,3);
    pks_ipsc_BC(pks_ipsc_BC==0)=[];
    ipsc_BC_on_BC_mean = mean(pks_ipsc_BC);
    ipsc_BC_on_BC_std = std(pks_ipsc_BC);
    ipsc_BC_on_BC = [ipsc_BC_on_BC_mean;ipsc_BC_on_BC_std];
    IPSC_BC_on_BC = [IPSC_BC_on_BC ipsc_BC_on_BC];
end 

IPSC_BC_on_BC_table = IPSC_BC_on_BC;
num = (1:15)';
IPSC_BC_on_BC_table = array2table(IPSC_BC_on_BC_table');
IPSC_BC_on_BC_table.num = num;
IPSC_BC_on_BC_table = [IPSC_BC_on_BC_table(:,end) IPSC_BC_on_BC_table(:,1) IPSC_BC_on_BC_table(:,2)];

IPSC_BC_on_BC_table.Properties.VariableNames = {'BC_Number', 'Mean_Peak', 'Standard_Deviation'};

subplot(3,1,1)
IPSC_BC_on_BC_mean = IPSC_BC_on_BC(1,:);
IPSC_BC_on_BC_std = IPSC_BC_on_BC(2,:);
x = linspace(0,14,length(IPSC_BC_on_BC_mean));
scatter(x,IPSC_BC_on_BC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_BC_on_BC_mean=0.1;
text(x+dx, IPSC_BC_on_BC_mean+dIPSC_BC_on_BC_mean, c);
xlabel('Individual BCs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BC_on_BC_mean,IPSC_BC_on_BC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BCs onto BCs','FontSize',15,'FontWeight','bold')

%% IPSCs only from BC on AAC

IPSC_BC_on_AAC = [];
ipsc_BC_on_AAC = [];
for i = 1:1:15
    pks_ipsc_BC = allcellsAAC{i}(:,3);
    pks_ipsc_BC(pks_ipsc_BC==0)=[];
    ipsc_BC_on_AAC_mean = mean(pks_ipsc_BC);
    ipsc_BC_on_AAC_std = std(pks_ipsc_BC);
    ipsc_BC_on_AAC = [ipsc_BC_on_AAC_mean;ipsc_BC_on_AAC_std];
    IPSC_BC_on_AAC = [IPSC_BC_on_AAC ipsc_BC_on_AAC];
end 

IPSC_BC_on_AAC_table = IPSC_BC_on_AAC;
num = (1:15)';
IPSC_BC_on_AAC_table = array2table(IPSC_BC_on_AAC_table');
IPSC_BC_on_AAC_table.num = num;
IPSC_BC_on_AAC_table = [IPSC_BC_on_AAC_table(:,end) IPSC_BC_on_AAC_table(:,1) IPSC_BC_on_AAC_table(:,2)];

IPSC_BC_on_AAC_table.Properties.VariableNames = {'AAC_Number', 'Mean_Peak', 'Standard_Deviation'};

subplot(3,1,2)
IPSC_BC_on_AAC_mean = IPSC_BC_on_AAC(1,:);
IPSC_BC_on_AAC_std = IPSC_BC_on_AAC(2,:);
x = linspace(0,15,length(IPSC_BC_on_AAC_mean));
scatter(x,IPSC_BC_on_AAC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_BC_on_AAC_mean=0.1;
text(x+dx, IPSC_BC_on_AAC_mean+dIPSC_BC_on_AAC_mean, c);
xlabel('Individual AACs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BC_on_AAC_mean,IPSC_BC_on_AAC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BCs onto AACs','FontSize',15,'FontWeight','bold')

%% IPSCs only from BC on BiC

IPSC_BC_on_BiC = [];
ipsc_BC_on_BiC = [];
for i = 1:1:15
    pks_ipsc_BC = allcellsBiC{i}(:,3);
    pks_ipsc_BC(pks_ipsc_BC==0)=[];
    ipsc_BC_on_BiC_mean = mean(pks_ipsc_BC);
    ipsc_BC_on_BiC_std = std(pks_ipsc_BC);
    ipsc_BC_on_BiC = [ipsc_BC_on_BiC_mean;ipsc_BC_on_BiC_std];
    IPSC_BC_on_BiC = [IPSC_BC_on_BiC ipsc_BC_on_BiC];
end 

IPSC_BC_on_BiC_table = IPSC_BC_on_BiC;
num = (1:15)';
IPSC_BC_on_BiC_table = array2table(IPSC_BC_on_BiC_table');
IPSC_BC_on_BiC_table.num = num;
IPSC_BC_on_BiC_table = [IPSC_BC_on_BiC_table(:,end) IPSC_BC_on_BiC_table(:,1) IPSC_BC_on_BiC_table(:,2)];

IPSC_BC_on_BiC_table.Properties.VariableNames = {'BiC_Number', 'Mean_Peak', 'Standard_Deviation'};

subplot(3,1,3)
IPSC_BC_on_BiC_mean = IPSC_BC_on_BiC(1,:);
IPSC_BC_on_BiC_std = IPSC_BC_on_BiC(2,:);
x = linspace(0,14,length(IPSC_BC_on_BiC_mean));
scatter(x,IPSC_BC_on_BiC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_BC_on_BiC_mean=0.1;
text(x+dx, IPSC_BC_on_BiC_mean+dIPSC_BC_on_BiC_mean, c);
xlabel('Individual BiCs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BC_on_BiC_mean,IPSC_BC_on_BiC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BCs onto BiCs','FontSize',15,'FontWeight','bold')

%% TRY DIFFERENT COMBINATIONS - IPSCs Only From BC and BiC onto BC, BiC and AAC gathered 
% Sum all ipsc currents
all_ipsc_on_BC_BiC_combo1 = [];
all_epsc_on_BC_BiC_combo1 = [];
all_ipsc_on_BC_BiC_AAC_combo1 = [];
all_ipsc_on_BC_AAC_combo1 = [];
all_epsc_on_BC_BiC_AAC_combo1 = [];
all_epsc_on_BC_AAC_combo1 = [];

for i = 1:1:15
    if i == 15
        t = 1;
    elseif 15 > i 
        t = i +1;
    end 
        tot_cur_ipsc_on_BC_combo1 =  current_BiC_on_BC(:,i) + current_BC_on_BC(:,i);
        tot_cur_ipsc_on_BiC_combo1 =  current_BiC_on_BiC(:,t) + current_BC_on_BiC(:,t);
        tot_cur_ipsc_on_AAC_combo1 =  current_BiC_on_AAC(:,t) + current_BC_on_AAC(:,t);
        
        tot_cur_ipsc_on_BC_BiC_combo1 = tot_cur_ipsc_on_BC_combo1 + tot_cur_ipsc_on_BiC_combo1;
        tot_cur_ipsc_on_BC_BiC_AAC_combo1 = tot_cur_ipsc_on_BC_combo1 + tot_cur_ipsc_on_AAC_combo1 + tot_cur_ipsc_on_BiC_combo1;
        tot_cur_ipsc_on_BC_AAC_combo1 = tot_cur_ipsc_on_BC_combo1 + tot_cur_ipsc_on_AAC_combo1;
    
        all_ipsc_on_BC_BiC_combo1 = [all_ipsc_on_BC_BiC_combo1 tot_cur_ipsc_on_BC_BiC_combo1];
        all_ipsc_on_BC_BiC_AAC_combo1 = [all_ipsc_on_BC_BiC_AAC_combo1 tot_cur_ipsc_on_BC_BiC_AAC_combo1];
        all_ipsc_on_BC_AAC_combo1 = [all_ipsc_on_BC_AAC_combo1 tot_cur_ipsc_on_BC_AAC_combo1];
    
     % EPSC
        tot_cur_epsc_on_BC_BiC_AAC_combo1 =  current_PYR_on_BC(:,i) + current_PYR_on_BiC(:,t) + current_PYR_on_AAC(:,t);
        all_epsc_on_BC_BiC_AAC_combo1 = [all_epsc_on_BC_BiC_AAC_combo1 tot_cur_epsc_on_BC_BiC_AAC_combo1];
    
        tot_cur_epsc_on_BC_AAC_combo1 =  current_PYR_on_BC(:,i) + current_PYR_on_AAC(:,t);
        all_epsc_on_BC_AAC_combo1 = [all_epsc_on_BC_AAC_combo1 tot_cur_epsc_on_BC_AAC_combo1];
    
        tot_cur_epsc_on_BC_BiC_combo1 =  current_PYR_on_BC(:,i) + current_PYR_on_BiC(:,t);
        all_epsc_on_BC_BiC_combo1 = [all_epsc_on_BC_BiC_combo1 tot_cur_epsc_on_BC_BiC_combo1];
end 

%% Find the peaks of the summed IPSCs from BC and BiC onto BC, BiC and AAC - and EPSCs too

peaks_all_PV_on_BC_BiC_combo1 = [];
f1 = figure;
for k = 1:1:15
    figure(f1);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Peak Detection on IPSCs from BCs, BiCs on BCs, BiC'];
    subplot(5,3,k);
    [pks, locs] = findpeaks(all_ipsc_on_BC_BiC_combo1(:,k),'MinPeakDistance',3000); % peak detection
    findpeaks(all_ipsc_on_BC_BiC_combo1(:,k),'MinPeakDistance',3000);
    hold on; 
    title (['BC, BiC Number #' num2str(k)])
    xlabel('Time (1/40 ms)')
    ylabel('IPSCs')    
    temp_cur = all_ipsc_on_BC_BiC_combo1(:,k);
    allrows = (1:40000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    peaks_all = temp_cur;
    peaks_all_PV_on_BC_BiC_combo1 = [peaks_all_PV_on_BC_BiC_combo1 peaks_all];
end


peaks_all_PV_on_BC_AAC_combo1 = [];
f2 = figure;
for k = 1:1:15
    figure(f2);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Peak Detection on IPSCs from BCs, AACs on BCs, AACs'];
    subplot(5,3,k);
    [pks, locs] = findpeaks(all_ipsc_on_BC_AAC_combo1(:,k),'MinPeakDistance',3000); % peak detection
    findpeaks(all_ipsc_on_BC_AAC_combo1(:,k),'MinPeakDistance',3000);
    hold on; 
    title (['BC, AAC Number #' num2str(k)])
    xlabel('Time (1/40 ms)')
    ylabel('IPSCs')     
    temp_cur = all_ipsc_on_BC_AAC_combo1(:,k);
    allrows = (1:40000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    peaks_all = temp_cur;
    peaks_all_PV_on_BC_AAC_combo1 = [peaks_all_PV_on_BC_AAC_combo1 peaks_all];
end

% Find the peaks of the summed ipsc currents
peaks_all_PV_on_BC_BiC_AAC_combo1 = [];
f3 = figure;
for k = 1:1:15
    figure(f3);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Peak Detection on IPSCs from BCs, BiCs, AACs onto BCs, BiC, AACs'];
    subplot(5,3,k);   
    [pks, locs] = findpeaks(all_ipsc_on_BC_BiC_AAC_combo1(:,k),'MinPeakDistance',3000); % peak detection
    findpeaks(all_ipsc_on_BC_BiC_AAC_combo1(:,k),'MinPeakDistance',3000);
    hold on; 
    title (['BC, BiC, AAC Number #' num2str(k)])
    xlabel('Time (1/40 ms)')
    ylabel('IPSCs from All Inhibitory Cells')   
    temp_cur = all_ipsc_on_BC_BiC_AAC_combo1(:,k);
    allrows = (1:40000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    peaks_all = temp_cur;
    peaks_all_PV_on_BC_BiC_AAC_combo1 = [peaks_all_PV_on_BC_BiC_AAC_combo1 peaks_all];
end

peaks_PYR_on_BC_BiC_AAC_combo1 = [];
f4 = figure;
for k = 1:1:15
    figure(f4);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Peak Detection on EPSCs onto BCs, BiC, AACs'];
    subplot(5,3,k);    
    [pks, locs] = findpeaks(-all_epsc_on_BC_BiC_AAC_combo1(:,k),'MinPeakDistance',3000); % peak detection
    findpeaks(-all_epsc_on_BC_BiC_AAC_combo1(:,k),'MinPeakDistance',3000);
    hold on; 
    title (['BC, BiC, AAC Number #' num2str(k)])
    xlabel('Time (1/40 ms)')
    ylabel('EPSCs')
    temp_cur = all_epsc_on_BC_BiC_AAC_combo1(:,k);
    allrows = (1:40000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    peaks_all = temp_cur;
    peaks_PYR_on_BC_BiC_AAC_combo1 = [peaks_PYR_on_BC_BiC_AAC_combo1 peaks_all];
end

peaks_PYR_on_BC_AAC_combo1 = [];
f5 = figure;
for k = 1:1:15
    figure(f5);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Peak Detection on EPSCs onto BCs, AACs'];
    subplot(5,3,k);    
    [pks, locs] = findpeaks(-all_epsc_on_BC_AAC_combo1(:,k),'MinPeakDistance',3000); % peak detection
    findpeaks(-all_epsc_on_BC_AAC_combo1(:,k),'MinPeakDistance',3000);
    hold on; 
    title (['BC, AAC Number #' num2str(k)])
    xlabel('Time (1/40 ms)')
    ylabel('EPSCs')
    temp_cur = all_epsc_on_BC_AAC_combo1(:,k);
    allrows = (1:40000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    peaks_all = temp_cur;
    peaks_PYR_on_BC_AAC_combo1 = [peaks_PYR_on_BC_AAC_combo1 peaks_all];
end

peaks_PYR_on_BC_BiC_combo1 = [];
f6 = figure;
for k = 1:1:15
    figure(f6);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Peak Detection on EPSCs onto BCs, BiCs'];
    subplot(5,3,k);
    [pks, locs] = findpeaks(-all_epsc_on_BC_BiC_combo1(:,k),'MinPeakDistance',3000); % peak detection
    findpeaks(-all_epsc_on_BC_BiC_combo1(:,k),'MinPeakDistance',3000);
    hold on; 
    title (['BC, BiC Number #' num2str(k)])
    xlabel('Time (1/40 ms)')
    ylabel('EPSCs ')
    temp_cur = all_epsc_on_BC_BiC_combo1(:,k);
    allrows = (1:40000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    peaks_all = temp_cur;
    peaks_PYR_on_BC_BiC_combo1 = [peaks_PYR_on_BC_BiC_combo1 peaks_all];
end
%% IPSCs from BC and BiC onto BC and BiC 

IPSC_all_on_BC_BiC_combo1 = [];
ipsc_all_on_BC_BiC_combo1 = [];
for i = 1:1:15
    pks_ipsc_all = peaks_all_PV_on_BC_BiC_combo1(:,i);
    pks_ipsc_all(pks_ipsc_all == 0) = [];
    ipsc_all_mean = mean(pks_ipsc_all);
    ipsc_all_std = std(pks_ipsc_all);
    ipsc_all_on_BC_BiC_combo1 = [ipsc_all_mean;ipsc_all_std];
    IPSC_all_on_BC_BiC_combo1 = [IPSC_all_on_BC_BiC_combo1 ipsc_all_on_BC_BiC_combo1];
end 

IPSC_all_on_BC_BiC_table = IPSC_all_on_BC_BiC_combo1;
num = (1:15)';
IPSC_all_on_BC_BiC_table = array2table(IPSC_all_on_BC_BiC_table');
IPSC_all_on_BC_BiC_table.num = num;
IPSC_all_on_BC_BiC_table = [IPSC_all_on_BC_BiC_table(:,end) IPSC_all_on_BC_BiC_table(:,1) IPSC_all_on_BC_BiC_table(:,2)];

IPSC_all_on_BC_BiC_table.Properties.VariableNames = {'BC_and_BiC_Number', 'Mean_Peak', 'Standard_Deviation'};

IPSC_all_on_BC_BiC_mean = IPSC_all_on_BC_BiC_combo1(1,:);
IPSC_all_on_BC_BiC_std = IPSC_all_on_BC_BiC_combo1(2,:);
x = linspace(0,15,length(IPSC_all_on_BC_BiC_mean));
figure
scatter(x,IPSC_all_on_BC_BiC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_all_on_BC_BiC_mean=0.1;
text(x+dx, IPSC_all_on_BC_BiC_mean+dIPSC_all_on_BC_BiC_mean, c);
xlabel('Each Point Includes 1 BC and 1 BiC ','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_all_on_BC_BiC_mean,IPSC_all_on_BC_BiC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BC, BiC onto BC, BiC','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSC from BC, BiC onto BC, BiC
fig = uitable('Data',IPSC_all_on_BC_BiC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'1 BC, 1 BiC Number #','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]); %% Display the table as a figure

%% IPSCs from BC and BiC onto BC, BiC and AAC - graph and table 

IPSC_all_on_BC_BiC_AAC_combo1 = [];
ipsc_all_on_BC_BiC_AAC_combo1 = [];
for i = 1:1:15
    pks_ipsc_all = peaks_all_PV_on_BC_BiC_AAC_combo1(:,i);
    pks_ipsc_all(pks_ipsc_all == 0) = [];
    ipsc_all_mean = mean(pks_ipsc_all);
    ipsc_all_std = std(pks_ipsc_all);
    ipsc_all_on_BC_BiC_AAC_combo1 = [ipsc_all_mean;ipsc_all_std];
    IPSC_all_on_BC_BiC_AAC_combo1 = [IPSC_all_on_BC_BiC_AAC_combo1 ipsc_all_on_BC_BiC_AAC_combo1];
end 

IPSC_all_on_BC_BiC_AAC_table = IPSC_all_on_BC_BiC_AAC_combo1;
num = (1:15)';
IPSC_all_on_BC_BiC_AAC_table = array2table(IPSC_all_on_BC_BiC_AAC_table');
IPSC_all_on_BC_BiC_AAC_table.num = num;
IPSC_all_on_BC_BiC_AAC_table = [IPSC_all_on_BC_BiC_AAC_table(:,end) IPSC_all_on_BC_BiC_AAC_table(:,1) IPSC_all_on_BC_BiC_AAC_table(:,2)];

IPSC_all_on_BC_BiC_AAC_table.Properties.VariableNames = {'BC_BiC_AAC_Number', 'Mean_Peak', 'Standard_Deviation'};

IPSC_all_on_BC_BiC_AAC_mean = IPSC_all_on_BC_BiC_AAC_combo1(1,:);
IPSC_all_on_BC_BiC_AAC_std = IPSC_all_on_BC_BiC_AAC_combo1(2,:);
x = linspace(0,14,length(IPSC_all_on_BC_BiC_AAC_mean));
figure
scatter(x,IPSC_all_on_BC_BiC_AAC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_all_on_BC_BiC_AAC_mean=0.1;
text(x+dx, IPSC_all_on_BC_BiC_AAC_mean+dIPSC_all_on_BC_BiC_AAC_mean, c);
xlabel('Each Point Includes 1 BC, 1 BiC, 1 AAC','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_all_on_BC_BiC_AAC_mean,IPSC_all_on_BC_BiC_AAC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BC, BiC and BC, BiC, AAC','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSC from BC, BiC onto BC, BiC, AAC

fig = uitable('Data',IPSC_all_on_BC_BiC_AAC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'1 BC, 1 BiC, 1 AAC Number #','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);%% Display the table as a figure


%% IPSCs from BC and AAC onto BC, AAC - graph and table 

IPSC_all_on_BC_AAC_combo1 = [];
ipsc_all_on_BC_AAC_combo1 = [];
for i = 1:1:15
    pks_ipsc_all = peaks_all_PV_on_BC_AAC_combo1(:,i);
    pks_ipsc_all(pks_ipsc_all == 0) = [];
    ipsc_all_mean = mean(pks_ipsc_all);
    ipsc_all_std = std(pks_ipsc_all);
    ipsc_all_on_BC_AAC_combo1 = [ipsc_all_mean;ipsc_all_std];
    IPSC_all_on_BC_AAC_combo1 = [IPSC_all_on_BC_AAC_combo1 ipsc_all_on_BC_AAC_combo1];
end 

IPSC_all_on_BC_AAC_table = IPSC_all_on_BC_AAC_combo1;
num = (1:15)';
IPSC_all_on_BC_AAC_table = array2table(IPSC_all_on_BC_AAC_table');
IPSC_all_on_BC_AAC_table.num = num;
IPSC_all_on_BC_AAC_table = [IPSC_all_on_BC_AAC_table(:,end) IPSC_all_on_BC_AAC_table(:,1) IPSC_all_on_BC_AAC_table(:,2)];

IPSC_all_on_BC_AAC_table.Properties.VariableNames = {'BC_AAC_Number', 'Mean_Peak', 'Standard_Deviation'};


IPSC_all_on_BC_AAC_mean = IPSC_all_on_BC_AAC_combo1(1,:);
IPSC_all_on_BC_AAC_std = IPSC_all_on_BC_AAC_combo1(2,:);
x = linspace(0,14,length(IPSC_all_on_BC_AAC_mean));
figure
scatter(x,IPSC_all_on_BC_AAC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_all_on_BC_AAC_mean=0.1;
text(x+dx, IPSC_all_on_BC_AAC_mean+dIPSC_all_on_BC_AAC_mean, c);
xlabel('Each Point Includes 1 BC, 1 AAC','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_all_on_BC_AAC_mean,IPSC_all_on_BC_AAC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BC, AAC onto BC, AAC','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSC from BC, AAC onto BC, AAC

fig = uitable('Data',IPSC_all_on_BC_AAC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'1 BC, 1 AAC Number #','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);%% Display the table as a figure

%% Find Mean EPSC and SD onto BC, AAC - graph and table 

EPSC_on_BC_AAC_combo1 = [];
epsc_on_BC_AAC_combo1 = [];
for i = 1:1:15
    pks_epsc_all = peaks_PYR_on_BC_AAC_combo1(:,i);
    pks_epsc_all(pks_epsc_all == 0) = [];
    epsc_all_mean = mean(pks_epsc_all);
    epsc_all_std = std(pks_epsc_all);
    epsc_on_BC_AAC_combo1 = [epsc_all_mean;epsc_all_std];
    EPSC_on_BC_AAC_combo1 = [EPSC_on_BC_AAC_combo1 epsc_on_BC_AAC_combo1];
end 

EPSC_on_BC_AAC_table = EPSC_on_BC_AAC_combo1;
num = (1:15)';
EPSC_on_BC_AAC_table = array2table(EPSC_on_BC_AAC_table');
EPSC_on_BC_AAC_table.num = num;
EPSC_on_BC_AAC_table = [EPSC_on_BC_AAC_table(:,end) EPSC_on_BC_AAC_table(:,1) EPSC_on_BC_AAC_table(:,2)];

EPSC_on_BC_AAC_table.Properties.VariableNames = {'BC_AAC_Number', 'Mean_Peak', 'Standard_Deviation'};

EPSC_on_BC_AAC_mean = EPSC_on_BC_AAC_combo1(1,:);
EPSC_on_BC_AAC_std = EPSC_on_BC_AAC_combo1(2,:);
x = linspace(0,14,length(EPSC_on_BC_AAC_mean));
figure
scatter(x,EPSC_on_BC_AAC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dEPSC_on_BC_AAC_mean=0.1;
text(x+dx, EPSC_on_BC_AAC_mean+dEPSC_on_BC_AAC_mean, c);
xlabel('Each Point Includes 1 BC, 1 AAC','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak EPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,EPSC_on_BC_AAC_mean,EPSC_on_BC_AAC_std,'b','LineStyle','none')
title('Mean Peak EPSC onto BC, AAC','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of EPSC on BC, AAC
fig = uitable('Data',EPSC_on_BC_AAC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'1 BC, 1 AAC Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% Find Mean EPSC and SD onto BC, BiC, AAC - graph and table 

EPSC_on_BC_BiC_AAC_combo1 = [];
epsc_on_BC_BiC_AAC_combo1 = [];
for i = 1:1:15
    pks_epsc_all = peaks_PYR_on_BC_BiC_AAC_combo1(:,i);
    pks_epsc_all(pks_epsc_all == 0) = [];
    epsc_all_mean = mean(pks_epsc_all);
    epsc_all_std = std(pks_epsc_all);
    epsc_on_BC_BiC_AAC_combo1 = [epsc_all_mean;epsc_all_std];
    EPSC_on_BC_BiC_AAC_combo1 = [EPSC_on_BC_BiC_AAC_combo1 epsc_on_BC_BiC_AAC_combo1];
end 

EPSC_on_BC_BiC_AAC_table = EPSC_on_BC_BiC_AAC_combo1;
num = (1:15)';
EPSC_on_BC_BiC_AAC_table = array2table(EPSC_on_BC_BiC_AAC_table');
EPSC_on_BC_BiC_AAC_table.num = num;
EPSC_on_BC_BiC_AAC_table = [EPSC_on_BC_BiC_AAC_table(:,end) EPSC_on_BC_BiC_AAC_table(:,1) EPSC_on_BC_BiC_AAC_table(:,2)];

EPSC_on_BC_BiC_AAC_table.Properties.VariableNames = {'BC_BiC_AAC_Number', 'Mean_Peak', 'Standard_Deviation'};

EPSC_on_BC_BiC_AAC_mean = EPSC_on_BC_BiC_AAC_combo1(1,:);
EPSC_on_BC_BiC_AAC_std = EPSC_on_BC_BiC_AAC_combo1(2,:);
x = linspace(0,14,length(EPSC_on_BC_AAC_mean));
figure
scatter(x,EPSC_on_BC_BiC_AAC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dEPSC_on_BC_BiC_AAC_mean=0.1;
text(x+dx, EPSC_on_BC_BiC_AAC_mean+dEPSC_on_BC_BiC_AAC_mean, c);
xlabel('Each Point Includes 1 BC, 1 BiC, 1 AAC','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak EPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,EPSC_on_BC_BiC_AAC_mean,EPSC_on_BC_BiC_AAC_std,'b','LineStyle','none')
title('Mean Peak EPSC onto BC, BiC, AAC','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of EPSC on BC, BiC, AAC
fig = uitable('Data',EPSC_on_BC_BiC_AAC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'1 BC, 1 BiC, 1 AAC Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% Find Mean EPSC and SD onto BC, BiC - graph and table 

EPSC_on_BC_BiC_combo1 = [];
epsc_on_BC_BiC_combo1 = [];
for i = 1:1:15
    pks_epsc_all = peaks_PYR_on_BC_BiC_combo1(:,i);
    pks_epsc_all(pks_epsc_all == 0) = [];
    epsc_all_mean = mean(pks_epsc_all);
    epsc_all_std = std(pks_epsc_all);
    epsc_on_BC_BiC_combo1 = [epsc_all_mean;epsc_all_std];
    EPSC_on_BC_BiC_combo1 = [EPSC_on_BC_BiC_combo1 epsc_on_BC_BiC_combo1];
end 

EPSC_on_BC_BiC_table = EPSC_on_BC_BiC_combo1;
num = (1:15)';
EPSC_on_BC_BiC_table = array2table(EPSC_on_BC_BiC_table');
EPSC_on_BC_BiC_table.num = num;
EPSC_on_BC_BiC_table = [EPSC_on_BC_BiC_table(:,end) EPSC_on_BC_BiC_table(:,1) EPSC_on_BC_BiC_table(:,2)];

EPSC_on_BC_BiC_table.Properties.VariableNames = {'BC_BiC_Number', 'Mean_Peak', 'Standard_Deviation'};

EPSC_on_BC_BiC_mean = EPSC_on_BC_BiC_combo1(1,:);
EPSC_on_BC_BiC_std = EPSC_on_BC_BiC_combo1(2,:);
x = linspace(0,14,length(EPSC_on_BC_BiC_mean));
figure
scatter(x,EPSC_on_BC_BiC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dEPSC_on_BC_BiC_mean=0.1;
text(x+dx, EPSC_on_BC_BiC_mean+dEPSC_on_BC_BiC_mean, c);
xlabel('Each Point Includes 1 BC, 1 BiC','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak EPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,EPSC_on_BC_BiC_mean,EPSC_on_BC_BiC_std,'b','LineStyle','none')
title('Mean Peak EPSC onto BC, BiC','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of EPSC on BC, BiC
fig = uitable('Data',EPSC_on_BC_BiC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'1 BC, 1 BiC Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% Excitatory/Inhibitory Ratios on BC, AAC population and BC, BiC, AAC population, and BC, BiC - Prep

IPSC_BC_BiC_AAC_on_BC_BiC_AAC_combo1 = IPSC_all_on_BC_BiC_AAC_combo1;
IPSC_BC_AAC_on_BC_AAC_combo1 = IPSC_all_on_BC_AAC_combo1;
IPSC_BC_BiC_on_BC_BiC_combo1 = IPSC_all_on_BC_BiC_combo1;


%% Excitatory/Inhibitory Ratios on BC, AAC population and BC, BiC, AAC population 
Ratios_BC_AAC = [];
Ratios_BC_BiC_AAC = [];
Ratios_BC_BiC =[];

E_I_BC_AAC = abs(EPSC_on_BC_AAC_combo1(1,:)./IPSC_BC_AAC_on_BC_AAC_combo1(1,:))';
E_I_BC_BiC_AAC = abs(EPSC_on_BC_BiC_AAC_combo1(1,:)./IPSC_BC_BiC_AAC_on_BC_BiC_AAC_combo1(1,:))';
E_I_BC_BiC = abs(EPSC_on_BC_BiC_combo1(1,:)./IPSC_BC_BiC_on_BC_BiC_combo1(1,:))';

cells = 1:15;
Ratios_BC_AAC = [Ratios_BC_AAC cells' E_I_BC_AAC];
Ratios_BC_AAC = array2table(Ratios_BC_AAC);
Ratios_BC_AAC.Properties.VariableNames = {'BC_AAC_number' 'Ratio_BC_AAC'};

Ratios_BC_BiC_AAC = [Ratios_BC_BiC_AAC cells' E_I_BC_BiC_AAC];
Ratios_BC_BiC_AAC = array2table(Ratios_BC_BiC_AAC);
Ratios_BC_BiC_AAC.Properties.VariableNames = {'BC_BiC_AAC_number' 'Ratio_BC_BiC_AAC'};

Ratios_BC_BiC = [Ratios_BC_BiC cells' E_I_BC_BiC];
Ratios_BC_BiC = array2table(Ratios_BC_BiC);
Ratios_BC_BiC.Properties.VariableNames = {'BC_BiC_number' 'Ratio_BC_BiC'};

%% E/I Ratio - BC, AAC to BC, AAC
fig = uitable('Data',Ratios_BC_AAC{:,:},...
    'RowName',[],...
    'ColumnName',{'1 BC, 1 AAC Number','BC, AAC to BC, AAC'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);%%

%% E/I Ratio -  BC,BiC, AAC to BC, BiC, AAC
fig = uitable('Data',Ratios_BC_BiC_AAC{:,:},...
    'RowName',[],...
    'ColumnName',{'1 BC, 1 BiC, 1 AAC Number','BC, BiC, AAC to BC, BiC, AAC'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% E/I Ratio - BC, BiC to BC, BiC
fig = uitable('Data',Ratios_BC_BiC{:,:},...
    'RowName',[],...
    'ColumnName',{'1 BC, 1 BiC Number','BC, BiC to BC, BiC'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);%%

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
