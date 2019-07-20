%% Calculate All Excitatory/Inhibitory Ratios onto BCs and AACs
%  Melisa Gumus
%  May 2018

%% Load Data From Netclamp Results
clear all
close all
clc

g = fullfile('~/','Documents',...
    'labs', 'skinnerlab','networkclamp_results','Sc10_SimD4000_DegSt048_40nodes_conn900_gpyrpyr30_gecca3pyr0p16_G9',...
    'BC',{...
    'BC_33281';...
    'BC_33333';...
    'BC_33554';...
    'BC_33775'...
    },{...
    'mytrace_33281_syns.dat';...
    'mytrace_33333_syns.dat';...
    'mytrace_33554_syns.dat';...
    'mytrace_33775_syns.dat'...
    });

f = fullfile('~/','Documents',...
    'labs', 'skinnerlab','networkclamp_results','Sc10_SimD4000_DegSt048_40nodes_conn900_gpyrpyr30_gecca3pyr0p16_G9',...
    'AAC',{...
    'AAC_24';...
    'AAC_60';...
    'AAC_102';...
    'AAC_111'...
    },{...
    'mytrace_24_syns.dat';...
    'mytrace_60_syns.dat';...
    'mytrace_102_syns.dat';...
    'mytrace_111_syns.dat'...
    });

h = fullfile('~/','Documents',...
    'labs', 'skinnerlab','networkclamp_results','Sc10_SimD4000_DegSt048_40nodes_conn900_gpyrpyr30_gecca3pyr0p16_G9',...
    'BiC',{...
    'BiC_147';...
    'BiC_157';...
    'BiC_232';...
    'BiC_337'...
    },{...
    'mytrace_147_syns.dat';...
    'mytrace_157_syns.dat';...
    'mytrace_232_syns.dat';...
    'mytrace_337_syns.dat'...
    });


%% Write Data on Matrix
alldataBC = [];
for m = 1:1:4
    temp_data = readtable(g{m},'Delimiter','\t');
    temp_data = table2array(temp_data);
    alldataBC = [alldataBC temp_data];
end

dataBC = mat2cell(alldataBC, 80000, ...
    [13, 13, 13, 13]);

alldataAAC = [];
for m = 1:1:4
    temp_data = readtable(f{m},'Delimiter','\t');
    temp_data = table2array(temp_data);
    alldataAAC = [alldataAAC temp_data];
end

dataAAC = mat2cell(alldataAAC, 80000, ...
    [13, 13, 13, 13]);

alldataBiC = [];
for m = 1:1:4
    temp_data = readtable(h{m},'Delimiter','\t');
    temp_data = table2array(temp_data);
    alldataBiC = [alldataBiC temp_data];
end

dataBiC = mat2cell(alldataBiC, 80000, ...
    [13, 13, 13, 13]);

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

for m = 1:4  % number of cells
    for k = 2:13  % number of input
        if k == 3
            temp_current_BiC = dataBC{m}(:,k);
            [pks, locs] = findpeaks(dataBC{m}(:,k),'MinPeakDistance',3000); % peak detection
            temp_BiC = dataBC{m}(:,k);
            allrows = (1:80000)';
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
            allrows = (1:80000)';
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
            allrows = (1:80000)';
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

allcellsBC = mat2cell(M_on_BC, 80000, ...
    [3, 3, 3, 3]);

%% Creates a big table consists of inputs from BiC, PYR, and BC... onto AAC

M = [];
current_BC_on_AAC = [];
current_PYR_on_AAC = [];
current_BiC_on_AAC = [];
for m = 1:4  % number of cells
    for k = 2:13  % number of input
        if k == 3
            temp_current_BiC = dataAAC{m}(:,k);
            [pks, locs] = findpeaks(dataAAC{m}(:,k),'MinPeakDistance',3000); % peak detection
            temp_BiC = dataAAC{m}(:,k);
            allrows = (1:80000)';
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
            allrows = (1:80000)';
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
            allrows = (1:80000)';
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

allcellsAAC = mat2cell(M, 80000, ...
    [3, 3, 3, 3]);

%% Creates a big table consists of inputs from BiC, PYR, and BC... onto BiC

M_BiC = [];
current_BC_on_BiC = [];
current_PYR_on_BiC = [];
current_BiC_on_BiC = [];
for m = 1:4  % number of cells
    for k = 2:13  % number of input
        if k == 3
            temp_current_BiC = dataBiC{m}(:,k);
            [pks, locs] = findpeaks(dataBiC{m}(:,k),'MinPeakDistance',3000); % peak detection
            temp_BiC = dataBiC{m}(:,k);
            allrows = (1:80000)';
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
            allrows = (1:80000)';
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
            allrows = (1:80000)';
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

allcellsBiC = mat2cell(M_BiC, 80000, ...
    [3, 3, 3, 3]);

%% TRY DIFFERENT COMBINATIONS - IPSCs Only From BC and AAC onto BC and AAC
% Sum all ipsc currents
all_ipsc_on_BC_AAC_combo1 = [];
all_epsc_on_BC_AAC_combo1 = [];

for i = 1:1:4
    for t = 1:1:4

        tot_cur_ipsc_on_BC_combo1 = current_BC_on_BC(:,i);
        tot_cur_ipsc_on_AAC_combo1 = current_BC_on_AAC(:,t);

        tot_cur_ipsc_on_BC_AAC_combo1 = tot_cur_ipsc_on_BC_combo1 + tot_cur_ipsc_on_AAC_combo1;
        all_ipsc_on_BC_AAC_combo1 = [all_ipsc_on_BC_AAC_combo1 tot_cur_ipsc_on_BC_AAC_combo1];

        tot_cur_epsc_on_BC_AAC_combo1 =  current_PYR_on_BC(:,i) + current_PYR_on_AAC(:,t);
        all_epsc_on_BC_AAC_combo1 = [all_epsc_on_BC_AAC_combo1 tot_cur_epsc_on_BC_AAC_combo1];
    end
end
%% Find the peaks of the summed IPSCs from BC and AAC onto BC, AAC - and EPSCs too

peaks_all_PV_on_BC_AAC_combo1 = [];
f1 = figure;
for k = 1:1:16 % all combinations (4*4)
    [pks, locs] = findpeaks(all_ipsc_on_BC_AAC_combo1(:,k),'MinPeakDistance',3000); % peak detection
    findpeaks(all_ipsc_on_BC_AAC_combo1(:,k),'MinPeakDistance',3000);

    temp_cur = all_ipsc_on_BC_AAC_combo1(:,k);
    allrows = (1:80000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    peaks_all = temp_cur;
    peaks_all_PV_on_BC_AAC_combo1 = [peaks_all_PV_on_BC_AAC_combo1 peaks_all];
end

%%
peaks_PYR_on_BC_AAC_combo1 = [];
for k = 1:1:16
    [pks, locs] = findpeaks(-all_epsc_on_BC_AAC_combo1(:,k),'MinPeakDistance',3000); % peak detection
    findpeaks(-all_epsc_on_BC_AAC_combo1(:,k),'MinPeakDistance',3000);

    temp_cur = all_epsc_on_BC_AAC_combo1(:,k);
    allrows = (1:80000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    peaks_all = temp_cur;
    peaks_PYR_on_BC_AAC_combo1 = [peaks_PYR_on_BC_AAC_combo1 peaks_all];
end

%% IPSCs from BC and AAC onto BC and AAC

IPSC_all_on_BC_AAC_combo1 = [];
ipsc_all_on_BC_AAC_combo1 = [];
for i = 1:1:16
    pks_ipsc_all = peaks_all_PV_on_BC_AAC_combo1(:,i);
    pks_ipsc_all(pks_ipsc_all == 0) = [];
    ipsc_all_mean = mean(pks_ipsc_all);
    ipsc_all_std = std(pks_ipsc_all);
    ipsc_all_on_BC_AAC_combo1 = [ipsc_all_mean;ipsc_all_std];
    IPSC_all_on_BC_AAC_combo1 = [IPSC_all_on_BC_AAC_combo1 ipsc_all_on_BC_AAC_combo1];
end

%% Find Mean EPSC and SD onto BC, AAC

EPSC_on_BC_AAC_combo1 = [];
epsc_on_BC_AAC_combo1 = [];
for i = 1:1:16
    pks_epsc_all = peaks_PYR_on_BC_AAC_combo1(:,i);
    pks_epsc_all(pks_epsc_all == 0) = [];
    epsc_all_mean = mean(pks_epsc_all);
    epsc_all_std = std(pks_epsc_all);
    epsc_on_BC_AAC_combo1 = [epsc_all_mean;epsc_all_std];
    EPSC_on_BC_AAC_combo1 = [EPSC_on_BC_AAC_combo1 epsc_on_BC_AAC_combo1];
end

%% Excitatory/Inhibitory Ratios on BC, AAC population and BC, AAC population

IPSC_BC_AAC_on_BC_AAC_combo1 = IPSC_all_on_BC_AAC_combo1;

Ratios_BC_AAC =[];
temp_ratio_BC_AAC = [];

E_I_BC_AAC = abs(EPSC_on_BC_AAC_combo1(1,:)./IPSC_BC_AAC_on_BC_AAC_combo1(1,:))';

temp_ratio_BC_AAC = [temp_ratio_BC_AAC E_I_BC_AAC];
Ratios_BC_AAC = array2table(reshape(temp_ratio_BC_AAC,[4,4]));

%% E/I Ratio - BC, AAC to BC, AAC
fig = uitable('Data',Ratios_BC_AAC{:,:},...
    'RowName',[],...
    'ColumnName',[],...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);
%%

less_than_1 = sum(temp_ratio_BC_AAC >1)/16;
