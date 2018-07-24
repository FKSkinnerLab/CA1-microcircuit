% Load Data 
clear all
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
    'pvbasket_336536_1000';...
    'pvbasket_336674_1000';...
    'pvbasket_337088_1000';...
    'pvbasket_337364_1000';...
    'pvbasket_337640_1000';...
    'pvbasket_338192_1000'...
    },{...
    'mytrace_332810_syns.dat';...
    'mytrace_333500_syns.dat';...
    'mytrace_333776_syns.dat';...
    'mytrace_334466_syns.dat';...
    'mytrace_335018_syns.dat';...
    'mytrace_335432_syns.dat';...
    'mytrace_335846_syns.dat';...
    'mytrace_336260_syns.dat';...
    'mytrace_336536_syns.dat';...
    'mytrace_336674_syns.dat';...
    'mytrace_337088_syns.dat';...
    'mytrace_337364_syns.dat';...
    'mytrace_337640_syns.dat';...
    'mytrace_338192_syns.dat'...
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
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load Data 
clear all
clc 
 
g = fullfile('C:\','Users','Melisa', ...
    'Desktop','Netclamp','BC','pvbasket',{...
    'pvbasket_332810_1000';...
    'pvbasket_333500_1000';...
    'pvbasket_333776_1000';...
    'pvbasket_334466_1000';...
    'pvbasket_335018_1000';...
    'pvbasket_335432_1000';...
    'pvbasket_335846_1000';...
    'pvbasket_336260_1000';...
    'pvbasket_336536_1000';...
    'pvbasket_336674_1000';...
    'pvbasket_337088_1000';...
    'pvbasket_337364_1000';...
    'pvbasket_337640_1000';...
    'pvbasket_338192_1000'...
    },{...
    'mytrace_332810_syns.dat';...
    'mytrace_333500_syns.dat';...
    'mytrace_333776_syns.dat';...
    'mytrace_334466_syns.dat';...
    'mytrace_335018_syns.dat';...
    'mytrace_335432_syns.dat';...
    'mytrace_335846_syns.dat';...
    'mytrace_336260_syns.dat';...
    'mytrace_336536_syns.dat';...
    'mytrace_336674_syns.dat';...
    'mytrace_337088_syns.dat';...
    'mytrace_337364_syns.dat';...
    'mytrace_337640_syns.dat';...
    'mytrace_338192_syns.dat'...
    });

f = fullfile('C:\','Users','Melisa', ...
    'Desktop','Netclamp','AAC','AAC',{...
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

h = fullfile('C:\','Users','Melisa', ...
    'Desktop','Netclamp','BiC','BiC',{...
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
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cd C:\Users\Melisa\Documents\MATLAB\CA1_SimTracker\pyr\
%files = dir('pyr*1000\mytrace*syns.dat');
%pyr_names = f.names;

alldataBC = [];
for m = 1:1:14
    temp_data = readtable(g{m},'Delimiter','\t');
    temp_data = table2array(temp_data);
    alldataBC = [alldataBC temp_data];
end

dataBC = mat2cell(alldataBC, 40000, ...
    [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]);


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

%% Creates a big table consists of inputs from BiC PYR and BC in order %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


M_on_BC = [];
current_BC_on_BC = [];
current_PYR_on_BC = [];
current_BiC_on_BC = [];
current_cck = [];
current_ivy = [];
current_ngf = [];
current_olm = [];
current_sca = []

for m = 1:14  % number of cells 
    for k = 2:12  % number of input 
        if k == 3
            temp_current_BiC = dataBC{m}(:,k);
            [pks, locs] = findpeaks(dataBC{m}(:,k),'MinPeakDistance',4000); % peak detection
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
            [pks, locs] = findpeaks(-dataBC{m}(:,k),'MinPeakDistance',4000); % peak detection
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
            [pks, locs] = findpeaks(dataBC{m}(:,k),'MinPeakDistance',4000); % peak detection
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
    current_PYR_on_BC = [current_PYR_on_BC temp_current_PYR];
    current_BiC_on_BC = [current_BiC_on_BC temp_current_BiC];
    current_BC_on_BC = [current_BC_on_BC temp_current_BC];
    current_cck = [current_cck temp_current_cck];
    current_ivy = [current_ivy temp_current_ivy];
    current_ngf = [current_ngf temp_current_ngf];
    current_olm = [current_olm temp_current_olm];
    current_sca = [current_sca temp_current_sca];
    
    M_on_BC = [M_on_BC BiC_on_BC PYR_on_BC BC_on_BC];
    
end 

allcellsBC = mat2cell(M_on_BC, 40000, ...
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]); 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = [];
current_BC_on_AAC = [];
current_PYR_on_AAC = [];
current_BiC_on_AAC = [];
for m = 1:15  % number of cells 
    for k = 2:12  % number of input 
        if k == 3
            temp_current_BiC = dataAAC{m}(:,k);
            [pks, locs] = findpeaks(dataAAC{m}(:,k),'MinPeakDistance',4000); % peak detection
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
            [pks, locs] = findpeaks(-dataAAC{m}(:,k),'MinPeakDistance',4000); % peak detection
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
            [pks, locs] = findpeaks(dataAAC{m}(:,k),'MinPeakDistance',4000); % peak detection
            %findpeaks(data{m}(:,k),M
            temp_BC = dataAAC{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_BC(element,:) = 0;
            end 
            BC_on_AAC = temp_BC;
        end 
    end
    current_PYR_on_AAC = [current_PYR_on_AAC temp_current_PYR];
    current_BiC_on_AAC = [current_BiC_on_AAC temp_current_BiC];
    current_BC_on_AAC = [current_BC_on_AAC temp_current_BC];
    M = [M BiC_on_AAC PYR_on_AAC BC_on_AAC];
    current_cck = [current_cck temp_current_cck];
    current_ivy = [current_ivy temp_current_ivy];
    current_ngf = [current_ngf temp_current_ngf];
    current_olm = [current_olm temp_current_olm];
    current_sca = [current_sca temp_current_sca];
end 

allcellsAAC = mat2cell(M, 40000, ...
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]); 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M_BiC = [];
current_BC_on_BiC = [];
current_PYR_on_BiC = [];
current_BiC_on_BiC = [];
for m = 1:15  % number of cells 
    for k = 2:12  % number of input 
        if k == 3
            temp_current_BiC = dataBiC{m}(:,k);
            [pks, locs] = findpeaks(dataBiC{m}(:,k),'MinPeakDistance',4000); % peak detection
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
            [pks, locs] = findpeaks(-dataBiC{m}(:,k),'MinPeakDistance',4000); % peak detection
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
            [pks, locs] = findpeaks(dataBiC{m}(:,k),'MinPeakDistance',4000); % peak detection
            %findpeaks(data{m}(:,k),M
            temp_BC = dataBiC{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_BC(element,:) = 0;
            end 
            BC_on_BiC = temp_BC;
        end 
    end
    current_PYR_on_BiC = [current_PYR_on_BiC temp_current_PYR];
    current_BiC_on_BiC = [current_BiC_on_BiC temp_current_BiC];
    current_BC_on_BiC = [current_BC_on_BiC temp_current_BC];
    M_BiC = [M_BiC BiC_on_BiC PYR_on_BiC BC_on_BiC];
end 

allcellsBiC = mat2cell(M_BiC, 40000, ...
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]); 

%% Graph EPSCs on BC

figure 
for i = 1:1:14
    temp = allcellsBC{i}(:,2);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(7,2,i)
    histfit(temp)
    hold on
    title (['Pyr' num2str(i)]) % name each graph with the corresponding pyr #
end 
suptitle('EPSC Distribution on PYR Cells')

%% Graph EPSCs on AAC
figure 
for i = 1:1:15
    temp = allcellsAAC{i}(:,2);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(5,3,i)
    histfit(temp)
    hold on
    title (['AAC' num2str(i)]) % name each graph with the corresponding pyr #
end 
suptitle('EPSC Distribution on AAC')

%% Graph EPSCs on BiC

figure 
for i = 1:1:15
    temp = allcellsBiC{i}(:,2);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(5,3,i)
    histfit(temp)
    hold on
    title (['BiC' num2str(i)]) % name each graph with the corresponding pyr #
end 
suptitle('EPSC Distribution on BiC')

%% Graph IPSCs from BiC onto BC

figure 
for i = 1:1:14
    temp = allcellsBC{i}(:,1);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(7,2,i)  
    histfit(temp,50)
    xlim([-2 2])
    hold on 
    title (['Pyr' num2str(i)]) % name each graph with the corresponding pyr #
end 
suptitle('IPSC Distribution from BiC onto PYR')

%% Graph IPSCs from BC onto BC

figure 
for i = 1:1:14
    temp = allcellsBC{i}(:,3);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(7,2,i)  
    histfit(temp,50)
    xlim([-2 2])
    hold on 
    title (['Pyr' num2str(i)]) % name each graph with the corresponding pyr #
end 
suptitle('IPSC Distribution from BC onto PYR')


%% EPSCs from PYR on BC

EPSC_BC = [];
epsc_BC = [];
for i = 1:1:14 % number of PYR cells
    pks_epsc = allcellsBC{i}(:,2);
    pks_epsc(pks_epsc==0)=[];
    epsc_mean = mean(pks_epsc);
    epsc_std = std(pks_epsc);
    epsc_BC = [epsc_mean;epsc_std];
    EPSC_BC = [EPSC_BC epsc_BC];
end 

EPSC_BC = array2table(EPSC_BC);
EPSC_BC.Properties.VariableNames = {'BC1'...
    'BC2' 'BC3' 'BC4' 'BC5' 'BC6'...
    'BC7' 'BC8' 'BC9' 'BC10' 'BC11'...
    'BC12' 'BC13' 'BC14' 'BC15'};

subplot(3,1,1)
EPSC_mean = table2array(EPSC_BC(1,:));
EPSC_std = table2array(EPSC_BC(2,:));
x = linspace(0,14,length(EPSC_mean));
scatter(x,EPSC_mean,'black','filled');
xlabel('Individual Basket Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak EPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,EPSC_mean,EPSC_std,'b','LineStyle','none')
title('Mean Peak EPSC onto BC','FontSize',15,'FontWeight','bold')

%% EPSCs from PYR on AAC

EPSC_on_AAC = [];
epsc_on_AAC = [];
for i = 1:1:15 % number of PYR cells
    pks_epsc = allcellsAAC{i}(:,2);
    pks_epsc(pks_epsc==0)=[];
    epsc_on_AAC_mean = mean(pks_epsc);
    epsc_on_AAC_std = std(pks_epsc);
    epsc_on_AAC = [epsc_on_AAC_mean;epsc_on_AAC_std];
    EPSC_on_AAC = [EPSC_on_AAC epsc_on_AAC];
end 

EPSC_on_AAC = array2table(EPSC_on_AAC);
EPSC_on_AAC.Properties.VariableNames = {'AAC1'...
    'AAC2' 'AAC3' 'AAC4' 'AAC5' 'AAC6'...
    'AAC7' 'AAC8' 'AAC9' 'AAC10' 'AAC11'...
    'AAC12' 'AAC13' 'AAC14' 'AAC15'};

subplot(3,1,1)
EPSC_mean = table2array(EPSC_on_AAC(1,:));
EPSC_std = table2array(EPSC_on_AAC(2,:));
x = linspace(0,15,length(EPSC_mean));
scatter(x,EPSC_mean,'black','filled');
xlabel('Individual AACs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak EPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,EPSC_mean,EPSC_std,'b','LineStyle','none')
title('Mean Peak EPSC onto AAC','FontSize',15,'FontWeight','bold')

%% EPSCs from PYR on BiC

EPSC_on_BiC = [];
epsc_on_BiC = [];
for i = 1:1:15 % number of PYR cells
    pks_epsc = allcellsBiC{i}(:,2);
    pks_epsc(pks_epsc==0)=[];
    epsc_on_BiC_mean = mean(pks_epsc);
    epsc_on_BiC_std = std(pks_epsc);
    epsc_on_BiC = [epsc_on_BiC_mean;epsc_on_BiC_std];
    EPSC_on_BiC = [EPSC_on_BiC epsc_on_BiC];
end 

EPSC_on_BiC = array2table(EPSC_on_BiC);
EPSC_on_BiC.Properties.VariableNames = {'BiC1'...
    'BiC2' 'BiC3' 'BiC4' 'BiC5' 'BiC6'...
    'BiC7' 'BiC8' 'BiC9' 'BiC10' 'BiC11'...
    'BiC12' 'BiC13' 'BiC14' 'BiC15'};

subplot(3,1,1)
EPSC_on_BiC_mean = table2array(EPSC_on_BiC(1,:));
EPSC_on_BiC_std = table2array(EPSC_on_BiC(2,:));
x = linspace(0,15,length(EPSC_on_BiC_mean));
scatter(x,EPSC_on_BiC_mean,'black','filled');
xlabel('Individual BiCs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak EPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,EPSC_on_BiC_mean,EPSC_on_BiC_std,'b','LineStyle','none')
title('Mean Peak EPSC onto BiC','FontSize',15,'FontWeight','bold')


%% IPSCs only from BiC on BC

IPSC_BiC_on_BC = [];
ipsc_BiC_on_BC = [];
for i = 1:1:14 % number of PYR cells
    pks_ipsc_BiC = allcellsBC{i}(:,1);
    pks_ipsc_BiC(pks_ipsc_BiC==0)=[];
    ipsc_BiC_mean = mean(pks_ipsc_BiC);
    ipsc_BiC_std = std(pks_ipsc_BiC);
    ipsc_BiC_on_BC = [ipsc_BiC_mean;ipsc_BiC_std];
    IPSC_BiC_on_BC = [IPSC_BiC_on_BC ipsc_BiC_on_BC];
end 

IPSC_BiC_on_BC = array2table(IPSC_BiC_on_BC);
IPSC_BiC_on_BC.Properties.VariableNames = {'BC1'...
    'BC2' 'BC3' 'BC4' 'BC5' 'BC6'...
    'BC7' 'BC8' 'BC9' 'BC10' 'BC11'...
    'BC12' 'BC13' 'BC14' 'BC15'};


subplot(3,1,2)
IPSC_BiC_on_BC_mean = table2array(IPSC_BiC_on_BC(1,:));
IPSC_BiC_on_BC_std = table2array(IPSC_BiC_on_BC(2,:));
x = linspace(0,14,length(IPSC_BiC_on_BC_mean));
scatter(x,IPSC_BiC_on_BC_mean,'black','filled');
xlabel('Individual Basket Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC from BiC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BiC_on_BC_mean,IPSC_BiC_on_BC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BiC onto BC','FontSize',15,'FontWeight','bold')
%% IPSCs only from BiC on AAC

IPSC_BiC_on_AAC = [];
ipsc_BiC_on_AAC = [];
for i = 1:1:15 % number of PYR cells
    pks_ipsc_BiC = allcellsAAC{i}(:,1);
    pks_ipsc_BiC(pks_ipsc_BiC==0)=[];
    ipsc_BiC_on_AAC_mean = mean(pks_ipsc_BiC);
    ipsc_BiC_on_AAC_std = std(pks_ipsc_BiC);
    ipsc_BiC_on_AAC = [ipsc_BiC_on_AAC_mean;ipsc_BiC_on_AAC_std];
    IPSC_BiC_on_AAC = [IPSC_BiC_on_AAC ipsc_BiC_on_AAC];
end 

IPSC_BiC_on_AAC = array2table(IPSC_BiC_on_AAC);
IPSC_BiC_on_AAC.Properties.VariableNames = {'AAC1'...
    'AAC2' 'AAC3' 'AAC4' 'AAC5' 'AAC6'...
    'AAC7' 'AAC8' 'AAC9' 'AAC10' 'AAC11'...
    'AAC12' 'AAC13' 'AAC14' 'AAC15'};

subplot(3,1,2)
IPSC_BiC_on_AAC_mean = table2array(IPSC_BiC_on_AAC(1,:));
IPSC_BiC_on_AAC_std = table2array(IPSC_BiC_on_AAC(2,:));
x = linspace(0,15,length(IPSC_BiC_on_AAC_mean));
scatter(x,IPSC_BiC_on_AAC_mean,'black','filled');
xlabel('Individual AACs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC from BiC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BiC_on_AAC_mean,IPSC_BiC_on_AAC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BiC onto AAC','FontSize',15,'FontWeight','bold')
%% IPSCs only from BiC on BiC

IPSC_BiC_on_BiC = [];
ipsc_BiC_on_BiC = [];
for i = 1:1:15 % number of PYR cells
    pks_ipsc_BiC = allcellsBiC{i}(:,1);
    pks_ipsc_BiC(pks_ipsc_BiC==0)=[];
    ipsc_BiC_on_BiC_mean = mean(pks_ipsc_BiC);
    ipsc_BiC_on_BiC_std = std(pks_ipsc_BiC);
    ipsc_BiC_on_BiC = [ipsc_BiC_on_BiC_mean;ipsc_BiC_on_BiC_std];
    IPSC_BiC_on_BiC = [IPSC_BiC_on_BiC ipsc_BiC_on_BiC];
end 

IPSC_BiC_on_BiC = array2table(IPSC_BiC_on_BiC);
IPSC_BiC_on_BiC.Properties.VariableNames = {'BiC1'...
    'BiC2' 'BiC3' 'BiC4' 'BiC5' 'BiC6'...
    'BiC7' 'BiC8' 'BiC9' 'BiC10' 'BiC11'...
    'BiC12' 'BiC13' 'BiC14' 'BiC15'};

subplot(3,1,2)
IPSC_BiC_on_BiC_mean = table2array(IPSC_BiC_on_BiC(1,:));
IPSC_BiC_on_BiC_std = table2array(IPSC_BiC_on_BiC(2,:));
x = linspace(0,15,length(IPSC_BiC_on_BiC_mean));
scatter(x,IPSC_BiC_on_BiC_mean,'black','filled');
xlabel('Individual BiCs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC from BiC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BiC_on_BiC_mean,IPSC_BiC_on_BiC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BiC onto BiC','FontSize',15,'FontWeight','bold')

%% IPSCs only from BC on BC

IPSC_BC_on_BC = [];
ipsc_BC_on_BC = [];
for i = 1:1:14 % number of PYR cells
    pks_ipsc_BC = allcellsBC{i}(:,3);
    pks_ipsc_BC(pks_ipsc_BC==0)=[];
    ipsc_BC_on_BC_mean = mean(pks_ipsc_BC);
    ipsc_BC_on_BC_std = std(pks_ipsc_BC);
    ipsc_BC_on_BC = [ipsc_BC_on_BC_mean;ipsc_BC_on_BC_std];
    IPSC_BC_on_BC = [IPSC_BC_on_BC ipsc_BC_on_BC];
end 

IPSC_BC_on_BC = array2table(IPSC_BC_on_BC);
IPSC_BC_on_BC.Properties.VariableNames = {'BC1'...
    'BC2' 'BC3' 'BC4' 'BC5' 'BC6'...
    'BC7' 'BC8' 'BC9' 'BC10' 'BC11'...
    'BC12' 'BC13' 'BC14' 'BC15'};


subplot(3,1,3)
IPSC_BC_mean = table2array(IPSC_BC_on_BC(1,:));
IPSC_BC_std = table2array(IPSC_BC_on_BC(2,:));
x = linspace(0,14,length(IPSC_BC_mean));
scatter(x,IPSC_BC_mean,'black','filled');
xlabel('Individual Basket Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC from BC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BC_mean,IPSC_BC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BC onto BC','FontSize',15,'FontWeight','bold')

%% IPSCs only from BC on AAC

IPSC_BC_on_AAC = [];
ipsc_BC_on_AAC = [];
for i = 1:1:15 % number of PYR cells
    pks_ipsc_BC = allcellsAAC{i}(:,3);
    pks_ipsc_BC(pks_ipsc_BC==0)=[];
    ipsc_BC_on_AAC_mean = mean(pks_ipsc_BC);
    ipsc_BC_on_AAC_std = std(pks_ipsc_BC);
    ipsc_BC_on_AAC = [ipsc_BC_on_AAC_mean;ipsc_BC_on_AAC_std];
    IPSC_BC_on_AAC = [IPSC_BC_on_AAC ipsc_BC_on_AAC];
end 

IPSC_BC_on_AAC = array2table(IPSC_BC_on_AAC);
IPSC_BC_on_AAC.Properties.VariableNames = {'AAC1'...
    'AAC2' 'AAC3' 'AAC4' 'AAC5' 'AAC6'...
    'AAC7' 'AAC8' 'AAC9' 'AAC10' 'AAC11'...
    'AAC12' 'AAC13' 'AAC14' 'AAC15'};

subplot(3,1,3)
IPSC_BC_on_AAC_mean = table2array(IPSC_BC_on_AAC(1,:));
IPSC_BC_on_AAC_std = table2array(IPSC_BC_on_AAC(2,:));
x = linspace(0,15,length(IPSC_BC_on_AAC_mean));
scatter(x,IPSC_BC_on_AAC_mean,'black','filled');
xlabel('Individual AACs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC from BC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BC_on_AAC_mean,IPSC_BC_on_AAC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BC onto AAC','FontSize',15,'FontWeight','bold')

%% IPSCs only from BC on BiC

IPSC_BC_on_BiC = [];
ipsc_BC_on_BiC = [];
for i = 1:1:15 % number of PYR cells
    pks_ipsc_BC = allcellsBiC{i}(:,3);
    pks_ipsc_BC(pks_ipsc_BC==0)=[];
    ipsc_BC_on_BiC_mean = mean(pks_ipsc_BC);
    ipsc_BC_on_BiC_std = std(pks_ipsc_BC);
    ipsc_BC_on_BiC = [ipsc_BC_on_BiC_mean;ipsc_BC_on_BiC_std];
    IPSC_BC_on_BiC = [IPSC_BC_on_BiC ipsc_BC_on_BiC];
end 

IPSC_BC_on_BiC = array2table(IPSC_BC_on_BiC);
IPSC_BC_on_BiC.Properties.VariableNames = {'BiC1'...
    'BiC2' 'BiC3' 'BiC4' 'BiC5' 'BiC6'...
    'BiC7' 'BiC8' 'BiC9' 'BiC10' 'BiC11'...
    'BiC12' 'BiC13' 'BiC14' 'BiC15'};

subplot(3,1,3)
IPSC_BC_on_BiC_mean = table2array(IPSC_BC_on_BiC(1,:));
IPSC_BC_on_BiC_std = table2array(IPSC_BC_on_BiC(2,:));
x = linspace(0,15,length(IPSC_BC_on_BiC_mean));
scatter(x,IPSC_BC_on_BiC_mean,'black','filled');
xlabel('Individual BiCs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC from BiC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BC_on_BiC_mean,IPSC_BC_on_BiC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BC onto BiC','FontSize',15,'FontWeight','bold')

%% BC and BiC ipsc current onto BC, BiC and AAC gathered %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sum all ipsc currents

all_ipsc_on_BC_BiC = [];
all_epsc_on_BC_BiC = [];
all_ipsc_on_BC_BiC_AAC = [];
all_ipsc_on_BC_AAC = [];
all_epsc_on_BC_BiC_AAC = [];
all_epsc_on_BC_AAC = [];

for i = 1:1:14
    tot_cur_ipsc_on_BC =  current_BiC_on_BC(:,i) + current_BC_on_BC(:,i);
    tot_cur_ipsc_on_BiC =  current_BiC_on_BiC(:,i) + current_BC_on_BiC(:,i);
    tot_cur_ipsc_on_AAC =  current_BiC_on_AAC(:,i) + current_BC_on_AAC(:,i);
    
    tot_cur_ipsc_on_BC_BiC = tot_cur_ipsc_on_BC + tot_cur_ipsc_on_BiC;
    tot_cur_ipsc_on_BC_BiC_AAC = tot_cur_ipsc_on_AAC + tot_cur_ipsc_on_AAC;
    tot_cur_ipsc_on_BC_AAC = tot_cur_ipsc_on_AAC + tot_cur_ipsc_on_AAC;
    
    all_ipsc_on_BC_BiC = [all_ipsc_on_BC_BiC tot_cur_ipsc_on_BC_BiC];
    all_ipsc_on_BC_BiC_AAC = [all_ipsc_on_BC_BiC_AAC tot_cur_ipsc_on_BC_BiC_AAC];
    all_ipsc_on_BC_AAC = [all_ipsc_on_BC_AAC tot_cur_ipsc_on_BC_AAC];
    
    % EPSC
    tot_cur_epsc_on_BC_BiC_AAC =  current_PYR_on_BC(:,i) + current_PYR_on_BiC(:,i) + current_PYR_on_AAC(:,i);
    all_epsc_on_BC_BiC_AAC = [all_epsc_on_BC_BiC_AAC tot_cur_epsc_on_BC_BiC_AAC];
    
    tot_cur_epsc_on_BC_AAC =  current_PYR_on_BC(:,i) + current_PYR_on_AAC(:,i);
    all_epsc_on_BC_AAC = [all_epsc_on_BC_AAC tot_cur_epsc_on_BC_AAC];
end

%% Find the peaks of the summed ipsc currents %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

peaks_all_PV_on_BC_BiC = [];
for k = 1:1:14
    [pks, locs] = findpeaks(all_ipsc_on_BC_BiC(:,k),'MinPeakDistance',4000); % peak detection
    temp_cur = all_ipsc_on_BC_BiC(:,k);
    allrows = (1:40000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    peaks_all = temp_cur;
    peaks_all_PV_on_BC_BiC = [peaks_all_PV_on_BC_BiC peaks_all];
end


peaks_all_PV_on_BC_AAC = [];
for k = 1:1:14
    [pks, locs] = findpeaks(all_ipsc_on_BC_AAC(:,k),'MinPeakDistance',4000); % peak detection
    temp_cur = all_ipsc_on_BC_BiC(:,k);
    allrows = (1:40000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    peaks_all = temp_cur;
    peaks_all_PV_on_BC_AAC = [peaks_all_PV_on_BC_AAC peaks_all];
end

% Find the peaks of the summed ipsc currents
peaks_all_PV_on_BC_BiC_AAC = [];
for k = 1:1:14
    [pks, locs] = findpeaks(all_ipsc_on_BC_BiC_AAC(:,k),'MinPeakDistance',4000); % peak detection
    temp_cur = all_ipsc_on_BC_BiC_AAC(:,k);
    allrows = (1:40000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    peaks_all = temp_cur;
    peaks_all_PV_on_BC_BiC_AAC = [peaks_all_PV_on_BC_BiC_AAC peaks_all];
end

peaks_PYR_on_BC_BiC_AAC = [];
for k = 1:1:14
    [pks, locs] = findpeaks(all_epsc_on_BC_BiC_AAC(:,k),'MinPeakDistance',4000); % peak detection
    temp_cur = all_epsc_on_BC_BiC_AAC(:,k);
    allrows = (1:40000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    peaks_all = temp_cur;
    peaks_PYR_on_BC_BiC_AAC = [peaks_PYR_on_BC_BiC_AAC peaks_all];
end

peaks_PYR_on_BC_AAC = [];
for k = 1:1:14
    [pks, locs] = findpeaks(all_epsc_on_BC_AAC(:,k),'MinPeakDistance',4000); % peak detection
    temp_cur = all_epsc_on_BC_AAC(:,k);
    allrows = (1:40000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    peaks_all = temp_cur;
    peaks_PYR_on_BC_AAC = [peaks_PYR_on_BC_AAC peaks_all];
end

%% BC and BiC ipsc currents together onto BC, BiC - graph and table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IPSC_all_on_BC_BiC = [];
ipsc_all_on_BC_BiC = [];
for i = 1:1:14 % number of PYR cells
    pks_ipsc_all = peaks_all_PV_on_BC_BiC(:,i);
    pks_ipsc_all(pks_ipsc_all == 0) = [];
    ipsc_all_mean = mean(pks_ipsc_all);
    ipsc_all_std = std(pks_ipsc_all);
    ipsc_all_on_BC_BiC = [ipsc_all_mean;ipsc_all_std];
    IPSC_all_on_BC_BiC = [IPSC_all_on_BC_BiC ipsc_all_on_BC_BiC];
end 

IPSC_all_on_BC_BiC = array2table(IPSC_all_on_BC_BiC);
IPSC_all_on_BC_BiC.Properties.VariableNames = {'BC_BiC1'...
    'BC_BiC2' 'BC_BiC3' 'BC_BiC4' 'BC_BiC5' 'BC_BiC6'...
    'BC_BiC7' 'BC_BiC8' 'BC_BiC9' 'BC_BiC10' 'BC_BiC11'...
    'BC_BiC12' 'BC_BiC13' 'BC_BiC14'};

IPSC_all_on_BC_BiC_mean = table2array(IPSC_all_on_BC_BiC(1,:));
IPSC_all_on_BC_BiC_std = table2array(IPSC_all_on_BC_BiC(2,:));
x = linspace(0,14,length(IPSC_all_on_BC_BiC_mean));
figure
scatter(x,IPSC_all_on_BC_BiC_mean,'black','filled');
xlabel('Each Point Includes BC, BiC ','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_all_on_BC_BiC_mean,IPSC_all_on_BC_BiC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BC, BiC onto BC, BiC','FontSize',15,'FontWeight','bold')

%% BC and BiC ipsc currents together onto BC, BiC and AAC - graph and table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


IPSC_all_on_BC_BiC_AAC = [];
ipsc_all_on_BC_BiC_AAC = [];
for i = 1:1:14 % number of PYR cells
    pks_ipsc_all = peaks_all_PV_on_BC_BiC_AAC(:,i);
    pks_ipsc_all(pks_ipsc_all == 0) = [];
    ipsc_all_mean = mean(pks_ipsc_all);
    ipsc_all_std = std(pks_ipsc_all);
    ipsc_all_on_BC_BiC_AAC = [ipsc_all_mean;ipsc_all_std];
    IPSC_all_on_BC_BiC_AAC = [IPSC_all_on_BC_BiC_AAC ipsc_all_on_BC_BiC_AAC];
end 

IPSC_all_on_BC_BiC_AAC = array2table(IPSC_all_on_BC_BiC_AAC);
IPSC_all_on_BC_BiC_AAC.Properties.VariableNames = {'BC_BiC_AAC1'...
    'BC_BiC_AAC2' 'BC_BiC_AAC3' 'BC_BiC_AAC4' 'BC_BiC_AAC5' 'BC_BiC_AAC6'...
    'BC_BiC_AAC7' 'BC_BiC_AAC8' 'BC_BiC_AAC9' 'BC_BiC_AAC10' 'BC_BiC_AAC11'...
    'BC_BiC_AAC12' 'BC_BiC_AAC13' 'BC_BiC_AAC14'};

IPSC_all_on_BC_BiC_AAC_mean = table2array(IPSC_all_on_BC_BiC_AAC(1,:));
IPSC_all_on_BC_BiC_AAC_std = table2array(IPSC_all_on_BC_BiC_AAC(2,:));
x = linspace(0,14,length(IPSC_all_on_BC_BiC_AAC_mean));
figure
scatter(x,IPSC_all_on_BC_BiC_AAC_mean,'black','filled');
xlabel('Each Point Includes BC, BiC, AAC','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_all_on_BC_BiC_AAC_mean,IPSC_all_on_BC_BiC_AAC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BC, BiC, AAC and BC, BiC, AAC','FontSize',15,'FontWeight','bold')

%% BC and BiC ipsc currents together onto BC, AAC - graph and table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
IPSC_all_on_BC_AAC = [];
ipsc_all_on_BC_AAC = [];
for i = 1:1:14 % number of PYR cells
    pks_ipsc_all = peaks_all_PV_on_BC_AAC(:,i);
    pks_ipsc_all(pks_ipsc_all == 0) = [];
    ipsc_all_mean = mean(pks_ipsc_all);
    ipsc_all_std = std(pks_ipsc_all);
    ipsc_all_on_BC_AAC = [ipsc_all_mean;ipsc_all_std];
    IPSC_all_on_BC_AAC = [IPSC_all_on_BC_AAC ipsc_all_on_BC_AAC];
end 

IPSC_all_on_BC_AAC = array2table(IPSC_all_on_BC_AAC);
IPSC_all_on_BC_AAC.Properties.VariableNames = {'BC_AAC1'...
    'BC_AAC2' 'BC_AAC3' 'BC_AAC4' 'BC_AAC5' 'BC_AAC6'...
    'BC_AAC7' 'BC_AAC8' 'BC_AAC9' 'BC_AAC10' 'BC_AAC11'...
    'BC_AAC12' 'BC_AAC13' 'BC_AAC14'};

IPSC_all_on_BC_AAC_mean = table2array(IPSC_all_on_BC_AAC(1,:));
IPSC_all_on_BC_AAC_std = table2array(IPSC_all_on_BC_AAC(2,:));
x = linspace(0,14,length(IPSC_all_on_BC_AAC_mean));
figure
scatter(x,IPSC_all_on_BC_AAC_mean,'black','filled');
xlabel('Each Point Includes BC, AAC','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_all_on_BC_AAC_mean,IPSC_all_on_BC_AAC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BC, AAC onto BC, AAC','FontSize',15,'FontWeight','bold')

%% EPSC currents together onto BC, AAC - graph and table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EPSC_on_BC_AAC = [];
epsc_on_BC_AAC = [];
for i = 1:1:14 % number of PYR cells
    pks_epsc_all = peaks_PYR_on_BC_AAC(:,i);
    pks_epsc_all(pks_epsc_all == 0) = [];
    epsc_all_mean = mean(pks_epsc_all);
    epsc_all_std = std(pks_epsc_all);
    epsc_on_BC_AAC = [epsc_all_mean;ipsc_all_std];
    EPSC_on_BC_AAC = [EPSC_on_BC_AAC epsc_on_BC_AAC];
end 

EPSC_on_BC_AAC = array2table(EPSC_on_BC_AAC);
EPSC_on_BC_AAC.Properties.VariableNames = {'BC_AAC1'...
    'BC_AAC2' 'BC_AAC3' 'BC_AAC4' 'BC_AAC5' 'BC_AAC6'...
    'BC_AAC7' 'BC_AAC8' 'BC_AAC9' 'BC_AAC10' 'BC_AAC11'...
    'BC_AAC12' 'BC_AAC13' 'BC_AAC14'};

EPSC_on_BC_AAC_mean = table2array(EPSC_on_BC_AAC(1,:));
EPSC_on_BC_AAC_std = table2array(EPSC_on_BC_AAC(2,:));
x = linspace(0,14,length(EPSC_on_BC_AAC_mean));
figure
scatter(x,EPSC_on_BC_AAC_mean,'black','filled');
xlabel('Each Point Includes BC, AAC','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak EPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,EPSC_on_BC_AAC_mean,EPSC_on_BC_AAC_std,'b','LineStyle','none')
title('Mean Peak EPSC onto BC, BiC','FontSize',15,'FontWeight','bold')

%% EPSC currents together onto BC, BiC, AAC - graph and table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EPSC_on_BC_BiC_AAC = [];
epsc_on_BC_BiC_AAC = [];
for i = 1:1:14 % number of PYR cells
    pks_epsc_all = peaks_PYR_on_BC_BiC_AAC(:,i);
    pks_epsc_all(pks_epsc_all == 0) = [];
    epsc_all_mean = mean(pks_epsc_all);
    epsc_all_std = std(pks_epsc_all);
    epsc_on_BC_BiC_AAC = [epsc_all_mean;ipsc_all_std];
    EPSC_on_BC_BiC_AAC = [EPSC_on_BC_BiC_AAC epsc_on_BC_BiC_AAC];
end 

EPSC_on_BC_BiC_AAC = array2table(EPSC_on_BC_BiC_AAC);
EPSC_on_BC_BiC_AAC.Properties.VariableNames = {'BC_BiC_AAC1'...
    'BC_BiC_AAC2' 'BC_BiC_AAC3' 'BC_BiC_AAC4' 'BC_BiC_AAC5' 'BC_BiC_AAC6'...
    'BC_BiC_AAC7' 'BC_BiC_AAC8' 'BC_BiC_AAC9' 'BC_BiC_AAC10' 'BC_BiC_AAC11'...
    'BC_BiC_AAC12' 'BC_BiC_AAC13' 'BC_BiC_AAC14'};

EPSC_on_BC_BiC_AAC_mean = table2array(EPSC_on_BC_BiC_AAC(1,:));
EPSC_on_BC_BiC_AAC_std = table2array(EPSC_on_BC_BiC_AAC(2,:));
x = linspace(0,14,length(EPSC_on_BC_AAC_mean));
figure
scatter(x,EPSC_on_BC_BiC_AAC_mean,'black','filled');
xlabel('Each Point Includes BC, BiC, AAC','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak EPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,EPSC_on_BC_BiC_AAC_mean,EPSC_on_BC_BiC_AAC_std,'b','LineStyle','none')
title('Mean Peak EPSC onto BC, BiC, AAC','FontSize',15,'FontWeight','bold')

%% E/I Ratios on BC, BiC Cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IPSC_BC_BiC_AAC_on_BC_BiC_AAC = table2array(IPSC_all_on_BC_BiC_AAC);
IPSC_BC_AAC_on_BC_AAC = table2array(IPSC_all_on_BC_AAC);

EPSC_on_BC_AAC = table2array(EPSC_on_BC_AAC);
EPSC_on_BC_BiC_AAC = table2array(EPSC_on_BC_BiC_AAC);

%IPSC_BC_on_BC = table2array(IPSC_BC_on_BC);
%IPSC_all_on_BC= table2array(IPSC_all_on_BC);  % all refers to BC and BiC together
%EPSC_BC = table2array(EPSC_BC);
%%
Ratios_BC_AAC = [];
Ratios_BC_BiC_AAC = [];
E_I_BC_AAC = abs(EPSC_on_BC_AAC(1,:)./IPSC_BC_AAC_on_BC_AAC(1,:))';
E_I_BC_BiC_AAC = abs(EPSC_on_BC_BiC_AAC(1,:)./IPSC_BC_BiC_AAC_on_BC_BiC_AAC(1,:))';

cells = 1:14;
Ratios_BC_AAC = [Ratios_BC_AAC cells' E_I_BC_AAC];
Ratios_BC_AAC = array2table(Ratios_BC_AAC);
Ratios_BC_AAC.Properties.VariableNames = {'cell_number' 'Ratio_BC_AAC'};

Ratios_BC_BiC_AAC = [Ratios_BC_BiC_AAC cells' E_I_BC_BiC_AAC];
Ratios_BC_BiC_AAC = array2table(Ratios_BC_BiC_AAC);
Ratios_BC_BiC_AAC.Properties.VariableNames = {'cell_number' 'Ratio_BC_BiC_AAC'};

%%
uitable('Data',Ratios_BC_AAC{:,:},'ColumnName',Ratios_BC_AAC.Properties.VariableNames,...
    'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

%%
uitable('Data',Ratios_BC_BiC_AAC{:,:},'ColumnName',Ratios_BC_BiC_AAC.Properties.VariableNames,...
    'Units', 'Normalized', 'Position',[0, 0, 1, 1]);


%% E/I Ratios on PYR Cells

IPSC_BiC_on_BC = table2array(IPSC_BiC_on_BC);
IPSC_BC_on_BC = table2array(IPSC_BC_on_BC);
IPSC_all_on_BC= table2array(IPSC_all_on_BC);  % all refers to BC and BiC together
EPSC_BC = table2array(EPSC_BC);

Ratios_BC = [];
E_I_BC = abs(EPSC_BC(1,:)./IPSC_BC_on_BC(1,:))';
E_I_BiC = abs(EPSC_BC(1,:)./IPSC_BiC_on_BC(1,:))';
E_I_all = abs(EPSC_BC(1,:)./IPSC_all_on_BC(1,:))'; 
 
bc = 1:14;
Ratios_BC = [Ratios_BC bc' E_I_BC E_I_BiC E_I_all];
Ratios_BC = array2table(Ratios_BC);
 
Ratios_BC.Properties.VariableNames = {'BC_no' 'Ratio_BC_BC'...
    'Ratio_BiC_BC' 'Ratio_All_BC'};

%% Voltage
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
    'pvbasket_336536_1000';...
    'pvbasket_336674_1000';...
    'pvbasket_337088_1000';...
    'pvbasket_337364_1000';...
    'pvbasket_337640_1000';...
    'pvbasket_338192_1000'...
    },{...
    'othertrace_332810_soma.dat';...
    'othertrace_333500_soma.dat';...
    'othertrace_333776_soma.dat';...
    'othertrace_334466_soma.dat';...
    'othertrace_335018_soma.dat';...
    'othertrace_335432_soma.dat';...
    'othertrace_335846_soma.dat';...
    'othertrace_336260_soma.dat';...
    'othertrace_336536_soma.dat';...
    'othertrace_336674_soma.dat';...
    'othertrace_337088_soma.dat';...
    'othertrace_337364_soma.dat';...
    'othertrace_337640_soma.dat';...
    'othertrace_338192_soma.dat'...
    });

%%
allvol = [];
for m = 1:1:14
    temp_data = readtable(g{m},'Delimiter','\t');
    temp_data = table2array(temp_data);
    allvol = [allvol temp_data];
end

voldata = mat2cell(allvol, 40000, ...
    [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10]);
%%
time = voldata{1}(:,1); % same for all cells
vol = [];
for k =1:14
    tempvol = voldata{k}(:,2);
    vol = [vol tempvol];
    subplot(7,2,k)
    plot(time,tempvol,'b')
    hold on
    title (['BC' num2str(k)]) % name each graph with the corresponding pyr #
end 
suptitle('Voltage Recordings from BC')


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
