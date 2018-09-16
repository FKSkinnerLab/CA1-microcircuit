%% Calculate Excitatory/Inhibitory Ratios onto BiCs
%  Melisa Gumus
%  May 2018 

%% Load Data From Netclamp Results
clear all
close all
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

%%
clear all
close all
clc 
 
f = fullfile('/Users','macklabadmin','Documents', ...
    'other','BiC',{...
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
alldata = [];
for m = 1:1:15
    temp_data = readtable(f{m},'Delimiter','\t');
    temp_data = table2array(temp_data);
    alldata = [alldata temp_data];
end

data = mat2cell(alldata, 40000, ...
    [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]);

%% Creates a big table consists of inputs from AAC, BiC, PYR, and BC... onto BiC 

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

%% Graph of EPSCs on BiC

figure 
for i = 1:1:15
    temp = allcells{i}(:,2);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(5,3,i)
    p = histfit(temp);
    set(p(1),'facecolor',[0.1 0.6 0.6])
    set(p(2),'color','k')
    hold on
    title (['BiC Number #' num2str(i)]) 
    xlabel('EPSC')
    ylabel('Number of EPSCs')
end 
suptitle('Distribution of EPSCs on BiCs')

%% Graph IPSCs from BiC onto BiC
figure 
for i = 1:1:15
    temp = allcells{i}(:,1);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(5,3,i)
    p = histfit(temp);
    set(p(1),'facecolor',[0.1 0.6 0.6])
    set(p(2),'color','k')
    hold on 
    title (['BiC Number #' num2str(i)]) 
    xlabel('IPSCs from BiC')
    ylabel('Number of IPSCs')
end 
suptitle('Distribution of IPSCs from BiCs onto BiCs')

%% Graph of IPSCs from BC onto BiC
figure 
for i = 1:1:15
    temp = allcells{i}(:,3);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(5,3,i)
    p = histfit(temp);
    set(p(1),'facecolor',[0.1 0.6 0.6])
    set(p(2),'color','k')
    hold on 
    title (['BiC Number #' num2str(i)]) 
    xlabel('IPSCs from BCs')
    ylabel('Number of IPSCs')
end 
suptitle('Distribution of IPSCs from BCs onto BiCs')

%% EPSCs from PYR onto BiC
EPSC = [];
epsc = [];
for i = 1:1:15 
    pks_epsc = allcells{i}(:,2);
    pks_epsc(pks_epsc==0)=[];
    epsc_mean = mean(pks_epsc);
    epsc_std = std(pks_epsc);
    epsc = [epsc_mean;epsc_std];
    EPSC = [EPSC epsc];
end 

EPSC_table = EPSC;
num = (1:15)';
EPSC_table = array2table(EPSC_table');
EPSC_table.num = num;
EPSC_table = [EPSC_table(:,end) EPSC_table(:,1) EPSC_table(:,2)];

EPSC_table.Properties.VariableNames = {'BiC_Number', 'Mean_Peak', 'Standard_Deviation'};


subplot(3,1,1)
EPSC_mean = EPSC(1,:);
EPSC_std = EPSC(2,:);
x = linspace(0,15,length(EPSC_mean));
scatter(x,EPSC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dEPSC_mean=0.1;
text(x+dx, EPSC_mean+dEPSC_mean, c);
xlabel('Individual BiCs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak EPSCs','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,EPSC_mean,EPSC_std,'b','LineStyle','none')
title('Mean Peak EPSCs onto BiCs','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of EPSC onto BiCs
fig = uitable('Data',EPSC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'BiC Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% IPSCs only from BiC onto BiCs
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

IPSC_BiC_table = IPSC_BiC;
num = (1:15)';
IPSC_BiC_table = array2table(IPSC_BiC_table');
IPSC_BiC_table.num = num;
IPSC_BiC_table = [IPSC_BiC_table(:,end) IPSC_BiC_table(:,1) IPSC_BiC_table(:,2)];

IPSC_BiC_table.Properties.VariableNames = {'BiC_Number', 'Mean_Peak', 'Standard_Deviation'};

subplot(3,1,2)
IPSC_BiC_mean = IPSC_BiC(1,:);
IPSC_BiC_std = IPSC_BiC(2,:);
x = linspace(0,15,length(IPSC_BiC_mean));
scatter(x,IPSC_BiC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_BiC_mean=0.1;
text(x+dx, IPSC_BiC_mean+dIPSC_BiC_mean, c);
xlabel('Individual BiCs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSCs from BiCs','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BiC_mean,IPSC_BiC_std,'b','LineStyle','none')
title('Mean Peak IPSCs from BiCs onto BiCs','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSC from BiCs onto BiCs
fig = uitable('Data',IPSC_BiC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'BiC Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% IPSCs only from BC onto BiCs
IPSC_BC = [];
ipsc_BC = [];
for i = 1:1:15 
    pks_ipsc_BC = allcells{i}(:,3);
    pks_ipsc_BC(pks_ipsc_BC==0)=[];
    ipsc_BC_mean = mean(pks_ipsc_BC);
    ipsc_BC_std = std(pks_ipsc_BC);
    ipsc_BC = [ipsc_BC_mean;ipsc_BC_std];
    IPSC_BC = [IPSC_BC ipsc_BC];
end 

IPSC_BC_table = IPSC_BC;
num = (1:15)';
IPSC_BC_table = array2table(IPSC_BC_table');
IPSC_BC_table.num = num;
IPSC_BC_table = [IPSC_BC_table(:,end) IPSC_BC_table(:,1) IPSC_BC_table(:,2)];

IPSC_BC_table.Properties.VariableNames = {'BiC_Number', 'Mean_Peak', 'Standard_Deviation'};

subplot(3,1,3)
IPSC_BC_mean = IPSC_BC(1,:);
IPSC_BC_std = IPSC_BC(2,:);
x = linspace(0,15,length(IPSC_BC_mean));
scatter(x,IPSC_BC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_BC_mean=0.1;
text(x+dx, IPSC_BC_mean+dIPSC_BC_mean, c);
xlabel('Individual BiCs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSCs from BiC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BC_mean,IPSC_BC_std,'b','LineStyle','none')
title('Mean Peak IPSCs from BCs onto BiCs','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSC from BC onto BiC
fig = uitable('Data',IPSC_BC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'BiC Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% Calculate Mean IPSCs Gathered from BiC and BC onto BiC
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

%% Mean IPSCs Gathered from BC and BiC onto BiC - Table and Graph

IPSC_all = [];
ipsc_all = [];
for i = 1:1:15 
    pks_ipsc_all = peaks_all_PV(:,i);
    pks_ipsc_all(pks_ipsc_all == 0) = [];
    ipsc_all_mean = mean(pks_ipsc_all);
    ipsc_all_std = std(pks_ipsc_all);
    ipsc_all = [ipsc_all_mean;ipsc_all_std];
    IPSC_all = [IPSC_all ipsc_all];
end 

IPSC_all_table = IPSC_all;
num = (1:15)';
IPSC_all_table = array2table(IPSC_all_table');
IPSC_all_table.num = num;
IPSC_all_table = [IPSC_all_table(:,end) IPSC_all_table(:,1) IPSC_all_table(:,2)];

IPSC_all_table.Properties.VariableNames = {'BiC_Number', 'Mean_Peak', 'Standard_Deviation'};

IPSC_all_mean = IPSC_all(1,:);
IPSC_all_std = IPSC_all(2,:);
x = linspace(0,15,length(IPSC_all_mean));
figure
scatter(x,IPSC_all_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_all_mean=0.1;
text(x+dx, IPSC_all_mean+dIPSC_all_mean, c);
xlabel('Individual BiCs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSCs','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_all_mean,IPSC_all_std,'b','LineStyle','none')
title('Mean Peak IPSCs from BCs and BiCs onto BiCs','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSCs from BiC, and BC onto PYR Cells
fig = uitable('Data',IPSC_all_table{:,:},...
    'RowName',[],...
    'ColumnName',{'BiC Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% All IPSCs from All inhibitory Neurons onto BiC - graph and table

IPSC_all_together = [];
ipsc_all_together = [];
for i = 1:1:15
    pks_ipsc_all_together = peaks_all_PV_together(:,i);
    pks_ipsc_all_together(pks_ipsc_all_together == 0) = [];
    ipsc_all_mean_together = mean(pks_ipsc_all_together);
    ipsc_all_std_together = std(pks_ipsc_all_together);
    ipsc_all_together = [ipsc_all_mean_together;ipsc_all_std_together];
    IPSC_all_together = [IPSC_all_together ipsc_all_together];
end 

IPSC_all_together_table = IPSC_all_together;
num = (1:15)';
IPSC_all_together_table = array2table(IPSC_all_together_table');
IPSC_all_together_table.num = num;
IPSC_all_together_table = [IPSC_all_together_table(:,end) IPSC_all_together_table(:,1) IPSC_all_together_table(:,2)];

IPSC_all_together_table.Properties.VariableNames = {'BiC_Number', 'Mean_Peak', 'Standard_Deviation'};

IPSC_all_together_mean = IPSC_all_together(1,:);
IPSC_all_std_together = IPSC_all_together(2,:);
x = linspace(1,15,length(IPSC_all_together_mean));
figure
scatter(x,IPSC_all_together_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_all_together_mean=0.1;
text(x+dx, IPSC_all_together_mean+dIPSC_all_together_mean, c);
xlabel('Individual BiCs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_all_together_mean,IPSC_all_std_together,'b','LineStyle','none')
title('Mean Peak IPSCs from All Inhibitory Cells onto BiCs','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSCs from All Inhibitory Neurons onto BiC
fig = uitable('Data',IPSC_all_together_table{:,:},...
    'RowName',[],...
    'ColumnName',{'BiC Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% Excitatory/Inhibitory Ratios on BiC
Ratios_BiC = [];
E_I_BC = abs(EPSC(1,:)./IPSC_BC(1,:))';
E_I_BiC = abs(EPSC(1,:)./IPSC_BiC(1,:))';
E_I_all = abs(EPSC(1,:)./IPSC_all(1,:))'; 
E_I_all_together = abs(EPSC(1,:)./IPSC_all_together(1,:))';

%% E/I Ratio - Table 
bic = 1:15;
Ratios_BiC = [Ratios_BiC bic' E_I_BC E_I_BiC E_I_all E_I_all_together];
Ratios_BiC = array2table(Ratios_BiC);
 
Ratios_BiC.Properties.VariableNames = {'BiC_no' 'Ratio_BC_on_BiC'...
    'Ratio_BiC_on_BiC' 'Ratio_BC_BiC_on_BiC' 'All_ipsc_onto_BiC'};

%% Display the table as a figure

uitable('Data',Ratios_BiC{:,:},...
    'RowName', [],...
    'ColumnName',{'BiC Number',...
    'BC to BiC',...
    'BiC to BiC',...
    'BC, BiC to BiC',...
    'All Inhibitory Neurons to BiC'},...
    'Units', 'Normalized',...
    'Position',[0, 0, 1, 1]);

%%






