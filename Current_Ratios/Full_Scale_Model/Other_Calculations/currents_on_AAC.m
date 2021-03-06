%% Calculate Excitatory/Inhibitory Ratios onto AACs
%  Melisa Gumus
%  May 2018 

%% Load Data From Netclamp Results
clear all
close all
clc
 
f = fullfile('/Users','macklabadmin','Documents', ...
    'other','AAC',{...
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

%%
clear all
close all
clc
 
f = fullfile('/Users','macklabadmin','Documents','other','aac',{...
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

%% Write Data on Matrix
alldata = [];
for m = 1:1:15
    temp_data = readtable(f{m},'Delimiter','\t');
    temp_data = table2array(temp_data);
    alldata = [alldata temp_data];
end

data = mat2cell(alldata, 40000, ...
    [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]);

%% Creates a big table consists of inputs from BiC, PYR, and BC... onto AAC

M = [];
current_BC = [];
current_PYR = [];
current_BiC = [];
current_cck = [];
current_ivy = [];
current_ngf = [];
current_olm = [];
current_sca = [];
figure1 = figure;
figure2 = figure;
figure3 = figure;
for m = 1:15  % number of cells 
    for k = 2:12  % number of input 
        if k == 3
            temp_current_BiC = data{m}(:,k);
            figure(figure1);
            t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
            t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
            t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
            t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
            t.String = ['Peak Detection on IPSCs from BiCs onto AACs'];
            subplot(5,3,m);
            [pks, locs] = findpeaks(data{m}(:,k),'MinPeakDistance',3000); % peak detection
            findpeaks(data{m}(:,k),'MinPeakDistance',3000);
            hold on; 
            title (['AAC Number #' num2str(m)])
            xlabel('Time (1/40 ms)')
            ylabel('IPSC')
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
            figure(figure2);
            t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
            t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
            t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
            t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
            t.String = ['Peak Detection on EPSCs onto AACs'];
            subplot(5,3,m);
            [pks, locs] = findpeaks(-data{m}(:,k),'MinPeakDistance',3000); % peak detection
            findpeaks(-data{m}(:,k),'MinPeakDistance',3000);
            hold on; 
            title (['AAC Number #' num2str(m)])
            xlabel('Time (1/40 ms)')
            ylabel('EPSC')
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
            figure(figure3);
            t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
            t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
            t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
            t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
            t.String = ['Peak Detection on IPSCs from BCs onto AACs'];
            subplot(5,3,m);
            [pks, locs] = findpeaks(data{m}(:,k),'MinPeakDistance',3000); % peak detection
            findpeaks(data{m}(:,k),'MinPeakDistance',3000);
            hold on;
            title (['AAC Number #' num2str(m)])
            xlabel('Time (1/40 ms)')
            ylabel('IPSC')
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

%% Graph of EPSCs from PYR onto AAC
figure 
for i = 1:1:15
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Distribution of EPSCs on AACs'];
    temp = allcells{i}(:,2);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(5,3,i)
    p = histfit(temp);
    set(p(1),'facecolor',[0.1 0.6 0.6]);
    set(p(2),'color','k')
    hold on
    title (['AAC Number #' num2str(i)])
    xlabel('EPSC')
    ylabel('Number of EPSCs')
end 
suptitle('Distribution of EPSCs on AAC')

%% Graph of IPSCs from BiC onto AAC
figure 
for i = 1:1:15
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Distribution of IPSCs from BiCs on AACs'];
    temp = allcells{i}(:,1);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(5,3,i)
    p = histfit(temp);
    set(p(1),'facecolor',[0.1 0.6 0.6])
    set(p(2),'color','k')
    hold on 
    title (['AAC Number #' num2str(i)])
    xlabel('IPSCs from AAC')
    ylabel('Number of IPSCs')
end 
suptitle('Distribution of IPSCs from BiC onto AAC')

%% Graph of IPSCs from BC onto AAC
figure 
for i = 1:1:15
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Distribution of IPSCs from BCs on AACs'];
    temp = allcells{i}(:,3);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(5,3,i)
    p = histfit(temp);
    set(p(1),'facecolor',[0.1 0.6 0.6])
    set(p(2),'color','k')
    hold on 
    title (['AAC Number #' num2str(i)])
    xlabel('IPSCs from BC')
    ylabel('Number of IPSCs')
end 
suptitle('Distribution of IPSCs from BC onto AAC')

%% Find Mean EPSCs and SD onto AAC
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

EPSC_table.Properties.VariableNames = {'AAC_Number', 'Mean_Peak', 'Standard_Deviation'};

subplot(3,1,1)
EPSC_mean = EPSC(1,:);
EPSC_std = EPSC(2,:);
x = linspace(0,15,length(EPSC_mean));
scatter(x,EPSC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dEPSC_mean=0.1;
text(x+dx, EPSC_mean+dEPSC_mean, c);
xlabel('Individual AACs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak EPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,EPSC_mean,EPSC_std,'b','LineStyle','none')
title('Mean Peak EPSC onto AACs','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of EPSC on AAC
fig = uitable('Data',EPSC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'AAC Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% IPSCs only from BiC onto AAC
IPSC_BiC = [];
ipsc_BiC = [];
for i = 1:1:15 
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

IPSC_BiC_table.Properties.VariableNames = {'AAC_Number', 'Mean_Peak', 'Standard_Deviation'};

subplot(3,1,2)
IPSC_BiC_mean = IPSC_BiC(1,:);
IPSC_BiC_std = IPSC_BiC(2,:);
x = linspace(0,15,length(IPSC_BiC_mean));
scatter(x,IPSC_BiC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_BiC_mean=0.1;
text(x+dx, IPSC_BiC_mean+dIPSC_BiC_mean, c);
xlabel('Individual AACs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC from BiC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BiC_mean,IPSC_BiC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BiCs onto AACs','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSC from BiC onto AAC
fig = uitable('Data',IPSC_BiC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'AAC Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% IPSCs only from BC onto AAC
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

IPSC_BC_table.Properties.VariableNames = {'BC_Number', 'Mean_Peak', 'Standard_Deviation'};

subplot(3,1,3)
IPSC_BC_mean = IPSC_BC(1,:);
IPSC_BC_std = IPSC_BC(2,:);
x = linspace(0,15,length(IPSC_BC_mean));
scatter(x,IPSC_BC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_BC_mean=0.1;
text(x+dx, IPSC_BC_mean+dIPSC_BC_mean, c);
xlabel('Individual AACs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC from BC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BC_mean,IPSC_BC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BCs onto AACs','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSC from BC onto AACs
fig = uitable('Data',IPSC_BC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'AAC Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% IPSCs Only from BC and BiC onto AAC gathered
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

% Find the peaks of the summed ipsc currents - BC and BiC
peaks_all_PV = [];
f1 = figure;
for k = 1:1:15
    figure(f1);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Peak Detection on IPSCs from BCs, BiCs on AACs'];
    subplot(5,3,k);
    [pks, locs] = findpeaks(all_ipsc(:,k),'MinPeakDistance',3000); % peak detection
    findpeaks(all_ipsc(:,k),'MinPeakDistance',3000);
    hold on; 
    title (['AAC Number #' num2str(k)])
    xlabel('Time (1/40 ms)')
    ylabel('IPSC')
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
f2 = figure;
for k = 1:1:15
    figure(f2);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Peak Detection on IPSCs from All Inhibitory Cells onto AACs'];
    subplot(5,3,k);
    [pks, locs] = findpeaks(all_ipsc_together(:,k),'MinPeakDistance',4000); % peak detection
    findpeaks(all_ipsc_together(:,k),'MinPeakDistance',4000);
    hold on; 
    title (['AAC Number #' num2str(k)])
    xlabel('Time (1/40 ms)')
    ylabel('IPSC')
    temp_cur_together = all_ipsc_together(:,k);
    allrows = (1:40000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur_together(element,:) = 0;
    end
    peaks_all_together = temp_cur_together;
    peaks_all_PV_together = [peaks_all_PV_together peaks_all_together];
end

%% IPSCs from BC and BiC onto AAC - graph and table
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

IPSC_all_table = IPSC_all;
num = (1:15)';
IPSC_all_table = array2table(IPSC_all_table');
IPSC_all_table.num = num;
IPSC_all_table = [IPSC_all_table(:,end) IPSC_all_table(:,1) IPSC_all_table(:,2)];

IPSC_all_table.Properties.VariableNames = {'AAC_Number', 'Mean_Peak', 'Standard_Deviation'};

IPSC_all_mean = IPSC_all(1,:);
IPSC_all_std = IPSC_all(2,:);
x = linspace(0,15,length(IPSC_all_mean));
figure
scatter(x,IPSC_all_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_all_mean=0.1;
text(x+dx, IPSC_all_mean+dIPSC_all_mean, c);
xlabel('Individual AACs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_all_mean,IPSC_all_std,'b','LineStyle','none')
title('Mean Peak IPSC from BCs and BiCs onto AACs','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSCs Only from BiC and BC onto AACs
fig = uitable('Data',IPSC_all_table{:,:},...
    'RowName',[],...
    'ColumnName',{'AAC Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% All IPSCs from All Inhibitory Neurons onto AAC - graph and table

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

IPSC_all_together_table.Properties.VariableNames = {'AAC_Number', 'Mean_Peak', 'Standard_Deviation'};

 
IPSC_all_together_mean = IPSC_all_together(1,:);
IPSC_all_std_together = IPSC_all_together(2,:);
x = linspace(0,14,length(IPSC_all_together_mean));
figure
scatter(x,IPSC_all_together_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_all_together_mean=0.1;
text(x+dx, IPSC_all_together_mean+dIPSC_all_together_mean, c);
xlabel('Individual AACs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSCs','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_all_together_mean,IPSC_all_std_together,'b','LineStyle','none')
title('Mean Peak IPSCs from All Inhibitory Cells onto AACs','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSCs from All Inhibitory Neurons onto AAC
fig = uitable('Data',IPSC_all_together_table{:,:},...
    'RowName',[],...
    'ColumnName',{'AAC Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% Excitatory/Inhibitory Ratios on PYR Cells
Ratios_AAC = [];
E_I_BC = abs(EPSC(1,:)./IPSC_BC(1,:))';
E_I_BiC = abs(EPSC(1,:)./IPSC_BiC(1,:))';
E_I_all = abs(EPSC(1,:)./IPSC_all(1,:))'; 
E_I_all_together = abs(EPSC(1,:)./IPSC_all_together(1,:))';

%% E/I Ratio - Table 
aac = 1:15;
Ratios_AAC = [Ratios_AAC aac' E_I_BC E_I_BiC E_I_all E_I_all_together];
Ratios_AAC = array2table(Ratios_AAC);
 
Ratios_AAC.Properties.VariableNames = {'AAC_no' 'Ratio_BC_on_AAC'...
    'Ratio_BiC_on_AAC' 'Ratio_BC_BiC_on_AAC' 'All_ipsc_onto_AAC'};

%% Display the E/I table as a figure

uitable('Data',Ratios_AAC{:,:},...
    'RowName', [],...
    'ColumnName',{'AAC Number',...
    'BC to AAC',...
    'BiC to AAC',...
    'BC, BiC to AAC',...
    'All Inhibitory Neurons to AAC'},...
    'Units', 'Normalized',...
    'Position',[0, 0, 1, 1]);

%% Voltage

g = fullfile('/home','melisagumus','Documents', ...
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
    'mytrace_0_soma.dat';...
    'mytrace_36_soma.dat';...
    'mytrace_180_soma.dat';...
    'mytrace_288_soma.dat';...
    'mytrace_360_soma.dat';...
    'mytrace_468_soma.dat';...
    'mytrace_576_soma.dat';...
    'mytrace_720_soma.dat';...
    'mytrace_828_soma.dat';...
    'mytrace_900_soma.dat';...
    'mytrace_1008_soma.dat';...
    'mytrace_1152_soma.dat';...
    'mytrace_1224_soma.dat';...
    'mytrace_1332_soma.dat';...
    'mytrace_1404_soma.dat'...
    });

%%

allvol= [];
for m = 1:15
    temp_vol = readtable(g{m},'Delimiter','\t');
    temp_vol = table2array(temp_vol);
    allvol = [allvol temp_vol(:,2)];
end

vol = mat2cell(allvol, 40000,...
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]);

%% Graph of Voltage of AAC
figure 
for i = 1:1:15
    temp = vol{i};
    %temp(temp == 0) = [];   %get rid of zeros
    subplot(5,3,i)
    plot(temp);
    hold on
    title (['AAC Number #' num2str(i)]) 
    xlabel('Time')
    ylabel('Voltage')
end 

%%







