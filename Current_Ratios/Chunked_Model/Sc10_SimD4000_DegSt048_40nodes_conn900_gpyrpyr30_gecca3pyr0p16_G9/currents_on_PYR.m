%% Calculate Excitatory/Inhibitory Ratios onto PYRs
%  Melisa Gumus
%  May 2018


%% Load Data From Netclamp Results
clear all
close all
clc

f = fullfile('~/','Documents',...
    'labs', 'skinnerlab','networkclamp_results','Sc10_SimD4000_DegSt048_40nodes_conn900_gpyrpyr30_gecca3pyr0p16_G9',...
    'PYR',{...
    'PYR_9911';...
    'PYR_11467';...
    'PYR_20025';...
    'PYR_30917'...
    },{...
    'mytrace_9911_syns.dat';...
    'mytrace_11467_syns.dat';...
    'mytrace_20025_syns.dat';...
    'mytrace_30917_syns.dat'...
    });

%% Write Data on Matrix
alldata = [];
for m = 1:1:4
    temp_data = readtable(f{m},'Delimiter','\t');
    temp_data = table2array(temp_data);
    alldata = [alldata temp_data];
end

data = mat2cell(alldata, 80000, ...
    [13, 13, 13, 13]); %extra one is the noise cell

%% Creates a big table consists of inputs from AAC, BiC, PYR, and BC... onto PYR

M = [];
current_AAC = [];
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
figure4 = figure;
for m = 1:4  % number of cells
    for k = 2:13  % number of input - excluding first column (time)
        if k ==2
            temp_current_AAC = data{m}(:,k);
            figure(figure1);
            t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character.
            t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
            t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
            t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
            t.String = ['Peak Detection on IPSCs from AAC onto PYR'];
            subplot(2,2,m);
            [pks, locs] = findpeaks(data{m}(:,k),'MinPeakDistance',3000); % peak detection with interval based threshold
            findpeaks(data{m}(:,k),'MinPeakDistance',3000);
            hold on;
            title (['PYR Number #' num2str(m)])
            xlabel('Time (1/40 ms)')
            ylabel('IPSCs from AACs')
            temp_AAC = data{m}(:,k);
            allrows = (1:80000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_AAC(element,:) = 0;
            end
            AAC = temp_AAC;
        elseif k == 3
            temp_current_BiC = data{m}(:,k);
            figure(figure2);
            t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character.
            t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
            t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
            t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
            t.String = ['Peak Detection on IPSCs from BiC onto PYR'];
            subplot(2,2,m);
            [pks, locs] = findpeaks(data{m}(:,k),'MinPeakDistance',3000); % peak detection
            findpeaks(data{m}(:,k),'MinPeakDistance',3000);
            hold on;
            title (['PYR Number #' num2str(m)])
            xlabel('Time (1/40 ms)')
            ylabel('IPSCs from BiCs')
            temp_BiC = data{m}(:,k);
            allrows = (1:80000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_BiC(element,:) = 0;
            end
            BiC = temp_BiC;
        elseif k == 8
            temp_current_PYR = data{m}(:,k);
            figure(figure3);
            t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character.
            t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
            t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
            t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
            t.String = ['Peak Detection on EPSCs onto PYR'];
            subplot(2,2,m);
            [pks, locs] = findpeaks(-data{m}(:,k),'MinPeakDistance',3000); % peak detection
            findpeaks(-data{m}(:,k),'MinPeakDistance',3000);
            hold on;
            title (['PYR Number #' num2str(m)])
            xlabel('Time (1/40 ms)')
            ylabel('EPSCs from PYRs')
            temp_PYR = data{m}(:,k);
            allrows = (1:80000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_PYR(element,:) = 0;
            end
            PYR = temp_PYR;
        elseif k == 9
            temp_current_BC = data{m}(:,k);
            figure(figure4);
            t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character.
            t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
            t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
            t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
            t.String = ['Peak Detection on IPSCs from BC onto PYR'];
            subplot(2,2,m);
            [pks, locs] = findpeaks(data{m}(:,k),'MinPeakDistance',3000); % peak detection
            findpeaks(data{m}(:,k),'MinPeakDistance',3000);
            hold on;
            title (['PYR Number #' num2str(m)])
            xlabel('Time (1/40 ms)')
            ylabel('IPSCs from BCs')
            temp_BC = data{m}(:,k);
            allrows = (1:80000)';
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
    current_AAC = [current_AAC temp_current_AAC];
    current_PYR = [current_PYR temp_current_PYR];
    current_BiC = [current_BiC temp_current_BiC];
    current_BC = [current_BC temp_current_BC];
    current_cck = [current_cck temp_current_cck];
    current_ivy = [current_ivy temp_current_ivy];
    current_ngf = [current_ngf temp_current_ngf];
    current_olm = [current_olm temp_current_olm];
    current_sca = [current_sca temp_current_sca];
    M = [M AAC BiC PYR BC];
end

allcells = mat2cell(M, 80000, ...
    [4, 4, 4, 4]);

%% Graph of EPSCs onto PYR
f1 = figure;
for i = 1:1:4
    figure(f1);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character.
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Distribution of EPSCs on PYRs'];
    temp = allcells{i}(:,3);
    temp(temp == 0) = [];   %get rid of zeros
    s1 = subplot(2,2,i);
    p = histfit(temp);
    set(p(1),'facecolor',[0.1 0.6 0.6]);
    set(p(2),'color','k');
    hold on;
    title (['Pyramidal Cell #' num2str(i)]); % name each graph with the corresponding pyr #
    xlabel('EPSC');
    ylabel('Number of EPSCs');
end


%% Graph of IPSCs from AAC onto PYR
f2 = figure;
for i = 1:1:4
    figure(f2);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character.
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Distribution of IPSCs from AACs onto PYRs'];
    temp = allcells{i}(:,1);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(2,2,i)
    p = histfit(temp);
    set(p(1),'facecolor',[0.1 0.6 0.6])
    set(p(2),'color','k')
    hold on
    title (['Pyramidal Cell #' num2str(i)]) % name each graph with the corresponding pyr #
    xlabel('IPSCs')
    ylabel('Number of IPSCs')
end

%suptitle('Distribution of IPSCs from AAC onto PYR Cells')

%% Graph of IPSCs from BiC onto PYR
f3 = figure;
for i = 1:1:4
    figure(f3);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character.
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Distribution of IPSCs from BiCs onto PYRs'];
    temp = allcells{i}(:,2);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(2,2,i)
    p = histfit(temp);
    set(p(1),'facecolor',[0.1 0.6 0.6])
    set(p(2),'color','k')
    hold on
    title (['Pyramidal Cell #' num2str(i)]) % name each graph with the corresponding pyr #
    xlabel('IPSCs')
    ylabel('Number of IPSCs')
end
%suptitle('Distribution of IPSCs from BiC onto PYR Cells')

%% Graph of IPSCs from BC onto PYR
f4 = figure;
for i = 1:1:4
    figure(f4);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character.
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Distribution of IPSCs from BCs onto PYRs'];
    temp = allcells{i}(:,4);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(2,2,i)
    p = histfit(temp);
    set(p(1),'facecolor',[0.1 0.6 0.6])
    set(p(2),'color','k')
    hold on
    title (['Pyramidal Cell #' num2str(i)]) % name each graph with the corresponding pyr #
    xlabel('IPSCs')
    ylabel('Number of IPSCs')
end
%suptitle('Distribution of IPSCs from BC onto PYR Cells')

%% Find Mean EPSCs and SD onto PYR
EPSC = [];
epsc = [];
for i = 1:1:4 % number of PYR cells
    pks_epsc = allcells{i}(:,3);
    pks_epsc(pks_epsc==0)=[];
    epsc_mean = mean(pks_epsc);
    epsc_std = std(pks_epsc);
    epsc = [epsc_mean;epsc_std];
    EPSC = [EPSC epsc];
end

EPSC_table = EPSC;
num = (1:4)';
EPSC_table = array2table(EPSC_table');
EPSC_table.num = num;
EPSC_table = [EPSC_table(:,end) EPSC_table(:,1) EPSC_table(:,2)];

EPSC_table.Properties.VariableNames = {'Pyramidal_Cell_Number', 'Mean_Peak', 'Standard_Deviation'};

subplot(2,2,1)
EPSC_mean = EPSC(1,:);
EPSC_std = EPSC(2,:);
x = linspace(1,4,length(EPSC_mean));
scatter(x,EPSC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:4]'; b =num2str(a); c=cellstr(b);
dx=0.1; dEPSC_mean=0.1;
text(x+dx, EPSC_mean+dEPSC_mean, c);
xlabel('Individual PYR Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak EPSCs','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,EPSC_mean,EPSC_std,'b','LineStyle','none')
title('Mean Peak EPSCs onto PYR Cells','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of EPSC on PYR Cells
fig = uitable('Data',EPSC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'Pyramidal Cell Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

% txt_title = uicontrol('Style','text',...
%     'Units','Normalized',...
%     'String','pyr',...
%     'Position',[0.1, 0.9, 0.47, 1]);


%% IPSCs only from AAC onto PYR
IPSC_AAC = [];
ipsc_AAC = [];
for i = 1:1:4 % number of PYR cells
    pks_ipsc_AAC = allcells{i}(:,1);
    pks_ipsc_AAC(pks_ipsc_AAC==0)=[];
    ipsc_AAC_mean = mean(pks_ipsc_AAC);
    ipsc_AAC_std = std(pks_ipsc_AAC);
    ipsc_AAC = [ipsc_AAC_mean;ipsc_AAC_std];
    IPSC_AAC = [IPSC_AAC ipsc_AAC];
end

IPSC_AAC_table = IPSC_AAC;
num = (1:4)';
IPSC_AAC_table = array2table(IPSC_AAC_table');
IPSC_AAC_table.num = num;
IPSC_AAC_table = [IPSC_AAC_table(:,end) IPSC_AAC_table(:,1) IPSC_AAC_table(:,2)];

IPSC_AAC_table.Properties.VariableNames = {'Pyramidal_Cell_Number', 'Mean_Peak', 'Standard_Deviation'};

subplot(2,2,2)
IPSC_AAC_mean = IPSC_AAC(1,:);
IPSC_AAC_std = IPSC_AAC(2,:);
x = linspace(0,4,length(IPSC_AAC_mean));
scatter(x,IPSC_AAC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:4]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_AAC_mean=0.1;
text(x+dx, IPSC_AAC_mean+dIPSC_AAC_mean, c);
xlabel('Individual PYR Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC from AAC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_AAC_mean,IPSC_AAC_std,'b','LineStyle','none')
title('Mean Peak IPSC from AAC onto PYR Cells','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSC from AAC onto PYR Cells
fig = uitable('Data',IPSC_AAC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'Pyramidal Cell Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% IPSCs only from BiC onto PYR
IPSC_BiC = [];
ipsc_BiC = [];
for i = 1:1:4 % number of PYR cells
    pks_ipsc_BiC = allcells{i}(:,2);
    pks_ipsc_BiC(pks_ipsc_BiC==0)=[];
    ipsc_BiC_mean = mean(pks_ipsc_BiC);
    ipsc_BiC_std = std(pks_ipsc_BiC);
    ipsc_BiC = [ipsc_BiC_mean;ipsc_BiC_std];
    IPSC_BiC = [IPSC_BiC ipsc_BiC];
end

IPSC_BiC_table = IPSC_BiC;
num = (1:4)';
IPSC_BiC_table = array2table(IPSC_BiC_table');
IPSC_BiC_table.num = num;
IPSC_BiC_table = [IPSC_BiC_table(:,end) IPSC_BiC_table(:,1) IPSC_BiC_table(:,2)];

IPSC_BiC_table.Properties.VariableNames = {'Pyramidal_Cell_Number', 'Mean_Peak', 'Standard_Deviation'};

subplot(2,2,3)
IPSC_BiC_mean = IPSC_BiC(1,:);
IPSC_BiC_std = IPSC_BiC(2,:);
x = linspace(0,4,length(IPSC_BiC_mean));
scatter(x,IPSC_BiC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:4]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_BiC_mean=0.1;
text(x+dx, IPSC_BiC_mean+dIPSC_BiC_mean, c);
xlabel('Individual PYR Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC from BiC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BiC_mean,IPSC_BiC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BiC onto PYR Cells','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSC from BiC onto PYR Cells
fig = uitable('Data',IPSC_BiC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'Pyramidal Cell Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% IPSCs only from BC onto PYR
IPSC_BC = [];
ipsc_BC = [];
for i = 1:1:4 % number of PYR cells
    pks_ipsc_BC = allcells{i}(:,4);
    pks_ipsc_BC(pks_ipsc_BC==0)=[];
    ipsc_BC_mean = mean(pks_ipsc_BC);
    ipsc_BC_std = std(pks_ipsc_BC);
    ipsc_BC = [ipsc_BC_mean;ipsc_BC_std];
    IPSC_BC = [IPSC_BC ipsc_BC];
end

IPSC_BC_table = IPSC_BC;
num = (1:4)';
IPSC_BC_table = array2table(IPSC_BC_table');
IPSC_BC_table.num = num;
IPSC_BC_table = [IPSC_BC_table(:,end) IPSC_BC_table(:,1) IPSC_BC_table(:,2)];

IPSC_BC_table.Properties.VariableNames = {'Pyramidal_Cell_Number', 'Mean_Peak', 'Standard_Deviation'};

subplot(2,2,4)
IPSC_BC_mean = IPSC_BC(1,:);
IPSC_BC_std = IPSC_BC(2,:);
x = linspace(0,4,length(IPSC_BC_mean));
scatter(x,IPSC_BC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:4]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_BC_mean=0.1;
text(x+dx, IPSC_BC_mean+dIPSC_BC_mean, c);
xlabel('Individual PYR Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC from BC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BC_mean,IPSC_BC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BC onto PYR Cells','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSC from BC onto PYR Cells
fig = uitable('Data',IPSC_BC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'Pyramidal Cell Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% Calculate Mean IPSCs Gathered from AAC, BiC and BC onto PYR

% Sum all IPSC currents
all_ipsc = [];
all_ipsc_together = [];

for i = 1:1:4
    tot_cur_ipsc = current_AAC(:,i) + current_BiC(:,i) + current_BC(:,i);
    tot_ipsc_together = current_BiC(:,i)...
        +current_AAC(:,i)...
        +current_BC(:,i)...
        +current_cck(:,i)...
        +current_ivy(:,i)...
        +current_ngf(:,i)...
        +current_olm(:,i)...
        +current_sca(:,i);

    all_ipsc = [all_ipsc tot_cur_ipsc];
    all_ipsc_together = [all_ipsc_together tot_ipsc_together];
end

% Find the peaks of the summed ipsc currents - from BC, BiC, AAC
peaks_all_PV = [];
f1 = figure;
for k = 1:1:4
    figure(f1);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character.
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Peak Detection on IPSCs from BCs, BiCs, AACs onto PYRs'];
    subplot(2,2,k);
    [pks, locs] = findpeaks(all_ipsc(:,k),'MinPeakDistance',3000); % peak detection
    findpeaks(all_ipsc(:,k),'MinPeakDistance',3000);
    hold on;
    title (['PYR Number #' num2str(k)])
    xlabel('Time (1/40 ms)')
    ylabel('IPSC')
    temp_cur = all_ipsc(:,k);
    allrows = (1:80000)';
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
for k = 1:1:4
    figure(f2);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character.
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Peak Detection on IPSCs from All inhibitory Cells onto PYRs'];
    subplot(2,2,k);
    [pks, locs] = findpeaks(all_ipsc_together(:,k),'MinPeakDistance',3000); % peak detection
    findpeaks(all_ipsc_together(:,k),'MinPeakDistance',3000);
    hold on;
    title (['PYR Number #' num2str(k)])
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

%% Mean IPSCs Gathered from AAC, BiC and BC onto PYR - Table and Graph

IPSC_all = [];
ipsc_all = [];
for i = 1:1:4 % number of PYR cells
    pks_ipsc_all = peaks_all_PV(:,i);
    pks_ipsc_all(pks_ipsc_all == 0) = [];
    ipsc_all_mean = mean(pks_ipsc_all);
    ipsc_all_std = std(pks_ipsc_all);
    ipsc_all = [ipsc_all_mean;ipsc_all_std];
    IPSC_all = [IPSC_all ipsc_all];
end

IPSC_all_table = IPSC_all;
num = (1:4)';
IPSC_all_table = array2table(IPSC_all_table');
IPSC_all_table.num = num;
IPSC_all_table = [IPSC_all_table(:,end) IPSC_all_table(:,1) IPSC_all_table(:,2)];

IPSC_all_table.Properties.VariableNames = {'Pyramidal_Cell_Number', 'Mean_Peak', 'Standard_Deviation'};

IPSC_all_mean = IPSC_all(1,:);
IPSC_all_std = IPSC_all(2,:);
x = linspace(0,4,length(IPSC_all_mean));
figure
scatter(x,IPSC_all_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:4]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_all_mean=0.1;
text(x+dx, IPSC_all_mean+dIPSC_all_mean, c);
xlabel('Individual PYR Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_all_mean,IPSC_all_std,'b','LineStyle','none')
title('Mean Peak IPSC from BC, BiC and AAC onto PYR Cells','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSCs from AAC, BiC, and BC onto PYR Cells
fig = uitable('Data',IPSC_all_table{:,:},...
    'RowName',[],...
    'ColumnName',{'Pyramidal Cell Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% All IPSCs from All Inhibitory Neurons onto PYR - graph and table

IPSC_all_together = [];
ipsc_all_together = [];
for i = 1:1:4% number of PYR cells
    pks_ipsc_all_together = peaks_all_PV_together(:,i);
    pks_ipsc_all_together(pks_ipsc_all_together == 0) = [];
    ipsc_all_mean_together = mean(pks_ipsc_all_together);
    ipsc_all_std_together = std(pks_ipsc_all_together);
    ipsc_all_together = [ipsc_all_mean_together;ipsc_all_std_together];
    IPSC_all_together = [IPSC_all_together ipsc_all_together];
end

IPSC_all_together_table = IPSC_all_together;
num = (1:4)';
IPSC_all_together_table = array2table(IPSC_all_together_table');
IPSC_all_together_table.num = num;
IPSC_all_together_table = [IPSC_all_together_table(:,end) IPSC_all_together_table(:,1) IPSC_all_together_table(:,2)];

IPSC_all_together_table.Properties.VariableNames = {'Pyramidal_Cell_Number', 'Mean_Peak', 'Standard_Deviation'};


IPSC_all_together_mean = IPSC_all_together(1,:);
IPSC_all_std_together = IPSC_all_together(2,:);
x = linspace(0,4,length(IPSC_all_together_mean));
figure
scatter(x,IPSC_all_together_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:4]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_all_together_mean=0.1;
text(x+dx, IPSC_all_together_mean+dIPSC_all_together_mean, c);
xlabel('Individual PYR Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_all_together_mean,IPSC_all_std_together,'b','LineStyle','none')
title('Mean Peak IPSC from All Inhibitory Cells onto PYRs','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSCs from All Inhibitory Neurons onto PYR Cells
fig = uitable('Data',IPSC_all_together_table{:,:},...
    'RowName',[],...
    'ColumnName',{'Pyramidal Cell Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% IPSCs Only from AAC and BC onto PYR gathered
% Sum all ipsc currents
AAC_BC_ipsc = [];
for i = 1:1:4
    AAC_BC_cur_ipsc = current_AAC(:,i) + current_BC(:,i);
    AAC_BC_ipsc = [AAC_BC_ipsc AAC_BC_cur_ipsc];
end

% Find the peaks of the summed ipsc currents
peaks_AAC_BC = [];
f1 = figure;
for k = 1:1:4
    figure(f1);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character.
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Peak Detection on IPSCs from AACs, BCs onto PYRs'];
    subplot(2,2,k);
    [pks, locs] = findpeaks(AAC_BC_ipsc(:,k),'MinPeakDistance',3000); % peak detection
    findpeaks(AAC_BC_ipsc(:,k),'MinPeakDistance',3000);
    hold on;
    title (['PYR Number #' num2str(k)])
    xlabel('Time (1/40 ms)')
    ylabel('IPSC')
    temp_cur = AAC_BC_ipsc(:,k);
    allrows = (1:80000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    pk_AAC_BC = temp_cur;
    peaks_AAC_BC = [peaks_AAC_BC pk_AAC_BC];
end

%% IPSC from AAC and BC onto PYR - graph and table
IPSC_AAC_BC = [];
ipsc_AAC_BC = [];
for i = 1:1:4 % number of PYR cells
    pks_ipsc_AAC_BC = peaks_AAC_BC(:,i);
    pks_ipsc_AAC_BC(pks_ipsc_AAC_BC == 0) = [];
    ipsc_AAC_BC_mean = mean(pks_ipsc_AAC_BC);
    ipsc_AAC_BC_std = std(pks_ipsc_AAC_BC);
    ipsc_AAC_BC = [ipsc_AAC_BC_mean;ipsc_AAC_BC_std];
    IPSC_AAC_BC = [IPSC_AAC_BC ipsc_AAC_BC];
end

IPSC_AAC_BC_table = IPSC_AAC_BC;
num = (1:4)';
IPSC_AAC_BC_table = array2table(IPSC_AAC_BC_table');
IPSC_AAC_BC_table.num = num;
IPSC_AAC_BC_table = [IPSC_AAC_BC_table(:,end) IPSC_AAC_BC_table(:,1) IPSC_AAC_BC_table(:,2)];

IPSC_AAC_BC_table.Properties.VariableNames = {'Pyramidal_Cell_Number', 'Mean_Peak', 'Standard_Deviation'};

IPSC_AAC_BC_mean = IPSC_AAC_BC(1,:);
IPSC_AAC_BC_std = IPSC_AAC_BC(2,:);
x = linspace(1,4,length(IPSC_AAC_BC_mean));
figure
scatter(x(1,:),IPSC_AAC_BC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:4]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_AAC_BC_mean=0.1;
text(x+dx, IPSC_AAC_BC_mean+dIPSC_AAC_BC_mean, c);
xlabel('Individual PYR Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_all_mean,IPSC_all_std,'b','LineStyle','none')
title('Mean Peak IPSC from BC and AAC onto PYR Cells','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSCs Only from AAC and BC onto PYR Cells
fig = uitable('Data',IPSC_AAC_BC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'Pyramidal Cell Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% IPSCs Only from BiC and BC onto PYR gathered
% Sum all ipsc currents
BiC_BC_ipsc = [];
for i = 1:1:4
    BiC_BC_cur_ipsc = current_BiC(:,i) + current_BC(:,i);
    BiC_BC_ipsc = [BiC_BC_ipsc BiC_BC_cur_ipsc];
end

% Find the peaks of the summed ipsc currents
peaks_BiC_BC = [];
f1 = figure;
for k = 1:1:4
    figure(f1);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character.
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Peak Detection on IPSCs from BiCs, BCs onto PYRs'];
    subplot(2,2,k);
    [pks, locs] = findpeaks(BiC_BC_ipsc(:,k),'MinPeakDistance',3000); % peak detection
    findpeaks(BiC_BC_ipsc(:,k),'MinPeakDistance',3000);
    hold on;
    title (['PYR Number #' num2str(k)])
    xlabel('Time (1/40 ms)')
    ylabel('IPSC')
    temp_cur = BiC_BC_ipsc(:,k);
    allrows = (1:80000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    pk_BiC_BC = temp_cur;
    peaks_BiC_BC = [peaks_BiC_BC pk_BiC_BC];
end

%% IPSC from BiC and BC onto PYR - graph and table
IPSC_BiC_BC = [];
ipsc_BiC_BC = [];
for i = 1:1:4 % number of PYR cells
    pks_ipsc_BiC_BC = peaks_BiC_BC(:,i);
    pks_ipsc_BiC_BC(pks_ipsc_BiC_BC == 0) = [];
    ipsc_BiC_BC_mean = mean(pks_ipsc_BiC_BC);
    ipsc_BiC_BC_std = std(pks_ipsc_BiC_BC);
    ipsc_BiC_BC = [ipsc_BiC_BC_mean;ipsc_BiC_BC_std];
    IPSC_BiC_BC = [IPSC_BiC_BC ipsc_BiC_BC];
end

IPSC_BiC_BC_table = IPSC_BiC_BC;
num = (1:4)';
IPSC_BiC_BC_table = array2table(IPSC_BiC_BC_table');
IPSC_BiC_BC_table.num = num;
IPSC_BiC_BC_table = [IPSC_BiC_BC_table(:,end) IPSC_BiC_BC_table(:,1) IPSC_BiC_BC_table(:,2)];

IPSC_BiC_BC_table.Properties.VariableNames = {'Pyramidal_Cell_Number', 'Mean_Peak', 'Standard_Deviation'};

IPSC_BiC_BC_mean = IPSC_BiC_BC(1,:);
IPSC_BiC_BC_std = IPSC_BiC_BC(2,:);
x = linspace(1,4,length(IPSC_BiC_BC_mean));
figure
scatter(x(1,:),IPSC_BiC_BC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:4]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_BiC_BC_mean=0.1;
text(x+dx, IPSC_BiC_BC_mean+dIPSC_BiC_BC_mean, c);
xlabel('Individual PYR Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_all_mean,IPSC_all_std,'b','LineStyle','none')
title('Mean Peak IPSC from BC and BiC onto PYR Cells','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSCs Only from BiC and BC onto PYR Cells
fig = uitable('Data',IPSC_BiC_BC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'Pyramidal Cell Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% Excitatory/Inhibitory Ratios on PYR Cells
Ratios_PYR = [];
E_I_AAC = abs(EPSC(1,:)./IPSC_AAC(1,:))';
E_I_BC = abs(EPSC(1,:)./IPSC_BC(1,:))';
E_I_BiC = abs(EPSC(1,:)./IPSC_BiC(1,:))';
E_I_AAC_BC = abs(EPSC(1,:)./IPSC_AAC_BC(1,:))';
E_I_BiC_BC = abs(EPSC(1,:)./IPSC_BiC_BC(1,:))';

E_I_all = abs(EPSC(1,:)./IPSC_all(1,:))';
E_I_all_together = abs(EPSC(1,:)./IPSC_all_together(1,:))';

%% E/I Ratio - Table
pyr = 1:4;
Ratios_PYR = [pyr' E_I_AAC E_I_BC E_I_BiC E_I_AAC_BC E_I_BiC_BC E_I_all E_I_all_together];
Ratios_PYR = array2table(Ratios_PYR);

Ratios_PYR.Properties.VariableNames = {'pyr_no' 'Ratio_AAC_on_PYR'...
    'Ratio_BC_on_PYR' 'Ratio_BiC_on_PYR' ...
    'Ratio_AAC_BC_on_PYR'...
    'Ratio_BiC_BC_on_PYR'...
    'Ratio_AAC_BiC_BC_on_PYR' 'All_ipsc_onto_PYR'};

%% Display the E/I table as a figure

uitable('Data',Ratios_PYR{:,:},...
    'RowName', [],...
    'ColumnName',{'Pyramidal Cell Number',...
    'AAC to PYR',...
    'BC to PYR',...
    'BiC to PYR',...
    'AAC, BC onto PYR',...
    'BiC, BC onto PYR',...
    'AAC, BiC and BC to PYR',...
    'All Inhibitory Neurons to PYR'},...
    'Units', 'Normalized',...
    'Position',[0, 0, 1, 1]);

%%

g = fullfile('~/','Documents',...
    'labs', 'skinnerlab','networkclamp_results','Sc10_SimD4000_DegSt048_40nodes_conn900_gpyrpyr30_gecca3pyr0p16_G9',...
    'PYR',{...
    'PYR_9911';...
    'PYR_11467';...
    'PYR_20025';...
    'PYR_30917'...
    },{...
    'mytrace_9911_soma.dat';...
    'mytrace_11467_soma.dat';...
    'mytrace_20025_soma.dat';...
    'mytrace_30917_soma.dat'...
    });


%%

allvol= [];
for m = 1:4
    temp_vol = readtable(g{m},'Delimiter','\t');
    temp_vol = table2array(temp_vol);
    allvol = [allvol temp_vol(:,2)];
end

vol = mat2cell(allvol, 80000,...
    [1, 1, 1, 1]);

%% Graph of Voltage of BC
figure
for i = 1:1:4
    temp = vol{i};
    %temp(temp == 0) = [];   %get rid of zeros
    subplot(2,2,i)
    plot(temp);
    hold on
    title (['PYR Number #' num2str(i)])
    xlabel('Time')
    ylabel('Voltage')
end

%%
