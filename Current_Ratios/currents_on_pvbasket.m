% Load Data 
clear all
clc 

%datafile = '/home/melisagumus/Documents/MATLAB/CA1_SimTracker/pvbasket/pyr_36884_1000/mytrace_36884_syns.dat';

%ls('/home/melisagumus/Documents/MATLAB/CA1_SimTracker/pvbasket/pvbasket*1000');

f = fullfile('/home','melisagumus','Documents', ...
    'MATLAB','CA1_SimTracker','pvbasket',{...
    'pvbasket_332810_1000';...
    'pvbasket_333500_1000';...
    'pvbasket_333776_1000';...
    'pvbasket_334466_1000';...
    'pvbasket_335018_1000';...
    'pvbasket_335432_1000';...
    'pvbasket_335846_1000';...
    'pvbasket_336260_1000'...
    },{...
    'mytrace_332810_syns.dat';...
    'mytrace_333500_syns.dat';...
    'mytrace_333776_syns.dat';...
    'mytrace_334466_syns.dat';...
    'mytrace_335018_syns.dat';...
    'mytrace_335432_syns.dat';...
    'mytrace_335846_syns.dat';...
    'mytrace_336260_syns.dat'...
    });


%% 
%cd C:\Users\Melisa\Documents\MATLAB\CA1_SimTracker\pvbasket\
%files = dir('pyr*1000\mytrace*syns.dat');
%pyr_names = f.names;

alldata = [];
for m = 1:1:8   % 8 pvbasket in total
    temp_data = readtable(f{m},'Delimiter','\t');
    temp_data = table2array(temp_data);
    alldata = [alldata temp_data];
end

data = mat2cell(alldata, 40000, ...      % 4000 rows for 1000 ms 
    [12, 12, 12, 12, 12, 12, 12, 12]);   % 12 column for each pvbasket

%% Creates a big table consists of inputs from BiC PYR and BC in order

M = [];
test = [];
for m = 1:8  % number of pvbasket cells 
    for k = 2:12  % number of input 
        [pks, locs] = findpeaks(abs(data{m}(:,k))); % peak detection 
        temp = data{m}(:,k);
        allrows = (1:40000)';
        notpeak = setdiff(allrows,locs);
        for t = 1:1:numel(notpeak)
            element = notpeak(t,:);
            temp(element,:) = 0;  
        end
        if k == 3
            BiC = abs(temp);
        elseif k == 8
            PYR = abs(temp);
        elseif k == 9 
            BC = abs(temp);       
        end         
    end 
    M = [M BiC PYR BC];
end 

allcells = mat2cell(M, 40000, ...  % 4000 rows for each pvbasket
    [3, 3, 3, 3, 3, 3, 3, 3]);     % 3 columns for each cell 

%% Graph EPSCs on BC

figure 
for i = 1:1:8
    temp = allcells{i}(:,2);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(5,2,i)
    histfit(temp)
    hold on
    title (['Pyr' num2str(i)]) % name each graph with the corresponding pyr #
end 
suptitle('EPSC Distribution on BC')

%% Graph IPSCs from AAC onto BC
% 
% figure 
% for i = 1:1:8
%     temp = allcells{i}(:,1);
%     temp(temp == 0) = [];   %get rid of zeros
%     subplot(5,2,i)  
%     histfit(temp,100)
%     xlim([-2 2])
%     hold on 
%     title (['Pyr' num2str(i)]) % name each graph with the corresponding pyr #
% end 
% suptitle('IPSC Distribution from AAC onto PV Basket')

%% Graph IPSCs from BiC onto BC
figure 
for i = 1:1:8
    temp = allcells{i}(:,1);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(5,2,i)  
    histfit(temp,50)
    %xlim([-5 5])
    hold on 
    title (['Pyr' num2str(i)]) % name each graph with the corresponding pyr #
end 
suptitle('IPSC Distribution from BiC onto BC')

%% Graph IPSCs from BC onto BC

figure 
for i = 1:1:8
    temp = allcells{i}(:,3);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(5,2,i)  
    histfit(temp,50)
    xlim([-50 50])
    hold on 
    title (['Pyr' num2str(i)]) % name each graph with the corresponding pyr #
end 
suptitle('IPSC Distribution from BC onto BC')


%% EPSCs from PYR
EPSC = [];
for i = 1:1:8 % number of PYR cells
    epsc = mean(allcells{i}(:,2));
    EPSC = [EPSC epsc];
end 

EPSC = array2table(EPSC);
EPSC.Properties.VariableNames = {'pyr1'...
    'pyr2' 'pyr3' 'pyr4' 'pyr5' 'pyr6' 'pyr7' 'pyr8'};

%% IPSCs only from AAC
% IPSC_AAC = [];
% for i = 1:1:10 % number of PYR cells
%     ipsc_aac = mean(allcells{i}(:,1));
%     IPSC_AAC = [IPSC_AAC ipsc_aac];
% end 
% 
% IPSC_AAC = array2table(IPSC_AAC);
% IPSC_AAC.Properties.VariableNames = {'pyr1'...
%     'pyr2' 'pyr3' 'pyr4' 'pyr5' 'pyr6' 'pyr7' 'pyr8' 'pyr9' 'pyr10'};

%% IPSCs only from BiC
IPSC_BiC = [];
for i = 1:1:8 % number of PYR cells
    ipsc_bic = mean(allcells{i}(:,1));
    IPSC_BiC = [IPSC_BiC ipsc_bic];
end 

IPSC_BiC = array2table(IPSC_BiC);
IPSC_BiC.Properties.VariableNames = {'pyr1'...
    'pyr2' 'pyr3' 'pyr4' 'pyr5' 'pyr6' 'pyr7' 'pyr8'};

%% IPSCs only from BC
IPSC_BC = [];
for i = 1:1:8 % number of PYR cells
    ipsc_bc = mean(allcells{i}(:,3));
    IPSC_BC = [IPSC_BC ipsc_bc];
end 

IPSC_BC = array2table(IPSC_BC);
IPSC_BC.Properties.VariableNames = {'pyr1'...
    'pyr2' 'pyr3' 'pyr4' 'pyr5' 'pyr6' 'pyr7' 'pyr8'};

%% all current onto PYR gathered



%% Raster Plot for one neuron
rasterfile = '/home/melisagumus/Documents/MATLAB/CA1_SimTracker/pvbasket/pvbasket_332810_1000/spikeraster.dat';

raster = readtable(rasterfile);
num_spikes = length(raster.VarName2);
var2 = (raster.Var2)';

plot([var2; var2], [ones(1,num_spikes); ones(1, num_spikes) + 1], 'b');
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
