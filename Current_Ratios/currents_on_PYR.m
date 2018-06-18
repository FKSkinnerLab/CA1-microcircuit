% Load Data 
clear all
clc 

% datafile = '/home/melisagumus/Documents/MATLAB/CA1_SimTracker/pyr/pyr_36884_1000/mytrace_36884_syns.dat';
% 
% ls('/home/melisagumus/Documents/MATLAB/CA1_SimTracker/pyr/pyr*1000');
 
f = fullfile('/home','melisagumus','Documents', ...
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
    'mytrace_36884_syns.dat';...
    'mytrace_68032_syns.dat';...
    'mytrace_83606_syns.dat';...
    'mytrace_99180_syns.dat';...
    'mytrace_106967_syns.dat';...
    'mytrace_114754_syns.dat';...
    'mytrace_200411_syns.dat';...
    'mytrace_254920_syns.dat';...
    'mytrace_286068_syns.dat';...
    'mytrace_301642_syns.dat'...
    });


%% 
%cd C:\Users\Melisa\Documents\MATLAB\CA1_SimTracker\pyr\
%files = dir('pyr*1000\mytrace*syns.dat');
%pyr_names = f.names;

alldata = [];
for m = 1:1:10
    temp_data = readtable(f{m},'Delimiter','\t');
    temp_data = table2array(temp_data);
    alldata = [alldata temp_data];
end

data = mat2cell(alldata, 40000, ...
    [12, 12, 12, 12, 12, 12, 12, 12, 12, 12]);

%% Creates a big table consists of inputs from AAC BiC PYR and BC in order

M = [];
test = [];
for m = 1:10  % number of cells 
    for k = 2:12  % number of input 
        [pks, locs] = findpeaks(abs(data{m}(:,k))); % peak detection 
        temp = data{m}(:,k);
        allrows = (1:40000)';
        notpeak = setdiff(allrows,locs);
        for t = 1:1:numel(notpeak)
            element = notpeak(t,:);
            temp(element,:) = 0;  
        end 
        if k == 2
            AAC = abs(temp);
        elseif k == 3
            BiC = abs(temp);
        elseif k == 8
            PYR_cell = abs(temp);
        elseif k == 9 
            BC = abs(temp);       
        end         
    end 
    M = [M AAC BiC PYR_cell BC];
end 

allcells = mat2cell(M, 40000, ...
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4]); 

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
for i = 1:1:10 % number of PYR cells
    epsc = mean(allcells{i}(:,3));
    EPSC = [EPSC epsc];
end 

EPSC = array2table(EPSC);
EPSC.Properties.VariableNames = {'pyr1'...
    'pyr2' 'pyr3' 'pyr4' 'pyr5' 'pyr6' 'pyr7' 'pyr8' 'pyr9' 'pyr10'};

%% IPSCs only from AAC
IPSC_AAC = [];
for i = 1:1:10 % number of PYR cells
    ipsc_aac = mean(allcells{i}(:,1));
    IPSC_AAC = [IPSC_AAC ipsc_aac];
end 

IPSC_AAC = array2table(IPSC_AAC);
IPSC_AAC.Properties.VariableNames = {'pyr1'...
    'pyr2' 'pyr3' 'pyr4' 'pyr5' 'pyr6' 'pyr7' 'pyr8' 'pyr9' 'pyr10'};

%% IPSCs only from BiC
IPSC_BiC = [];
for i = 1:1:10 % number of PYR cells
    ipsc_bic = mean(allcells{i}(:,2));
    IPSC_BiC = [IPSC_BiC ipsc_bic];
end 

IPSC_BiC = array2table(IPSC_BiC);
IPSC_BiC.Properties.VariableNames = {'pyr1'...
    'pyr2' 'pyr3' 'pyr4' 'pyr5' 'pyr6' 'pyr7' 'pyr8' 'pyr9' 'pyr10'};

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
