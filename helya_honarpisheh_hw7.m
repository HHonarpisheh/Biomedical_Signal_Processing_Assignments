close all;
clear all;

% Load the spike train data from 'spike.mat' (assuming it contains two spike trains)
load('spike.mat');

% Assuming spike1 and spike2 are your spike trains (each of size 40000x1)
spike1 = x(:, 1);
spike2 = x(:, 2);

% Define the number of repetitions (intervals) and the number of trials in each repetition
numRepetitions = 10;
numTrialsPerRepetition = 4000;

% Initialize arrays to store spike counts
spikeCounts1 = zeros(numRepetitions, 1);
spikeCounts2 = zeros(numRepetitions, 1);

% Divide the spike data into repetitions and count the spikes in each repetition
for rep = 1:numRepetitions
    startIndex = (rep - 1) * numTrialsPerRepetition + 1;
    endIndex = rep * numTrialsPerRepetition;
    spikeCounts1(rep) = sum(spike1(startIndex:endIndex));
    spikeCounts2(rep) = sum(spike2(startIndex:endIndex));
end

% Calculate the mean and variance of spike counts for each spike train
meanSpikeCount1 = mean(spikeCounts1);
varSpikeCount1 = var(spikeCounts1);

meanSpikeCount2 = mean(spikeCounts2);
varSpikeCount2 = var(spikeCounts2);

% Calculate the Fano factors for spike1 and spike2
FanoFactor1 = varSpikeCount1 / meanSpikeCount1
FanoFactor2 = varSpikeCount2 / meanSpikeCount2

% Display the Fano factors and assess if the data appears Poisson-distributed
if abs(FanoFactor1 - 1) < 0.02
    disp('Spike train 1 appears to be Poisson-distributed.');
else
    disp('Spike train 1 does not appear to be Poisson-distributed.');
end

if abs(FanoFactor2 - 1) < 0.02
    disp('Spike train 2 appears to be Poisson-distributed.');
else
    disp('Spike train 2 does not appear to be Poisson-distributed.');
end


lambda = 5; 
numSamples = 1000; % Number of samples
% Generate Poisson-distributed samples
poisson_samples = zeros(numSamples, 1);

for i = 1:numSamples
    % Generate a Poisson sample
    k = poisson_sample(lambda);
    poisson_samples(i) = k;
end

% Calculate the Fano factor
F = var(poisson_samples) / mean(poisson_samples);

% Display Fano factor
disp(['Fano factor: ' num2str(F)]);

% Custom Poisson sample generator
function sample = poisson_sample(lambda)
    sample = 0;
    p = exp(-lambda);
    F = p;
    u = rand();
    
    while u > F
        sample = sample + 1;
        p = p * lambda / sample;
        F = F + p;
    end
end
