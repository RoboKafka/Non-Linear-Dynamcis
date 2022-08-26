%
% Author: Cristiano Martinelli
%
% Code: single frequency sine generator for NI9263 

%% Create Data Acquisition Session
% Create a session for the specified vendor.

s = daq.createSession('ni');
%% Set Session Properties
% Set properties that are not using default values.

s.Rate = 10000;
%% Add Channels to Session
% Add channels and set channel properties.

addAnalogOutputChannel(s,'cDAQ2Mod4','ao0','Voltage');
%% Define Test Signal
% Create a test sine wave signal of specified peak-to-peak amplitude for each 
% channel.

amplitudePeakToPeak_ch1 = 1;

sineFrequency = 13.65; % Hz
totalDuration = 300; % seconds

outputSignal(:,1) = createSine(amplitudePeakToPeak_ch1/2, ...
                        sineFrequency, s.Rate, 'bipolar', totalDuration);
outputSignal(end+1,:) = 0;
%% Queue Signal Data
% Make signal data available to session for generation.

queueOutputData(s,outputSignal);
%% Generate Signal
% Start foreground generation

startForeground(s);
%% Clean Up
% Clear the session and channels.

clear s outputSignal
%% Create Test Signal
% Helper function for creating test sine wave signal.

function sine = createSine(amplitude, frequency, sampleRate, type, duration)

sampleRatePerCycle = floor(sampleRate/frequency);
period = 1/frequency;
s = period/sampleRatePerCycle;
t = (0 : s : period-s)';

if strcmpi(type, 'bipolar')
    y = amplitude*sin(2*pi*frequency*t);
elseif strcmpi(type, 'unipolar')
    y = amplitude*sin(2*pi*frequency*t) + amplitude;
end

numCycles = round(frequency*duration);
sine = repmat(y, numCycles, 1);
end