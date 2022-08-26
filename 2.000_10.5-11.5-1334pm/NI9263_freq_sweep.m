% Author: Cristiano Martinelli
%
% Code: generating a frequency sweep
%
% date: 08/07/2022

% controlling NI 9263 DAQ

% NI9263
daq.getDevices()
s = daq.createSession('ni');
s.Rate = 10000;

% Frequence Sweep
frq = 9.5:0.05:12;
%frq = flip(frq);

% Position of changes
pos_chg = zeros(size(frq));
% Time for each frequency
delta_time = zeros(size(frq));
% Time in seconds for each frequency block
Tp = 10;
% amplitude of sine wave
Amp = 1; % V

for ii = 1:size(frq,2)
    sampleRatePerCycle = floor(s.Rate/frq(ii));
    % Period
    T = 1/frq(ii);
    % Number of periods
    nP = ceil(Tp/T);
    % Initial time
    dt = T/sampleRatePerCycle;
    tr = 0:dt:T*nP;
    % Delta time for the frequency block
    delta_time(ii) = tr(end)-tr(1);
    % Sinusoidal frequency
    ys = sin(frq(ii)*2*pi*tr)*Amp; %[V]
    % Concatenating the vector
    if ii == 1
        t = tr;
        y = ys;
        pos_chg(ii) = size(t,2);
    else
        t = [t,tr(2:end)+t(end)];
        y = [y,ys(2:end)];
        pos_chg(ii) = size(t,2);
    end
end

%
figure()
plot(t,y,'-b','linewidth',1)
hold on
xline(t(pos_chg),'-k')
grid on
box on
xlabel('Time [s]')
ylabel('Amplitude [V]')
title('Signal in time')

figure()
plot(diff(t))
ylim([0 0.002])
grid on
box on
xlabel('Points')
ylabel('Delta t [s]')
title('Difference in time')


figure()
plot(delta_time)
grid on
box on
xlabel('Points')
ylabel('Time lenght [s]')
title('Time lenght for each block of signal')



s.addAnalogOutputChannel('cDAQ1Mod1', 'ao0', 'voltage'); %in R2014b 'Voltage' !!!
s.queueOutputData (y');
% start the signal
s.startForeground;
% end of signal
s.removeChannel(1)

