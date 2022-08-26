%% aquiring all data
close all
clear all
clc
% forward sweep
A = load('10p5_f.mat');
A(:,:,2) = load('10p6_f.mat');
A(:,:,3) =load('10p7_f.mat');
A(:,:,4) = load('10p8_f.mat');
A(:,:,5) = load('10p9_f.mat');
A(:,:,6) = load('11p0_f.mat');
A(:,:,7) = load('11p1_f.mat');
A(:,:,8) = load('11p2_f.mat');
A(:,:,9) = load('11p3_f.mat');
A(:,:,10) = load('11p4_f.mat');
A(:,:,11) = load('11p5_f.mat');

% Backward sweep
B = load('10p5_b.mat');
B(:,:,2) = load('10p6_b.mat');
B(:,:,3) =load('10p7_b.mat');
B(:,:,4) = load('10p8_b.mat');
B(:,:,5) = load('10p9_b.mat');
B(:,:,6) = load('11p0_b.mat');
B(:,:,7) = load('11p1_b.mat');
B(:,:,8) = load('11p2_b.mat');
B(:,:,9) = load('11p3_b.mat');
B(:,:,10) = load('11p4_b.mat');
B(:,:,11) = load('11p5_b.mat');

% initiate amplitude matrix
amplitude_f = [];
ampltiude_b = [];
freq_sweep = (10.5:0.1:11.5);

%% sensitivity and initiating arguments
g = 9.80665;     % m/s^2 (gravity acceleration)

% Accelerometer Data
sens_44 = 10.88; % mV/g  (accelerometer sensitivity 44)
v044 = 0;
sens_46 = 10.71; % mV/g  (accelerometer sensitivity 46)
v046 = 0;
for i = [1:11]
    % data extraction
    t_f = cell2mat(struct2cell(A(:,:,i)));
    t_f = t_f(:,1);% time [s]
    dt_f = t_f(2)-t_f(1); % delta t [s]
    fs_f = 1/dt_f; % Sampling Frequency [Hz]
    
    t_b = cell2mat(struct2cell(B(:,:,i)));
    t_b = t_b(:,1);% time [s]
    dt_b = t_b(2)-t_b(1); % delta t [s]
    fs_b = 1/dt_b; % Sampling Frequency [Hz]

    % acc 46
    volt_acc46_f = cell2mat(struct2cell(A(:,:,i)));% accelerometrt voltage [V]
    volt_acc46_f = volt_acc46_f(:,3);
    aA46_f = ((volt_acc46_f-v046/1000)/(sens_46/1000))*g; % acceleromenter acceleration [m/s^2]

    
    volt_acc46_b = cell2mat(struct2cell(B(:,:,i)));% accelerometrt voltage [V]
    volt_acc46_b = volt_acc46_b(:,3);
    aA46_b = ((volt_acc46_b-v046/1000)/(sens_46/1000))*g; % acceleromenter acceleration [m/s^2]

    %acc 44
    volt_acc44_f = cell2mat(struct2cell(A(:,:,i)));% accelerometrt voltage [V]
    volt_acc44_f = volt_acc44_f(:,2);
    aA44_f = ((volt_acc44_f-v044/1000)/(sens_44/1000))*g; % acceleromenter acceleration [m/s^2]

    volt_acc44_b = cell2mat(struct2cell(B(:,:,i)));% accelerometrt voltage [V]
    volt_acc44_b = volt_acc44_b(:,2);
    aA44_b = ((volt_acc44_b-v044/1000)/(sens_44/1000))*g; % acceleromenter acceleration [m/s^2]

    %% Filter
    fc_h_f = 3/(fs_f/2);  % Cut off Frequency (high)1
    fc_l_f = 100/(fs_f/2);  % Cut off Frequency (low)
    order = 6; % 6th Order Filter
    [b1_f,a1_f] = butter(order,fc_h_f,'high');
    [b2_f,a2_f] = butter(order,fc_l_f,'low');
    
    fc_h_b= 3/(fs_b/2);  % Cut off Frequency (high)1
    fc_l_b = 100/(fs_b/2);  % Cut off Frequency (low)
    order = 6; % 6th Order Filter
    [b1_b,a1_b] = butter(order,fc_h_b,'high');
    [b2_b,a2_b] = butter(order,fc_l_b,'low');


    %% Accelerometer

    % Acceleration accelerometer
    aAf44_f = filtfilt(b1_f,a1_f,aA44_f);
    aAf44_f = filtfilt(b2_f,a2_f,aAf44_f);

    aAf46_f = filtfilt(b1_f,a1_f,aA46_f);
    aAf46_f = filtfilt(b2_f,a2_f,aAf46_f);
    
    aAf44_b = filtfilt(b1_b,a1_b,aA44_b);
    aAf44_b = filtfilt(b2_b,a2_b,aAf44_b);

    aAf46_b = filtfilt(b1_b,a1_b,aA46_b);
    aAf46_b = filtfilt(b2_b,a2_b,aAf46_b);

    % % subplot(3,2,1)
    % % plot(t,aA44,t,aAf44,'linewidth',2);
    % % hold on
    %figure()
    %plot(t,aA46,t,aAf46,'linewidth',2);
    % xlim([3 4])
    % grid on
    % box on
    % xlabel('Time [s]')
    % ylabel('Acceleration [m/s^2]')
    % legend('Acc44-raw','Acc44-filtered','Acc46-raw','Acc46-filtered')
    % title('Acceleration accelerometer')
    %% Velocity from ACC
    vA46_f = cumtrapz(t_f,aAf46_f);
    vA46_f = detrend(vA46_f);
    vAf46_f = filtfilt(b1_f,a1_f,vA46_f);
    vAf46_f = filtfilt(b2_f,a2_f,vAf46_f);

    vA44_f = cumtrapz(t_f,aAf44_f);
    vA44_f = detrend(vA44_f);
    vAf44_f = filtfilt(b1_f,a1_f,vA44_f);
    vAf44_f= filtfilt(b2_f,a2_f,vAf44_f);
    % back sweep 
    vA46_b = cumtrapz(t_b,aAf46_b);
    vA46_b = detrend(vA46_b);
    vAf46_b = filtfilt(b1_b,a1_b,vA46_b);
    vAf46_b = filtfilt(b2_b,a2_b,vAf46_b);

    vA44_b = cumtrapz(t_b,aAf44_b);
    vA44_b = detrend(vA44_b);
    vAf44_b = filtfilt(b1_b,a1_b,vA44_b);
    vAf44_b= filtfilt(b2_b,a2_b,vAf44_b);

    %figure()
    %plot(t,vA46,t,vAf46,'linewidth',2);
    % xlim([3 4])
    % grid on
    % box on
    % xlabel('Time [s]')
    % ylabel('Velocity [m/s]')
    %% Displacement from velocity
    dA46_f = cumtrapz(t_f,vAf46_f);
    dA46_f = detrend(dA46_f);
    dAf46_f = filtfilt(b1_f,a1_f,dA46_f);
    dAf46_f = filtfilt(b2_f,a2_f,dAf46_f);

    dA44_f = cumtrapz(t_f,vAf44_f);
    dA44_f = detrend(dA44_f);
    dAf44_f = filtfilt(b1_f,a1_f,dA44_f);
    dAf44_f = filtfilt(b2_f,a2_f,dAf44_f);
    
    dA46_b = cumtrapz(t_b,vAf46_b);
    dA46_b = detrend(dA46_b);
    dAf46_b = filtfilt(b1_b,a1_b,dA46_b);
    dAf46_b = filtfilt(b2_b,a2_b,dAf46_b);

    dA44_b = cumtrapz(t_b,vAf44_b);
    dA44_b = detrend(dA44_b);
    dAf44_b = filtfilt(b1_b,a1_b,dA44_b);
    dAf44_b = filtfilt(b2_b,a2_b,dAf44_b);

    %figure()
    %plot(t,dA46,t,dAf46,'linewidth',2);
    % xlim([3 4])
    % grid on
    % box on
    % xlabel('Time [s]')
    % ylabel('Displacement [m]')
    % legend('Acc46-raw','Acc46-filtered')
    % title('Displacement accelerometer')


    %% Orbits
    %figure()
    l_f = length(aA46_f);
    start_f = ceil(1/3*l_f);
    endt_f = ceil(2/3*l_f);
    
    l_b = length(aA46_b);
    start_b = ceil(1/3*l_b);
    endt_b = ceil(2/3*l_b);
    %plot(dAf46(start:endt),vAf46(start:endt),'-b','linewidth',1);
    % grid on
    % %box on
    % xlabel('Displacement [m]')
    % ylabel('Velocity [m/2]')
    % legend('Acc46')
    % title('Orbit accelerometer forced Mass')

    %% FRF
    TF_f= islocalmax(dAf44_f);
    TT_f = islocalmin(dAf44_f);
    
    TF_b= islocalmax(dAf44_b);
    TT_b = islocalmin(dAf44_b);
    
    amplitude_f(i) = (mean(dAf44_f(TF_f))+abs(mean(dAf44_f(TT_f))))/2;
    amplitude_b(i) = (mean(dAf44_b(TF_b))+abs(mean(dAf44_b(TT_b))))/2;
end
plot(freq_sweep , amplitude_f,'-*b','LineWidth',2)
hold on 
plot(freq_sweep , amplitude_b,'-*r','LineWidth',2)
legend ('forward','backward')
grid on
%%


%figure()
%plot(t,dAf46,t(TF),dAf46(TF),t(TT),dAf46(TT))

% sf = [10:0.1:12.5];
% amp = [ 0.0049 0.0049 0.0050...
%     0.0050 0.0051 0.0051...
%     0.0052 0.0052 0.0053...
%     0.0053 0.0054 0.0054...
%     0.0054 0.0055 0.0056...
%     0.0056 0.0057 0.0057 ...
%     0.0058 0.0058 0.0059 ...
%     0.0059 0.0059 0.0060 0.0061 0.0061];
% p = polyfit(sf,amp,10);
% y_fit = polyval(p,sf);
% plot(sf,amp,'ro',sf,y_fit)

