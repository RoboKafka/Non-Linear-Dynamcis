%% aquiring all data
close all
clear all
clc

% acc 44 = tip
% acc 46 = base
% forward sweep

% A = load('10.5_0.100_f.mat');
% A(:,:,2) = load('10.5_0.200_f.mat');
% A(:,:,3) =load('10.5_0.300_f.mat');
% A(:,:,4) = load('10.5_0.400_f.mat');
% A(:,:,5) = load('10.5_0.500_f.mat');
% A(:,:,6) = load('10.5_0.600_f.mat');
% A(:,:,7) = load('10.5_0.700_f.mat');
% A(:,:,8) = load('10.5_0.800_f.mat');
% A(:,:,9) = load('10.5_0.900_f.mat');
% A(:,:,10) = load('10.5_1.000_f.mat');

A = load('10p7_f_0.1.mat');
A(:,:,2) = load('10p7_f_0.2.mat');
A(:,:,3) =load('10p7_f_0.3.mat');
A(:,:,4) = load('10p7_f_0.4.mat');
A(:,:,5) = load('10p7_f_0.5.mat');
A(:,:,6) = load('10p7_f_0.6.mat');
A(:,:,7) = load('10p7_f_0.7.mat');
A(:,:,8) = load('10p7_f_0.8.mat');
A(:,:,9) = load('10p7_f_0.9.mat');
A(:,:,10) = load('10p7_f_1.0.mat');
A(:,:,11) = load('10p7_f_1.1.mat');
A(:,:,12) = load('10p7_f_1.2.mat');
A(:,:,13) =load('10p7_f_1.3.mat');
A(:,:,14) = load('10p7_f_1.4.mat');
A(:,:,15) = load('10p7_f_1.5.mat');
A(:,:,16) = load('10p7_f_1.6.mat');
A(:,:,17) = load('10p7_f_1.7.mat');
A(:,:,18) = load('10p7_f_1.8.mat');
A(:,:,19) = load('10p7_f_1.9.mat'); 
A(:,:,20) = load('10p7_f_2.0.mat');
%backward sweep
% A = load('C:\Users\rohit\Desktop\FRF_thess\model_1\0.600_9.7-11.5-1301pm\10p7_b.mat');
% A(:,:,2) = load('C:\Users\rohit\Desktop\FRF_thess\model_1\0.900_10.5-11.5-1334pm\10p7_b.mat');
% A(:,:,3) =load('C:\Users\rohit\Desktop\FRF_thess\model_1\1.200_10.5-11.5-1334pm\10p7_b.mat');
% A(:,:,4) = load('C:\Users\rohit\Desktop\FRF_thess\model_1\1.500_10.5-11.5-1334pm\10p7_b.mat');
% A(:,:,5) = load('C:\Users\rohit\Desktop\FRF_thess\model_1\1.600_10.5-11.5-1334pm\10p7_b.mat');
% A(:,:,6) = load('C:\Users\rohit\Desktop\FRF_thess\model_1\1.700_10.5-11.5-1334pm\10p7_b.mat');
% A(:,:,7) = load('C:\Users\rohit\Desktop\FRF_thess\model_1\1.800_10.5-11.5-1334pm\10p7_b.mat');
% A(:,:,8) = load('C:\Users\rohit\Desktop\FRF_thess\model_1\1.900_10.5-11.5-1334pm\10p7_b.mat');
% A(:,:,9) = load('C:\Users\rohit\Desktop\FRF_thess\model_1\2.000_10.5-11.5-1334pm\10p7_b.mat');
% A(:,:,10) = load('C:\Users\rohit\Desktop\FRF_thess\model_1\2.100_10.5-11.5-1334pm\10p7_b.mat');
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
XX_max = [];
XX_min = [];
FX_max = [];
FX_min = [];
XX_dot_max = [];
XX_dot_min = [];
FV_max = [];
FV_min = [];
for i = [1:10] % frequency sweep
    % data extraction
    t_f = cell2mat(struct2cell(A(:,:,i)));
    t_f = t_f(:,1);% time [s]
    dt_f = t_f(2)-t_f(1); % delta t [s]
    fs_f = 1/dt_f; % Sampling Frequency [Hz]
    
    % acc 46
    volt_acc46_f = cell2mat(struct2cell(A(:,:,i)));% accelerometrt voltage [V]
    volt_acc46_f = volt_acc46_f(:,3);
    aA46_f = ((volt_acc46_f-v046/1000)/(sens_46/1000))*g; % acceleromenter acceleration [m/s^2]

   
    %acc 44
    volt_acc44_f = cell2mat(struct2cell(A(:,:,i)));% accelerometrt voltage [V]
    volt_acc44_f = volt_acc44_f(:,2);
    aA44_f = ((volt_acc44_f-v044/1000)/(sens_44/1000))*g; % acceleromenter acceleration [m/s^2]
    
    %% Filter
    fc_h_f = 3/(fs_f/2);  % Cut off Frequency (high)1
    fc_l_f = 100/(fs_f/2);  % Cut off Frequency (low)
    order = 6; % 6th Order Filter
    [b1_f,a1_f] = butter(order,fc_h_f,'high');
    [b2_f,a2_f] = butter(order,fc_l_f,'low');
    
    %% Accelerometer

    % Acceleration accelerometer
    aAf44_f = filtfilt(b1_f,a1_f,aA44_f);
    aAf44_f = filtfilt(b2_f,a2_f,aAf44_f);

    aAf46_f = filtfilt(b1_f,a1_f,aA46_f);
    aAf46_f = filtfilt(b2_f,a2_f,aAf46_f)
   
    lacc = aAf46_f(islocalmin(aAf46_f));
    min_idx = lacc<0; % create logical index
    lacc = mean(lacc(min_idx));
    uacc = mean(aAf46_f(islocalmax(aAf46_f)));
    acc_amp(i) = (uacc + abs(lacc))/2;
    
    %% Velocity from ACC
    vA46_f = cumtrapz(t_f,aAf46_f);
    vA46_f = detrend(vA46_f);
    vAf46_f = filtfilt(b1_f,a1_f,vA46_f);
    vAf46_f = filtfilt(b2_f,a2_f,vAf46_f);

    vA44_f = cumtrapz(t_f,aAf44_f);
    vA44_f = detrend(vA44_f);
    vAf44_f = filtfilt(b1_f,a1_f,vA44_f);
    vAf44_f= filtfilt(b2_f,a2_f,vAf44_f);
   
    %% Displacement from velocity
    dA46_f = cumtrapz(t_f,vAf46_f);
    dA46_f = detrend(dA46_f);
    dAf46_f = filtfilt(b1_f,a1_f,dA46_f);
    dAf46_f = filtfilt(b2_f,a2_f,dAf46_f);
    dAf46_f(1:10)

    dA44_f = cumtrapz(t_f,vAf44_f);
    dA44_f = detrend(dA44_f);
    dAf44_f = filtfilt(b1_f,a1_f,dA44_f);
    dAf44_f = filtfilt(b2_f,a2_f,dAf44_f);
    

    %% Orbits
    %figure()
    l_f = length(aA44_f);
    start_f = ceil(1/4*l_f);
    endt_f = ceil(3/4*l_f);
    

    X = (dAf44_f(start_f:endt_f)-dAf46_f(start_f:endt_f));% absolute dis
    Y = (vAf44_f(start_f:endt_f)-vAf46_f(start_f:endt_f)); % absolute vel
    Z = -aAf44_f(start_f:endt_f)*0.008290;% relative acc
    
   
    %% surface section to plot force -displacment
    delta = 0.01;
    xx = zeros(length(X),1);
    fx = zeros(length(Z),1);
    for j = 1:length(Y)
        if abs(Y(j))<=delta
            xx(j,1)=X(j);
            fx(j,1)=Z(j);
        else
            xx(j,1) = 100;
            fx(j,1) = 100;
        end
    end
    xx(xx==100)=[];
    fx(fx==100)=[];
    
    XX_max=[XX_max;max(xx)];
    XX_min =[XX_min;min(xx)];
    FX_max= [FX_max;max(fx)];
    FX_min=[FX_min;min(fx)];
    
    figure(1)% sliced at velocity = 0
    
    plot(xx,fx,'*r','LineWidth',0.2)
    xlabel('tip displacement')
    ylabel('relative acceleration')
    title('displacement vs acceleration')
    hold on
    grid on
    
   
    %% surface section to plot force - Velocity
    delta_s = 0.00001;
    xx_dot = zeros(length(Y),1);
    fv = zeros(length(Z),1);
    
    for k = 1:length(X)
        if abs(X(k)) <= delta_s
            xx_dot(k,1) = Y(k);
            fv(k,1) = Z(k);   
             
        else
            xx_dot(k,1) = 999;
            fv(k,1) = 999;
        end
    end
    
    xx_dot(xx_dot == 999) = [];
    fv(fv == 999) = [];
   
    XX_dot_max=[XX_dot_max;max(xx_dot)];
    XX_dot_min =[XX_dot_min;min(xx_dot)];
    FV_max= [FV_max;max(fv)];
    FV_min=[FV_min;min(fv)];

    figure(2)
    subplot(2,1,2),plot(xx_dot,fv,'*r','LineWidth',0.2)
    hold on
    grid on
    xlabel('tip velocity')
    ylabel('relative acceleration')
    title('veloity vs acceleration')
    legend('Damping property')
    subplot(2,1,1),plot(xx,fx,'*b','LineWidth',0.2)
    xlabel('tip displacement')
    ylabel('relative acceleration')
    title('displacement vs acceleration')
    legend('stiffness property')
    hold on
    grid on
    
    %%
    figure(3) %acc vs displcement vs velocity
    plot3(X,Y,Z,'--b','linewidth',0.1);
    hold on
    grid on
    xlabel('tip displacement','FontSize', 12)
    ylabel('tip velocity','FontSize', 12)
    zlabel('relative acceleration','FontSize', 12)
    title('Restoring Force Surface')
    
    figure(4)% displacement vs velocity
    plot(X,Y,'--b','linewidth',0.1);
    hold on
    grid on
    ylabel('tip velocity' ,'FontSize', 20)
    xlabel('tip displacement', 'FontSize', 20)
    legend('experimental', 'FontSize', 20)
    title('Orbit')
    
%     TF_f= islocalmax(dAf46_f);
%     TT_f = islocalmin(dAf46_f);
    
end
%% polyfit for damping and stiffness
 XX = [flip(XX_max);XX_min];
 FX = [flip(FX_max);FX_min];
 XX_dot = [flip(XX_dot_max);XX_dot_min];
 FV = [flip(FV_max);FV_min];
 
 [p,pint] = polyfit(XX,FX,2);
 f = polyval(p,XX);
 T_1 = table(XX,FX,f,FX-f,'VariableNames',{'X','Y','Fit','FitError'})
 p_eq = poly2sym(p)

 [q,qint] = polyfit(XX_dot,FV,3);
 q_ = polyval(q,XX_dot);
 T_2 = table(XX_dot,FV,q_,FV-q_,'VariableNames',{'X','Y','Fit','FitError'})
 q_eq = poly2sym(q)


 figure(5)
 
 subplot(2,1,1),plot(XX,FX,'*b',XX,f,'-k','LineWidth',2)
 xlabel('Displacement','FontSize', 12)
 ylabel('Acceleration','FontSize', 12)
 title('stiffness')
 legend('aqq-data','polyfit order 1')
 grid on
 subplot(2,1,2),plot(XX_dot,FV,'*r',XX_dot,q_,'-k','LineWidth',2)
 xlabel('Velocity','FontSize', 12)
 ylabel('Acceleration','FontSize', 12)
 title('damping')
 legend('aqq-data','fitdata')
 legend('aqq-data','polyfit order 3')
 grid on
 
