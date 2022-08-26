%% Frequency response of a duffin oscillator
% This script integrated the ODE defined in the function duffin iteratively
% for different values of Omega to create the frequency response of a
% duffin oscillator with constant parameters
close all 
clear all
%% Parameter definition 

m=0.0087;    % Mass
zita = 0.0487;%damping ratio
wn = 11.41 *2*pi; %natural frequency
c=2*m*wn*zita;    % Damping constant 
k=44;     % Linear stiffness
alpha=0; % Nonlinear stiffness

    Omega=10.5*2*pi;  % Frequency of excitation 

    B= cell2mat(struct2cell(load('amplitude.mat'))) ;
    % B = [1.3259 2.0216];
    for i = [1:10] ;   % Amplitude of excitation 
        A = B(i);

        %% Parameter Storage 
        % All the parameters are stored in the variable PAR
        PAR=[m,c,k,alpha,A,Omega]; 

        %% State variable initialization 
        % Here we initialise the state variable y(1) and y(2). We hypothesize that
        % the system start from zero displacement and zero velocity

        INIT=[0 0];

        %% System integration 
        % We integrate the differential equation using ODE45

        [t,y]=ode45(@(t,y)duffin(t,y,PAR),[0,1000*2*pi/PAR(6)],INIT);

        %% Plot results 
        figure(1)
        dis = y(:,2);
        vel = y(:,1);
        acc = A*sin(Omega.*t);
        plot3(dis(end-1000:end),vel(end-1000:end),-acc(end-1000:end),'r','LineWidth',1)
        hold on
        xlabel('dis')
        ylabel('vel')
        legend('analytical plot')

    end
    

%% Frequency Response - Forward sweep
%
% we use a while loop to sweep the frequency from 4 $rad/s$ to 8 $rad/s$. We
% use the final state of the previous integration as initial condition for
% the next one
% figure %open a figure for the plots
OMf=[]; 
Df=[];
while PAR(6)<12*2*pi
    INIT=y(end,:);
    [t,y]=ode45(@(t,y)duffin(t,y,PAR),[0,500*2*pi/PAR(6)],INIT);
    mD=max(y(t>99*2*pi/PAR(6),2)); %find the max displacement in the last period
    figure(2)
    plot(PAR(6)/(2*pi),mD,'.b') %plot 
    xlabel('rad/sec')
    hold on % set the figure so that it does not delete the previous plots before drawing a new one
    drawnow % allow the figure to refresh at each iteration 
    OMf=[OMf;PAR(6)]; % store the frequencies in a growing vector
    Df=[Df;mD];       % store the maximum displacement in a growing vector
    PAR(6)=PAR(6)+0.1;% increment the frequency of a 0.1 rad/s
end
% %% Frequency Response - Backward sweep
% % We use a for loop to sweep the frequency backward from $8 rad/s$ to $4 rad/s$. We
% % use the final state of the previous integration as initial condition for
% % the next one
% OMb=[]; 
% Db=[];
% while PAR(6)>4
%     PAR(6)=PAR(6)-0.1;  %decrease the frequency or 0.1 rad/s
%     INIT=y(end,:);      %initialise the state variable to the last value reached in the previous integration 
%     [t,y]=ode45(@(t,y)duffin(t,y,PAR),[0,100*2*pi/PAR(6)],INIT);
%     mD=max(y(t>99*2*pi/PAR(6),2)); %find the max displacement in the last period
%     plot(PAR(6),mD,'.r')  %plot 
%     drawnow % allow the figure to refresh at each iteration 
%     OMb=[OMb;PAR(6)];% store the frequencies in a growing vector
%     Db=[Db;mD];% store the maximum displacement in a growing vector
% end
% %% 
% %We replot the results stored in the vectors to have continuous lines
% plot(OMf,Df,'b')  
% plot(OMb,Db,'r')
% xlabel('frequency')
% ylabel('displacement')