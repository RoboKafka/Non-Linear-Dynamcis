%% Frequency response of a duffin oscillator
% This script integrated the ODE defined in the function duffin iteratively
% for different values of Omega to create the frequency response of a
% duffin oscillator with constant parameters
close all 
clear all
clc

%% Parameter definition 

for A = 0.0004
    m=0.0087;      % Mass
    c=0.0595;    % Damping constant 
    k=42.6075;     % Linear stiffness
     % Nonlinear stiffness
        % Amplitude of excitation 
    Omega=10.7;  % Frequency of excitation 

    %% Parameter Storage 
    % All the parameters are stored in the variable PAR
    PAR=[m,c,k,A,Omega]; 

    %% State variable initialization 
    % Here we initialise the state variable y(1) and y(2). We hypothesize that
    % the system start from zero displacement and zero velocity

    INIT=[0 0];

    %% System integration 
    % We integrate the differential equation using ODE45

    [t,y]=ode45(@(t,y)duffin(t,y,PAR),[0,100*2*pi/PAR(5)],INIT);

    %% Plot results 
    
    figure
    plot(y(:,1),y(:,2))
    xlabel('velocity')
    ylabel('displacement')
    
end
