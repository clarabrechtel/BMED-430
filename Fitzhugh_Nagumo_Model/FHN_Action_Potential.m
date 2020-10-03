% file: FHN_Action_Potential.m
%	This is the Fitzhugh-Nagumo action potential modified to fit 
%   the dynamics of a cardiomyocyte.  
%
%   Inputs:  SimData, Stim (structures)
%   Outputs: Results (structure)
%
%	Author: James Eason
%	All rights reserved
%   Edited: Jonathon Stearns
%   Date: 9/1/2020
%
%   Edited: Clara Brechtel
%   Date: 9/23/2020

%% Define SimData and Stimulus structures 
clc
clear all
close all
% Define constants in a structure
SimData.epsilon = 0.0025;	% epsilon
SimData.beta = 0.5;         % beta
SimData.alpha = 10;         % alpha = threshold membrane voltage
SimData.Vp = 100;           % peak voltage
SimData.V_0 = 0;            % initial value of the membrane voltage, V
SimData.W_0 = 0;            % initial value of the recovery variable, W
SimData.endtime = 500;      % simulation duration (ms)
SimData.dt = 1;             % time step (ms)

% Stimulus
Stim.delay = 20;        % Stimulus begins at 0.2 ms
Stim.duration = 1;      % Stimulus lasts for 1 ms
Stim.magnitude = 50;    % Stimulus amplitude in mA/cm^3

Results = FHN_cardiac(SimData, Stim);
%% This function applies Fitzhugh-Nagumo action potential model 

function [Results] = FHN_cardiac(SimData, Stim)

% Assign initial conditions of V, W, and time
Results.V(1) = SimData.V_0;
Results.W(1) = SimData.W_0;
Results.time(1) = 0;

% calculate number of time steps required
steps = SimData.endtime/SimData.dt;      % Number of time steps to reach end 

% Begin time stepping 
for i = 1:steps
    Results.time(i+1) = i*SimData.dt;   
% define stimulus based on the time 
if (Results.time(i) >= Stim.delay) && (Results.time(i) <= (Stim.delay + Stim.duration))
    stim(i) = Stim.magnitude;
else 
    stim(i) = 0;
end

% calculate derivatives
dVdt(i) = (-Results.V(i)/SimData.Vp.^2)*(Results.V(i)-SimData.Vp)*(Results.V(i)-SimData.alpha)-Results.W(i)+stim(i);
dWdt(i) = SimData.epsilon*(SimData.beta*Results.V(i)-Results.W(i));
    
% Use Forward Euler method to calculate V and W at next time step
Results.V(i+1) = Results.V(i) + (SimData.dt*dVdt(i));
Results.W(i+1) = Results.W(i) + (SimData.dt*dWdt(i));

end   % end of time stepping loop

%% Find Action Potential Duration and print

j=1;
t_up=0;
t_down=0;
for j = 1:steps
   
    if (Results.V(j+1) > SimData.alpha) && (Results.V(j) <= SimData.alpha)
        t_up = Results.time(j+1);
    end
    if (Results.V(j) >= SimData.alpha) && (Results.V(j+1) < SimData.alpha)
        t_down = Results.time(j+1);
        APD = t_down - t_up
        break;
    end

j=j+1;
end

%% plot V and W vs time
hold off;
plot(Results.time,Results.V);
hold on;
plot(Results.time,Results.W);

% add title, legend, and axis labels
title(['AP Duration: ', num2str(APD), ' ms'], 'FontSize',16);  
legend('V','W');
xlabel('time(ms)', 'FontSize',14);
ylabel('transmembrane voltage (mV)', 'FontSize',14);

% Be sure that you return time, V, and W arrays 
% in the structure called Results

disp('Optimal Value for dt = 0.1 ms.');

end