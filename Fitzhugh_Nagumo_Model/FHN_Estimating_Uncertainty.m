% file: FHN_Estimating_Uncertainty.m
%	This is the Fitzhugh-Nagumo action potential modified to fit 
%   the dynamics of a cardiomyocyte. An evaluation for this model 
%   has been added using the Monte Carlo Model. 
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
%   Date: 10/1/2020

%% Define SimData and Stimulus structures 
clc
clear all
close all

% Define constants in a structure
SimData.Vp = 100;           % peak voltage
SimData.V_0 = 0;            % initial value of the membrane voltage, V
SimData.W_0 = 0;            % initial value of the recovery variable, W
SimData.endtime = 500;      % simulation duration (ms)
SimData.dt = 0.01;             % time step (ms)

% Stimulus
Stim.delay = 20;        % Stimulus begins at 0.2 ms
Stim.duration = 1;      % Stimulus lasts for 1 ms
Stim.magnitude = 50;    % Stimulus amplitude in mA/cm^3

%% This section implements Monte Carlo method 
% Define uncertain parameters 
alpha_mean = 10.00;
alpha_sigma = 0.93;
beta_mean = 0.500;
beta_sigma = 0.064;
epsilon_mean = 0.0025;
epsilon_sigma = 0.00007;

N = 1000;               % Number of Monte Carlo simulations 
% Make probably density functions for the 6 parameters with 2% sigma
alpha_norm_dist = makedist('Normal','mu',alpha_mean,'sigma',alpha_sigma);
beta_norm_dist = makedist('Normal','mu',beta_mean,'sigma',beta_sigma);
epsilon_norm_dist = makedist('Normal','mu',epsilon_mean,'sigma',epsilon_sigma);

% Create distributions for the sample measurements of the length and width
alpha = random(alpha_norm_dist,N,1);
beta = random(beta_norm_dist,N,1);
epsilon = random(epsilon_norm_dist,N,1);

% Running a Monte Carlo simulation for the area calculation

h = waitbar(0,'Running Monte Carlo simulation...');
tic;
for i=1:N
    if (~mod(i,10)) waitbar(i/N,h,sprintf('Running %d Monte Carlo simulations...', N));end
    SimData.alpha = alpha(i);
    SimData.beta = beta(i);
    SimData.epsilon = epsilon(i);
    Results = FHN_cardiac_2(SimData, Stim);
    apd(i) = Results.apd;
end
Runtime = toc;
close(h);

% Calculate mean, SD, and SE of apd 
apd = reshape(apd,1000,1);
apd_mean = mean(apd);
apd_SD = std(apd);
apd_SE = apd_SD/sqrt(N);

% Determine correlation between apd and three parameters 
apd_alpha = corr(alpha, apd);
apd_beta = corr(beta, apd);
apd_epsilon = corr(epsilon, apd);

% Find best correlated variable 
if (abs(apd_alpha) > abs(apd_beta)) && (abs(apd_alpha) > abs(apd_epsilon))
    rand_variable = alpha;
    CC = apd_alpha;
    var = 'Alpha';
elseif (abs(apd_beta) > abs(apd_alpha)) && (abs(apd_beta) > abs(apd_epsilon))
    rand_variable = beta;
    CC = apd_beta;
    var = 'Beta';
else 
    rand_variable = epsilon;
    CC = apd_epsilon;
    var = 'Epsilon';
end
        
% Output results and plot the histogram for calculated area
subplot(2,1,1);
histogram(apd,25);
title(sprintf('Estimated APD (N=%d)',N));
xlabel('Action Potential Duration (ms)');
ylabel('Occurrences');

subplot(2,1,2);
scatter(rand_variable, apd, 'filled');
hold on;
trendline = fit(rand_variable, apd, 'poly1');
plot(trendline, rand_variable, apd);
title(sprintf('Sensitivity Analysis: CC = %f', CC));
xlabel(sprintf('Correlation Coefficient: %s', var));
ylabel('Action Potential Duration (ms)');
legend({'APD','Trendline'});

fprintf('\n Loop Runtime = %d [seconds]', Runtime);
fprintf('\nSample size = %d', N);
fprintf('Mean APD = %6.4f cm^2, std dev = %6.4f cm^2', apd_mean, apd_SD);
fprintf('Standard error of the estimated mean = %6.4f cm^2',apd_SE);
fprintf('Correlation coefficient of alpha = %6.4f, beta = %6.4f, epsilon = %6.4f', apd_alpha, apd_beta, apd_epsilon);

%% This is the function that calculates the APD, magnitude, and "shape"

function [Results] = FHN_cardiac_2(SimData, Stim)

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

%% Calculate APD 

j=1;
t_up=0;
t_down=0;
for j = 1:steps
   
    if (Results.V(j+1) > SimData.alpha) && (Results.V(j) <= SimData.alpha)
        t_up = Results.time(j+1);
    end
    if (Results.V(j) >= SimData.alpha) && (Results.V(j+1) < SimData.alpha)
        t_down = Results.time(j+1);
        if (t_down ~= 0) && (t_up ~= 0)
            Results.apd = t_down - t_up;
        else
            Results.apd = 0;
        end 
        break;
    end
j=j+1;
end

end
