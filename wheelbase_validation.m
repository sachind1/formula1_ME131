% Wheelbase Length Validation
% ME 131 Lab 7
% March 2018
% Author: Charlott Vallon

clear 
clc

% Load validation experiment signal
exp = load('15_2.mat','sig'); 
sig = exp.sig;

%% Simulate system
% simulate_system simulates the system forward using the kinematic bicycle
% model, and compares the simulated x-y coordinates to the true
% experimental x-y coordinates (as determined by a GPS). 

% N = # of re-initializations of the simulation to the true states (as
% determined by GPS)

% lf = distance from center of front wheels to center of mass (in [m])
[x,y] = simulate_system(sig, N, lf);


