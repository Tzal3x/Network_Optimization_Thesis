%% Space 3D:--------------------------------------------------------
%{

MATLAB R2018a

%}
% Use clc commmand to clear command window
cd C:\Users\User\Documents\GitHub\Network_Optimization_Thesis\part1-satelites-stations
addpath C:\Users\User\Documents\GitHub\Network_Optimization_Thesis\part1-satelites-stations\m_map %adding m_map package
addpath C:\Users\User\Documents\GitHub\Network_Optimization_Thesis\part1-satelites-stations\thesis_matlab_library % adding thesis_matlab_library

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% - - - - - - - - - - - + + + + + + \ M E N U / + + + + + + - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Main parameters (if you want to experiment with program's parameters just change ONLY the following): - - -  

% NUMBER_OF_SATELLITES = 2; %5 integer, default 4
% NUMBER_OF_STATIONS = 1; %3 integer, default 2
% RANDOM_VELOCITIES = false; % boolean, default false
% INVERSE_VELOCITIES_SATEL = ones(1,NUMBER_OF_SATELLITES) * 3; % smaller value -> faster, it can be a vector of the desirable speeds [v1 v2 ... vn], where n == NUMBER_OF_SATELLITES
% INVERSE_VELOCITIES_STATIONS = ones(1,NUMBER_OF_STATIONS) * 80; % larger value -> slower, >> >> >> >> >> >> >> >> >> >> >> >> 
% STOP_AT_TIME = 30;% integer, declares when the time should be stopped
% THETA_PHI = [30  60];
% LINK_CAPACITY = 100; % WARNING! LINK_CAPACITY must be equal to ...
% PRINT_DETAILS = true; % true/false: Displays optimization problem's details (distance matrix, parameters (Aeq, beq, A, b, l, ...) etc)
% PRINT_MAIN_PARAMETERS = true;
% SHOW_TOPOLOGY = true;
% 
% ena = solve_part1(NUMBER_OF_SATELLITES, NUMBER_OF_STATIONS, RANDOM_VELOCITIES, INVERSE_VELOCITIES_SATEL,...
%                     INVERSE_VELOCITIES_STATIONS, STOP_AT_TIME, THETA_PHI, LINK_CAPACITY, PRINT_DETAILS, PRINT_MAIN_PARAMETERS,SHOW_TOPOLOGY);

 
NUMBER_OF_SATELLITES = 4; %4; %10; % integer, default 2
NUMBER_OF_STATIONS = 2; %2; %3; % integer, default 1
RANDOM_VELOCITIES = false; % boolean, default false
INVERSE_VELOCITIES_SATEL = ones(1,NUMBER_OF_SATELLITES) * 1; % smaller value -> faster, it can be a vector of the desirable speeds [v1 v2 ... vn], where n == NUMBER_OF_SATELLITES
INVERSE_VELOCITIES_STATIONS = ones(1,NUMBER_OF_STATIONS) * 800; %80=moving, 800=almost imovable % larger value -> slower, >> >> >> >> >> >> >> >> >> >> >> >> 
STOP_AT_TIME = 3; %20; % == EPOCHS integer, declares when the time should be stopped
% THETA_PHI = [0  20]; % 30, 60
LINK_CAPACITY = 100; % WARNING! LINK_CAPACITY must be equal to ...
COMMUNICATION_RANGE = 150;
PRINT_DETAILS = true; % true/false: Displays optimization problem's details (distance matrix, parameters (Aeq, beq, A, b, l, ...) etc)
PRINT_MAIN_PARAMETERS = true;
SHOW_TOPOLOGY = false; % 3D topology with spheres
THETA_PHI_ca = {}; % I have defined different inclinations of orbits for satellites and stations

% Demo 1: (same orbits) -------------------------------------------------
% for i = 1:(NUMBER_OF_SATELLITES + NUMBER_OF_STATIONS)
%     if i <= NUMBER_OF_SATELLITES
%        THETA_PHI_ca{i} = [-10 20]; % set satellite orbits. It's not necessary to be the same for everyone.
%     else
%        THETA_PHI_ca{i} = [0 10]; % set station "orbits". Necessary to be the same.
%     end
% end
% -----------------------------------------------------------------------

% Demo 2: (different orbits) --------------------------------------------
seeds = 1:NUMBER_OF_SATELLITES; % satel = 10, stations = 3.
for i = 1:(NUMBER_OF_SATELLITES + NUMBER_OF_STATIONS)
    if i <= NUMBER_OF_SATELLITES
       rng(seeds(i))
       random_factor1 = ceil(rand()*100-50); % [-50, 50]
       random_factor2 = ceil(rand()*100-50); % [-50, 50]
       THETA_PHI_ca{i} = [-random_factor1  random_factor2]; % set satellite orbits. It's not necessary to be the same for everyone.
    else
       THETA_PHI_ca{i} = [0 10]; % set station "orbits". Necessary to be the same.
    end
end
% ----------------------------------------------------------------------

% THETA_PHI_ca{1} = [-30 50];
% THETA_PHI_ca{2} = [-50 20];
% THETA_PHI_ca{3} = [-60 60];
% THETA_PHI_ca{4} = [-80 10];
%
% THETA_PHI_ca{5} = [0 10]; % station


dyo = solve_part2(NUMBER_OF_SATELLITES, NUMBER_OF_STATIONS, RANDOM_VELOCITIES, INVERSE_VELOCITIES_SATEL,...
                    INVERSE_VELOCITIES_STATIONS, STOP_AT_TIME, THETA_PHI_ca, LINK_CAPACITY, PRINT_DETAILS, PRINT_MAIN_PARAMETERS, SHOW_TOPOLOGY, COMMUNICATION_RANGE);

    
% Plot 2D world map and links
%{
 /!\ UNDER CONSTRUCTION ...

 THIS SCRIPT IS FOR THE VISUALIZATION OF THE SATELITE ORBITS ON A 2D MAP
 OF THE EARTH. BEFORE CONTINUING IT SHOULD BE CONSIDERED THAT THERE IS NO 
 COMPLETE MATCHING OF THE 3D ORBITS TO THE 2D MAPPING BECAUSE STATIONS 
 SHOULD BE REGARDED AS STATIONARY (ROTATION OF THE EARTH WILL BE IGNORED) 
 TO KEEP THE VELOCITY CALCULATIONS OF EACH SATELITE SIMPLE. 
 m_map PACKAGE SHOULD BE USEFUL FOR THIS KIND OF JOB.
%}
 


