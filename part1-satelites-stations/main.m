%% Space 3D:--------------------------------------------------------
%{

MATLAB R2018a

%}

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

 
NUMBER_OF_SATELLITES = 2; % integer, default 2
NUMBER_OF_STATIONS = 1; % integer, default 1
RANDOM_VELOCITIES = false; % boolean, default false
INVERSE_VELOCITIES_SATEL = ones(1,NUMBER_OF_SATELLITES) * 3; % smaller value -> faster, it can be a vector of the desirable speeds [v1 v2 ... vn], where n == NUMBER_OF_SATELLITES
INVERSE_VELOCITIES_STATIONS = ones(1,NUMBER_OF_STATIONS) * 80; % larger value -> slower, >> >> >> >> >> >> >> >> >> >> >> >> 
STOP_AT_TIME = 3;% == EPOCHS integer, declares when the time should be stopped
THETA_PHI = [30  60];
LINK_CAPACITY = 100; % WARNING! LINK_CAPACITY must be equal to ...
PRINT_DETAILS = true; % true/false: Displays optimization problem's details (distance matrix, parameters (Aeq, beq, A, b, l, ...) etc)
PRINT_MAIN_PARAMETERS = true;
SHOW_TOPOLOGY = false;

dyo = solve_part2(NUMBER_OF_SATELLITES, NUMBER_OF_STATIONS, RANDOM_VELOCITIES, INVERSE_VELOCITIES_SATEL,...
                    INVERSE_VELOCITIES_STATIONS, STOP_AT_TIME, THETA_PHI, LINK_CAPACITY, PRINT_DETAILS, PRINT_MAIN_PARAMETERS, SHOW_TOPOLOGY);



                
%% Plot 2D world map and links
%{
 /!\ UNDER CONSTRUCTION ...

 THIS SCRIPT IS FOR THE VISUALIZATION OF THE SATELITE ORBITS ON A 2D MAP
 OF THE EARTH. BEFORE CONTINUING IT SHOULD BE CONSIDERED THAT THERE IS NO 
 COMPLETE MATCHING OF THE 3D ORBITS TO THE 2D MAPPING BECAUSE STATIONS 
 SHOULD BE REGARDED AS STATIONARY (ROTATION OF THE EARTH WILL BE IGNORED) 
 TO KEEP THE VELOCITY CALCULATIONS OF EACH SATELITE SIMPLE. 
 m_map PACKAGE SHOULD BE USEFUL FOR THIS KIND OF JOB.
%}



