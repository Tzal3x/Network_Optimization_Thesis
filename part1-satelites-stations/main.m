%% Space 3D:--------------------------------------------------------
%{

MATLAB R2018a

%}
% Use clc commmand to clear command window
cd C:\Users\User\Documents\GitHub\Network_Optimization_Thesis\part1-satelites-stations
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

 
NUMBER_OF_SATELLITES = 2; %4; %10; % integer, default 2
NUMBER_OF_STATIONS = 1; %2; %3; % integer, default 1
RANDOM_VELOCITIES = false; % boolean, default false
% INVERSE_VELOCITIES_SATEL = ones(1,NUMBER_OF_SATELLITES) * 1 * 100; % smaller value -> faster, it can be a vector of the desirable speeds [v1 v2 ... vn], where n == NUMBER_OF_SATELLITES
INVERSE_VELOCITIES_SATEL = [-1, 1]* 100; % 1000
INVERSE_VELOCITIES_STATIONS = ones(1,NUMBER_OF_STATIONS) * 800 * 1000; %80=moving, 800=almost imovable % larger value -> slower, >> >> >> >> >> >> >> >> >> >> >> >> 
STOP_AT_TIME = 3; % == EPOCHS integer, declares when the time should be stopped
% THETA_PHI = [0  20]; % 30, 60
LINK_CAPACITY = 30; % WARNING! LINK_CAPACITY must be equal to ...
COMMUNICATION_RANGE = 150; % Default: 150; (oldest)
PRINT_DETAILS = false; % true/false: Displays optimization problem's details (distance matrix, parameters (Aeq, beq, A, b, l, ...) etc)
PRINT_MAIN_PARAMETERS = true;
SHOW_TOPOLOGY = false; % 3D topology with spheres
THETA_PHI_ca = {}; % I have defined different inclinations of orbits for satellites and stations
% SOLVER = "heuristic_1"; % "fmincon"/"linprog"/"heuristic_1"
SOLVER = "linprog";
% SOLVER = "fmincon";
INIT_POS = [-70 -20 30];

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
% seeds = 1:NUMBER_OF_SATELLITES; % satel = 10, stations = 3.
% disp('> THETA_PHI: __________________')
% for i = 1:(NUMBER_OF_SATELLITES + NUMBER_OF_STATIONS)
%     if i <= NUMBER_OF_SATELLITES
%        rng(seeds(i))
%        random_factor1 = ceil(rand()*100-50); % [-50, 50]
%        random_factor2 = ceil(rand()*100-50); % [-50, 50]
%        THETA_PHI_ca{i} = [-random_factor1  random_factor2]; % set satellite orbits. It's not necessary to be the same for everyone.
%     else
%        THETA_PHI_ca{i} = [0 10]; % set station "orbits". Necessary to be the same.
%     end
%     disp('node: '+string(i))
%     disp(THETA_PHI_ca{i});
% end
% disp('_____________________________')
% ----------------------------------------------------------------------

THETA_PHI_ca{1} = [10 20];
THETA_PHI_ca{2} = [10 20];
THETA_PHI_ca{3} = [0 -30];

% THETA_PHI_ca{3} = [-60 60];
% THETA_PHI_ca{4} = [-80 10];
% 
% THETA_PHI_ca{5} = [0 10]; % station
% THETA_PHI_ca{5} = [0 10]; % station

dyo = solve_part2(NUMBER_OF_SATELLITES, NUMBER_OF_STATIONS, RANDOM_VELOCITIES, INVERSE_VELOCITIES_SATEL,...
                    INVERSE_VELOCITIES_STATIONS, STOP_AT_TIME, THETA_PHI_ca, LINK_CAPACITY, PRINT_DETAILS,...
                    PRINT_MAIN_PARAMETERS, SHOW_TOPOLOGY, COMMUNICATION_RANGE, SOLVER, INIT_POS);

 