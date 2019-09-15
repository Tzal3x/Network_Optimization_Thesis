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

NUMBER_OF_SATELLITES = 3; %4; %10; % integer, default 2
NUMBER_OF_STATIONS = 1; %2; %3; % integer, default 1
RANDOM_VELOCITIES = false; % boolean, default false
INVERSE_VELOCITIES_SATEL = ones(1,NUMBER_OF_SATELLITES) * 1 * 100; % smaller value -> faster, it can be a vector of the desirable speeds [v1 v2 ... vn], where n == NUMBER_OF_SATELLITES
% INVERSE_VELOCITIES_SATEL = [-1, 1]* 100; % 1000
INVERSE_VELOCITIES_STATIONS = ones(1,NUMBER_OF_STATIONS) * 800 * 1000; %80=moving, 800=almost imovable % larger value -> slower, >> >> >> >> >> >> >> >> >> >> >> >> 
STOP_AT_TIME = 20; % == EPOCHS integer, declares when the time should be stopped
% THETA_PHI = [0  20]; % 30, 60
LINK_CAPACITY = 20; % WARNING! LINK_CAPACITY must be equal to ...
COMMUNICATION_RANGE = 150; % Default: 150; (oldest)
PRINT_DETAILS = false; % true/false: Displays optimization problem's details (distance matrix, parameter s (Aeq, beq, A, b, l, ...) etc)
PRINT_MAIN_PARAMETERS = false;
SHOW_TOPOLOGY = false; % 3D topology with spheres
THETA_PHI_ca = {}; % I have defined different inclinations of orbits for satellites and stations

% SOLVER = "heuristic_1_fastest"; % "fmincon"/"linprog"/"heuristic_1_fastest"/"heuristic_1_closest"

% INIT_POS = [-30 -30 30]; % simple example for solver comparison
INIT_POS = 1:(NUMBER_OF_SATELLITES + NUMBER_OF_STATIONS)*800; % GO TO CREATE_NODES TO ACTIVATE IT
 
%{
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
Simple version were optimal is better than heuristic:

THETA_PHI_ca{1} = [10 20]; % satellite 1
THETA_PHI_ca{2} = [10 20]; % satellite 2
THETA_PHI_ca{3} = [0 -30]; % station
INVERSE_VELOCITIES_SATEL = [-1, 1]* 100; % 1000
INVERSE_VELOCITIES_STATIONS = ones(1,NUMBER_OF_STATIONS) * 800 * 1000;
STOP_AT_TIME = 2;
INIT_POS = [-30 -30 30];

GO TO create_nodes AND UNCOMMENT THE OF THE init_pos ARGUMENT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%}
% THETA_PHI_ca{1} = [10 20];
% THETA_PHI_ca{2} = [10 20];
% THETA_PHI_ca{3} = [0 -30];

%%%%%%% Demo 1: (same orbits) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:(NUMBER_OF_SATELLITES + NUMBER_OF_STATIONS)
    if i <= NUMBER_OF_SATELLITES
       THETA_PHI_ca{i} = [-10 20]; % set satellite orbits. It's not necessary to be the same for everyone.
    else
       THETA_PHI_ca{i} = [0 10]; % set station "orbits". Necessary to be the same.
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SOLVER = "heuristic_1_closest"; 
% SOLVER = "linprog";
% SOLVER = "fmincon";


dyo = solve_part2(NUMBER_OF_SATELLITES, NUMBER_OF_STATIONS, RANDOM_VELOCITIES, INVERSE_VELOCITIES_SATEL,...
                    INVERSE_VELOCITIES_STATIONS, STOP_AT_TIME, THETA_PHI_ca, LINK_CAPACITY, PRINT_DETAILS,...
                    PRINT_MAIN_PARAMETERS, SHOW_TOPOLOGY, COMMUNICATION_RANGE, SOLVER, INIT_POS);


disp("__________~ Results ~__________")
disp("[abs(total_utility), delay, execution_time]:")
disp(dyo)                
                

%%
%%%%%%%% Demo 2: (different orbits) --------------------------------------------
seeds = 1:NUMBER_OF_SATELLITES; % satel = 10, stations = 3.
% disp('> THETA_PHI: __________________')
for i = 1:(NUMBER_OF_SATELLITES + NUMBER_OF_STATIONS)
    if i <= NUMBER_OF_SATELLITES
       rng(seeds(i))
       random_factor1 = ceil(rand()*100-50); % [-50, 50]
       random_factor2 = ceil(rand()*100-50); % [-50, 50]
       THETA_PHI_ca{i} = [-random_factor1  random_factor2]; % set satellite orbits. It's not necessary to be the same for everyone.
    else
       THETA_PHI_ca{i} = [0 10]; % set station "orbits". Necessary to be the same.
    end
%     disp('node: '+string(i))
%     disp(THETA_PHI_ca{i});
end
% disp('_____________________________')
%%%%%%% ----------------------------------------------------------------------
res_present_text = ["Utility score: ", "Delay: ","Execution time: "];
% dyo = solve_part2(NUMBER_OF_SATELLITES, NUMBER_OF_STATIONS, RANDOM_VELOCITIES, INVERSE_VELOCITIES_SATEL,...
%                     INVERSE_VELOCITIES_STATIONS, STOP_AT_TIME, THETA_PHI_ca, LINK_CAPACITY, PRINT_DETAILS,...
%                     PRINT_MAIN_PARAMETERS, SHOW_TOPOLOGY, COMMUNICATION_RANGE, SOLVER, INIT_POS);
% disp(SOLVER + " results:")
% disp(res_present_text+string(dyo));
% disp("paused..."); pause;


%%%%% AUTOMATIC EXPERIMENT GENERATOR: =====================================
disp('AUTOMATIC EXPERIMENT GENERATOR ACTVATED:')
LINK_CAPACITY = 20; % CONSTANT
COMMUNICATION_RANGE = 150; % CONSTANT
MAX_TOTAL_EPOCHS = 10:10:100; % total epochs senarios
fid1 = fopen('C:\Users\User\Desktop\thesis_solver_results\solver_utility_results.txt', 'w');
fclose(fid1);
fid1 = fopen('C:\Users\User\Desktop\thesis_solver_results\solver_utility_results.txt', 'a+');
fprintf(fid1, '%s,%s,%s,%s,%s\n', ["Total_epochs","Node Group","heuristic_1_closest","heuristic_1_fastest","linprog"]); % Node groups are (sats=4,stats=lg4=2)

fid2 = fopen('C:\Users\User\Desktop\thesis_solver_results\solver_delay_results.txt', 'w');
fclose(fid2);
fid2 = fopen('C:\Users\User\Desktop\thesis_solver_results\solver_delay_results.txt', 'a+');
fprintf(fid2, '%s,%s,%s,%s,%s\n', ["Total_epochs","Node Group","heuristic_1_closest","heuristic_1_fastest","linprog"]); % Node groups are (sats=4,stats=lg4=2)

fid3 = fopen('C:\Users\User\Desktop\thesis_solver_results\solver_exectime_results.txt', 'w');
fclose(fid3);
fid3 = fopen('C:\Users\User\Desktop\thesis_solver_results\solver_exectime_results.txt', 'a+');
fprintf(fid3, '%s,%s,%s,%s,%s\n', ["Total_epochs","Node Group","heuristic_1_closest","heuristic_1_fastest","linprog"]); % Node groups are (sats=4,stats=lg4=2)
SATELLITE_VECTOR = [4,16,32,64]; % number of satellites per epoch
STATION_VECTOR = ceil(log2(SATELLITE_VECTOR)); % number of stations per epoch = log2 of #sats
for epoch = MAX_TOTAL_EPOCHS
    for numnodes = 2:2%1:length(SATELLITE_VECTOR)
        NUMBER_OF_SATELLITES = SATELLITE_VECTOR(numnodes); 
        NUMBER_OF_STATIONS =  STATION_VECTOR(numnodes);
        STOP_AT_TIME = epoch;
        INVERSE_VELOCITIES_SATEL = ones(1,NUMBER_OF_SATELLITES) * 1 * 100;
        INVERSE_VELOCITIES_STATIONS = ones(1,NUMBER_OF_STATIONS) * 800 * 1000; 
        INIT_POS = 1:(NUMBER_OF_SATELLITES + NUMBER_OF_STATIONS) * 800; % GO TO CREATE_NODES TO ACTIVATE IT
        
        disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        disp('Total epochs: ' + string(epoch))
        disp('SATELLITES: ' + string(NUMBER_OF_SATELLITES))
        disp('STATIONS: ' + string(NUMBER_OF_STATIONS))
        
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
%             disp('node: '+string(i))
%             disp(THETA_PHI_ca{i});
        end
        % Get results and append them in file to process them in R and make a
        % descriptive analysis.
        SOLVERS = ["heuristic_1_closest","heuristic_1_fastest","linprog"];
        utility_results = null(1,1);
        delay_results2 = null(1,1);
        exectime_results3 = null(1,1);
        for i = 1:length(SOLVERS)
            SOLVER = SOLVERS(i);  
            disp('Running solver: '+SOLVER)
            dyo = solve_part2(NUMBER_OF_SATELLITES, NUMBER_OF_STATIONS, RANDOM_VELOCITIES, INVERSE_VELOCITIES_SATEL,...
                                INVERSE_VELOCITIES_STATIONS, STOP_AT_TIME, THETA_PHI_ca, LINK_CAPACITY, PRINT_DETAILS,...
                                PRINT_MAIN_PARAMETERS, SHOW_TOPOLOGY, COMMUNICATION_RANGE, SOLVER, INIT_POS);
            utility_results = [utility_results, dyo(1)]; % first element of dyo is the utility score
            delay_results2 = [delay_results2, dyo(2)]; % scond element of dyo is the delay score
            exectime_results3 = [exectime_results3, dyo(3)]; % third element of dyo is the execution time
        end
        utility_results = [epoch, numnodes, utility_results];
        delay_results2 = [epoch, numnodes, delay_results2];
        exectime_results3 = [epoch, numnodes, exectime_results3];
        fprintf(fid1, '%i,%i,%f,%f,%f\n', utility_results); %numnodes specifies the group of number of nodes (e.g. 2 = [sats=16, stats=4])
        fprintf(fid2, '%i,%i,%f,%f,%f\n', delay_results2); %numnodes specifies the group of number of nodes (e.g. 2 = [sats=16, stats=4])
        fprintf(fid3, '%i,%i,%f,%f,%f\n', exectime_results3); %numnodes specifies the group of number of nodes (e.g. 2 = [sats=16, stats=4])
    end
end
fclose(fid1); fclose(fid2); fclose(fid3);