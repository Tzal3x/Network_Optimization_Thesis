%% Space 3D:--------------------------------------------------------
%{

MATLAB R2018a

%}

cd C:\Users\User\Documents\GitHub\Network_Optimization_Thesis\part1-satelites-stations
addpath C:\Users\User\Documents\GitHub\Network_Optimization_Thesis\part1-satelites-stations\m_map %adding m_map package

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% - - - - - - - - - - - + + + + + + \ M E N U / + + + + + + - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Main parameters (if you want to experiment with program's parameters just change ONLY the following): - - -  
NUMBER_OF_SATELLITES = 1; %5 integer, default 2
NUMBER_OF_STATIONS = 2; %2 integer, default 1
RANDOM_VELOCITIES = false; % boolean, default false
INVERSE_VELOCITIES_SATEL = ones(1,NUMBER_OF_SATELLITES) * 3; % smaller value -> faster, it can be a vector of the wanted speeds [v1 v2 ... vn], where n == NUMBER_OF_SATELLITES
INVERSE_VELOCITIES_STATIONS = ones(1,NUMBER_OF_STATIONS) * 80; % larger value -> slower, >> >> >> >> >> >> >> >> >> >> >> >> 
STOP_AT_TIME = 30;% integer, declares when the time should be stopped
THETA_PHI = [30  60];
LINK_CAPACITY = 100; % WARING! LINK_CAPACITY must be equal to 
PRINT_DETAILS = true; % true/false: Displays optimization problem's details (distance matrix, parameters (Aeq, beq, A, b, l, ...) etc)
PRINT_MAIN_PARAMETERS = true;
SHOW_TOPOLOGY = true;

ena = kyria(NUMBER_OF_SATELLITES, NUMBER_OF_STATIONS, RANDOM_VELOCITIES, INVERSE_VELOCITIES_SATEL,...
                    INVERSE_VELOCITIES_STATIONS, STOP_AT_TIME, THETA_PHI, LINK_CAPACITY, PRINT_DETAILS, PRINT_MAIN_PARAMETERS,SHOW_TOPOLOGY);

function out = kyria(NUMBER_OF_SATELLITES, NUMBER_OF_STATIONS, RANDOM_VELOCITIES, INVERSE_VELOCITIES_SATEL,...
                    INVERSE_VELOCITIES_STATIONS, STOP_AT_TIME, THETA_PHI, LINK_CAPACITY, PRINT_DETAILS, PRINT_MAIN_PARAMETERS,SHOW_TOPOLOGY)
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    disp(' ')% line-break to seperate previous runs of the programm.
    disp(' ')
    disp(' ')
    disp(' ')
    disp("|========================================= S T A R T =============================================|")
    nodes = create_nodes(NUMBER_OF_SATELLITES, NUMBER_OF_STATIONS, INVERSE_VELOCITIES_SATEL, ... 
                         INVERSE_VELOCITIES_STATIONS, RANDOM_VELOCITIES, THETA_PHI); 
%     prompt = '[~I/O:] Display the program''s main parameters? 1/0 (one=yes, zero=no)...';
%     PRINT_MAIN_PARAMETERS = input(prompt); % Boolean; If true, displays each one of the above parameters


    if PRINT_MAIN_PARAMETERS
       disp("NUMBER_OF_SATELLITES: " + string(NUMBER_OF_SATELLITES))
       disp("NUMBER_OF_STATIONS: " + string(NUMBER_OF_STATIONS))
       disp("RANDOM_VELOCITIES: " + string(RANDOM_VELOCITIES))
       disp("INVERSE_VELOCITIES_SATEL: " + print_list(INVERSE_VELOCITIES_SATEL))
       disp("INVERSE_VELOCITIES_STATIONS: " + print_list(INVERSE_VELOCITIES_STATIONS))
       disp("STOP_AT_TIME: " + string(STOP_AT_TIME))
       disp("THETA_PHI: " + print_list(THETA_PHI))
       disp("LINK_CAPACITY: " + string(LINK_CAPACITY))
       disp("PRINT_DETAILS: " + string(PRINT_DETAILS))
    end

    % Setting up figure display options: - - - - - - - - - - - - - - - - -
    earth_radius = 200;
    if SHOW_TOPOLOGY
    figure('Name','3D Simulation of satelite orbits');
    title("satelites:"+string(NUMBER_OF_SATELLITES)+...
          ", stations:"+string(NUMBER_OF_STATIONS)+...
          ", random velocities:"+string(RANDOM_VELOCITIES)+...
          ", stop:"+string(STOP_AT_TIME));
        hold on; % keep plotting on the existing figure
        axis equal
        view(30,0)% setting azemuth and elevation angles of camera for nicer visualization 
    axis([-earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200]); %setting up figure size
    end    
    
    %%%% WARNING! 'times' must be always < lesser from linspace length of all satelites and stations (avoiding index out of bounds error)
    stop = STOP_AT_TIME; % 30
    times = 10000; % upper bound of iterations
    [x,y,z] = sphere();%get coordinates in order to create 3D spheres. A sphere can be the earth, a satelite or a station.   
    for i = 1:times
       %Coordinates ---------------------------
       coords = [];
       for j = 1:length(nodes)
          coords = [coords; nodes(j).lifetime_coordinates(:,i)' ];
       end
       
       if SHOW_TOPOLOGY
           %Visualize ---------------------------
           earth = surf( earth_radius*x, earth_radius*y, earth_radius*z );
           displays = [];
           for j = 1:length(nodes)
              if j <= NUMBER_OF_SATELLITES
                displays = [displays, surf(coords(j,1)+10*x,coords(j,2)+10*y,coords(j,3)+10*z)]; 
              else
                displays = [displays, surf(coords(j,1)+20*x,coords(j,2)+20*y,coords(j,3)+20*z)]; 
              end
           end
       end
           % Placing here the 'stop' break for the figure to stay active with all
           % objects:
           if i == stop
               break
           end
       if SHOW_TOPOLOGY
           %Delete (to create the animation illusion)
           pause(0.01) %WARNING: all 'delete' (for graphic objects) functions must be after the 'pause' function 
           for j = 1:length(displays)
              delete(displays(j))
           end
           delete(earth)
       end
    end
    % hold off %not necessary(?)
    disp('[~Report:] Time stopped! A network has been created.')

    % Creating objective function and constraints: ========================================================================================================================================================================================================================================================
    disp('Creating objective function and constraints...')
    n = NUMBER_OF_SATELLITES + NUMBER_OF_STATIONS;% number of total nodes
    DISTANCES = zeros(n,n); %(N+M)x(N+M) , N: #sats, M: #stats
    LINKS = zeros(n,n); %(N+M)x(N+M)^2
    communication_range = 150; % links exist at a distance smaller or equal of communication_range
    for i = 1:n
        for j = 1:n
            DISTANCES(i,j) = euclidean_dist(coords(i,:),coords(j,:)); 
            if euclidean_dist(coords(i,:),coords(j,:)) <= communication_range %200 is arbitrary
                if strcmp(nodes(i).name,'satelite') && strcmp(nodes(j).name,'satelite')
                    LINKS(i,j) = 1; %satelite to satelite link
                else
                    LINKS(i,j) = -1; %satelite to station link, '-1' since stations are sinks in graph terms
                end
            else
                LINKS(i,j) = 0;
            end
            if (i == j) % || ((i >= n-M) && j <= i ) % if self or (if station and lower trianglular matrix because "stations sink data only one way")
               LINKS(i,j)=0; % do not link the node to itself
            end
        end
    end
    num_of_links = length(find(LINKS~=0))/2;
    disp('Done!')
    %Each row concerns -> a satelite, last 2 rows -> stations

    if PRINT_DETAILS
        disp('|=================================================================================================|')
        disp('Distance matrix:------------------------------------------------------------------')
        disp(DISTANCES)
        disp('Link matrix:----------------------------------------------------------------------')
        disp(LINKS)%should I include the diagonal elements (selfs)?
        disp('-- [1] = satelite to satelite connection, [-1] = satelite to station -------------')
        disp('-- Number of links: '+string(num_of_links))
        disp('|=================================================================================================|')
    end

    % Constructing fmincon parameters (objective function, contraints: Aeq, beq, A, b) ========================================================================================================================================================================================================================================================
    %%%% Constructing the function
    disp('Constructing optimization parameters...')

    % Creating Aeq (Kirchhoff's law matrix)
    % xij = zeros(num_of_links,3); % xij[value][parent node 1][parent node 2]
    counter = 1;
    xij={}; % cell array (like lists in the R programming language)
    for i = 1:n
        for j = i:n
           if LINKS(i,j) ~= 0 
              xij{counter} = 1:3;
              xij{counter}(1) = LINKS(i,j);
              xij{counter}(2) = i;
              xij{counter}(3) = j;
              counter = counter + 1;
           end
        end
    end
    if PRINT_DETAILS
        disp('|[value][parent node 1][parent node 2]------|')
        for i = 1:(counter-1)
        %     disp(string(xij{i}(1)) + "|" + string(xij{i}(2)) + "|" + string(xij{i}(3)) + "|")
           disp(xij{i})
        end
        disp('') % break line
    end

    num_station_links = 0; % number of links that have a station as a parent
    for i = 1:length(xij)
        if xij{i}(2) > NUMBER_OF_SATELLITES || xij{i}(3) > NUMBER_OF_SATELLITES
            num_station_links = num_station_links + 1;
        end
    end

    A_kirchhoff = zeros(0,num_of_links*2+n); % Each row is a node. Each column until column[num of satellites] expresses which links belong to the node. Multiplying A_kirchhoff with optimization vector creates all the kirchhoff equalities that must hold.
    for i = 1:n
        temp = xijvec(i,xij,n,num_station_links,NUMBER_OF_SATELLITES);%,NUMBER_OF_STATIONS); % DEBUG , I WAS PREVIOUSLY USING xijvec()
        if isempty(A_kirchhoff)
            A_kirchhoff = [A_kirchhoff , temp]; %column bind
        else
            A_kirchhoff = [A_kirchhoff ; temp]; %row bind
        end
    end
    %A_kirchhoff = [A_kirchhoff ; [zeros(1,num_of_links*2), ones(1,NUMBER_OF_SATELLITES + NUMBER_OF_STATIONS)]]; % Adding sum(divergencies)==0 (NOT NECESSARY, IT SHOULD HOLD NEVERTHELESS)

    beq = zeros(1,length(A_kirchhoff(:,1))); % OTHER KIRCHOFF EQUALITY PART
    A = [];% A = diag(ones(1,length(A_kirchhoff(1,:))));
    b = [];% b = [zeros(1,num_of_links*2-num_station_links)+LINK_CAPACITY,zeros(1,NUMBER_OF_SATELLITES)+inf, zeros(1,NUMBER_OF_STATIONS)]; %zerosNUMBER_OF_SATELLITES+100 % CAPACITIES (UPPER BOUNDS)
    upper_bounds = [zeros(1,num_of_links*2-num_station_links)+LINK_CAPACITY, zeros(1,NUMBER_OF_SATELLITES+NUMBER_OF_STATIONS)+inf];
    % lower_bounds = [zeros(1,(num_of_links*2-num_station_links)) zeros(1,NUMBER_OF_SATELLITES)+lbsi (1:NUMBER_OF_STATIONS)*(-inf)]; % zeros for rates (and sources(?)), -inf for sinks
    lower_bounds = zeros(1,length(A_kirchhoff(1,:)));
    disp('Done!')
    % Prints: ------------------------------------------------------------
    if PRINT_DETAILS
        disp('|- - -------\ (Aeq) K I R C H H O F F /-------------- - -|')
        disp(A_kirchhoff)
        disp('|- - -------\ (beq) EQUALITY VECTOR /---------------- - -|')
        disp(beq)
    %     disp('|- - -------\ (A) INEQUALITY MATRIX /---------------- - -|')
    %     disp(A)
        disp('|- - -------\ (b) CAPACITIES (UPPER BOUNDS)/--------- - -|')
    %     disp(b)
        disp(upper_bounds)
        disp('|- - -------\ (lower_bounds) LOWER BOUNDS /---------- - -|')
        disp(lower_bounds)
    end
    % --------------------------------------------------------------------
    % Using fmincon! ========================================================================================================================================================================================================================================================
    %%%%--------- fmincon(fun,x0,A,b,Aeq,beq) % x0 is the initial point used by the optimizer
    objective_function = @(xs)[(10^(-3))*ones(1,num_of_links*2 - num_station_links), zeros(1,NUMBER_OF_SATELLITES),-ones(1,NUMBER_OF_STATIONS)]*xs'; %xs is the optimization vector. xs = [x1,x2,...,x(links_num),x(links_num+1),...,x(2*links_num),s1,s2,...,sn]
    x0 = zeros(1,(num_of_links*2-num_station_links+n));
    % x0 = [ones(1,(num_of_links*2-num_station_links + n - 3)) 1 2 -100]; % debug
    opt_results = fmincon(objective_function, x0, A, b', A_kirchhoff, beq', lower_bounds, upper_bounds); % zeros(size(1:(num_of_links*2+n))) -> every rate should be positive
    % for i = 1:length(opt_results) % Print in a 'nicer' format. Only prints non zero values
    %     if opt_results(i)>1
    %         if i<=num_of_links*2
    %             sprintf("x%i= %e",[i,opt_results(i)])
    %         else
    %             sprintf("s%i= %e",[i-num_of_links*2,opt_results(i)])
    %         end    
    %     end
    % end
    disp('---------- OPTIMIZATION RESULTS ----------------------------------------------------------------')
    disp(opt_results)
    disp('------------------------------------------------------------------------------------------------')
    out = opt_results;
end% end of main function

%% Plot 2D world map and links
%{
 THIS SCRIPT IS FOR THE VISUALIZATION OF THE SATELITE ORBITS ON A 2D MAP
 OF THE EARTH. BEFORE CONTINUING IT SHOULD BE CONSIDERED THAT THERE IS NO 
 COMPLETE MATCHING OF THE 3D ORBITS TO THE 2D MAPPING BECAUSE STATIONS 
 SHOULD BE REGARDED AS STATIONARY (ROTATION OF THE EARTH WILL BE IGNORED) 
 TO KEEP THE VELOCITY CALCULATIONS OF EACH SATELITE SIMPLE.
%}

%% FUNCTIONS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% - - - - - - - - - - - - - - - - - - - - F U N C T I O N S - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

function out = create_nodes(num_satelites, num_stations, sat_inverse_vel, stat_inverse_vel, random_factor, theta_phi)
% Creates a list of nodes in the form of [satelite_objects .... station_objects]
%
% - num_satelites: number of satelites (set default as 9)
% - num_stations: number of stations (set default as 2)
% - sat_inverse_vel: vector, inverse velocity of satelites (smaller values imply high
% velocity) use values around [3 3 3 ... 3]
% - stat_inverse_vel: vector, inverse velocity of satelites (every station has the
% same value) (smaller values imply high velocity) use values around [80 80 80 ... 80]
% - random_factor: boolean, if true it adds some random noise from uniform distribution in [-10,10] on the satelite velocities
% - theta_phi: arg_theta and arg_phi values. In order to be a circle they
% should be equal

    earth_radius = 200;
    rounds = 10;
    matr = [ ];
    inverse_velocities = [sat_inverse_vel, stat_inverse_vel]; % the bigger the value the slower the velocity
    
    if random_factor
        rng(79) %setting seed
        random_numbers_sats = (rand(1,num_satelites)*8+1)-(4-1);% (rand(1,num_satelites)*number*2)-number -> rand() returns values in (0,1) so
                                                                % I do this in order to turn it (the output range) into (-number,number)
        rands = [random_numbers_sats, zeros(1,num_stations)];% adding random factor only at satelite velocities
        
        inverse_velocities = inverse_velocities + rands;
    end
    inverse_velocities = inverse_velocities * 1000;
       
    % Remember: satelite3D(arg_theta,arg_phi, arg_alt, arg_init_pos, arg_periods, arg_vel, arg_name)
    % Constructing satelites & stations:
    initializer_difs = (0:(num_satelites+num_stations-1))*27; %used to initialiaze satelites at different positions in order to not overlap each other
    for ii = 1:(num_satelites + num_stations)
        if ii <= num_satelites 
            matr = [matr, satelite3D(theta_phi(1), theta_phi(2), earth_radius+50, 20 + initializer_difs(ii), rounds, inverse_velocities(ii), 'satelite')];
        else
            matr = [matr, satelite3D(theta_phi(1), theta_phi(2), earth_radius, 200 + initializer_difs(ii)*4, rounds, inverse_velocities(ii), 'station')];
        end
    end
    
    disp('[~Report:] Satellites and stations successfully created')
    out = matr;
end


function out = euclidean_dist(vec1, vec2)
% Definition: euclidean_dist(vec1, vec2)
% Calculates euclidean distance of two vectors
    out = sqrt(sum((vec1 - vec2) .^ 2));
end


function out = print_list(lista)
    temp = '';
    for i = 1:length(lista)
        temp = [temp  ' '  convertStringsToChars(string(lista(i)))];
    end
    out = temp;
end


function out = xijvec(i,x,num_nodes,num_station_links,NUMBER_OF_SATELLITES) 
% This function is used for the creation of A_kirchoff matrix.
% However, each time it is called, one row of kirchoff matrix is created
% (regarding link i). Therefore the matrix is created in the main script.
% The output is a vector of A_kirchhoff matrix's row (fianl_Aeq).
% INPUTS: - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% i = current position of outter iterator (used in main script at assemblying)
% x = xij where xij is a cell array. xij{i} contains value, parent_node_1, parent_node2 of link i 
% num_nodes = number of nodes
% num_sats = number of satelites

    temp = zeros(1,length(x));
    for it =1:length(x) % for every link, check if any edge is node i
        parent1 = x{it}(2);
        parent2 = x{it}(3);
        if parent1 == i  % if node i is parent 1 of the link
           temp(it) = 1;
        else %elif:
            if parent2 == i % if node i is parent 2 of the link
                temp(it) = -1;
            end
        end
    end
    % Final step: out == [Aeq1 Aeq2 si]
    temp3 = diag(ones(1,num_nodes));
    if i > NUMBER_OF_SATELLITES
       temp3 = -temp3; 
    end
    temp2 = temp(1:(length(temp)-num_station_links)); % se periptwsh pou den valoume tis sthles
    out = [temp, -temp2, -temp3(i,:)];
end
