%{
!> PART 2: Formulating a problem like the "problem 1" of the paper:
" Delay-Tolerant Network Utility Maximization, Toumpis et al."
*** [Self-note:] A 'cd C:/...' may be needed for the function to run smoothly *** 
%}

function out = solve_part2(NUMBER_OF_SATELLITES, NUMBER_OF_STATIONS, RANDOM_VELOCITIES, INVERSE_VELOCITIES_SATEL,...
                    INVERSE_VELOCITIES_STATIONS, STOP_AT_TIME, THETA_PHI, LINK_CAPACITY, PRINT_DETAILS,...
                    PRINT_MAIN_PARAMETERS, SHOW_TOPOLOGY, COMMUNICATION_RANGE, SOLVER, INIT_POS)
    % Solves the second version of the problem Delay-Tolerant Network
    % Utility Maximization Problem [Problem I: DTNUM]
    % Having added data buffers
    %
    % Input parameters:
    % NUMBER_OF_SATELLITES = integer, number of satellites 
    % NUMBER_OF_STATIONS = integer, number of stations 
    % RANDOM_VELOCITIES = boolean, true if a random factor should be added at the satellite velocities
    % INVERSE_VELOCITIES_SATEL = vector, each element is the inverse velocity of a satellite
    % INVERSE_VELOCITIES_STATIONS = vector, each element is the inverse
    % velocity of a station. Shouldl be the same for all, since it coincide
    % with earths angular speed around itself
    % STOP_AT_TIME = integer, time when the optimization should take place
    % (number of epochs)
    % THETA_PHI = vector, consisting of 2 elements [azimuth elevation]
    % LINK_CAPACITY = float, maximum link capacity
    % PRINT_DETAILS = boolean
    % COMMUNICATION_RANGE = maximum distance between nodes that a
    % connection can be established.
    % SOLVER = type of method used to solve DTNUM.'fmincon'/'linprog'/'heuristic_1'
    % Output:
    % A vector consisting of optimal values
    % ---------------------------------------------------------------------
    
    % Setting up figure display options: - - - - - - - - - - - - - - - - -
    earth_radius = 200;
    if SHOW_TOPOLOGY
    figure('Name','3D Simulation of satelite orbits');
    title("satelites:"+string(NUMBER_OF_SATELLITES)+...
          ", stations:"+string(NUMBER_OF_STATIONS)+...
          ", random velocities:"+string(RANDOM_VELOCITIES)+...
          ", total epochs:"+string(STOP_AT_TIME));
        hold on; % keep plotting on the existing figure
        axis equal
        view(30,0)% setting azemuth and elevation angles of camera for nicer visualization 
    axis([-earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200]); %setting up figure size
    end       

    disp("|========================================= INSIDE PART 2 =============================================|")
    nodes = create_nodes(NUMBER_OF_SATELLITES, NUMBER_OF_STATIONS, INVERSE_VELOCITIES_SATEL, ... 
                         INVERSE_VELOCITIES_STATIONS, RANDOM_VELOCITIES, THETA_PHI, INIT_POS); 
%     for node = 1:length(nodes)
%         temp = nodes(node).lifetime_coordinates; 
% %         disp(temp)
% %         pause
%         for j = 1:(length(temp(1,:))-1)
%             x = [temp(1,j), temp(1,j+1)];
%             y = [temp(2,j), temp(2,j+1)];
%             z = [temp(3,j), temp(3,j+1)];
%             line(x,y,z,'LineStyle','--')
%         end
%     end
    
    if PRINT_MAIN_PARAMETERS
       disp("NUMBER_OF_SATELLITES: " + string(NUMBER_OF_SATELLITES))
       disp("NUMBER_OF_STATIONS: " + string(NUMBER_OF_STATIONS))
       disp("RANDOM_VELOCITIES: " + string(RANDOM_VELOCITIES))
       disp("INVERSE_VELOCITIES_SATEL: " + print_list(INVERSE_VELOCITIES_SATEL))
       disp("INVERSE_VELOCITIES_STATIONS: " + print_list(INVERSE_VELOCITIES_STATIONS))
       disp("STOP_AT_TIME: " + string(STOP_AT_TIME))
%        disp("THETA_PHI: " + print_list(THETA_PHI))
       disp("LINK_CAPACITY: " + string(LINK_CAPACITY))
       disp("PRINT_DETAILS: " + string(PRINT_DETAILS))
    end
    
    n = NUMBER_OF_SATELLITES + NUMBER_OF_STATIONS; % number of total nodes (#sats + #stats)
    
    %%%% WARNING! 'times' must be always < lesser from linspace length of all satelites and stations (avoiding index out of bounds error)
    stop = STOP_AT_TIME; % epochs, keep it low like 30
    times = 10000; % upper bound of iterations
    [x,y,z] = sphere();%get coordinates in order to create 3D spheres. A sphere can be the earth, a satelite or a station.       
%     epoch_link_list = get_list_num_links(STOP_AT_TIME, nodes); % number of links per epoch 
    
%===================================================================================================================================================
%===================================================================================================================================================
%===================================================================================================================================================
%===================================================================================================================================================
%===================================================================================================================================================
    Aeqs_cell_array = {} ; % each element is the Aeq corresponding to an epoch.
    list_num_station_links = null(1,1); % number of station undirectional links per epoch (used later at upper bound calculation)
    xij_ca = {}; % xij cell array containing the xij vector of every epoch
    coords_ca = {}; % coords cell array containing the coords matrix of every epoch
    distances_ca = {}; % euclidean distances between entities cell array
    list_num_links = null(1,1); % number of undirectional links per epoch
    for epoch = 1:times
       disp(''); disp(''); % break 2 lines
       disp('> EPOCH ' + string(epoch) + ' ===============================================================================================|')
       %Coordinates ---------------------------
       coords = [];
       for j = 1:length(nodes) % [longitude][latitude][altitude]
          coords = [coords; ...
                    nodes(j).lifetime_coordinates(:,epoch)' ];
       end
       coords_ca{epoch} = coords;
       
       distances_ca{epoch} = create_DISTANCES(coords,n);
       LINKS = create_LINKS(coords, nodes, n, COMMUNICATION_RANGE);

       xij = create_flow_info(LINKS, n); % get xij (parent1 parent2 vector)
       xij_ca{epoch} = xij;
       LINKS = abs(LINKS);
       num_links = length(xij); % number of undirectional links
       list_num_links = [list_num_links, num_links];
       
       num_station_links = 0; % number of links that have a station as a parent
       for it = 1:length(xij)
           if xij{it}(2) > NUMBER_OF_SATELLITES || xij{it}(3) > NUMBER_OF_SATELLITES
               num_station_links = num_station_links + 1;
           end
       end
       list_num_station_links = [list_num_station_links, num_station_links];
       
       if PRINT_DETAILS
           disp('| Distance matrix '+string(epoch)+':____________________________________________________________________________________________|')
           disp(distances_ca{epoch}); % debug
           disp('|(abs) Link matrix '+string(epoch)+':____________________________________________________________________________________________|')
           disp(abs(LINKS))
           disp('| xij '+string(epoch)+':____________________________________________________________________________________________________|')
           disp('|[value][parent node 1][parent node 2]------|')
           for i = 1:length(xij)
               disp(xij{i})
           end
           
           disp('') % break line
       end
       
%        num_station_links = get_nstat_links(xij, NUMBER_OF_SATELLITES);   % number of satellite-station undirectional links
       temp_aeq_rows = {}; % contains rows of "unfinished" Aeq of epoch i
       for node = 1:n
           temp_aeq_rows{node} = xijvec3(node, xij, n, NUMBER_OF_SATELLITES);
       end
       Aeqs_cell_array{epoch} = assemble_Aeq(temp_aeq_rows);
       if PRINT_DETAILS
          disp('| Aeq of epoch '+string(epoch)+' (missing previous epoch buffers):____________________________________________________|')
          disp(Aeqs_cell_array{epoch})
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
           if epoch == stop
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
    disp('[~Report:] Final epoch reached!')
    Aeq = asMultAeqs(Aeqs_cell_array);
    if PRINT_DETAILS
         disp('| FINAL Aeq: ____________________________________________________________________________________________|')
         disp(Aeq)
    end
%===================================================================================================================================================
%===================================================================================================================================================
%===================================================================================================================================================
%===================================================================================================================================================
%===================================================================================================================================================
%===================================================================================================================================================
    
    % Constructing fmincon parameters (objective function, contraints: Aeq, beq, A, b) 
    %%%% Constructing the function
    beq = zeros(1,length(Aeq(:,1))); % OTHER KIRCHOFF EQUALITY PART
    A = [];
    b = [];
    
    % At upper_bounds vector, flows & buffers should have a bound. Divergencies not
    upper_bounds = create_upper_bounds(NUMBER_OF_SATELLITES, list_num_links, list_num_station_links, nodes, STOP_AT_TIME, LINK_CAPACITY);
    lower_bounds = zeros(1,length(Aeq(1,:))); % remains the same as in PART_1
    
    disp('[~Report:] Done!')
    if PRINT_DETAILS
        disp('|- - -------\ (beq) EQUALITY VECTOR /---------------- - -|')
        disp(beq)
        disp('|- - -------\ UPPER BOUNDS /--------- - -|')
        disp(upper_bounds)
        disp('|- - -------\ LOWER BOUNDS /---------- - -|')
        disp(lower_bounds)
    end
    
    % Using fmincon! ============================================================================================================================================================================================================================
    %%%%--------- fmincon(fun,x0,A,b,Aeq,beq) % x0 is the initial point used by the optimizer
    coefs = objective_function_coefs_2(list_num_links, list_num_station_links, n, STOP_AT_TIME, NUMBER_OF_SATELLITES);
    objective_function = @(xs)coefs*xs'; %xs is the optimization vector. xs = [x1,x2,...,x(links_num),x(links_num+1),...,x(2*links_num),s1,s2,...,sn]
%     objective_function2 = @(xs)coefs*log(xs+0.05)'; 
    
    if PRINT_DETAILS
       disp('|- - -------\ Objective Function Coefficients /--------- - -|')
       disp(coefs)
    end
    x0 = zeros(1,length(Aeq(1,:))); %initial point
    if PRINT_DETAILS
       disp('|- - -------\ Initial Point x0: /--------- - -|')
       disp(x0)
    end
    
    % Optimal results (using convex optimizasion)
    if strcmp(SOLVER,'fmincon')
        opt_results = fmincon(objective_function, x0, A, b, Aeq, beq', lower_bounds', upper_bounds'); % zeros(size(1:(num_of_links*2+n))) -> every rate should be positive
    elseif strcmp(SOLVER,'linprog')
    % Optimal results (using linear programming)
        opt_results = linprog(coefs, A, b, Aeq, beq', lower_bounds', upper_bounds'); % zeros(size(1:(num_of_links*2+n))) -> every rate should be positive
        opt_results = opt_results';
    elseif strcmp(SOLVER,'heuristic_1')
    % Heuristic results:
        opt_results = heuristic_1(distances_ca, nodes, xij_ca, LINK_CAPACITY, NUMBER_OF_SATELLITES, 30);
    end

    disp('---------- OPTIMIZATION RESULTS ----------------------------------------------------------------')
    disp(opt_results)
    disp('------------------------------------------------------------------------------------------------')
    out = opt_results;

    
    %%% Create Simple Graph (figure): _____________________________________
%     total_epochs = length(Aeqs_cell_array);
%     ncols = null(1,1); % number of columns of Aeq at epoch i
%     for i = 1:total_epochs
%        ncols = [ncols, length(Aeqs_cell_array{i}(1,:))];
%     end
%     for i = 1:total_epochs
%        if i ~= 1
%            temp_opt_results = opt_results((ncols(i-1) + 1):(ncols(i-1) + ncols(i) - n));
%        elseif i == 1
%            temp_opt_results = opt_results(1:(ncols(i) - n));
%        end
%        createGraph(NUMBER_OF_SATELLITES, xij_ca{i}, temp_opt_results, nodes, i)
%     end
%     

    % Seperating optimization results by epoch:
    opt_flows_ca = {}; % vectors, optimization results per epoch without buffers
    opt_buffers_ca = {}; % optimal buffers per epoch
    opt_divergencies_ca = {}; % optimal divergencies per epoch
    total_epochs = length(Aeqs_cell_array); 
    ncols = null(1,1); % number of columns of Aeq at epoch i
    for i = 1:total_epochs
       ncols = [ncols, length(Aeqs_cell_array{i}(1,:))]; % +n because previous buffers of Aeq of every epoch are missing
    end

    for i = 1:total_epochs
       if i ~= 1
           ncols_of_ALL_previous = sum(ncols(1:(i-1)));
           temp_opt_flows = opt_results((ncols_of_ALL_previous + 1):(ncols_of_ALL_previous + ncols(i) - 2*n));
           temp_divergencies = opt_results((ncols_of_ALL_previous + ncols(i) - 2*n + 1):(ncols_of_ALL_previous + ncols(i)-n));
           temp_buffers = opt_results((ncols_of_ALL_previous + ncols(i) - n + 1):(ncols_of_ALL_previous + ncols(i)));
       elseif i == 1
           temp_opt_flows = opt_results(1:(ncols(i) - 2*n));
           temp_divergencies = opt_results((ncols(i) - 2*n + 1):(ncols(i)-n));
           temp_buffers = opt_results((ncols(i) - n + 1):ncols(i));
       end
       opt_flows_ca{i} = temp_opt_flows;
       opt_divergencies_ca{i} = temp_divergencies;
       opt_buffers_ca{i} = temp_buffers;
    end
    
    
    %{
    Total utility gained by each solver is the sum of station divergencies.
    Be careful to change this if you chose a different utility maximization
    function (objective function) such as log(d1)+log(d2)+...+log(dn).
    %}
    total_utility = 0;
    for i = 1:length(opt_divergencies_ca)
        temp = length(opt_divergencies_ca{i});
        total_utility = total_utility + sum(opt_divergencies_ca{i}((NUMBER_OF_SATELLITES+1):temp));
    end
    disp('+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +')
    disp('+ + + '+ SOLVER + ' total utility (sum of station divs)(abs): ' + string(abs(total_utility)))
    disp('+ + +  Average total information received by stations per epoch: '+ string(calculate_delay_v1(opt_divergencies_ca, NUMBER_OF_SATELLITES)))
    disp('+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +')
    
    GraphMap(NUMBER_OF_STATIONS, NUMBER_OF_SATELLITES, coords_ca, opt_flows_ca, opt_buffers_ca, opt_divergencies_ca, xij_ca, nodes,'light',SOLVER); % using fmincon/linprog results
    
end% end of main function