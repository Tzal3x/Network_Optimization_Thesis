function results = heuristic_1(DISTANCE_MATRIX_CA, nodes, XIJ_CA, LINK_CAPACITY, NUMBER_OF_SATELLITES, BUFFER_BOUND, GREEDY_CRITERION)
% Definition: heuristic_1(DISTANCE_MATRIX_CA, INF_TR_VECTOR, NODES, XIJ_CA, CAPACITY_BOUNDS, COMMUNICATION_RANGE)
%  Heuristic (greedy-like) algorithm that aims to solve the satellite-station DTNUM
%  problem. i.e. the same we were trying to solve at 'solve_part2.m'
%  The Algorithm: --------------------------------------------------------
%  The algorithm contains epochs and epochs contain rounds.
%  An epoch is a single time unit that is followed by the corresponding
%  topology (a graph).
%  A round contains the below steps:
%  1) For every node, find the neighbor that is nearest
%  to a station.
%  2) Send it as much information as possible, then do the same to the second
%  nearest neighbor to a station until there can be no further information
%  transferred. 
%  3) Repeat until no node can transfer any data.
%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%  DISTANCE_MATRIX_CA = cell array of distance matrixes. each element is the
%  distance matrix of a given epoch.
% 
%  <not included!> INF_TR_VECTOR = (Information Transmission Vector) information that is generated by each node and is needed to be transmitted. 
%
%  nodes = vector of <node> objects.
%
%  XIJ_CA = cell arrayof cell arrays, each element contains a link's info about it's two
%  ends (reffered to as parents)
%
%  LINK_CAPACITY = vector containing the capacity bounds of links
%
%  COMMUNICATION_RANGE = the maximum distance between nodes that a
%  connection can be established
%
%  BUFFER_BOUND = the upper limit of buffer values
%
%  GREEDY_CRITERION = 'closest'/'fastest'. 'closest' is the satellite that
%  is closest to a station and 'fastest corresponds to the satellite that
%  moves faster to the direction of a station

%     disp('=======================================================================================================================================')
%     disp(' >>> Inside heuristic_1_'+GREEDY_CRITERION+'!!! ==========================================================================================================')    
%     disp('=======================================================================================================================================')
    NUMBER_OF_STATIONS = length(nodes) - NUMBER_OF_SATELLITES; 
    total_epochs = length(DISTANCE_MATRIX_CA); % get total epochs
    GENERATED_INFO = 10;
    % Results (each epochs results are an element of each cell array):
    flows_results_ca = {}; % contains each epoch's LINK_RATES_VECTOR
    buffer_results_ca = {};
    tic
    for epoch = 1:total_epochs % for every epoch
%         disp('-----------------------------------------------------------------------------------------------------------')
%         disp(' >> epoch '+ string(epoch) + '-----------------------------------------------------------------------------')
%         disp('-----------------------------------------------------------------------------------------------------------')
        current_distance_matrix = DISTANCE_MATRIX_CA{epoch}; % get distances between nodes
        
        current_xij = XIJ_CA{epoch}; % get current xij cell array
         %{
         Creating new xij vector. {1}[ 1 | 1 | 2 ]    [12;  
                                  {2}[ 1 | 2 | 4 ] to  24]
         Cell array to vector were 12 (twelve) means parent1=1 parent2=2.
         %}
        new_xij = null(1,1);
        unwanted = null(1,1); % deleting station to satellite links
        total_links = length(current_xij);
        for i = 1:(2*total_links)
            if i <= total_links
                new_xij = [new_xij, str2double(string(current_xij{i}(2))+string(current_xij{i}(3)))];
            elseif i > total_links
                ii = i - total_links;
                if current_xij{ii}(3) > NUMBER_OF_SATELLITES
                    unwanted = [unwanted, i];
                end
                new_xij = [new_xij, str2double(string(current_xij{ii}(3))+string(current_xij{ii}(2)))];
            end
        end
        new_xij(unwanted) = [];
        
        % These vectors play main roles:
        num_station_links = get_nstat_links(XIJ_CA{epoch}, NUMBER_OF_SATELLITES); % get number of satellite-to-station links
        INF_TR_VECTOR = [GENERATED_INFO * ones(1,NUMBER_OF_SATELLITES), zeros(1,NUMBER_OF_STATIONS)]; % creating ITV 
        if epoch > 1
            previous_buffer = buffer_results_ca{epoch-1};
            for iter = 1:length(previous_buffer) 
               exceeds = (GENERATED_INFO + previous_buffer(iter)) > BUFFER_BOUND; % boolean
                if exceeds % IF previous buffers plus generated info exceed buffer bounds...
                    INF_TR_VECTOR(iter) = GENERATED_INFO + previous_buffer(iter) - BUFFER_BOUND;
                else
                    INF_TR_VECTOR(iter) = GENERATED_INFO + previous_buffer(iter); % add remaining info of previous epoch to current
                end
            end
        end
        INF_TR_VECTOR_START = INF_TR_VECTOR; % INF_TR_VECTOR at the start of the epoch
        link_rates_vector = zeros(1,length(XIJ_CA{epoch})*2-num_station_links); % similar to optimal flows
        link_capacity_vector = LINK_CAPACITY * ones(1,length(link_rates_vector)); % available capacity of current round and epoch
        
        % Me ti tha isoutai to divergence sth telikh? % VERSION 2
%         divergence_results = zeros(1,length(nodes)); % osi pliroforia perase mesa to kathe node
%         divergence_results = [10*ones(1,NUMBER_OF_SATELLITES), zeros(1,NUMBER_OF_STATIONS)];
%         divergence_results = INF_TR_VECTOR_START; % VERSION 1 
        
        % Report messages: (to asure that everything goes as planned)
%         disp('> Information Transmission Vector(INF_TR_VECTOR):')
%         disp(INF_TR_VECTOR);
%         disp('link_rates_vector:')
%         disp(link_rates_vector)
%         disp('link_capacity_vector:')
%         disp(link_capacity_vector)
        
        sat_to_stat_matrix = null(3,NUMBER_OF_SATELLITES);% 3xNUMBER_OF_SATELLITES (dimensions) matrix containing the distance to the closest station of each satellite
        
        for node = 1:NUMBER_OF_SATELLITES
            if strcmp(GREEDY_CRITERION, 'closest')
                to_station_dist = min(current_distance_matrix(node,(NUMBER_OF_SATELLITES+1):length(nodes))); % find minimum distance between stations and other_node
                closest_station = find(current_distance_matrix(node,:) == to_station_dist); % CRITERION 1, FIND CLOSEST STATION
                sat_to_stat_matrix(1,node) = node;
                sat_to_stat_matrix(2,node) = closest_station;  
                sat_to_stat_matrix(3,node) = to_station_dist;
            elseif strcmp(GREEDY_CRITERION, 'fastest')
                if epoch > 1
                    speed_matrix = DISTANCE_MATRIX_CA{epoch-1} - DISTANCE_MATRIX_CA{epoch};
%                     disp('PREVIOUS DISTANCE MATRIX:');disp(DISTANCE_MATRIX_CA{epoch-1}); % debug
                else
                    speed_matrix = 0 - DISTANCE_MATRIX_CA{epoch}; % an na vrethei to -max = min twn apostasewn. alliws tha epelege th megalyterh apostash
                end
                    satel_to_station_speeds = speed_matrix(node,(NUMBER_OF_SATELLITES+1):length(nodes));
                    max_speed_to_station = max(satel_to_station_speeds); % me poia taxythta kateuthinetai grigorotera pros tous stathmous
                    which_station = find(speed_matrix(node,:) == max_speed_to_station); % se POION STATHMO katefthinetai pio grigora
                    sat_to_stat_matrix(1,node) = node;
                    sat_to_stat_matrix(2,node) = which_station;  
                    sat_to_stat_matrix(3,node) = max_speed_to_station;
%                     disp('CURRENT DISTANCE MATRIX:');disp(DISTANCE_MATRIX_CA{epoch});% debug
%                     disp('SPEED MATRIX:');disp(speed_matrix); % debug
%                     disp('sat_to_stat_matrix:');disp(sat_to_stat_matrix)
            end
%             pause % debug
        end   


        % "Neighbor" closest to station *increasing order*: (Could be a self, or the closest could not be a neighbor)
        if strcmp(GREEDY_CRITERION, 'closest')
            ordered_sat_to_stat =  sortrows(sat_to_stat_matrix',3)';% a little shady, MIGHT WANT TO INVESTIGATE WHILE DEBUGGING
        elseif strcmp(GREEDY_CRITERION, 'fastest')
            ordered_sat_to_stat =  -sortrows(-sat_to_stat_matrix',3)'; % decreasing order! The fastest should be first in row!
        end
%         disp('current_distance_matrix:');
%         disp(current_distance_matrix);
%         disp('ordered_sat_to_stat :');
%         disp(ordered_sat_to_stat );

        %{
         If minimum is a neighbor or self (a connection can be
         established) then use the link connecting current node and
         neighbor:
        %}
        next_round_bool = true;
        round = 0;
        forbidden_links = null(1,1); % here are registered links that have already been used
        while(next_round_bool) % Rounds of algorithm start here
            round = round + 1;
            next_round_bool = false; % if no more info can be transferred it will remain <false> and algorithm will terminate
%             disp('> round '+string(round)+': \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/')
            for sat = 1:NUMBER_OF_SATELLITES
                % Checking first if self is connected with a station:
                temp_osts = ordered_sat_to_stat(1:2,:); % temp_osts: temporary ordered_sat_to_stat
                temp2 = ordered_sat_to_stat(1,:) == sat; % delete that column to avoid duplicates
                temp_osts(:,temp2) = [];
                
                for k = (NUMBER_OF_SATELLITES+NUMBER_OF_STATIONS):-1:(NUMBER_OF_SATELLITES+1) % adding sat to every station pseudo-links
                    temp_osts = [[sat; k] ,temp_osts];
                end
                
                for i = 1:length(temp_osts(1,:)) % checking for neighbors
                    sat2 = temp_osts(1,i); % candidate                    
                    if sat == sat2 % this means that sat is close to a station and therefore he transmits the info directly
                        sat2 = temp_osts(2,i); % sat2 becomes the target station, else sat would transmit to itself
                    end
                    str2num_link = str2double(string(sat)+string(sat2)); % using this to spot which link will be used
                    which_link = find(new_xij == str2num_link);
                    if ~isempty(which_link)  % if they are neighbors there exists at least one conneting them
                        which_link = which_link(1); %to avoid duplicates, get unique
%                         disp('neighbors: yes')
                        if (INF_TR_VECTOR(sat) > 0) & (link_capacity_vector(which_link) > 0) & isempty(find(forbidden_links == str2num_link)) %#ok<AND2,EFIND> % IF INFO CAN BE TRANSMITTED
                                                           
                            % Registering used link:
                            link2 = str2double(string(sat2) + string(sat)); % opposite direction (because if one link is used, the other should be silent
                            forbidden_links = [forbidden_links, link2];
                                      
                            cap_minus_info = link_capacity_vector(which_link) - INF_TR_VECTOR(sat); % the absolute value of this is the information transmitted
                            
                            % Check buffer bound:
                            if epoch > 1 
                                if buffer_results_ca{epoch-1}(sat2) + abs(cap_minus_info) > BUFFER_BOUND % If previous buffer plus the incoming information exceeds buffer bound then cancel transmission
                                    continue
                                end
                            end
                            
                            previous_capacity = link_capacity_vector(which_link); % used to calculate link_rates_vector by finding the difference of current minus previous
                            if cap_minus_info >= 0 % capacity > info
                                INF_TR_VECTOR(sat2) = INF_TR_VECTOR(sat2) + INF_TR_VECTOR(sat) ;% info incoming to sat2, so sat2 needs to transfer more info
                                INF_TR_VECTOR(sat) =  0; % info gone
                                link_capacity_vector(which_link) = cap_minus_info;
                            else % info > capacity
                                INF_TR_VECTOR(sat2) = INF_TR_VECTOR(sat2) + previous_capacity;
                                INF_TR_VECTOR(sat) =  abs(cap_minus_info);
                                link_capacity_vector(which_link) = 0;
                            end
                            transmitted_info = previous_capacity - link_capacity_vector(which_link);
                            link_rates_vector(which_link) = link_rates_vector(which_link) + transmitted_info;
                            next_round_bool = true;
                            
%                             divergence_results(sat) = divergence_results(sat) - transmitted_info; % apostoleas, xanei pliroforia
%                             divergence_results(sat2) = divergence_results(sat2) + transmitted_info; % dektis, pairnei pliroforia

                              % Quality control report: ----------------
%                             disp(' ');
%                             disp('> Action:');
%                             disp('flow: '+string(sat)+' -> '+string(sat2)+': '+string(transmitted_info));
%                             disp('[~Report]: Main role vectors change to:');
%                             disp('INF_TR_VECTOR:'); disp(1:length(nodes)); disp(INF_TR_VECTOR);
%                             disp('link_capacity_vector:'); disp(new_xij); disp(link_capacity_vector);
%                             disp('link_rates_vector:'); disp(link_rates_vector);
                       elseif INF_TR_VECTOR(sat) == 0
%                              disp('Outs of info');
                             break; % if there is no remaining info o obe transmitted for now, then move to next satellite
                       elseif link_capacity_vector(which_link) > 0
%                              disp('Out of capacity!')
                       end % end if info can be transferred
                   else
%                       disp('neighbors: no')
                   end % end if neighbors
               end % end for closest to stations nodes
           end % end for every satellite
        end % end of rounds 
%         disp('> End of epoch results:')
%         disp('INF_TR_VECTOR:'); disp(INF_TR_VECTOR);
%         disp('link_capacity_vector:'); disp(link_capacity_vector);
%         disp('link_rates_vector:'); disp(link_rates_vector);
        
        buffer_results_ca{epoch} = [INF_TR_VECTOR(1:(length(INF_TR_VECTOR)-NUMBER_OF_STATIONS)), zeros(1,NUMBER_OF_STATIONS)]; % PETAW OTI INFO PERISSEPSE TWN STATIONS
        previous_buffer = 0;
        if epoch > 1
            previous_buffer = buffer_results_ca{epoch-1};
        end
        divergence_results =  INF_TR_VECTOR_START - INF_TR_VECTOR + buffer_results_ca{epoch} - previous_buffer; % DIVERGENCE VERSION 3, apla elysa ws pros to divergence doulevwntas panw sth gnwsth isothta
        flows_results_ca{epoch} = [link_rates_vector, divergence_results, buffer_results_ca{epoch}]; 
    end % end of epochs
%     disp('Paused after algorithm terminated...');
%     pause;
    
    final = null(1,1);
    for i = 1:length(flows_results_ca)
        final = [final, flows_results_ca{i}];
    end
%     disp('[~Report:] End of heuristic_1.')
%     disp('Final result:');
%     disp(final); 
%     disp('________________________________________________________________________________________________________________________');
    
    results_ca{1} = final; % final result
    results_ca{2} = toc; % time performance
    results = results_ca ;
end