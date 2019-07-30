function results = heuristic_1(DISTANCE_MATRIX_CA, INF_TR_VECTOR, nodes, XIJ_CA, LINK_CAPACITY, COMMUNICATION_RANGE, NUMBER_OF_SATELLITES)
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
%  INF_TR_VECTOR = information transmission vector. 
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

    disp(' >> Inside heuristic_1 !!! -----------------------------------------------')
    disp('> Information Transaction Vector:')
    disp(INF_TR_VECTOR);

    NUMBER_OF_STATIONS = length(nodes) - NUMBER_OF_SATELLITES;
    total_epochs = length(DISTANCE_MATRIX_CA);
    temporary_results_ca = {}; % contains each epoch's LINK_RATES_VECTOR
    for epoch = 1:total_epochs % for every epoch
        disp(' > epoch '+ string(epoch) + ' - - - - - - - - - - - - - - - - - - - - ')
        current_distance_matrix = DISTANCE_MATRIX_CA{epoch};
        current_xij = XIJ_CA{epoch};        
        
        link_rates_vector = zeros(1,length(XIJ_CA{epoch})*2-NUMBER_OF_STATIONS); % starts from zero and increases in every epoch
        link_capacity_vector = LINK_CAPACITY * ones(1,length(link_rates_vector)); % can be SPENDED
        
        % Report messages: (to asure that everything goes as planned)
%         disp('current_distance_matrix:')
%         disp(current_distance_matrix)
%         disp('current_xij:')
%         for x = 1:length(current_xij)
%             disp(current_xij{x})
%         end
%         disp('link_rates_vector:')
%         disp(link_rates_vector)
%         disp('link_capacity_vector:')
%         disp(link_capacity_vector)
        
        sat_to_stat_matrix = null(3,NUMBER_OF_SATELLITES);% 2xNUMBER_OF_SATELLITES matrix containing the distance to the closest station of each satellite
        for node = 1:NUMBER_OF_SATELLITES % for every satellite
            % Get nearest neighbors to a station (of current node):
            to_station_dist = min(current_distance_matrix(node,(NUMBER_OF_SATELLITES+1):length(nodes))); % find minimum distance between stations and other_node 
            closest_station = find(current_distance_matrix(node,:) == to_station_dist);% which station is closer to other_node
            sat_to_stat_matrix(1,node) = node;
            sat_to_stat_matrix(2,node) = closest_station;  
            sat_to_stat_matrix(3,node) = to_station_dist;
            % ^Interpretation: (Column 1: satellite 1 is closer to node [row1,column1] (i.e. station x) at a distance [row2,column 1])
        end% end of <for every node>
        disp('sat_to_stat_matrix:');
        disp(sat_to_stat_matrix);
        
        % "Neighbor" closest to station *increasing order*: (Could be a self, or the closest could not be a neighbor)
        ordered_sat_to_stat =  sortrows(sat_to_stat_matrix',3)';% a little shady, MIGHT WANT TO INVESTIGATE WHILE DEBUGGING
        disp('ordered_sat_to_stat :');
        disp(ordered_sat_to_stat );
        pause;%debug

        % If minimum is a neighbor or self (a connection can be
        % established) then use the link connecting current node and
        % neighbor:
        for i = 1:length(ordered_sat_to_stat)
            if ordered_sat_to_stat(i) <= COMMUNICATION_RANGE% if it is a neighbor
%                 link_rates_vector - U N D E R  C O N S T R U C T I O N
%                 link_capacity_vector - U N D E R  C O N S T R U C T I O N
            end
        end


        temporary_results_ca{epoch} = link_rates_vector;  
    end % end of epochs
    
    final = null(1,1);
    for i = 1:length(temporary_results_ca)
        final = [final, temporary_results_ca{i}];
    end
    results = link_rates_vector ;
end