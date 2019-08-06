function out = create_upper_bounds(NUMBER_OF_SATELLITES, list_num_links, list_num_station_links, nodes, total_epochs, LINK_CAPACITY)
% Creates upper_bound vector for problem of "part_2"
    bounds = null(1);
    num_of_nodes = length(nodes);
    NUMBER_OF_STATIONS = num_of_nodes - NUMBER_OF_SATELLITES;
    flow_capacity = LINK_CAPACITY;
%     buffer_capacity = 10; %old version
    buffer_capacity = inf; % heuristic does not have any buffer capacity
    for epoch = 1:total_epochs      
        num_flows = list_num_links(epoch)*2 - list_num_station_links(epoch);
        bounds = [ bounds, zeros(1, num_flows) + flow_capacity];% flows
        bounds = [ bounds, zeros(1, NUMBER_OF_SATELLITES) + 10]; % divergencies, (old version: inf), (heur1: 10)
        bounds = [ bounds, zeros(1, NUMBER_OF_STATIONS) + inf]; % divergencies
        
        bounds = [ bounds, zeros(1, num_of_nodes) + buffer_capacity]; % buffers 
    end
    out = bounds;   
end