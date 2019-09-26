function out = get_list_num_links(total_epochs, nodes)
% Syntax: get_list_num_links(epochs,nodes)
% Saves in a list the number of links existing in the network in every
% epoch.
% total_epochs: integer, total number of epochs
% nodes: list of satellite/station objects 
    n = length(nodes);
    stop = total_epochs; % epochs, keep it low like 30
    times = 10000; % upper bound of iterations
    list_num_links = [];
    for epoch = 1:times
       coords = []; % coordinates of each
       for j = 1:length(nodes)
          coords = [coords; nodes(j).lifetime_coordinates(:,epoch)' ];
       end
       
       LINKS = create_LINKS(coords, nodes, n);
       xij = create_flow_info(LINKS,n); 
       list_num_links = [list_num_links, length(xij)];
       if epoch == stop
           break
       end
    end % end for
    out =  list_num_links ;
end