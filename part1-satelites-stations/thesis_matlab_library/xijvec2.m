function out = xijvec2(epoch, node, total_epochs, xij, nodes, NUMBER_OF_SATELLITES)
    % xijvec version of part 2. Creates an Aeq row under the hypothesis
    % that the graph remains the same in every epoch.
    epoch = epoch-1;
    n = length(nodes);
    num_links = length(xij); % number of unidirectional links
    num_station_links = 0; % number of links that have a station as a parent
    for it = 1:length(xij)
        if xij{it}(2) > NUMBER_OF_SATELLITES || xij{it}(3) > NUMBER_OF_SATELLITES
            num_station_links = num_station_links + 1;
        end
    end
    step = length(xij)*2-num_station_links+n*2; % integer length of an Aeq where there is only one epoch
    final = zeros(1,(step)*total_epochs);
    for it =1:length(xij) % for every link, check if any edge is node i
        parent1 = xij{it}(2);
        parent2 = xij{it}(3);
        if parent1 == node  % if current node is parent 1 of the link
           final(step*epoch + it) = 1;
        elseif parent2 == node % if current node is parent 2 of the link
           final(step*epoch + it) = -1;
        end
    end
    for it = 1:(num_links-num_station_links) % second half
       final(step*epoch +  num_links + it) = -final(step*epoch + it);
    end
    % Divergencies:
    if node > NUMBER_OF_SATELLITES % if node is a station (sink)
       final(step*epoch + num_links*2 - num_station_links + node) = 1; 
    else
       final(step*epoch + num_links*2 - num_station_links + node) = -1;
    end
    % Buffers:
    if epoch+1 == 1 
        final(step*epoch + num_links*2 - num_station_links + n + node ) = -1;
    elseif epoch+1 >= 2
        final(step*epoch + num_links*2 - num_station_links + n + node ) = -1;
        final(step*(epoch-1) + num_links*2 - num_station_links + n + node) = 1;
    end    
    out = final;
end