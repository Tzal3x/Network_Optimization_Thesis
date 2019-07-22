function GraphMap(num_stations, num_satellites, coords_ca, opt_buffers_ca, opt_divergencies_ca, xij_ca, nodes)
    % Create a Graph projected on the World Map. Similar to GraphDists
 
    n = length(nodes);
    figure;
    worldmap('World')
    load coastlines %#ok<LOAD>
    plotm(coastlat,coastlon)

    total_epochs = length(coords_ca);
%     new_coords = {}; % new geographic/polar coordinates
    for epoch = 1:total_epochs
        xij =  xij_ca{epoch};
%         coords = coords_ca{epoch};
%         opt_results = opt_results_ca{epoch}(1:(length(opt_results_ca{epoch})-n)) ; % (leaving this here for later)
        opt_buffers = opt_buffers_ca{epoch}; % buffers of nodes of given epoch
        opt_divergencies = opt_divergencies_ca{epoch}; % divergencies of nodes of given epoch
        
        % Creating the edge flows (edge i starts from node X and ends to Y -parents-)
        outgoing1 = null(1,1); outgoing2 = null(1,1);
        incoming1 = null(1,1); incoming2 = null(1,1);
        unwanted = null(1,1);
        for i = 1:length(xij)
            outgoing1 = [outgoing1, xij{i}(2)];
            incoming1 = [incoming1, xij{i}(3)];
            outgoing2 = [outgoing2, xij{i}(3)];
            incoming2 = [incoming2, xij{i}(2)];
            if xij{i}(1) == -1
               unwanted = [unwanted, i];  
            end
        end
        
        % Remove station to satellite edges:
        outgoing2(unwanted) = [];
        incoming2(unwanted) = [];
        outgoing = [outgoing1, outgoing2];
        incoming = [incoming1, incoming2];

        % Transform coordinates from cartesian to polar coordinates
        temp = null(1,1);
        for j = 1:length(coords_ca{epoch}) % for every node
            x = coords_ca{epoch}(j,1);
            y = coords_ca{epoch}(j,2);
            z = coords_ca{epoch}(j,3);
            long = pi/2 - acos(z/sqrt(x^2+y^2+z^2)); % polar - latitude
            lat = atan2(y,x); % polar - longitude
            temp = [temp ; [lat, long]];
        end
        new_coords = rad2deg(temp);
        
        % coordinates: radians to(->)degrees
        animation_ca = {}; % saving figure objects to this cell array, 
        lat = new_coords(:,1);
        long = new_coords(:,2);
        len = length(long);

        % Draw nodes
        animation_ca{1} = scatterm(long(1:(len-num_stations)),lat(1:(len-num_stations)),'red','filled'); % satellite nodes
        animation_ca{2} = scatterm(long((len-num_stations+1):len),lat((len-num_stations+1):len),'green','filled'); % station nodes
        
        % Adding edges (connecting lines):
        lines_ca = {}; % saving line objects to be able to delete them for animation purposes
        edge_label_pos = {};
        for i = 1:length(outgoing)
            node1 = outgoing(i); 
            node2 = incoming(i);
            x1 = new_coords(node1,1); y1 = new_coords(node1,2);
            x2 = new_coords(node2,1); y2 = new_coords(node2,2);
            
            edge_label_pos{i} = [(x1+x2)/2, (y1+y2)/2]; % position of edge's label (optimal rates)
%             lines_ca{i} = line([x1, x2], [y1 y2],'Color',[0 0 1]); %lines are red
            [l,g] = gcwaypts(y1,x1,y2,x2,20);

            lines_ca{i} = geoshow(l,g,'displaytype','line','color','r',...
                          'markeredgecolor','r','markerfacecolor','r','marker','.');
        end
        
        pause(0.1)
        if epoch == total_epochs
            break
        end
        for j = 1:length(lines_ca)
            delete(lines_ca{j})
        end
        delete(animation_ca{1});
        delete(animation_ca{2});
    end


    
end