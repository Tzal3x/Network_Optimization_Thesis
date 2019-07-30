function GraphMap(num_stations, num_satellites, coords_ca, opt_results_ca, opt_buffers_ca, opt_divergencies_ca, xij_ca, nodes)
    % Create a Graph projected on the World Map. Similar to GraphDists
 
    n = length(nodes);
    figure;
    worldmap('World')
    load coastlines %#ok<LOAD>
    plotm(coastlat,coastlon,'Color',[0 0 1, 0.05]) % Transparent map for less visual complexity
    disp('[~Report:] Inside GraphMap: ========================================|')
    total_epochs = length(coords_ca);


    % Showing lifetime coordinates:
    for node = 1:num_satellites%length(nodes)
        node_coords = nodes(node).lifetime_coordinates';
        x = node_coords(:,1);
        y = node_coords(:,2);
        z = node_coords(:,3);
        long = pi/2 - acos(z./sqrt(x.^2+y.^2+z.^2)); % polar - latitude
        lat = atan2(y,x); % polar - longitude
        if node <= num_satellites
            color = [0 0 1,0.03];
        else
            color = [1 0 0,0.03];
        end            
        geoshow(rad2deg(long),rad2deg(lat),'LineStyle','--','displaytype','line','color',color,...
                      'markeredgecolor','r','marker','none');
        
%         geoshow(rad2deg(long), rad2deg(lat),'color',color,...
%                       'markeredgecolor','r','marker','.'); % KOULO!

    end

    for epoch = 1:total_epochs
        disp('>> Epoch ' + string(epoch) + '- - - - - - - - - - - - - - - - - -|')

        xij =  xij_ca{epoch};
        opt_results = opt_results_ca{epoch}(1:(length(opt_results_ca{epoch})-n)) ; % (leaving this here for later)
        opt_buffers = opt_buffers_ca{epoch}; % buffers of nodes of given epoch
        opt_divergencies = opt_divergencies_ca{epoch}; % divergencies of nodes of given epoch
        disp('> opt_results:'); disp(opt_results);
        disp('> opt_divergencies:'); disp(opt_divergencies);
        disp('> opt_buffers:'); disp(opt_buffers);

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
        animation_ca{1} = scatterm(long(1:(len-num_stations)),lat(1:(len-num_stations)),'blue','filled'); % satellite nodes
        animation_ca{2} = scatterm(long((len-num_stations+1):len),lat((len-num_stations+1):len),50,'red','filled'); % station nodes
        
        % Adding NODE info (divergence, and buffers at the end of current
        % and previous epoch. s(t) & s(t-1) respectively.
        node_texts_ca = {};
        for i = 1:n
            if i <= num_satellites
                onoma = "D"; % Satellite = (D)oryforos
                colorr = 'blue';
                AA = string(i);
            else
                onoma = "S"; % Station = (S)tathmos
                colorr = 'red';
                AA = string(i - num_satellites);
            end
            % Adding previous buffer
            if epoch > 1 
              prev_buf = opt_buffers_ca{epoch-1}(i);
            else
              prev_buf = 0;
            end
            d1 = string(round(opt_divergencies(i),2));
            b1 = string(round(opt_buffers(i),2));
            b2 = string(round(prev_buf,2));
            
            label = '  ' + onoma + AA + ':[' + d1 +', '+ b2 +', '+ b1 +']'; % text of node i
            disp('> labels:')
            disp(label)

            x = long(i);
            y = lat(i);
            node_texts_ca{i} = textm(x,y,label,'Color',colorr);
        end   
        
        % Adding edges (connecting lines):
        lines_ca = {}; % saving line objects to be able to delete them for animation purposes
        edge_label_pos = {};
        counter = 0; % counts number of lines_ca cell array's elements
        for i = 1:length(outgoing)
            node1 = outgoing(i); 
            node2 = incoming(i);
            x1 = new_coords(node1,1); y1 = new_coords(node1,2);
            x2 = new_coords(node2,1); y2 = new_coords(node2,2);
            
            [l,g] = gcwaypts(y1,x1,y2,x2,20);
            
            deviation = [linspace(0,2,20/2 + 1),linspace(2,0,20/2)];

            if outgoing(i) < incoming(i)
                temp1 = l+deviation'; temp2 = g+deviation';
                
                % Lines with arrows:
                first_piece = 1:(length(temp1)-3); % body of arrow
                second_piece = (length(temp1)-2):(length(temp1)-1); % arrow
                counter = counter + 1;
                lines_ca{counter} = geoshow(temp1(first_piece), temp2(first_piece),'displaytype','line','color',[.5 0 .5],...
                          'markeredgecolor','r','marker','none');
                
                counter = counter + 1;
                lines_ca{counter} = geoshow(temp1(second_piece), temp2(second_piece),'displaytype','line','color',[.5 0 .5],...
                          'markeredgecolor','r','marker','diamond');

                edge_label_pos{i} = [temp1(ceil(length(temp1)/2)), temp2(ceil(length(temp2)/2))];
            else
                temp1 = l-deviation'; temp2 = g-deviation';
                
                first_piece = 1:(length(temp1)-3); % body of arrow
                second_piece = (length(temp1)-2):(length(temp1)-1); % arrow
                
                counter = counter + 1;
                lines_ca{counter} = geoshow(temp1(first_piece),temp2(first_piece),'displaytype','line','color',[.5 0 .5],...
                          'markeredgecolor','r','marker','none');
                
                counter = counter + 1;
                lines_ca{counter} = geoshow(temp1(second_piece),temp2(second_piece),'displaytype','line','color',[.5 0 .5],...
                          'markeredgecolor','r','marker','diamond');
                

                edge_label_pos{i} = [temp1(ceil(length(temp1)/2)), temp2(ceil(length(temp2)/2))];
            end
        
        end
        
        
        % Adding edge flow rates: - - - - - - - - -
        edge_texts_ca = {};
        for i = 1:length(opt_results)
            pos_x = edge_label_pos{i}(1);
            pos_y= edge_label_pos{i}(2);
            edge_label = "["+string(outgoing(i))+"->"+string(incoming(i))+"]:"+string(round(opt_results(i),2));

            if outgoing(i) < incoming(i)
                edge_texts_ca{i} = textm(pos_x, pos_y, edge_label,'color',[.5 0 .5]);
            else
                edge_texts_ca{i} = textm(pos_x+2, pos_y, edge_label,'color',[.5 0 .5]);
            end
        end


        if epoch == total_epochs
            break
        end
        pause; % Delete current displays:
        for j = 1:length(lines_ca)
            delete(lines_ca{j})
        end
        for j = 1:length(edge_texts_ca)
            delete(edge_texts_ca{j})
        end
        for j = 1:length(node_texts_ca)
            delete(node_texts_ca{j})
        end
        delete(animation_ca{1}); % delete satellite nodes figures
        delete(animation_ca{2}); % delete station nodes figures
    end    
end