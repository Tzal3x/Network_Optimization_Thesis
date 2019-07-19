function GraphDists(num_satellites, opt_buffers_ca, opt_divergencies_ca, coords_ca, xij_ca, opt_results_ca, nodes, total_epochs)
% Creates a graph with the connecting nodes having the corresponding edge length     
    n = length(nodes);
    for epoch = 1:total_epochs
        
        % Get the current epoch's data:
        xij =  xij_ca{epoch};
        coords = coords_ca{epoch};
        opt_results = opt_results_ca{epoch}(1:(length(opt_results_ca{epoch})-n)) ; % (leaving this here for later)
        opt_buffers = opt_buffers_ca{epoch}; % buffers of nodes of given epoch
        opt_divergencies = opt_divergencies_ca{epoch}; % divergencies of nodes of given epoch
        % The above are used at node info (text that includes divergence
        % and buffers s(t) and s(t-1).
        
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
        % divs1 = [(length(nodes)+1):(length(nodes)+num_satellites),    (num_satellites+1):length(nodes)]; %source
        % divs2 = [1:num_satellites,   (length(nodes) + 1 + num_satellites):(length(nodes) * 2)]; %sink
        outgoing = [outgoing1, outgoing2]; %, divs1];
        incoming = [incoming1, incoming2]; %, divs2];
        
        % Adding nodes to graph:
        temp_scatter = scatter(coords(:,1),coords(:,2),'filled');
        earth_circle = viscircles([0 0], 200,'LineStyle','--','LineWidth',0.8,'Color',[1 0 0]); % WARNING! radius = 200 might be wrong!
        axis([-300, 300, -300, 300]);
        axis equal
        title('2D Graph of solved DTNUM problem')
        
        % Get node positions:
        node_positions_ca = {};
        for i = 1:length(coords(:,1))
            node_positions_ca{i} = [coords(i,1), coords(i,2)]; % node{i}[long, lat]
        end
        
        % Adding edges (connecting lines):
        lines_ca = {}; % saving line objects to be able to delete them for animation purposes
        edge_label_pos = {};
        for i = 1:length(outgoing)
            node1 = outgoing(i); 
            node2 = incoming(i);
            x1 = node_positions_ca{node1}(1); y1 = node_positions_ca{node1}(2);
            x2 = node_positions_ca{node2}(1); y2 = node_positions_ca{node2}(2);
            edge_label_pos{i} = [(x1+x2)/2, (y1+y2)/2];
            lines_ca{i} = line([x1, x2], [y1 y2],'Color',[0 0 1]); %lines are red
        end
         
        % Adding NODE info (divergence, and buffers at the end of current
        % and previous epoch. s(t) & s(t-1) respectively.
        node_texts_ca = {};
        for i = 1:n
            if i <= num_satellites
                onoma = "Satellite ";
                AA = string(i);
            else
                onoma = "Station ";
                AA = string(n - num_satellites);
            end
            % Adding previous buffer
            if epoch > 1 
              prev_buf = opt_buffers_ca{epoch-1}(i);
            else
              prev_buf = 0;
            end
            label = [onoma+ AA+ ":",...
                     "Divergence: "+string(opt_divergencies(i)),...
                     "s(t): " + string(opt_buffers(i)),...
                     "s(t-1): " + prev_buf]; % text of node i
                 
            x = node_positions_ca{i}(1)+6;
            y = node_positions_ca{i}(2)-3;
            node_texts_ca{i} = text(x,y,label,'Color','red');
        end        
        
        % Adding edge flow rates: - - - - - - - - -
        edge_texts_ca = {};
        for i = 1:length(opt_results)
            pos_x = edge_label_pos{i}(1);
            pos_y = edge_label_pos{i}(2);
            edge_label = [string(outgoing(i))+"->"+string(incoming(i)),string(opt_results(i))];
            if outgoing(i) < incoming(i)
                edge_texts_ca{i} = text(pos_x, pos_y, edge_label);
            else
                edge_texts_ca{i} = text(pos_x, pos_y+8, edge_label);
            end
        end
        
        % Using delete to make an animation ... 
        pause; %pause(0.8); %fps of animation
        if epoch == total_epochs
            break;
        end
        for i = 1:n
           delete(node_texts_ca{i});
        end
        for i = 1:length(lines_ca)
            delete(lines_ca{i});
        end
        delete(earth_circle);
        delete(temp_scatter);
    end
%     close all % close all figures
end