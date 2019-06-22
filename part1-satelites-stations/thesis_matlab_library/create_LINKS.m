function out = create_LINKS(coords, nodes, n)
    DISTANCES = zeros(n,n); %(N+M)x(N+M) , N: #sats, M: #stats
    LINKS = zeros(n,n); %(N+M)x(N+M)^2
    communication_range = 150; % links exist at a distance smaller or equal of communication_range
    for i = 1:n
        for j = 1:n
            DISTANCES(i,j) = euclidean_dist(coords(i,:),coords(j,:)); % calculating the distance between object's center point
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
    out = LINKS;
end