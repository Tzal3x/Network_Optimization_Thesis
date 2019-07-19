function out = create_DISTANCES(coords, n)
    dists = zeros(n);
    for i = 1:n
       for j = 1:n
           dists(i,j) = euclidean_dist(coords(i,:),coords(j,:));
       end
    end
    out = dists;
end