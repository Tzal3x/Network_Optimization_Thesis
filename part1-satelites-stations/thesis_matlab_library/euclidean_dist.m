function out = euclidean_dist(vec1, vec2)
% Definition: euclidean_dist(vec1, vec2)
% Calculates euclidean distance of two vectors
    out = sqrt(sum((vec1 - vec2) .^ 2));
end