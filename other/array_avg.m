function Dist = array_avg(Dists)
n = length(Dists);
minLength = min(cellfun(@length, Dists));
truncated_arrays = cellfun(@(x) x(1:minLength), Dists, 'UniformOutput', false);

result = zeros(1, minLength);
for i = 1:n
    result = result + truncated_arrays{i};
end

Dist = result/n;
end

