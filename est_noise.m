function [sigma, DataCube_clean] = est_noise(DataCube)
[m, n, p] = size(DataCube);
DataCube_clean = zeros(m, n, p);
data = reshape(DataCube, [p, m * n])';
sigma = zeros(p, 1);

for i = 1:p
    y = data(:, i);
    X = data;
    X(:, i) = [];
    y_clean = X * ((X' * X) \ (X' * y));
    DataCube_clean(:, :, i) = reshape(y_clean, m, n);
    sigma(i) = std(y_clean-y);
end
