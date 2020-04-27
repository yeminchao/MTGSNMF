function [DataCubeOut, U, V, C_value, t_elapsed] = MTGSNMF_denoising(DataCube, patch_size, stride, r, W, mu, lambda)

opts.maxit = 10000;
opts.tol = 1e-6;
opts.verbose = true;

[m, n, p] = size(DataCube);
patches = cell(p, 1);

for k = 1:p
    patches{k} = im2colstep(DataCube(:, :, k), patch_size, [stride, stride]);
end

[U, V, C_value, t_elapsed] = MTGSNMF(patches, W, r, mu, lambda, opts);

DataCubeOut = zeros(m, n, p);

countImg = countcover([m, n], patch_size, [stride, stride]);

for k = 1:p
    DataCubeOut(:, :, k) = col2imstep(U{k}*V, [m, n], patch_size, [stride, stride]) ./ countImg;
end

end
