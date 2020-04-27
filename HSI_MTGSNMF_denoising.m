function HSI_MTGSNMF_denoising(opts)

if ~exist('data', 'dir')
    mkdir('data');
end

load('DataCubeNoised.mat', 'DataCube');

ts = tic;

min_val = min(min(min(DataCube)));
max_val = max(max(max(DataCube)));
DataCube = (DataCube - min_val) / (max_val - min_val);

%%%%%%%%%%%%%%%%%%%%%  parameter setting  %%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('opts', 'var')
    opts = struct();
end

% patch size
if ~isfield(opts, 'patch_size')
    opts.patch_size = [7, 7];
end

% stride of overlapping patch sampling (the offset between two neighboring patches)
if ~isfield(opts, 'stride')
    opts.stride = 1;
end

% dictionary size
if ~isfield(opts, 'r')
    opts.r = round(4.0*prod(opts.patch_size));
end

% graph regularization parameter
if ~isfield(opts, 'mu')
    opts.mu = 100;
end

if exist(fullfile('data', sprintf('DataCubeOut_MTGSNMF-r=%d-mu=%g.mat', opts.r, opts.mu)), 'file')
    return;
end

% noise estimation, separate the estimated clean signal
[sigma, DataCube_clean] = est_noise(DataCube);

% graph weights
graph_data_file = fullfile('data', 'Graph.mat');
if exist(graph_data_file, 'file')
    load(graph_data_file, 'W');
else
    p = size(DataCube, 3);
    patches = cell(p, 1);
    for k = 1:p
        patches{k} = im2colstep(DataCube_clean(:, :, k), opts.patch_size, [opts.stride, opts.stride]);
    end

    patches = cat(1, patches{:})';
    num_patches = size(patches, 1);
    cluster_idx = sparse(kmeans(patches, round(sqrt(num_patches))));
    W = double(bsxfun(@eq, cluster_idx, cluster_idx'));
    save(graph_data_file, 'W');
    clear patches cluster_idx
end
clear DataCube_clean

% sparse regularization parameter
if ~isfield(opts, 'lambda')
    opts.lambda = sigma * sqrt(2*log(opts.r));
end

%%%%%%%%%%%%%%%%%%%  end parameter setting  %%%%%%%%%%%%%%%%%%%%

[DataCubeOut, U, ~, C_value, t_elapsed] = MTGSNMF_denoising(DataCube, opts.patch_size, opts.stride, opts.r, W, opts.mu, opts.lambda);

DataCubeOut = max(DataCubeOut, 0);
DataCubeOut = min(DataCubeOut, 1);
DataCubeOut = DataCubeOut * (max_val - min_val) + min_val;

t_total = toc(ts);

DataCube_clean = load('DataCube.mat', 'DataCube');
DataCube_clean = DataCube_clean.DataCube;

PSNR_val = psnr(DataCube_clean, DataCubeOut);
SSIM_val = ssim_index3d(DataCube_clean, DataCubeOut);
save(fullfile('data', sprintf('DataCubeOut_MTGSNMF-r=%d-mu=%g.mat', opts.r, opts.mu)), 'DataCubeOut', 'U', 'C_value', 't_elapsed', 'sigma', 't_total', 'opts', 'PSNR_val', 'SSIM_val');
