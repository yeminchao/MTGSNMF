function [U, V, C_value, t_elapsed] = MTGSNMF(X, W, r, mu, lambda, opts)

% error(nargchk(5, 6, nargin, 'struct'));

maxit = 10000;
tol = 1e-6;
verbose = true;

if exist('opts', 'var')
    if isfield(opts, 'maxit')
        maxit = opts.maxit;
    end
    if isfield(opts, 'tol')
        tol = opts.tol;
    end
    if isfield(opts, 'verbose')
        verbose = opts.verbose;
    end
end

if ~iscell(lambda)
    if numel(lambda) == 1
        temp = lambda;
        lambda = cell(size(X));
        lambda(:) = {temp};
    else
        lambda = num2cell(lambda);
    end
elseif numel(lambda) == 1
    temp = lambda{1};
    lambda = cell(size(X));
    lambda(:) = {temp};
end
if numel(lambda) ~= length(lambda) || length(lambda) ~= length(X)
    error('Size of lambda is invalid.');
end

if isa(X{1}, 'double')
    eps_val = 1e-14;
elseif isa(X{1}, 'single')
    eps_val = 1e-5;
else
    error('X invalid.');
end

if verbose
    tStart = tic;
end


K = length(X);
iscolumn_X = size(X, 2) == 1;
[n, m] = size(X{1});
X = cat(1, X{:});
class_type = class(X);
U = rand(n*K, r, class_type);
U = max(U, max(U(:))*eps_val);

V = rand(r, m, class_type);
V = max(V, max(V(:))*eps_val);

C = inf;

lambda = sum(cell2mat(lambda));

D = diag(sum(W, 2));
L = D - W;

for itrcount = 1:maxit
    UtX = U' * X;
    V = V .* (UtX + mu * V * W') ./ max(U'*U*V+mu*V*D'+lambda/2, max(UtX(:))*eps_val);
    V = max(V, max(V(:))*eps_val);
    XVt = X * V';
    U = U .* XVt ./ max(U*(V * V'), max(XVt(:))*eps_val);
    U = max(U, max(U(:))*eps_val);
    C_old = C;
    C = norm(X-U*V, 'fro')^2 + mu * sum(sum((V * L).*V)) + lambda * sum(sum(V));
    time_elapsed = toc(tStart);
    if verbose
        fprintf('Iteration %d.\tObjective function value C=%f\n', itrcount, C);
        fprintf('Elapsed time is %f seconds.\n', time_elapsed);
    end
    t_elapsed(itrcount) = time_elapsed;
    C_value(itrcount) = C;
    if C_old - C < tol * C && C_old >= C
        break;
    end
end

U = mat2cell(U, n*ones(1, K), r);
if ~iscolumn_X
    U = U';
end

end
