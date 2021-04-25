clear;
clc;

for r = 100:100:200
    for mu = [0, 1,10,100]
        opts = struct('r', r, 'mu', mu);
        HSI_MTGSNMF_denoising(opts);
    end
end
