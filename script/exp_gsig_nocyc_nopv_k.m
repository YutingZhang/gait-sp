try
    matlabpool open local.copy1
catch e
    disp('use default matlab pool configuration');
    matlabpool open
end

param_range.p = 0.08;
% param_range.p = [0.08 0.02:0.01:0.07 0.09:0.11];
% param_range.p = 0.18:0.01:0.32;
param_range.kernel_type = 0;
param_range.p_t = -1;
param_range.p_v = -1;
param_range.nearest_k = [1 5 20 40];
crossTest('gsig-no-cyc',param_range)
matlabpool close
