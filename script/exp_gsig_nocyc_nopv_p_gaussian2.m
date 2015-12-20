try
    matlabpool open local.copy1
catch e
    disp('use default matlab pool configuration');
    matlabpool open
end

% param_range.p = 0.03:0.02:0.11;
% param_range.p = 0.18:0.01:0.32;
param_range=[];
param_range.p = 0.41:0.01:0.5;
param_range.kernel_type = 1;
param_range.p_t = -1;
param_range.p_v = -1;
param_range.cycleNUM = inf;
crossTest('gsig-no-cyc',param_range)
matlabpool close
