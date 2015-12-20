try
    matlabpool open local.copy3
catch e
    disp('use default matlab pool configuration');
    matlabpool open
end
% param_range.p = 0.18:0.01:0.32;
% param_range.p = [ 0.05 0.02:0.01:0.04 0.06:0.01:0.12 ];
param_range.p = 0.08;
param_range.kernel_type = 0;
param_range.p_t = -1;
param_range.p_v = -1;

% param_range.cycleNUM = 1:10;
% crossTest('gsig-with-cyc',param_range);

param_range.cycleNUM = -1;
param_range.cycleNUM_SAMPLE   = 1:10;
param_range.cycleNUM_EXEMPLAR = inf;
crossTest('gsig-with-cyc',param_range);

% param_range.cycleNUM = -1;
% param_range.cycleNUM_SAMPLE   = inf;
% param_range.cycleNUM_EXEMPLAR = 1:10;
% crossTest('gsig-with-cyc',param_range);

matlabpool close
