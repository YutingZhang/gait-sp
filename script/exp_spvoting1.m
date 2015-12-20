try
    try
        matlabpool close
    catch e
    end
    matlabpool open local.copy1
catch e
    disp('use default matlab pool configuration');
    matlabpool open
end

param_range = [];
options = [];
% param_range.k = 1:3;
load score_param
options.CorruptedPerson = score_param.CorruptedPerson;

% crossTest('spvoting',param_range,options)

param_range.cycleNUM   = 1:10;
crossTest('spvoting',param_range,options)

% param_range.cycleNUM = -1;
% param_range.cycleNUM_SAMPLE   = 1:10;
% param_range.cycleNUM_EXEMPLAR = inf;
% crossTest('spvoting',param_range,options)
% 
% param_range.cycleNUM = -1;
% param_range.cycleNUM_SAMPLE   = inf;
% param_range.cycleNUM_EXEMPLAR = 1:10;
% crossTest('spvoting',param_range,options)

matlabpool close
