% try
%     try
%         matlabpool close
%     catch e
%     end
%     matlabpool open local
% catch e
%     disp('use default matlab pool configuration');
%     matlabpool open
% end

try
    matlabpool open
catch e
end

param_range = [];
options = [];

load score_param
options.CorruptedPerson = score_param.CorruptedPerson;

param_range.hasAuth = 0;
% param_range.wordNumBase = 8:2:12;
% param_range.chosenClusterNum = [1 2 4];
% param_range.wordNumBase = [10 15 20 30 40 50 ];

% param_range.wordNumBase = [10 15 30 40 50 ];
% param_range.chosenClusterNum  = 2;
% param_range.chosenExemplarNum = [10 20 30 40 50];

idxStartRange = 1;
param_range.wordNumBase = 100;
param_range.chosenClusterNum  = 1;
% param_range.chosenExemplarNum = [ 10 20 40 60 ];
% param_range.chosenExemplarNum = [ 4:2:8 12:2:16 ]; %10;
param_range.chosenExemplarNum = 8;
param_range.softmaxSigma      = 5;


param_range.cycleNUM = inf;
crossTest('spsrc',{ param_range idxStartRange },options);

% param_range.cycleNUM   = 1:10;
% crossTest('spvoting',param_range,options)

% param_range.cycleNUM = -1;
% param_range.cycleNUM_SAMPLE   = 1:10;
% param_range.cycleNUM_EXEMPLAR = inf;
% crossTest('spvoting',param_range,options)

% param_range.cycleNUM = -1;
% param_range.cycleNUM_SAMPLE   = inf;
% param_range.cycleNUM_EXEMPLAR = 1:10;
% crossTest('spvoting',param_range,options)

try
    matlabpool close
catch e
end