function crossTest( method, param_range, options, working_dir0 )

if ~exist( 'options', 'var' )
    options = [];
end

if ~exist( 'working_dir0', 'var' )
    working_dir0 = '';
end

method = lower(method);

if iscell( param_range )
    idxStartRange = param_range{2};
    param_range   = param_range{1};
else
    idxStartRange = 1;
end

ff_options = [];
param_init_fun = @empty_cell;
switch method 
    case 'spsrc'
        working_dir = 'spsrc';
        feat_fun  = @zgGSignature_WithCycNorm;
        score_fun = 'spsrc';
        param_range = fill_param4gsig_src(param_range);
    case 'spsrc2'
        working_dir = 'spsrc2';
        feat_fun  = @zgGSignature_WithCycNorm;
        score_fun = 'spsrc2';
        param_range = fill_param4gsig_src(param_range);
    case 'spsrc_nn'
        working_dir = 'spsrc_nn';
        feat_fun  = @zgGSignature_WithCycNorm;
        score_fun = 'spsrc_nn';
        param_range = fill_param4gsig_src(param_range);
    case 'gsig-kde-no-cyc'
        working_dir = 'gsig-no-cyc';
        feat_fun  = @zgGSignature_NoCycNorm;
        score_fun = @zgScore_GSig;
        param_range = fill_param4gsig_kde(param_range);
    case 'gsig-kde-with-cyc'
        working_dir = 'gsig-with-cyc';
        feat_fun  = @zgGSignature_WithCycNorm;
        score_fun = @zgScore_GSig;
        param_range = fill_param4gsig_kde(param_range);
    case 'cycle-matching'
        working_dir = 'cycle-matching';
        feat_fun  = @zgFeature_AllCycles;
        score_fun = @zgScore_AllCycles;
    case 'cycle-template-mean'
        working_dir = 'cycle-template-mean';
        feat_fun  = @zgFeature_CycleTemplate_Mean;
        score_fun = @zgScore_CycleTemplate;
        if ~isfield( param_range, 'l_order' )
            param_range.l_order = 2;
        end
    case 'spvoting'
        working_dir = 'spvoting';
        feat_fun  = @zgGSignature_WithCycNorm;
        score_fun = 'spvoting';
    case 'eigen-step'
        working_dir = 'eigen-step';
%         ff_options.useCyclicAllRecord = true;
        param_init_fun = @zgInitParam_EigenStep;
        feat_fun  = @zgFeature_EigenStep;
        score_fun = @zgScore_EigenStep;
        options.isAuthFold = true;
        options.PreNormFunName = 'neg';
    otherwise
        error( 'No such method' );
end

if ~isempty( working_dir0 )
    fprintf(1,'Notice: The working dir is override to `%s'' from the default `%s''\n', ...
        working_dir0, working_dir );
    working_dir = working_dir0;
end

wpath = [ '../results/' working_dir ];
if ~exist( wpath, 'dir' )
    mkdir( wpath );
end

feat_path = [ wpath '/feat.mat' ];
if ~exist( feat_path , 'file' )
    load gait.mat
    feat_fun_param = zgCall_InitParamFun( param_init_fun, DATA, options );
    FEAT = zgBatchExtractFeatures( DATA, feat_fun, ff_options, feat_fun_param );
    save(feat_path, 'FEAT', 'ff_options', 'feat_fun_param', '-v7.3');
    clear DATA FEAT
end

zgCrossValidationParam( feat_path, score_fun, {param_range idxStartRange}, wpath, options );

disp('OK');

end

function param_range = fill_param4gsig_kde( param_range )

if ( ~isfield( param_range, 'p_t' ) )
    param_range.p_t = -1;
end
if ( ~isfield( param_range, 'p_v' ) )
    param_range.p_v = -1;
end
if ( ~isfield( param_range, 'cycleNUM' ) )
    param_range.cycleNUM = inf;
end

end

function param_range = fill_param4gsig_src( param_range )

if ( ~isfield( param_range, 'cycleNUM' ) )
    param_range.cycleNUM = inf;
end

end


function C = empty_cell( varargin )
C = {};
end
