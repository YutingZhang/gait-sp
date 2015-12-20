function [SArr Info] = readEXPResult( method , param_range, working_dir0 )

load score_param

if ~exist( 'working_dir0', 'var' )
    working_dir0 = '';
end

%% methods

switch method 
    case 'spsrc'
        working_dir = 'spsrc';
        norm_method = 'none';
    case 'spsrc2'
        working_dir = 'spsrc2';
        norm_method = 'none';
    case 'spsrc_nn'
        working_dir = 'spsrc_nn';
        norm_method = 'none';
    case 'gsig-no-cyc'
        working_dir = 'gsig-no-cyc';
        norm_method = 'exp';
    case 'gsig-with-cyc'
        working_dir = 'gsig-with-cyc';
        norm_method = 'exp';
    case 'cycle-matching'
        working_dir = 'cycle-matching';
        norm_method = 'neg';
    case 'cycle-template-mean'
        working_dir = 'cycle-template-mean';
        norm_method = 'neg';
        if isempty(param_range)
            param_range.l_order = 2;
        end
    case 'spvoting'
        working_dir = 'spvoting';
        norm_method = 'none';
    case 'eigen-step'
        working_dir = 'eigen-step';
        norm_method = 'neg';
    otherwise
        error( 'No such method' );
end

if ~isempty( working_dir0 )
    working_dir = working_dir0;
end

statistic_param.PreNormFunName = norm_method;

%% seperated statistic
target_dir = ['../results/' working_dir ];
[SArr ParamInfo] = zgReadExpResults( target_dir, param_range, score_param, { statistic_param } );

%% combined statistics

Info = zgResultArrStatistic( SArr );

Info.ParamInfo = ParamInfo;

end
