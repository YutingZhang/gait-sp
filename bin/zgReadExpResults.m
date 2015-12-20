function [ SArr ParamInfo ] = zgReadExpResults( TargetDir, param_range, ...
    score_param, statistic_options )

%% read dir
EXP_FILES = dir( [ TargetDir '/*;.mat' ] );
fns = {EXP_FILES.name};

PArr = cell(length(fns),1);

for k=1:length(fns)
    PArr{k} = zgStr2Param( fns{k} );
end

%% parse param_range

pre = zgParamEmutatorInit( param_range );

ParamInfo_Element.range = [];
ParamInfo_Element.field = '';
ParamInfo = repmat( ParamInfo_Element, length(pre.fields), 1 );
for k=1:length(pre.fields)
    eval( [ 'ParamInfo(k).range = param_range.' pre.fields{k}  ';' ] );
    ParamInfo(k).field = pre.fields{k};
end

%% parse other score_param
if ~isfield( score_param, 'CorruptedPerson' )
    score_param.CorruptedPerson = [];
end

%% Iterate param_range

SArr = cell( [pre.total 1] );

[ paramc pre idx ] = zgParamEmutatorNext( pre );
while ~iscell(paramc)
    fprintf(1,'%s: ',zgParam2Str(paramc));
    % searching file
    selectedFileID = 0;
    for k=1:length(PArr)
        fpara = PArr{k};
        eqr = false;
        if isstruct(fpara)
            if ~isempty( paramc )
                eqr = true;
                for m = 1:length(pre.fields)
                    if ~isfield( fpara, pre.fields{m} )
                        eqr = false; break;
                    end
                    eval(['t = fpara.' pre.fields{m} ';']);
                    eval(['u = paramc.' pre.fields{m} ';']);
                    if abs( t - u ) > 1e-5  % eps
                        eqr = false; break;
                    end
                end
            end
        elseif isempty( fpara )
            if isempty( paramc )
                eqr = true;
            end
        end
        if eqr
            if isempty( paramc )
                selectedFileID = k;
            elseif selectedFileID<1
                selectedFileID = k;
            elseif length( fieldnames( PArr{k} ) )<length( fieldnames( PArr{selectedFileID} ) )
                selectedFileID = k;
            end
        end
    end

    if selectedFileID>0
        k = selectedFileID;
        fprintf(1,'%s: ', fns{k});
        clear( 'R','W','param','score_fun','options','origPersonNum' );
        load( [ TargetDir '/' fns{k}] );
        if ~exist( 'origPersonNum', 'var' )
            load gait_info
            origPersonNum.both   = DATA_BOTH_NUM;
            origPersonNum.single = DATA_SINGLE_NUM;
        end
        Corruputed_Test = [];
        if exist( 'options', 'var' )
            if isfield( options, 'CorruptedPerson' )
                Corruputed_Test = options.CorruptedPerson;
            end
        end
        [ R W PersonNum ] = zgEliminateCorrupted( R, W, origPersonNum, ...
            score_param.CorruptedPerson, Corruputed_Test );

        c_sub = num2cell(idx);
        s_ind = sub2ind( size(SArr), c_sub{:} );
        SArr{s_ind} = zgResultStatistic( R, W, PersonNum, statistic_options{:} );
        SArr{s_ind} = rmfield( SArr{s_ind}, 'IDENT_PD' );
        clear R W
        fprintf(1,'OK\n');
    else
        fprintf(1,'Warning - Empty\n');
    end
    
    [ paramc pre idx ] = zgParamEmutatorNext( pre );
end

% SArr = cell2mat( SArr );

end
