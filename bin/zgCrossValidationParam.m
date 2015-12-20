function zgCrossValidationParam( FEAT , score_fun , param_range , output_dir, options )
% R = zgCrossValidationS( FEAT , score_fun , param )
% score_fun: function pointer
%            or 'spvoting'

cvs_fun = @zgCrossValidationS;
if exist('options','var')
    if isfield( options, 'isAuthFold' )
        if options.isAuthFold
            cvs_fun = @zgCrossValidationS4authFold;
        end
    end
else
    options = [];
end

if ischar(score_fun)
    switch score_fun
        case 'spvoting'
            cvs_fun = @zgCrossValidationS4SPVoting;
        case 'spsrc'
            cvs_fun = @zgCrossValidationS4SPsrc;
        case 'spsrc2'
            cvs_fun = @zgCrossValidationS4SPsrc2;
        case 'spsrc_nn'
            cvs_fun = @zgCrossValidationS4SPsrc_nn;
        otherwise
            error('Unknown Score Function');
    end
end

clear('origPersonNum','featPersonNum');

if isfield( options, 'CorruptedPerson' )
    if ischar(FEAT)
        load( FEAT );
    end
    origPersonNum.both   = size(FEAT.BOTH,2);
    origPersonNum.single = size(FEAT.SINGLE,2);

    person_total_both0   = origPersonNum.both;
    person_total_single0 = origPersonNum.single;
    opn = person_total_both0 + person_total_single0;
    PIDX = true( opn, 1 );
    PIDX(options.CorruptedPerson) = false;
    FEAT.BOTH   = FEAT.BOTH( :, PIDX(1:person_total_both0), : );
    FEAT.SINGLE = FEAT.SINGLE( :, PIDX((person_total_both0+1):end), : );
    fprintf(1,'Notice: options.CorruptedPerson is applied in CrossValidation.\n');
end

if iscell( param_range )
    startIdx4Param = param_range{2};
    param_range = param_range{1};
else
    startIdx4Param = 1;
end

if isfield( options, 'ChosenRecord' )
    FEAT.BOTH   = FEAT.BOTH  (options.ChosenRecord,:,:);
    FEAT.SINGLE = FEAT.SINGLE(options.ChosenRecord,:);
    fprintf(1,'Notice: options.ChosenRecord is applied in CrossValidation.\n');
end

p_e = zgParamEmutatorInit( param_range );

[param p_e] = zgParamEmutatorNext( p_e );
idxParam = 0;
while ~iscell( param )
    ps = zgParam2Str( param );
    idxParam = idxParam + 1;
    if idxParam>=startIdx4Param
        NowTime = now;
        disp( [ 'No. ' int2str(idxParam) ': Start at ' datestr(clock,'yyyy-mm-dd HH:MM:SS') ] );
        fprintf(1,'%s: ',ps);
        [R W featPersonNum] = feval( cvs_fun, FEAT, score_fun, param, options );
        if ~exist( 'origPersonNum','var' )
            origPersonNum = featPersonNum;
        end
        save([output_dir '/' ps '.mat'],'R','W','param','score_fun','options','origPersonNum','-v7.3');
        clear R W
        disp( [ 'Elapsed time: ' datestr(now-NowTime, 'dd HH:MM:SS') ] ) ;
    end
    [param p_e] = zgParamEmutatorNext( p_e );
end

end
