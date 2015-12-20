function feat_fun_param = zgCall_InitParamFun( fun_handle, DATA, options )

if isfield( options, 'CorruptedPerson' )

    origPersonNum.both  = size(DATA.BOTH,2);
    origPersonNum.single = size(DATA.SINGLE,2);

    person_total_both0   = origPersonNum.both;
    person_total_single0 = origPersonNum.single;
    opn = person_total_both0 + person_total_single0;
    PIDX = true( opn, 1 );
    PIDX(options.CorruptedPerson) = false;
    DATA.BOTH   = DATA.BOTH( :, PIDX(1:person_total_both0), : );
    DATA.SINGLE = DATA.SINGLE( :, PIDX((person_total_both0+1):end), : );
    fprintf(1,'Notice: options.CorruptedPerson is applied in feature extraction parameter initialization.\n');
end

feat_fun_param = fun_handle( DATA, options );

end

