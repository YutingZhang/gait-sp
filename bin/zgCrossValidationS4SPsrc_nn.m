function [R W paramPersonNum] = zgCrossValidationS4SPsrc_nn( FEAT , score_fun , param, options )

    options.srcType = 3;
    [R W paramPersonNum] = zgCrossValidationS4SPsrc( FEAT , score_fun , param, options );

end
