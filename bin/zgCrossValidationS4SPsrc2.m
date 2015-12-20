function [R W paramPersonNum] = zgCrossValidationS4SPsrc2( FEAT , score_fun , param, options )

    options.srcType = 2;
    [R W paramPersonNum] = zgCrossValidationS4SPsrc( FEAT , score_fun , param, options );

end
