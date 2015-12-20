function D = zgGSignature_NoCycNorm( RECORD )

%% Use default parameters
D = zgGSignature_Base( RECORD );

%% Processing
D = zgGSignature_Base2NoCycNorm( D, RECORD );

end
