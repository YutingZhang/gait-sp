function D = zgGSignature_WithCycNorm( RECORD )

%% Call non-normalized version
D = zgGSignature_Base( RECORD );

%% Processing
D = zgGSignature_Base2WithCycNorm( D, RECORD );

end
