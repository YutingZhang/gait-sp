function D = zgFeature_AllCycles( RECORD, NormalizedLen )

if ~exist('NormalizedLen','var')
    NormalizedLen = 100;
end

CN = size( RECORD.data ,3 );
D=cell(1,CN);

cycles = RECORD.cycles{5};  % defautly use channel 5
% cycleLEN = diff(cycles);

[ cycST cycEN ] = zgCycleEnds( cycles , RECORD.splitters );

data = shiftdim( sqrt(sum(RECORD.data.*RECORD.data)), 1 );
for c=1:CN
    ds = data(:,c);
    A  = zeros( NormalizedLen, cycEN-cycST );
    % split & normalize cycles
    for k=cycST:(cycEN-1)
        st = cycles(k);
        en = cycles(k+1);
        X = (0:1/(en-st):1)*(NormalizedLen-1)+1;
        Y = ds( st:en );
        A(:,k-cycST+1) = interp1( X , Y , 1:NormalizedLen );
    end
    
    data(:,c);
    D{c} = A;
end

end
