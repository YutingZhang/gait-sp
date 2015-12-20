function D = zgGSignature_Base( RECORD )

%% Use default parameters
dspHalfSize = 10;
Smax = 3;
sigma0 = 1.2;
sigmaN = 0.2;
layer_per_oct = 1;


%% Processing
CN = size( RECORD.data ,3 );
D=cell(1,CN);

data = shiftdim( sqrt(sum(RECORD.data.*RECORD.data)), 1 );

for c=1:CN
    [ A.descriptor A.loc A.scale ]= zgSIFT1D( data(:,c), dspHalfSize, Smax, sigma0, sigmaN, layer_per_oct, 0.008/layer_per_oct/sigma0 );
    D{c} = A;
end

end
