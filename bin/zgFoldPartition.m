function [elementN foldRange] = zgFoldPartition( N, foldN )

elementN  = zeros( foldN, 1 );

n = N;
for k=1:foldN
    elementN(k) = ceil(n/(foldN-k+1));
    n = n-elementN(k);
end
foldRange = [ 1 ; (cumsum(elementN)+1) ];
foldRange = foldRange(1:end-1);
foldRange = [ foldRange , (foldRange+elementN)-1 ];

end
