load ../results/spsrc/feat.mat

subjList = 1:10;

idCh = 5;

numClusters = 9;
numCylcePerSubj = 15;
numSubj = length(subjList);

CM = vivid(numClusters,[0.1 0.45]);

D = cell(numSubj,1);
L = cell(numSubj,1);
C = cell(numSubj,1);
I = cell(numSubj,1);
S = cell(numSubj,1);

for k=1:numSubj
    S_k = FEAT.BOTH{1,subjList(k),1}{idCh}.scale;
%     validIdx = and(S_k<2,S_k>0.4);
    validIdx = true(size(S_k));
    S{k} = FEAT.BOTH{1,subjList(k),1}{idCh}.scale(validIdx);
    D{k} = FEAT.BOTH{1,subjList(k),1}{idCh}.descriptor(validIdx,:);
    L{k} = FEAT.BOTH{1,subjList(k),1}{idCh}.loc(validIdx);
    C{k} = FEAT.BOTH{1,subjList(k),1}{idCh}.cycleID(validIdx);
    I{k} = repmat( k, length(L{k}), 1 );
end
D = cell2mat( D );
L = cell2mat( L );
C = cell2mat( C );
I = cell2mat( I );
S = cell2mat( S );


IDX = kmeans( [D], numClusters, 'start', 'cluster' );



figure(1)
clf
for k=1:numClusters
    sidx = (IDX==k);
    plot( L(sidx), I(sidx)*numCylcePerSubj-C(sidx), 'd', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', CM(k,:), 'MarkerSize', 4 );
    hold on
end

for p=1:numSubj-1
    h = line( [0 1], [1 1]*p*numCylcePerSubj );
    set( h, 'Color', [1 1 1]*0.7, 'LineWidth', 0.5 );
end

set(gca, 'XLim', [0 1], 'XTick', 0:0.1:1 );
set(gca, 'YLim', [0 numCylcePerSubj*numSubj], ...
    'YTick', (1:numSubj)*numCylcePerSubj-numCylcePerSubj/2, ...
    'YTickLabel', arrayfun( @(x) {int2str(x)}, numSubj:-1:1 ) );

xlabel( 'In-cycle relative location' );
ylabel( 'Subject ID' );

hold off

figure(2)
clf
for k=1:numClusters
    sidx = (IDX==k);
    L_k = L(sidx);
%     [~,ds,xmesh]=kde( [L_k;L_k(L_k>0.7)-1;L_k(L_k<0.3)+1], 2^14, -1, 2 );
%     [ds,xmesh] = ksdensity( [L_k;L_k(L_k>0.7)-1;L_k(L_k<0.3)+1] );
    [ds,xmesh] = ksdensity( [L_k;L_k-1;L_k+1], 0:1e-4:1, 'width', 0.05 );
    ds = ds.*length(L_k);
    validIdx = and( xmesh>=0 , xmesh<=1 );
%     plot( xmesh, ds, 'Color', CM(k,:) );
    plot( xmesh(validIdx), ds(validIdx), 'Color', CM(k,:) );
    hold on;
end

ylabel('Density');
xlabel('In-cycle relative location');
yrange=get(gca,'YLim');
line( [0 1], [0 0], [0 0], 'Color', 'black' )
line( [0 0], yrange(2), [0 0], 'Color', 'black' );

set(gca, 'XLim', [0 1], 'XTick', 0:0.1:1 );

hold off


