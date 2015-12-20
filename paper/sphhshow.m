function sphhshow

if ~exist('DATA','var')
    load ../script/gait.mat
    load ../results/spsrc/feat.mat
end

idCh = 5;
% idCycles = [3 4 5];
idCycles = [4 5 6];

idxChosenSP = [ ... 
    4 8 15 
    4 7 16
    4 7 18 ];

halignChosenSP = { ...
    'center', 'center', 'right'
    'center', 'center', 'right'
    'center', 'center', 'right'
    };

valignChosenSP = { ...
    'bottom', 'top', 'bottom'
    'bottom', 'top', 'bottom'
    'bottom', 'top', 'bottom'
    };

hoffsetChosenSP = [ ...
    0, 0.01, -0.015
    0, 0.01, -0.015
    0, 0.01, -0.015
    ];
voffsetChosenSP = [ ...
    0.1, -0.05, 0.05
    0.1, -0.05, 0.05
    0.1, -0.05, 0.05
    ];


numSelectedSP = 3;

dataId={ 1 60 1 };

% baseDesSize = 0.2;
baseDesSize = 0.04;

d = DATA.BOTH{dataId{:}};
f = FEAT.BOTH{dataId{:}};
c = d.cycles{5};

N = length(idCycles);

clf

xrange = zeros(N,1);

for k=1:N
    cidx = (f{idCh}.baseCycleID+idCycles(k)-1);
    rangeCycle = [c(cidx) c(cidx+1)];
    dataCycle = d.data(:,rangeCycle(1):rangeCycle(2),idCh);
    dataCycle = sqrt( sum( dataCycle.*dataCycle, 1 ) );
    
    idxFeature = ( f{idCh}.cycleID == cidx );
    loc        = f{idCh}.loc(idxFeature);
    scale      = f{idCh}.scale(idxFeature);
    descriptor = f{idCh}.descriptor(idxFeature,:);
    loc        = loc*( rangeCycle(2)-rangeCycle(1) )*0.01;
    x          = ((1:length(dataCycle))-1)*0.01;
    
    [ loc, idxSP ] = sort(loc);
    scale = scale(idxSP);
    descriptor = descriptor(idxSP,:);
    
    spY = interp1( x, dataCycle, loc );
    
    xrange(k) = max(x);
    
    subplot( N, 2, k*2-1 );
    plot( x, dataCycle, 'b-' );
    hold on
    for u=1:length(loc)
        ds = baseDesSize*scale(u);
        line( loc(u)+[-ds ds]/2, spY(u)+[0 0], [0 0], 'Color', 'green' );
    end
    plot( loc, spY, 'or' );
    plot( loc, spY, '.g' );
    
    for u=1:numSelectedSP
        text( loc(idxChosenSP(k,u))+hoffsetChosenSP(k,u), ...
            spY(idxChosenSP(k,u))+voffsetChosenSP(k,u), ['\omega_' int2str(u)], ...
            'VerticalAlignment', valignChosenSP{k,u}, 'HorizontalAlignment', halignChosenSP{k,u} );
    end
    
    hold off
    
    subplot( N, 2, k*2 );
        hold on

    set(gca, 'XLim', [0 numSelectedSP], ... 
        'XTick', 0:1/(size(descriptor,2)-1):numSelectedSP, ...
        'box', 'off' );
    
    desXTickS = arrayfun( @(t) { int2str_interval(t,5,1) } , 1:size(descriptor,2) );
    desXTick = [ desXTickS(1:end-1) , ...
        repmat([ {[desXTickS{end} '/' desXTickS{1}]}, desXTickS(2:end-1) ],1,numSelectedSP-1), ...
        desXTickS(end) ];
    set(gca, 'XTickLabel', desXTick );
    
    ylimDes = 0.5;
    set(gca, 'YAxisLocation', 'right', 'YLim', [0 ylimDes], ...
        'YTick', 0:ylimDes/2:ylimDes );
    
    for u=1:numSelectedSP
        des = descriptor(idxChosenSP(k,u),:);
        des = des./norm(des);
        plot( (u-1)+(0:1/(length(des)-1):1), des, 'b.-', 'MarkerSize', 6 );
        text( u-0.5, ylimDes, ['{\bfs}(\omega_' int2str(u) ')'], ...
            'VerticalAlignment', 'top', ...
            'HorizontalAlignment', 'center' );
    end
    
    for u=1:numSelectedSP
        line( [u u]-1, [0 ylimDes*0.8], [0 0], 'Color', [1 1 1]*0.7, 'LineStyle', '--' );
    end
    hold off
    
end

mxrange = max(xrange);

for k=1:N
    subplot( N, 2, k*2-1 );
    set( gca, 'box', 'off', ...
        'xlim', [0 mxrange] );
    xlabel( 'Time / s' );
    ylabel( 'Acceleration / g' );
    if k~=N
        title( ['Step cycle ' int2str(k) ] );
    else
        title( 'Step cycle N' );
    end
    
    subplot( N, 2, k*2 );
    xlabel( 'Dimension index' );
    ylabel( 'Normalized strength' );
    if k~=N
        title( ['SPs on Step cycle ' int2str(k) ] );
    else
        title( 'SPs on from Step cycle N' );
    end

end

end

function s = int2str_interval(a,n,m)

if mod(a,n)~=m
    s = '';
else
    s = int2str(a);
end

end
