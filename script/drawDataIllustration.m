
idxData = 1;

A = squeeze( sqrt( sum( DATA.BOTH{idxData}.data.*DATA.BOTH{idxData}.data,1 ) ) );

L = {'Right wrist', 'Left upper arm', 'Right side of pelvis', 'Left thigh', 'the right ankle' };

C = DATA.BOTH{idxData}.cycles{5};
S = DATA.BOTH{idxData}.splitters;

XL = 1400;

clf
for k=1:5
    subplot( 5, 1, k );
   
    P = get(gca, 'Position');
    P(4) = P(4)*0.75;
    set( gca, 'Position', P );
   
    set( gca, 'box', 'off', ...
        'YLim', [0 4], 'YTick', 0:4, 'YTickLabel', {' '}, ...
        'XLim', [0 XL], 'XTick', 0:100:XL, 'XTickLabel', {' '} );
    
    xlabel( ' ' );
    ylabel( ' ' );
    title ( ' ' );



    
    hold on
    fill( [S(1) S(1) 0  0 ], [0 4 4 0], [1 1 1]*0.9, 'EdgeColor', 'none' );
    fill( [S(2) S(2) XL XL], [0 4 4 0], [1 1 1]*0.9, 'EdgeColor', 'none' );
    hold off

    
    axes
    
    set( gca, 'Position', P );
    
    set( gca, 'box', 'off', ...
        'YLim', [0 4], 'YTick', 0:4, ...
        'XLim', [0 XL], 'XTick', 0:100:XL);
    set( gca, 'Color', 'none' );
    
    xlabel( 'time / ms' );
    ylabel( 'acc. / g' );
    title ( L{k} );

    set( gca, 'Position', P );
 
    hold on
    
    plot( A(:,k), '-k', 'LineWidth', 0.5 );
    
    
    for r = 1:length(C)
        line( [C(r) C(r)], [0 4], [0 0], 'LineWidth', 0.5, ...
            'Color', 'blue', 'LineStyle', '--' );
    end
    
    if k==5
        plot( C, A(C,5), 'xr', 'LineWidth', 1 );
    end

    
    
    hold off
end
