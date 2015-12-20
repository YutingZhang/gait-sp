function savepdf( filename )

print( '-depsc2', [ filename '.eps' ] );
system( [ 'epstopdf ' filename '.eps' ] );
system( [ 'pdfcrop ' filename '.pdf' ] );
delete( [ filename '.pdf' ] );
movefile( [ filename '-crop.pdf' ], [ filename '.pdf' ] );

end
