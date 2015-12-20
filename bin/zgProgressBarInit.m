function pb = zgProgressBarInit( N )
pb.current = 0;
pb.printed = 0;
pb.total   = N;
fprintf(1,[repmat('=',1,100) '\n' ]);
end

