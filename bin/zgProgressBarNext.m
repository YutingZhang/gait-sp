function pb = zgProgressBarNext( pb0, inc )
pb = pb0;
pb.current = inc + pb0.current;
pb.printed = floor(pb.current/pb.total*100);
fprintf(1, repmat('>',1,pb.printed-pb0.printed) )
if pb.current == pb.total
    fprintf(1, '\n' )
end
end
