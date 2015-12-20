function [cycST cycEN] = zgCycleEnds( cycles, splitters )

extent_portion = 0.2;

Splitters = splitters;
SIDX = quantiz(Splitters,cycles);
cycleLEN = diff(cycles);

if SIDX(1)<1
    cycST = 1;
elseif SIDX(1)<length(cycles)
    cycST = SIDX(1);
    if (Splitters(1)-cycles(SIDX(1)))>cycleLEN(SIDX(1))*extent_portion   % shrink
        cycST = SIDX(1)+1;
    end
end

cycEN = SIDX(2);
if cycEN<length(cycles)
    if (cycles(SIDX(2)+1)-Splitters(2))>cycleLEN(SIDX(2))*extent_portion   % extend
        cycEN = SIDX(2)+1;
    end
end
cycST = min(cycST, cycEN);

end