function  alphan = mandelZeros(f, rangeIn, nm)
    rangeDim = 10;
    alphan = [];
    count = 0;
    diff = 10;
    toll = 0.0001;
    while diff>toll && count<nm
        alphaTmp = findzeros(f, [rangeIn, rangeIn+rangeDim]);
        alphan = [alphan double(alphaTmp)];
        count = length(alphan);
        rangeIn = rangeIn+rangeDim;
        diff = abs(alphan(end)-alphan(end-1)-pi);
    end
   
    if count<nm
        alphan = [alphan(1:end-1) alphan(end):pi:alphan(end)+(nm-count)*pi];
    end
end
