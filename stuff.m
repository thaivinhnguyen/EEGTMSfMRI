function x=stuff(x)
    % dont ask
    if isnan(x(1)), x(1)=nanmean(x); end;
    for i=2:length(x),
        if isnan(x(i)), x(i)=x(i-1); end
    end
    x = medfilt1(x-x(1),10);

