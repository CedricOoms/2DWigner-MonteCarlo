function out=chi_sqrd(N,M)
    E_i=N/M;
    x=rand([1 N]); %N random numbers between 0 and 1
    H=histogram(x,M);
    y=H.BinCounts; %y_i
    out=sum(((y-E_i).^2)./E_i); %Chi_sqrd
end