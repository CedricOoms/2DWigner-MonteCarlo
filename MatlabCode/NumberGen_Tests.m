%Tests on MATLAB's default random number generator (Mersenne twister
%algorithm)
clear all
N=10^6; %Amount of generated numbers
M=500; %Number of bins (BinWidth=1/M)
N_runs=10^3; %Number of times to calculate chi^2


chi_sqrdd=[]; %Counts how often chi_sqrd > M-2/3
for i=1:N_runs
    chi_sqrdd(end+1)=chi_sqrd(N,M);
end
average_chi_sqrd=sum(chi_sqrdd)/length(chi_sqrdd);
ref=M-2/3;
ref2=M;

diff=abs(average_chi_sqrd-ref);
fprintf('Difference between average value of chi^2 and M-2/3 is %d \n',diff);
diff2=abs(average_chi_sqrd-ref2);
fprintf('Difference between average value of chi^2 and M is %d \n',diff2);

%Auto-correlation test:
x=rand([1 N]);
autocorr(x)

