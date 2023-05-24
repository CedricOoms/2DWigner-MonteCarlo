%Script to plot mean displacements as a function of temperature
clear
N=30;
T_start=0.01;
T_end=0.025;
dT=0.0001;
T=T_start:dT:T_end;
average_dR_T=[]; %Average mean displacement wrt equilibrium positions as a function of temperature
for T_i=T
    dR=monte_carlo2D('Coulomb',5,N,1000,T_i);
    average_dR=sum(dR)/length(dR); %Average of the mean displacement wrt equilibrium position over all 'equivalent' simulation runs (same temperature)
    average_dR_T(end+1)=average_dR;
end
close all;
figure();
plot(T,average_dR_T);
xlabel('T [T_0]')
ylabel('<dR> [r_0]')
yline(0.5)
