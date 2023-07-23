%Runs MC simulation of for a 2D Wigner crystal
%its = number of different initial configs
% N = number of particles
% d_max = max random displacement
% steps = number of MC steps
% X, Y = size of starting domain (rectangle) X x Y (units of r_0)
% T = temperature (units of T_0)
function Monte_Carlo_Heating(its,N,mc_steps,T_start,T_fin)
    dT = 10^-3;
    E_GS_it=[]; %Vector that stores end energy of the crystallized state (T=0K) every MC simulation (aka iteration)
    close all;
    mean_dR=[]; %Vector holding mean displacement for every simulation run
    for i=1:its
        d_max=.2; %initial max displacement
        %Initial positions (random):
        X=2;
        Y=2;
        x=-X+2*X.*rand([1 N]);
        y=-Y+2*Y.*rand([1 N]);
        config_init=[x ; y]; %Positions of all N particles
        %Start Monte-Carlo to obtain T=0K config:
        [config_T0,E_T0,~]=MC_Routine('Coulomb',N,mc_steps,d_max,config_init,0,0); %Configuration at T=0K
        E_GS_it(end+1)=E_T0;
        %Plot for result:
        f=figure();
        %Start Monte-Carlo for heated up sample, starting from the
        %crystallized configuration at T=0K.
        %Temperature is incrementally increased and after each increment the system is allowed to equilibrialize for 10^5 MC steps:
        config_T = config_T0;
        T = T_start;
        while T <= T_fin
            %Equilibrialize at the new temperature for 10^5 MC steps:
            T = T + dT;
            mc_steps_eq = 10^5;
            [config_T_eq,E_T,~]=MC_Routine('Coulomb',N,mc_steps_eq,d_max,config_T,T,1);
        end
        
        %T_fin has been reached and the system has equilibrialized at that temperature.
        %Particle trajectories are plotted for 10^3 MC steps at T = T_fin:
        mc_steps_T=10^3;
        [config_T,E_T,~]=MC_Routine(int_pot,N,mc_steps_T,d_max,config_T_eq,T,1);
	for k=1:N
	    r_k=[];
	    for j=0:(mc_steps_T-1)
	        index_step=k+j*N;
		r_k=[r_k config_T(:,index_step)];
	    end
	    p_traj=plot(r_k(1,:),r_k(2,:),'blue'); %Particle trajectories
	    hold on;
	end
	
        p_T=scatter(config_T(1,(mc_steps_T*N+1):end),config_T(2,(mc_steps_T*N+1):end),25,'red','displayname','Starting configuration (T=0)'); %Position at end of MC sim with applied heating
        hold on;
        p_T0=scatter(config_T0(1,:),config_T0(2,:),25,'filled','black','displayname','End configuration'); %Starting config at T=0K (before heating) 
        hold off;
        legend([p_traj p_T0 p_T],'Trajectories','Starting configuration (T=0)','End configuration')
        title_text=append(sprintf('T=%.3f $T_0$, applied for %d MC steps.',[T mc_steps_T]));
        title(title_text,'interpreter','latex')
        box on %apply edges around each plot 'box'
        axis image %Keep aspectratio (keep ratio when rescaling plot)
        saveas(f,append('Configuration ', num2str(i)));
        close(f);
    end
    ff=figure();
    scatter(1:1:length(E_GS_it),round(E_GS_it,4),'filled','d');
    title(['Energy per obtained ground state configuration for N = ',num2str(N)]);
    xlabel('Configuration');
    ylabel('E/N');
    curtick = get(gca, 'xTick');
    xticks(unique(round(curtick)));
    saveas(ff,append('EnergyPlot_N',num2str(N)));
end
