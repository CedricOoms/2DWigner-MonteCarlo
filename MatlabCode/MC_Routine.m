%Monte Carlo simulation
%If track = 1 -> track particle motion each MC step and calculate Lindemann index (useful for heating)
%E_out = E/N for the final configuration (N=number of particles)
%int_pot (str) = interaction potential -> 'Coulomb' or 'LJ'
function [config_out,E_out,mean_displacement]=MC_Routine(int_pot,N,mc_steps,d_max,pos_in,T,track)
    dr_mean=zeros(1,N);
    updater_i=0;
    if track == 1
        pos_track=pos_in;
    end
    pos=pos_in;
    E=Energy_total(pos_in,int_pot);
    for i=1:mc_steps
        updater_i=updater_i+1;
        accepted=0;
        %Loop over all particles:
        for j=1:N
            %New position for particle j
            pos2=pos;
            if T==0
                [pos2(:,j),~]=new_pos(pos(:,j),d_max);
            else
                [pos2(:,j),dr]=new_pos(pos(:,j),d_max);
            end
            %Energy with this new position:
            E_new=Energy_total(pos2,int_pot);
            %If lower, then keep new position:
            dE=E_new-E; %Energy difference, if negative the new position yields a lower energy configuration
            if dE < 0
                pos=pos2;
                E=E_new;
                accepted=accepted+1;
            elseif (T ~= 0) && (dr < exp(-dE/T))
                pos=pos2;
                E=E_new;
                accepted=accepted+1;
            end
        end
        if updater_i==100
            acceptance_rate=accepted/100;
            dmax=d_max_updater(d_max,acceptance_rate,.2);
            updater_i=0;
            accepted=0;
            acceptance_rate=0;
        end
        if track==1
            pos_track=[pos_track pos]; %Keep configuration for every MC step, in order to track particle movement
            dR=pos_in-pos; 
            dr_vec=[]; %Vector containing distance between new position and equilibrium position
            for k=1:length(dR(1,:))
                dr_vec(end+1)=sum(dR(:,k).^2);
            end
            dr_mean=dr_mean+dr_vec;
        end
    end
    %Output:
    E_out=Energy_total(pos,int_pot)/N;
    dr_mean=dr_mean./mc_steps;
    mean_displacement=sum(dr_mean)/N; %Mean distance of particles wrt their equilibrium positions, averaged over all particles
    if track==1
        config_out=pos_track; %Contains config at each MC step (2 x mc_step*N) array
    else
        config_out=pos; %Just final configuration of MC simulation
    end
end