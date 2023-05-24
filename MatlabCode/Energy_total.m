function out=Energy_total(pos,int_pot)
    N=length(pos(1,:)); %Number of particles
    E=0;
    if strcmp(int_pot,'Coulomb') %Coulomb interaction
        %Loop over all N particles:
        for i=1:N
            r_i=pos(:,i);
            E=E+sum(r_i.^2); %Confining potential
            if i~=N
                %Interaction energy:
                for j=i+1:N
                    r_j=pos(:,j);
                    dr=r_i-r_j; %Vector
                    d=sqrt(sum(dr.^2)); %Distance between particles i and j (|dr|)
                    E=E+(d^-1);
                end
            end
        end
    elseif strcmp(int_pot,'LJ') %LJ interaction potential
        %Loop over all N particles:
        for i=1:N
            r_i=pos(:,i);
            E=E+sum(r_i.^2); %Confining potential
            if i~=N
                %Interaction energy:
                for j=i+1:N
                    r_j=pos(:,j);
                    dr=r_i-r_j; %Vector
                    d=sqrt(sum(dr.^2)); %Distance between particles i and j (|dr|)
                    E=E+(d^(-12)+d^(-6));
                end
            end
        end
    else
        sprintf('Invalid option for the interaction potential: choose Coulomb or LJ')
    end
    out=E;
end