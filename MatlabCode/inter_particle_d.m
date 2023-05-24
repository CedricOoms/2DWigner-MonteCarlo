%Calculates the interparticle distance r_{ij} between all particles
%out is a matrix r_{ij}
function out=inter_particle_d(pos)
    N=length(pos(1,:)); %Total number of particles in the system
    r_ij=zeros(N,N);
    for i=1:N
        r_i=pos(:,i);
        if i~=N
            for j=i+1:N
                r_j=pos(:,j);
                dr=r_i-r_j; %Vector
                d=sqrt(sum(dr.^2)); %Distance between particles i and j (|dr|)
                r_ij(i,j)=d; 
                r_ij(j,i)=d; %Because of the symmetry of the matrix
            end
        end
    end
    out=r_ij;
end