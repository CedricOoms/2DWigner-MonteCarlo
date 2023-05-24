function out=lindemann_index(cell_in,MC_steps,N)
    r_ij=zeros(N,N);
    r_ij_sqrd=zeros(N,N);
    local_index=[]; %Vector containing local Lindemann indices q_i
    %Time average:
    for i=1:MC_steps
        r_ij=r_ij+cell_in{1,i};
        r_ij_sqrd=r_ij_sqrd+(r_ij.^2);
    end
    r_ij_avg=r_ij./MC_steps; %Time-average of r_ij (<r_ij>)
    r_ij_avg_sqrd=r_ij_avg.^2; %<r_ij>^2
    r_ij_sqrd_avg=r_ij_sqrd./MC_steps; %Time-average of r_ij^2 (<r_ij^2>)
    %Calculation of local lindemann indices q_i:
    for i=1:N
        r_ij_avg_n=r_ij_avg(i,:); r_ij_avg_n(i)=[]; %Delete r_ii
        r_ij_avg_sqrd_n=r_ij_avg_sqrd(i,:); r_ij_avg_sqrd_n(i)=[]; %Delete r_ii
        r_ij_sqrd_avg_n=r_ij_sqrd_avg(i,:); r_ij_sqrd_avg_n(i)=[]; %Delete r_ii
        local_index(end+1)=((N-1)^(-1))*sum(sqrt(r_ij_sqrd_avg_n-r_ij_avg_sqrd_n)./r_ij_avg_n);
    end
    Lindemann_system=sum(local_index)/N; %Global Lindemann index (system average of Local Lindemann indices)
    out=Lindemann_system;
end