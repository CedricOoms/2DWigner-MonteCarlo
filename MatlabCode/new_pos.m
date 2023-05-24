%Generates new random position for particle with position 'pos' of one
%particle, d_max is the max allowed shift in coordinates along an axis
function [out,dr]=new_pos(pos,d_max)
    r_new=[];
    dr_v=[]; %Displacement vector
    for i=1:length(pos)
        rnd=-1+2*rand(); %Random number between -1 and 1
        r_new(end+1)=pos(i)+rnd*d_max; %New i^th coordinate of the particle
        dr_v(end+1)=rnd*d_max;
    end
    out=r_new';
    dr=sqrt(sum(dr_v.^2)); %Length of the displacement
end