%Checks the acceptance rate, and changes d_max accordingly
function out=d_max_updater(d_max,acceptance_rate,proc)
    if acceptance_rate < 0.5
        out=(1+proc)*d_max; %Increase of d_max by proc%
    elseif acceptance_rate > 0.5
        out=(1-proc)*d_max;
    else
        out=d_max;
    end
end