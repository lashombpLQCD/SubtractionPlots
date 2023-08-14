%% Defining Base-10 Log-Averaging Function 

% Takes a container v(1:ndegs,1:nconfig) containing the variance estimates 
% as a function of the polynomial degree and configurations and returns the
% base-10 log-average of the variance estimates 
function log_ave_var = log10mean(v, ndegs, nconfigs)
    for i = 1:ndegs 
        temp = 0; 
        for j = 1:nconfigs 
            temp = temp + log10( v(i,j) );
        end 
        temp = temp / nconfigs; 
        log_ave_var(i,1) = temp; 
    end 
end 