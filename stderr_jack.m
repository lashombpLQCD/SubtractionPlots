%% Defining Jackknife Standard Error Estimate Function 

% Takes jackknife estimates and computes the jackknife standard error
% estimate for the estimator 
function [error_jack,ave_v_J] = stderr_jack(v_J, ndegs, nconfigs) 
    ave_v_J(:,1) = mean(v_J,2); % Average jackknife estimate 
    for i = 1:ndegs
        sq_diff = 0; 
        for j = 1:nconfigs 
            sq_diff = sq_diff + (v_J(i,j) -  ave_v_J(i,1))^2;
        end
        error_jack(i,1) = sqrt( (nconfigs - 1)/(nconfigs) * sq_diff ); % Jackknife Error as a function of degree
    end
end 