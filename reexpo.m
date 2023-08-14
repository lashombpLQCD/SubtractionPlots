%% Function that performs base-10 re-exponentiation and calculates top and 
% bottom for error bars of re-exponentiated value 

function [re_expo_rel_var_poly, top_poly, bottom_poly, re_expo_rel_var_hfpoly, top_hfpoly, bottom_hfpoly, re_expo_rel_var_hfes, top_hfes, bottom_hfes] = reexpo(log_rel_var_poly, log_rel_var_poly_err, log_rel_var_hfpoly, log_rel_var_hfpoly_err, log_rel_var_hfes, log_rel_var_hfes_err)

neigH = size(log_rel_var_hfpoly,1); 

% In case only POLY and HFPOLY are given and not HFES 
% Will need to be reworked to allow for more in the future (possibly with
% switch statements) 
if (nargin < 6) 
    log_rel_var_hfes = zeros(neigH,1,1);
    log_rel_var_hfes_err = zeros(neigH,1,1);
end
    

% POLY
top_poly = 10.^(log_rel_var_poly + log_rel_var_poly_err) - 10.^(log_rel_var_poly); 
bottom_poly = 10.^(log_rel_var_poly) - 10.^(log_rel_var_poly - log_rel_var_poly_err);
re_expo_rel_var_poly = 10.^(log_rel_var_poly); 

% HFPOLY 
for ieig = 1:neigH
    top_hfpoly(ieig,:,1) = 10.^(log_rel_var_hfpoly(ieig,:,1) + log_rel_var_hfpoly_err(ieig,:,1)) - 10.^(log_rel_var_hfpoly(ieig,:,1)); 
    bottom_hfpoly(ieig,:,1) = 10.^(log_rel_var_hfpoly(ieig,:,1)) - 10.^(log_rel_var_hfpoly(ieig,:,1) - log_rel_var_hfpoly_err(ieig,:,1));
    re_expo_rel_var_hfpoly(ieig,:,1) = 10.^(log_rel_var_hfpoly(ieig,:,1)); 
end 

% HFES 
for ieig = 1:neigH
    top_hfes(ieig,:,1) = 10.^(log_rel_var_hfes(ieig,:,1) + log_rel_var_hfes_err(ieig,:,1)) - 10.^(log_rel_var_hfes(ieig,:,1)); 
    bottom_hfes(ieig,:,1) = 10.^(log_rel_var_hfes(ieig,:,1)) - 10.^(log_rel_var_hfes(ieig,:,1) - log_rel_var_hfes_err(ieig,:,1));
    re_expo_rel_var_hfes(ieig,:,1) = 10.^(log_rel_var_hfes(ieig,:,1)); 
end 