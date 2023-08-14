%% Base-10 Log of Log-averaged Relative Variance and Jackknife Error 
% Function that takes STANDARD ERROR from NS, POLY, and HFPOLY and returns 
% the base-10 logged relative variance, and the jackknife error for the 
% base-10 logged relative variance. (i.e. the exponent of the relative 
% variance and the error of the exponent). 
% The variances of each method are averaged using base-10 log-averaging and
% then re-exponentiated to compute the relative variance. Then the base-10
% log is performed to get the exponent of this relative variance

% If you feed it variance instead of standard error, then comment out the
% lines, 
% ns_scalar(ideg,i) = ( ns_scalar(ideg,i) * sqrt(nrhs) ).^2;, 
% etc. as they are calculating variance from the standard error


function [log_rel_var_poly, log_rel_var_hfpoly, log_rel_var_hfes, log_rel_var_poly_err,log_rel_var_hfpoly_err,log_rel_var_hfes_err] = log_rel_var_calc(ns_scalar,poly_scalar,hfpoly_scalar,nrhs,hfes_scalar)

    ndeg = size(poly_scalar,1);
    nconfigs = size(poly_scalar,2); 
    neigH = size(hfpoly_scalar,1); 
    
    % In case you only use NS, POLY, and HFPOLY and don't include HFES. 
    % Needs to be reworked if you want more in the future. -PL 
    if (nargin < 5) 
        hfes_scalar = zeros(neigH,ndeg,nconfigs); 
    end 
    
    
    % Inputs are Standard Error from Matlab program 
    ns_scalar_SE = ns_scalar; 
    poly_scalar_SE = poly_scalar; 
    hfpoly_scalar_SE = hfpoly_scalar;
    hfes_scalar_SE = hfes_scalar; 
    
    
    
    for ideg = 1:ndeg
        %% NS, POLY, HFPOLY, & HFES Variance Estimates 
        %%%%% Computing Variance Estimates from Standard Error Estimates
        for i = 1:nconfigs
            
            % NS
            ns_scalar(ideg,i) = ( ns_scalar_SE(ideg,i) * sqrt(nrhs) ).^2;
            
            % POLY 
            poly_scalar(ideg,i) = ( poly_scalar_SE(ideg,i) * sqrt(nrhs) ).^2;
            
            % HFPOLY 
            for ieig = 1:neigH
                hfpoly_scalar(ieig,ideg,i) = ( hfpoly_scalar_SE(ieig,ideg,i) * sqrt(nrhs) ).^2; 
            end 
            
            % HFES
            for ieig = 1:neigH
                hfes_scalar(ieig,ideg,i) = ( hfes_scalar_SE(ieig,ideg,i) * sqrt(nrhs) ).^2; 
            end 
        end
        
        
        %% For NS 
        %%%%% Perform log-averaging for NS variance 
        log_ave_ns_var(1,1)   = log10mean(ns_scalar(ideg,:),1,nconfigs); 
        
        
        %% For POLY
        %%%%% Perform log-averaging for POLY variance 
        log_ave_poly_var(ideg,1) = log10mean(poly_scalar(ideg,:),1,nconfigs); 


        %% For HFPOLY  
        %%%%% Perform log-averaging for HFPOLY variance 
        for ieig = 1:neigH
            log_ave_hfpoly_var(ieig,ideg,1) = log10mean(hfpoly_scalar(ieig,ideg,:),1,nconfigs); 
        end

        
        %% For HFES  
        %%%%% Perform log-averaging for HFES variance 
        for ieig = 1:neigH
            log_ave_hfes_var(ieig,ideg,1) = log10mean(hfes_scalar(ieig,ideg,:),1,nconfigs); 
        end

        %% Relative Log-averaged Variances

        % POLY 
        rel_var_poly(ideg,1) = 10.^(log_ave_poly_var(ideg,1)) / 10.^(log_ave_ns_var(1,1)); 

        % HFPOLY 
        for ieig = 1:neigH
            rel_var_hfpoly(ieig,ideg,1) = 10.^(log_ave_hfpoly_var(ieig,ideg,1)) / 10.^(log_ave_ns_var(1,1)); 
        end

        % HFES 
        for ieig = 1:neigH
            rel_var_hfes(ieig,ideg,1) = 10.^(log_ave_hfes_var(ieig,ideg,1)) / 10.^(log_ave_ns_var(1,1)); 
        end

        %% Exponent of Relative Log-averaged Variances 

        % POLY 
        log_rel_var_poly(ideg,1) = log10(rel_var_poly(ideg,1));
        
        % HFPOLY 
        for ieig = 1:neigH
            log_rel_var_hfpoly(ieig,ideg,1) = log10(rel_var_hfpoly(ieig,ideg,1));     
        end
        
        % HFES 
        for ieig = 1:neigH
            log_rel_var_hfes(ieig,ideg,1) = log10(rel_var_hfes(ieig,ideg,1));     
        end    

    end 




    %% %%%%% JACKKNIFE ERROR ANALYSIS %%%%% 


    %% For POLY 
    %%%%% Perform log-averaging for POLY and NS variance  
    for j = 1:nconfigs
        log_ave_poly_var_J(:,j) = log10mean(poly_scalar(:,[1:j-1, j+1:nconfigs]),ndeg,nconfigs-1); 
        log_ave_ns_var_J(1,j)   = log10mean(ns_scalar(:,[1:j-1, j+1:nconfigs]),1,nconfigs-1); 
    end 


    %% For HFPOLY  
    %%%%% Perform log-averaging for HFPOLY variance 
    for j = 1:nconfigs 
        for ieig = 1:neigH
            temp(1:ndeg,1:nconfigs) = 0; 
            temp(1:ndeg,1:nconfigs-1) = hfpoly_scalar(ieig,:,[1:j-1, j+1:nconfigs]); 
            log_ave_hfpoly_var_J(ieig,:,j) = log10mean(temp(:,:),ndeg,nconfigs-1); 
        end
    end   
    
    
    %% For HFES 
    %%%%% Perform log-averaging for HFPOLY variance 
    for j = 1:nconfigs 
        for ieig = 1:neigH
            temp(1:ndeg,1:nconfigs) = 0; 
            temp(1:ndeg,1:nconfigs-1) = hfes_scalar(ieig,:,[1:j-1, j+1:nconfigs]); 
            log_ave_hfes_var_J(ieig,:,j) = log10mean(temp(:,:),ndeg,nconfigs-1); 
        end
    end   

    %% Relative Log-averaged Variances

    for j = 1:nconfigs
        
        % POLY 
        rel_var_poly_J(:,j) = 10.^(log_ave_poly_var_J(:,j)) / 10.^(log_ave_ns_var_J(1,j)); 

        % HFPOLY 
        for ieig = 1:neigH
            rel_var_hfpoly_J(ieig,:,j) = 10.^(log_ave_hfpoly_var_J(ieig,:,j)) / 10.^(log_ave_ns_var_J(1,j)); 
        end
        
        % HFES 
        for ieig = 1:neigH
            rel_var_hfes_J(ieig,:,j) = 10.^(log_ave_hfes_var_J(ieig,:,j)) / 10.^(log_ave_ns_var_J(1,j)); 
        end
    end 

    %% Exponent of Relative Log-averaged Variances 

    for j = 1:nconfigs
        
        % POLY 
        log_rel_var_poly_J(:,j) = log10(rel_var_poly_J(:,j));
        
        % HFPOLY 
        for ieig = 1:neigH
            log_rel_var_hfpoly_J(ieig,:,j) = log10(rel_var_hfpoly_J(ieig,:,j));     
        end    
        
        % HFES
        for ieig = 1:neigH
            log_rel_var_hfes_J(ieig,:,j) = log10(rel_var_hfes_J(ieig,:,j));     
        end    
    end 


    %% Jackknife error estimates 

    % Rel POLY 
    [log_rel_var_poly_err(:,1),~] = stderr_jack(log_rel_var_poly_J(:,:),ndeg,nconfigs); 

    % Rel HFPOLY 
    for ieig = 1:neigH
        temp(1:ndeg,1:nconfigs) = 0; 
        temp(1:ndeg,1:nconfigs) = log_rel_var_hfpoly_J(ieig,:,:); 
       [log_rel_var_hfpoly_err(ieig,:,1),~] = stderr_jack(temp(:,:),ndeg,nconfigs); 
    end  
    
    % Rel HFES
    for ieig = 1:neigH
        temp(1:ndeg,1:nconfigs) = 0; 
        temp(1:ndeg,1:nconfigs) = log_rel_var_hfes_J(ieig,:,:); 
       [log_rel_var_hfes_err(ieig,:,1),~] = stderr_jack(temp(:,:),ndeg,nconfigs); 
    end  
end
