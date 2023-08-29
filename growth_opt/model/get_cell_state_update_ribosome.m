function [newcellstate] = get_cell_state_update_ribosome(cell_state, par, mutant, snf1_vals, jgy_vals, mgl)

%Experiment run time - run for a long time to get steady-state conditions
t_batch_start = 0; 
t_batch_final = 300; 

%Set solver options and run the ode solver
Rel_tol  = 1.0E-03; 
Abs_tol  = 1.0E-06; 

options  = odeset('RelTol',Rel_tol, ...
                  'AbsTol',Abs_tol, ...
                  'NonNegative',(1:length(cell_state)));

if strcmp(mutant,'none') || ...
   strcmp(mutant,'const_gl') || ...
   strcmp(mutant,'const_eh') || ...
   strcmp(mutant,'const_jgy') || ...
   strcmp(mutant,'const_jgy_gl') || ...
   strcmp(mutant,'hxt') 

    mutant_type = 'const_gl_eh_aaex_cell'; % constant glucose, ethanol, constant cells 
    
elseif strcmp(mutant,'const_snf1_gl') ||...
       strcmp(mutant,'const_snf1') ||...
       strcmp(mutant,'const_snf1_eh')
    
    mutant_type = 'const_snf1_gl_eh_aaex_cell'; % constant snf1, glucose, ethanol, constant cells 

elseif strcmp(mutant,'const_tor_gl') || ...
       strcmp(mutant,'const_tor') ||...
       strcmp(mutant,'const_tor_eh')
    
    mutant_type = 'const_tor_gl_eh_aaex_cell'; % constant snf1, glucose, ethanol, constant cells 

end

[~,y_batch] = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant_type, snf1_vals, jgy_vals, mgl), [t_batch_start,t_batch_final], cell_state, options); 

% new cell state is the steady-state condition
newcellstate = y_batch(end,:);  
end
