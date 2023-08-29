function [beta, alpha, rib, tRNA, eIF_a_s, eIF_a_tau] = module_gene_expression_update_ribosome(r0, met, prot, snf1, tor, par, mutant)

% useful functions
hill_transc_s_tau = @(s, t, ep, xi_s, xi_t, w_s, w_t, theta_s, theta_t)(ep + (xi_s .* ((s/w_s).^theta_s)) ...
                                                                           + (xi_t .* ((t/w_t).^theta_t))) ...
                                                                          ./ (1 + ((s/w_s).^theta_s) + ((t/w_t).^theta_t));                                                                                                                                 
                                                                          
hill_transc_s = @(s, ep, xi_s, w_s, theta_s)(ep + (xi_s .* ((s/w_s).^theta_s)) ...
                                         ./ (1 + ((s/w_s).^theta_s)));   
                                     
%% metabolites 
 
aa_in = met(5); 
ae    = met(6);

%% transcription factors 

alpha.r  = par.a0_r + (1 - par.a0_r)*hill_transc_s_tau(snf1, tor, par.ep_r, ...
                                                       par.xi_r_s, par.xi_r_tau, ...
                                                       par.w_r_s, par.w_r_tau, ...
                                                       par.theta_r_s, par.theta_r_tau);
alpha.z  = 1;
alpha.gy = 1;
alpha.fe = 1;
alpha.gn = par.a0_gn + (1 - par.a0_gn)*hill_transc_s(snf1, par.ep_gn, par.xi_gn_s, par.w_gn_s, par.theta_gn_s);
alpha.mt = par.a0_mt + (1 - par.a0_mt)*hill_transc_s(snf1, par.ep_mt, par.xi_mt_s, par.w_mt_s, par.theta_mt_s);                         
                               
alpha.as = par.a0_as + (1 - par.a0_as)*hill_transc_s_tau(snf1, tor, par.ep_as, ...
                                                         par.xi_as_s, par.xi_as_tau, ...
                                                         par.w_as_s, par.w_as_tau, ...
                                                         par.theta_as_s, par.theta_as_tau);  
                                         
alpha.at = 1; 

%% tRNA 
tRNA.t0 = par.mt  * r0;                                                                         % total
tRNA.tc = tRNA.t0 * hill(aa_in, par.K_tc_aa, par.n_tc_aa) * hill(ae, par.K_tc_ae, par.n_tc_ae); % charged
tRNA.tu = tRNA.t0 - tRNA.tc;                                                                    % uncharged

%% ribosome partitioning

% term for tRNA binding to ribosomes 
trna_tc = tRNA.tc./par.K_rib_tc;  
trna_tu = tRNA.tu./par.K_rib_tu;
trna_contri_tc = trna_tc./(1 + trna_tc + trna_tu); 
trna_contri_tu = trna_tu./(1 + trna_tc + trna_tu); 

eIF_a_s   = hill(par.D_rib_s, snf1, par.theta_rib_s);
eIF_a_tau = hill(tor, par.D_rib_tau, par.theta_rib_tau);

f_to_at   = par.k_tj .* table2array(struct2table(alpha))' .* trna_contri_tc .* hill(ae, par.K_rib_ae, par.n_rib_ae) * eIF_a_s * eIF_a_tau;
f_to_as   = par.k_tj .* table2array(struct2table(alpha))' .* trna_contri_tu .* hill(ae, par.K_rib_ae, par.n_rib_ae) * eIF_a_s * eIF_a_tau;

rib.rf    = r0./(1 + sum(f_to_at) + sum(f_to_as));
rib.rat_j = rib.rf .* f_to_at; 
rib.rat   = sum(rib.rat_j);
rib.ras_j = rib.rf .* f_to_as;
rib.ras   = sum(rib.ras_j);

fat = sum(f_to_at)/(1 + sum(f_to_at));

%% protein synthesis rate

beta = ((par.k_e .* rib.rat_j)./(par.l)); 

prot_syn_rate = sum(par.l .* beta); % aa units 
g_rate = prot_syn_rate./par.pro_den;

rib.others = [trna_tc trna_tu ...
              trna_contri_tc trna_contri_tu...
              hill(ae, par.K_rib_ae, par.n_rib_ae) sum(table2array(struct2table(alpha))')...
              eIF_a_s eIF_a_tau...
              rib.rf sum(f_to_at) rib.rat sum(f_to_as) rib.ras prot_syn_rate g_rate fat ...
              tRNA.t0  aa_in hill(aa_in, par.K_tc_aa, par.n_tc_aa)]; 
          
end 