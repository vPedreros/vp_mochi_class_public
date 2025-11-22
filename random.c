  /************************/
  /* For use with CONCEPT */
  /************************/

  /**************************/
  /* ^For use with CONCEPT^ */
  /**************************/
  //vp: Create perturbations for effDE Horndeski fluid
  // double delta_rho_smg = (1./3.) * (-pow(H,2)*M2*(kin + bra)*res*ppw->pvecmetric[ppw->index_mt_x_prime_smg] - 
  //                                   (3.*H*ppw->pvecback[pba->index_bg_H_prime]*M2*(2.-bra) + 
  //                                   (kin + 3.*bra)*a*pow(H,3)*M2 + 9.*a*H*(rho_tot+p_tot) +
  //                                    bra*k2*H*M2/a)*res*ppw->pvecmetric[ppw->index_mt_x_smg]);

  // double delta_p_smg = H*M2/3./a*(run*ppw->pvecmetric[ppw->index_mt_h_prime]/3. + 
  //                                               bra*res*ppw->pvecmetric[ppw->index_mt_x_prime_prime_smg] + 
  //                                               k2*run*res*ppw->pvecmetric[ppw->index_mt_x_smg]) + 
  //                                               (rho_tot+p_tot)*(res*ppw->pvecmetric[ppw->index_mt_x_prime_smg] +a*H*res*ppw->pvecmetric[ppw->index_mt_x_smg]) +
  //                                               a*H*pvecback_derivs[pba->index_bg_p_tot_wo_smg]*res*ppw->pvecmetric[ppw->index_mt_x_smg]+
  //                                               (2.*ppw->pvecback[pba->index_bg_H_prime]/a+pow(H,2)*pvecback_derivs[pba->index_bg_braiding_smg] + 4.*pow(H,2)*bra+pow(H,2)*bra*run*M2/3.*res*ppw->pvecmetric[ppw->index_mt_x_prime_smg]) +
  //                                               (6.*H*ppw->pvecback[pba->index_bg_H_prime] + 2.*ppw->pvecback[pba->index_bg_H_prime_prime] + a*pow(H,3)*pvecback_derivs[pba->index_bg_braiding_smg] + a*pow(H,3)*bra*(4.+run)+H*ppw->pvecback[pba->index_bg_H_prime]*(run+bra))*res*ppw->pvecmetric[ppw->index_mt_x_smg];

  // double delta_rho_plus_p_theta_smg = k2*(H/a*res*ppw->pvecmetric[ppw->index_mt_x_prime_smg]*M2/3. + M2/3.*(2.*ppw->pvecback[pba->index_bg_H_prime]/a+pow(H,2)*bra)*res*ppw->pvecmetric[ppw->index_mt_x_smg] + (rho_tot+p_tot)*res*ppw->pvecmetric[ppw->index_mt_x_prime_smg]);

  // double delta_rho_plus_p_shear_smg = H*M2*run/9./a*(ppw->pvecmetric[ppw->index_mt_h_prime]+6.*ppw->pvecmetric[ppw->index_mt_eta_prime]-2.*k2*res*ppw->pvecmetric[ppw->index_mt_x_smg]);

  double a_min = 1.7e-02;
  if (a<a_min){
    ppw->pvecmetric[ppw->index_mt_delta_smg] = 0;
    ppw->pvecmetric[ppw->index_mt_theta_smg] = 0;
    ppw->pvecmetric[ppw->index_mt_shear_smg] = 0;
  }

  if (a>a_min){
    double delta_rho_smg = ppw->delta_rho * (1-M2)/M2 +3*H/2/a*(beh*k2*ppw->pvecmetric[ppw->index_mt_eta]/a/H
                            - c14*res*ppw->pvecmetric[ppw->index_mt_x_prime_smg]
                            - res/a/H*(c15*k2 + c16*pow(a*H,2))*ppw->pvecmetric[ppw->index_mt_x_smg]
                            +bra*ppw->pvecmetric[ppw->index_mt_h_prime]/4);

    double delta_rho_plus_p_theta_smg_ = 3*pow(a,2)/2/k2*(
                                        ppw->rho_plus_p_theta*(1-M2)/M2
                                        - 2*k2/3/a*res*c0*H*ppw->pvecmetric[ppw->index_mt_x_smg]
                                        -k2/3/a2*res*cB*ppw->pvecmetric[ppw->index_mt_x_prime_smg]
                                      );
    double delta_rho_plus_p_shear_smg_ = -ppw->rho_plus_p_shear*(1-M2)/M2 - 2*k2*ten*ppw->pvecmetric[ppw->index_mt_eta]/9/a2
                                      +2*k2/9/a*H*run*ppw->pvecmetric[ppw->index_mt_alpha]
                                      +2*k2/9/a2*c8*res*ppw->pvecmetric[ppw->index_mt_x_smg]
                                      -2*k2/9/a2/a/H*cH*res*ppw->pvecmetric[ppw->index_mt_x_prime_smg];
                                      
    ppw->pvecmetric[ppw->index_mt_delta_smg] += delta_rho_smg_/ppw->pvecback[pba->index_bg_rho_smg];
    ppw->pvecmetric[ppw->index_mt_theta_smg] += delta_rho_plus_p_theta_smg_/(ppw->pvecback[pba->index_bg_rho_smg] + ppw->pvecback[pba->index_bg_p_smg]);
    ppw->pvecmetric[ppw->index_mt_shear_smg] += delta_rho_plus_p_shear_smg_/(ppw->pvecback[pba->index_bg_rho_smg] + ppw->pvecback[pba->index_bg_p_smg]);
  }

  /**************************/
  /* ^For use with CONCEPT^ */
  /**************************/