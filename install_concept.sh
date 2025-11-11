###########################
# Install and patch CLASS #
###########################
if [ "${class_install}" == "True" ] && [ "${class_installed}" == "False" ]; then
    init_install "CLASS"
    # Move the content of the current directory (the files of CLASS,
    # except .gitignore) to the CLASS installation directory.
    mkdir -p "${class_dir}"
    mv ./* "${class_dir}/"
    cd "${class_dir}"
    # Below we will do a lot of patching on the CLASS source code.
    # To aid us, we define the following functions.
    patch_class() {
        local filename="$1"
        local linenr="$2"
        action="$3"
        belonging="$4"
        new_lines="$5"
        wrap_in_comments="$6"
        # Determine file type ("c" or "python")
        extension="${filename##*.}"
        filetype=""
        if     [ "${extension}" == "c"   ] \
            || [ "${extension}" == "h"   ] \
            || [ "${extension}" == "cpp" ]; then
            filetype="c"
        elif   [ "${extension}" == "py"  ] \
            || [ "${extension}" == "pyx" ] \
            || [ "${extension}" == "pxd" ]; then
            filetype="python"
        fi
        # Add comments around inserted lines
        if [ "${wrap_in_comments}" != "False" ]; then
            if [ "${filetype}" == "c" ]; then
                new_lines=(""                      \
                    "/************************/"   \
                    "/* For use with CONCEPT */"   \
                    "/************************/"   \
                    "${new_lines[@]}"              \
                    "/**************************/" \
                    "/* ^For use with CONCEPT^ */" \
                    "/**************************/" \
                    ""                             \
                )
            elif [ "${filetype}" == "python" ]; then
                new_lines=(""                    \
                    "########################"   \
                    "# For use with CONCEPT #"   \
                    "########################"   \
                    "${new_lines[@]}"            \
                    "##########################" \
                    "# ^For use with CONCEPT^ #" \
                    "##########################" \
                    ""                           \
                )
            fi
        fi
        # Find indentation at linenr
        indentation=""
        if [ "${belonging}" == "below" ]; then
            n_lines=$(wc -l "${filename}" | awk '{print $1}')
            ((n_lines_down = n_lines - linenr + 1)) || :
            content="$(tail -n ${n_lines_down} "${filename}")"
        elif [ "${belonging}" == "above" ]; then
            ((n_lines_down = linenr - 1)) || :
            content="$(head -n ${n_lines_down} "${filename}" | tac)"
        fi
        local IFS=''
        while read -r line; do
            if [ -n "${line// }" ]; then
                line_unindented="$(echo "${line}" | awk '{gsub(/^ +/,"")} {print $0}')"
                ((indentation_size = ${#line} - ${#line_unindented})) || :
                if [ "${filetype}" == "c" ] && [ "${line_unindented:0:1}" == "}" ]; then
                    ((indentation_size += 2))
                fi
                indentation="$(printf "%${indentation_size}s")"
                break
            fi
        done <<< "${content}"
        # Construct string of indented lines from the new_lines array
        new_line_nr=0
        for new_line in "${new_lines[@]}"; do
            if [ -n "${new_line}" ]; then
                indentation_use="${indentation}"
            else
                indentation_use=""
            fi
            if [ ${new_line_nr} -eq 0 ]; then
                new_lines_str="${indentation_use}${new_line}"
                if [[ "${new_lines_str}" == " "* ]] || [[ "${new_lines_str}" == "\\n"* ]]; then
                    new_lines_str="\\${new_lines_str}"
                fi
            else
                if [ -n "${new_lines_str}" ]; then
                    new_lines_str="${new_lines_str}\\n${indentation_use}${new_line}"
                else
                    new_lines_str="\\\\n${indentation_use}${new_line}"
                fi
            fi
            ((new_line_nr += 1))
        done
        # Insert the new lines in the file
        if [ "${action}" == "insert" ]; then
            sed -i "${linenr}i${new_lines_str}" "${filename}"
        elif [ "${action}" == "replace" ]; then
            sed -i "${linenr}d" "${filename}"
            sed -i "${linenr}i${new_lines_str}" "${filename}"
        fi
    }
    redefine_class() {
        local filename="$1"
        name="$2"
        value="$3"
        sed -i "s/^\(\s*#define\s\+${name}\s\+\)\(\S*\)\(.*\)$/\1${value}\3/" "${filename}"
    }
    echo "Patching CLASS"
    # Change the values of some preprocessing directives
    # in header files, allowing larger inputs and outputs.
    # Careful though, as values too large will not fit on the stack.
    redefine_class "${class_dir}/include/common.h"        _MAXTITLESTRINGLENGTH_   100000  # 10⁵
    redefine_class "${class_dir}/include/parser.h"        _LINE_LENGTH_MAX_         10000  # 10⁴
    redefine_class "${class_dir}/include/parser.h"        _ARGUMENT_LENGTH_MAX_     10000  # 10⁴
    redefine_class "${class_dir}/include/perturbations.h" _MAX_NUMBER_OF_K_FILES_  100000  # 10⁵
    # As only the (relatively) late-time evolution is needed from CLASS,
    # we hard-code the perturbations output to only be printed here.
    a_min="3e-4"
    pattern=' +a *= *pvecback'
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" \
        "${class_dir}/source/perturbations.c" | head -n 1)
    ((linenr += 1))
    new_lines=(                                  \
        "/* Only return output at late times */" \
        "double a_min = ${a_min};"               \
        "if (a < a_min)"                         \
        "  return _SUCCESS_;"                    \
    )
    patch_class "${class_dir}/source/perturbations.c" ${linenr} "insert" "below" "${new_lines[@]}"
    # When using the Runge-Kutta evolver, the derivatives are only
    # computed at the beginning of each time step, which is not
    # necessarily precise enough. Here we remove the derivative
    # computation in tools/evolver_rkck.c entirely and place it
    # in the perturb_print_variables function of
    # source/perturbations.c instead. This also ensures that
    # e.g. the ppw struct has been updated correctly
    # when it is time to print the perturbation results.
    pattern='( *x1 *== *x_ini *)'
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" \
        "${class_dir}/tools/evolver_rkck.c" | head -n 1)
    new_lines=(                                          \
        "/* derivs will be called in print_variables */" \
        "if (0 == 1) {  /* (x1 == x_ini) { */"           \
    )
    patch_class "${class_dir}/tools/evolver_rkck.c" ${linenr} "replace" "below" "${new_lines[@]}"
    pattern='double *\\* *dataptr *;'
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" \
        "${class_dir}/source/perturbations.c" | head -n 1)
    ((linenr += 1))
    new_lines=(                                                                  \
        "/**"                                                                    \
        " * Compute perturbation derivatives. This also ensures that the"        \
        " * ppw (and other) structs are up-to-date. This is important"           \
        " * when using the Runge-Kutta evolver, as this is otherwise"            \
        " * not taken care off correctly."                                       \
        " */"                                                                    \
        "class_call("                                                            \
        "  perturb_derivs(tau, y, dy, parameters_and_workspace, error_message)," \
        "  error_message,"                                                       \
        "  error_message);"                                                      \
    )
    patch_class "${class_dir}/source/perturbations.c" ${linenr} "insert" "below" "${new_lines[@]}"
    # Include ncdm Psi0[q] in perturbation output
    include_ncdm_Psi0="False"
    if [ "${include_ncdm_Psi0}" == "True" ]; then
        pattern=' *char +tmp *\\[ *40 *\\] *;'
        linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" \
            "${class_dir}/source/perturbations.c" | head -n 1)
        new_lines=(           \
            "char tmp[1024];" \
            "int index_q;"    \
        )
        patch_class "${class_dir}/source/perturbations.c" ${linenr} \
            "replace" "below" "${new_lines[@]}"
        pattern='sprintf *\\( *tmp *, *\"cs2_ncdm\\[%d\\]\" *, *n_ncdm *\\) *;'
        linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" \
            "${class_dir}/source/perturbations.c" | head -n 1)
        ((linenr += 2))
        new_lines=(                                                                   \
            "/* Include ncdm Psi0[q] in perturbation output */"                       \
            "for (index_q=0; index_q<pba->q_size_ncdm[n_ncdm]; index_q++) {"          \
            "  sprintf(tmp,\"Psi0[%d](%.16f)\",n_ncdm,pba->q_ncdm[n_ncdm][index_q]);" \
            "  class_store_columntitle(ppt->scalar_titles,tmp,_TRUE_);"               \
            "}"                                                                       \
        )
        patch_class "${class_dir}/source/perturbations.c" ${linenr} \
            "insert" "below" "${new_lines[@]}"
        pattern='class_store_double *\\( *dataptr *, \
*delta_p_over_delta_rho_ncdm *\\[ *n_ncdm *\\] *,'
        linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" \
            "${class_dir}/source/perturbations.c" | head -n 1)
        ((linenr += 1))
        new_lines=(                                                            \
            "/* Include ncdm Psi0[q] in perturbation output */"                \
            "if (ppw->approx[ppw->index_ap_ncdmfa] == (int)ncdmfa_on) {"       \
            "  for (index_q=0; index_q<pba->q_size_ncdm[n_ncdm]; index_q++) {" \
            "    class_store_double(dataptr, 0.0, _TRUE_, storeidx);"          \
            "  }"                                                              \
            "}"                                                                \
            "else {"                                                           \
            "  idx = ppw->pv->index_pt_psi0_ncdm1;"                            \
            "  for (index_q=0; index_q<pba->q_size_ncdm[n_ncdm]; index_q++) {" \
            "    class_store_double(dataptr, y[idx], _TRUE_, storeidx);"       \
            "    /* Jump to next momentum bin */"                              \
            "    idx += (ppw->pv->l_max_ncdm[n_ncdm]+1);"                      \
            "  }"                                                              \
            "}"                                                                \
        )
        patch_class "${class_dir}/source/perturbations.c" ${linenr} \
            "insert" "below" "${new_lines[@]}"
    fi
    # Correctly implement the fld pressure perturbation,
    # both with and without PPF.
    pattern='double *rho_plus_p_theta_fld *;'
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" \
        "${class_dir}/include/perturbations.h" | head -n 1)
    ((linenr += 1))
    new_lines=(                                                                               \
"double delta_p_fld;  /**< pressure perturbation of fluid, very non-trivial in PPF scheme */" \
    )
    patch_class "${class_dir}/include/perturbations.h" ${linenr} \
        "insert" "above" "${new_lines[@]}"
    pattern='/\\* *fluid *contribution *\\*/'
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" \
        "${class_dir}/source/perturbations.c" | head -n 1)
    ((linenr += 1))
    new_lines=(                                                                                  \
        "/**"                                                                                    \
        " * Count up total pressure and conformal time derivative of pressure,"                  \
        " * excluding the fld species. These are used for the PPF formalism of fld."             \
        " */"                                                                                    \
        "double p_tot = 0.;"                                                                     \
        "double p_tot_prime = 0.;"                                                               \
        "if (pba->has_fld == _TRUE_ && pba->use_ppf == _TRUE_) {"                                \
        "  /* Photons */"                                                                        \
        "  p_tot += 1./3.*ppw->pvecback[pba->index_bg_rho_g];"                                   \
        "  p_tot_prime += -3.*a_prime_over_a*(1. + 1./3.)*1./3."                                 \
        "    *ppw->pvecback[pba->index_bg_rho_g];"                                               \
        "  /* Baryons have no pressure */"                                                       \
        "  /* Ultra relativistic species */"                                                     \
        "  if (pba->has_ur == _TRUE_) {"                                                         \
        "    p_tot += 1./3.*ppw->pvecback[pba->index_bg_rho_ur];"                                \
        "    p_tot_prime += -3.*a_prime_over_a*(1. + 1./3.)*1./3."                               \
        "      *ppw->pvecback[pba->index_bg_rho_ur];"                                            \
        "  }"                                                                                    \
        "  /* Cold dark matter has no pressure */"                                               \
        "  /* Non-cold dark matter */"                                                           \
        "  if (pba->has_ncdm == _TRUE_) {"                                                       \
        "    for(n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++) {"                                  \
        "      p_tot += ppw->pvecback[pba->index_bg_p_ncdm1 + n_ncdm];"                          \
        "      p_tot_prime += -a_prime_over_a*(5.*ppw->pvecback[pba->index_bg_p_ncdm1 + n_ncdm]" \
        "        - ppw->pvecback[pba->index_bg_pseudo_p_ncdm1 + n_ncdm]);"                       \
        "    }"                                                                                  \
        "  }"                                                                                    \
        "  /* Decaying cold dark matter has no pressure */"                                      \
        "  /* Decay radiation */"                                                                \
        "  if (pba->has_dr == _TRUE_) {"                                                         \
        "    p_tot += 1./3.*ppw->pvecback[pba->index_bg_rho_dr];"                                \
        "    p_tot_prime += -3.*a_prime_over_a*(1. + 1./3.)*1./3."                               \
        "      *ppw->pvecback[pba->index_bg_rho_dr]"                                             \
        "      + 1./3.*a*pba->Gamma_dcdm*ppw->pvecback[pba->index_bg_rho_dcdm];"                 \
        "  }"                                                                                    \
        "  /* Importantly, we skip the dark energy fluid */"                                     \
        "  /* Scalar field */"                                                                   \
        "  if (pba->has_scf == _TRUE_) {"                                                        \
        "    p_tot += ppw->pvecback[pba->index_bg_p_scf];"                                       \
        "    p_tot_prime += -a_prime_over_a/(a*a)*ppw->pvecback[pba->index_bg_phi_prime_scf]"    \
        "      *ppw->pvecback[pba->index_bg_phi_prime_scf]"                                      \
        "      - 2./3.*ppw->pvecback[pba->index_bg_dV_scf]"                                      \
        "        *ppw->pvecback[pba->index_bg_phi_prime_scf];"                                   \
        "  }"                                                                                    \
        "  /* Lambda has constant pressure */"                                                   \
        "}"                                                                                      \
    )
    patch_class "${class_dir}/source/perturbations.c" ${linenr} \
        "insert" "above" "${new_lines[@]}"
    for n in $(grep -n 'class_call *( *background_w_fld *(' "${class_dir}/source/perturbations.c" \
        | awk '{print $1}'); do
        n="${n//:}"
        if [ ${n} -gt ${linenr} ]; then
            ((linenr = n + 1))
            new_lines=(                                                 \
                "double w_prime_fld = dw_over_da_fld*a_prime_over_a*a;" \
            )
            patch_class "${class_dir}/source/perturbations.c" ${linenr} \
                "insert" "above" "${new_lines[@]}"
            break
        fi
    done
    for n in $(grep -n 'ppw *-> *rho_plus_p_theta_fld *=' "${class_dir}/source/perturbations.c" \
        | awk '{print $1}'); do
        n="${n//:}"
        if [ ${n} -gt ${linenr} ]; then
            ((linenr = n + 1))
            new_lines=(                                                                          \
                "/* Pressure perturbation of fld without PPF */"                                 \
                "double ca2_fld = w_fld - w_prime_fld/(3.*a_prime_over_a*(1. + w_fld));"         \
                "ppw->delta_p_fld = pba->cs2_fld*ppw->delta_rho_fld"                             \
                "  + (pba->cs2_fld - ca2_fld)*(3.*a_prime_over_a*ppw->rho_plus_p_theta_fld/k2);" \
            )
            patch_class "${class_dir}/source/perturbations.c" ${linenr} \
                "insert" "above" "${new_lines[@]}"
            break
        fi
    done
    pattern='s2sq *= ppw *->'
    linenr_1=$(awk "\$0 ~ \"${pattern}\" {print NR}" \
        "${class_dir}/source/perturbations.c" | head -n 1)
    ((linenr_1 += 1))
    pattern='ppw *-> *S_fld *='
    linenr_2=$(awk "\$0 ~ \"${pattern}\" {print NR}" \
        "${class_dir}/source/perturbations.c" | head -n 1)
    ((linenr_2 -= 1)) || :
    sed -i "${linenr_1},${linenr_2}d" "${class_dir}/source/perturbations.c"
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" \
        "${class_dir}/source/perturbations.c" | head -n 1)
    new_lines=(                                                                                  \
        "double alpha_prime, X, Y, Z, X_prime, Y_prime, Z_prime;"                                \
        "double rho_plus_p_theta_fld_prime, metric_euler;"                                       \
        "double rho_t, rho_t_prime, p_t, p_t_prime, rho_fld, rho_fld_prime, p_fld, p_fld_prime;" \
        "double H, H_prime;"                                                                     \
        "double theta_t,theta_t_prime, S, S_prime;"                                              \
        "if (ppt->gauge == synchronous) {"                                                       \
        "  alpha = (y[ppw->pv->index_pt_eta] + 1.5*a2/k2/s2sq*(ppw->delta_rho"                   \
        "    + 3.*a_prime_over_a/k2*ppw->rho_plus_p_theta)"                                      \
        "    - y[ppw->pv->index_pt_Gamma_fld])/a_prime_over_a;"                                  \
        "  alpha_prime = -2.*a_prime_over_a*alpha + y[ppw->pv->index_pt_eta]"                    \
        "    - 4.5*(a2/k2)*ppw->rho_plus_p_shear;"                                               \
        "  metric_euler = 0.;"                                                                   \
        "} else {"                                                                               \
        "  alpha = 0.;"                                                                          \
        "  alpha_prime = 0.;"                                                                    \
        "  metric_euler = k2*y[ppw->pv->index_pt_phi] - 4.5*a2*ppw->rho_plus_p_shear;"           \
        "}"                                                                                      \
    )
    patch_class "${class_dir}/source/perturbations.c" ${linenr} \
        "insert" "above" "${new_lines[@]}"
    pattern='ppw->delta_rho_fld *='
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" \
        "${class_dir}/source/perturbations.c" | tail -n 1)
    ((linenr += 1))
    new_lines=(                                                                                 \
        "rho_t = rho_plus_p_tot - p_tot;"                                                       \
        "p_t = p_tot;"                                                                          \
        "rho_t_prime = -3.*a_prime_over_a*(rho_t + p_t);"                                       \
        "p_t_prime = p_tot_prime;"                                                              \
        "rho_fld = ppw->pvecback[pba->index_bg_rho_fld];"                                       \
        "p_fld = w_fld*rho_fld;"                                                                \
        "rho_fld_prime = -3.*a_prime_over_a*(rho_fld + p_fld);"                                 \
        "p_fld_prime = w_prime_fld*rho_fld - 3.*a_prime_over_a*(1. + w_fld)*p_fld;"             \
        ""                                                                                      \
        "H = ppw->pvecback[pba->index_bg_H];"                                                   \
        "H_prime = ppw->pvecback[pba->index_bg_H_prime];"                                       \
        "X = c_gamma_k_H_square;"                                                               \
        "X_prime = -2.*X*(a_prime_over_a + H_prime/H);"                                         \
        "Y = 4.5*a2/k2/s2sq*(rho_t + p_t);"                                                     \
        "Y_prime = Y*(2.*a_prime_over_a + (rho_t_prime + p_t_prime)/(rho_t + p_t));"            \
        "Z = 2./3.*k2*H/a;"                                                                     \
        "Z_prime = Z*(H_prime/H - a_prime_over_a);"                                             \
        ""                                                                                      \
        "theta_t = ppw->rho_plus_p_theta/rho_plus_p_tot;"                                       \
        "theta_t_prime = -a_prime_over_a*theta_t + (-p_t_prime*theta_t + k2*ppw->delta_p"       \
        "  - k2*ppw->rho_plus_p_shear)/rho_plus_p_tot+metric_euler;"                            \
        ""                                                                                      \
        "S = ppw->S_fld;"                                                                       \
        "S_prime = -Z_prime/Z*S + 1./Z*(rho_fld_prime + p_fld_prime)*(theta_t + k2*alpha)"      \
        "  + 1./Z*(rho_fld + p_fld)*(theta_t_prime + k2*alpha_prime);"                          \
        "rho_plus_p_theta_fld_prime = Z_prime*(S - 1./(1. + Y)*(S/(1. + 1./X)"                  \
        "  + y[ppw->pv->index_pt_Gamma_fld]*X))"                                                \
        "  + Z*(S_prime + Y_prime/(1. + Y*Y + 2.*Y)*(S/(1. + 1./X)"                             \
        "    + y[ppw->pv->index_pt_Gamma_fld]*X)"                                               \
        "    - 1./(1. + Y)*(S_prime/(1. + 1./X) + S*X_prime/(1. + X*X + 2.*X)"                  \
        "      + ppw->Gamma_prime_fld*X + y[ppw->pv->index_pt_Gamma_fld]*X_prime))"             \
        "  - k2*alpha_prime*(rho_fld + p_fld) - k2*alpha*(rho_fld_prime + p_fld_prime);"        \
        ""                                                                                      \
        "ppw->delta_p_fld = (rho_plus_p_theta_fld_prime"                                        \
        "  + 4.*a_prime_over_a*ppw->rho_plus_p_theta_fld - (rho_fld + p_fld)*metric_euler)/k2;" \
    )
    patch_class "${class_dir}/source/perturbations.c" ${linenr} \
        "insert" "above" "${new_lines[@]}"
    pattern='ppw *-> *delta_p *\\+= *pba *-> *cs2_fld *\\* *ppw *-> *delta_rho_fld *;'
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" \
        "${class_dir}/source/perturbations.c" | head -n 1)
    new_lines=(                             \
        "ppw->delta_p += ppw->delta_p_fld;" \
    )
    patch_class "${class_dir}/source/perturbations.c" ${linenr} \
        "replace" "above" "${new_lines[@]}"
    # For the PPF scheme, include a maximum value for (c_Γ*k*a/H)² above
    # which Γ = Γ' = 0, for the sake of numerical stability. In CAMB
    # such a maximum value is present as well and is set to 30. I have
    # found that a value of 10³ or even 10⁴ ensures stability as well,
    # while perturbing the fld δ_fld, θ_fld and δp_fld solutions
    # significantly less. Instead of a discontinuities cutoff as
    # in CAMB, we implement a smooth transition to zero, as defects in
    # δ_fld, θ_fld and δp_fld have been observed otherwise.
    pattern='c_gamma_k_H_square *='
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" \
        "${class_dir}/source/perturbations.c" | head -n 1)
    sed -i "${linenr}d" "${class_dir}/source/perturbations.c"
    pattern='s2sq *= ppw *->'
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" \
        "${class_dir}/source/perturbations.c" | head -n 1)
    ((linenr += 1))
    new_lines=(                                                                                   \
        "/**"                                                                                     \
        " * The computation of Gamma_fld and Gamma_prime_fld becomes unstable"                    \
        " * at large c_Gamma*k/H. To stabilise the system we set these to zero"                   \
        " * at some large c_Gamma*k/(aH)."                                                        \
        " * As to not introduce discontinuities, we have a smooth transition"                     \
        " * phase between the untouched values and completely nullified values."                  \
        " * This transition is given the shape of an error function in"                           \
        " * log(c_Gamma*k/(aH)) space. The parameters c_gamma_k_H_square_max_{0|1}"               \
        " * specify the borders of the transition."                                               \
        " * Here we nullify/shrink Gamma_fld only."                                               \
        " */"                                                                                     \
        "double Gamma_fld, Gamma_weight, Gamma_weight_steepness;"                                 \
        "double c_gamma_k_H_square_max_0, c_gamma_k_H_square_max_1;"                              \
        "c_gamma_k_H_square_max_0 = 1e+3;"                                                        \
        "c_gamma_k_H_square_max_1 = 1e+4;"                                                        \
        "c_gamma_k_H_square = pow(pba->c_gamma_over_c_fld*k/a_prime_over_a, 2)*pba->cs2_fld;"     \
        "if (c_gamma_k_H_square > c_gamma_k_H_square_max_1){"                                     \
        "    Gamma_fld = 0.;"                                                                     \
        "} else {"                                                                                \
        "  Gamma_fld = y[ppw->pv->index_pt_Gamma_fld];"                                           \
        "  if (c_gamma_k_H_square > c_gamma_k_H_square_max_0){"                                   \
        "    Gamma_weight_steepness = 5.; /* 5 results in double precision perfect transition */" \
        "    Gamma_weight = 0.5*(erf(Gamma_weight_steepness*("                                    \
        "      0.5*(log(c_gamma_k_H_square_max_0) + log(c_gamma_k_H_square_max_1))"               \
        "      - log(c_gamma_k_H_square)"                                                         \
        "    )) + 1.);"                                                                           \
        "    Gamma_fld *= Gamma_weight;"                                                          \
        "  }"                                                                                     \
        "}"                                                                                       \
    )
    patch_class "${class_dir}/source/perturbations.c" ${linenr} \
        "insert" "above" "${new_lines[@]}"
    pattern='double *alpha_prime *, *X *, *Y *, *Z'
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" \
        "${class_dir}/source/perturbations.c" | head -n 1)
    replacements_todo=7
    replacements=0
    for n in $(grep -n 'y *\[ *ppw *-> *pv *-> *index_pt_Gamma_fld *\]' \
        "${class_dir}/source/perturbations.c" | awk '{print $1}'); do
        n="${n//:}"
        if [ ${n} -ge ${linenr} ]; then
            sed -i "${n}s/y *\[ *ppw *-> *pv *-> *index_pt_Gamma_fld *\]/Gamma_fld/" \
                "${class_dir}/source/perturbations.c"
            ((replacements += 1))
            if [ ${replacements} -eq ${replacements_todo} ]; then
                break
            fi
        fi
    done
    pattern='ppw->Gamma_prime_fld *='
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" \
        "${class_dir}/source/perturbations.c" | head -n 1)
    new_lines=(                                                                         \
        "/* Nullify/shrink Gamma_prime_fld as done for Gamma_fld above */"              \
        "if (c_gamma_k_H_square > c_gamma_k_H_square_max_1){"                           \
        "    ppw->Gamma_prime_fld = 0.;"                                                \
        "} else {"                                                                      \
        "  ppw->Gamma_prime_fld = a_prime_over_a*(ppw->S_fld/(1. + c_gamma_k_H_square)" \
        "    - (1. + c_gamma_k_H_square)*Gamma_fld);"                                   \
        "  if (c_gamma_k_H_square > c_gamma_k_H_square_max_0){"                         \
        "      ppw->Gamma_prime_fld *= Gamma_weight;"                                   \
        "  }"                                                                           \
        "}"                                                                             \
    )
    patch_class "${class_dir}/source/perturbations.c" ${linenr} \
        "replace" "above" "${new_lines[@]}"
    # Include fld in perturbation output
    new_class_perturbation_linenr() {
        if [ "${1}" == "perturb_prepare_output" ]; then
            pattern='ppt *-> *number_of_scalar_titles *='
            linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" \
                "${class_dir}/source/perturbations.c" | head -n 1)
            ((linenr -= 1)) || :
        elif [ "${1}" == "perturb_print_variables" ]; then
            pattern='class_store_double *\\( *dataptr *, *theta_scf *, *pba *-> *has_scf *, \
*storeidx *\\)*;'
            linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" \
                "${class_dir}/source/perturbations.c" | head -n 1)
            n_lines_total=$(wc -l "${class_dir}/source/perturbations.c" | awk '{print $1}')
            ((n_lines_down = n_lines_total - linenr)) || :
            local IFS=''
            while read -r line; do
                ((linenr += 1))
                if [[ "${line}" == "  }" ]]; then
                    break
                fi
            done <<< "$(tail -n ${n_lines_down} "${class_dir}/source/perturbations.c")"
        fi
        echo ${linenr}
    }
    linenr=$(new_class_perturbation_linenr "perturb_prepare_output")
    new_lines=(                                                                     \
        "/* Include fld in perturbation output */"                                  \
        "class_store_columntitle(ppt->scalar_titles, \"delta_fld\", pba->has_fld);" \
        "class_store_columntitle(ppt->scalar_titles, \"theta_fld\", pba->has_fld);" \
        "/**"                                                                       \
        " * We choose to store cs2_fld = delta_p_fld/delta_rho_fld rather than"     \
        " * simply delta_p_fld itself, as is done for massive neutrinos."           \
        " */"                                                                       \
        "class_store_columntitle(ppt->scalar_titles, \"cs2_fld\", pba->has_fld);"   \
    )
    patch_class "${class_dir}/source/perturbations.c" ${linenr} "insert" "below" "${new_lines[@]}"
    linenr=$(new_class_perturbation_linenr "perturb_print_variables")
    new_lines=(                                                                             \
        "/* Include fld in perturbation output */"                                          \
        "double w_fld, dw_over_da_fld, integral_fld, theta_fld;"                            \
        "if (pba->has_fld) {"                                                               \
        "  class_call(background_w_fld(pba, a, &w_fld, &dw_over_da_fld, &integral_fld),"    \
        "    pba->error_message, ppt->error_message);"                                      \
        "  class_store_double(dataptr, ppw->delta_rho_fld/pvecback[pba->index_bg_rho_fld]," \
        "    pba->has_fld, storeidx);"                                                      \
        "  /* For w_fld = -1 (Lambda), we have theta = 0 */"                                \
        "  if (w_fld == -1.) {"                                                             \
        "    theta_fld = 0.;"                                                               \
        "  }"                                                                               \
        "  else {"                                                                          \
        "    theta_fld = ppw->rho_plus_p_theta_fld/"                                        \
        "      ((1. + w_fld)*pvecback[pba->index_bg_rho_fld]);"                             \
        "  }"                                                                               \
        "  class_store_double(dataptr, theta_fld, pba->has_fld, storeidx);"                 \
        "  /**"                                                                             \
        "   * We choose to store cs2_fld = delta_p_fld/delta_rho_fld rather than"           \
        "   * simply delta_p_fld itself, as is done for massive neutrinos."                 \
        "   *"                                                                              \
        "   */"                                                                             \
        "  class_store_double(dataptr,"                                                     \
        "    ppw->delta_p_fld/ppw->delta_rho_fld, pba->has_fld, storeidx);"                 \
        "}"                                                                                 \
    )
    patch_class "${class_dir}/source/perturbations.c" ${linenr} "insert" "below" "${new_lines[@]}"
    # Use proper initial conditions for the growth factor D(a) and rate
    # f(a), implement second-order versions D2(a) and f2(a) as well as
    # third-order versions D3a(a), D3b(a), D3c(a), and f3a(a), f3b(a),
    # f3c(a), and take dcdm and (conditionally) ncdm into account.
    filename="${class_dir}/python/cclassy.pxd"
    pattern='int +index_bg_f'
    new_lines=(                                  \
        "# Second-order growth factor and rate"  \
        "int index_bg_D2"                        \
        "int index_bg_f2"                        \
        "# Third-order growth factors and rates" \
        "int index_bg_D3a"                       \
        "int index_bg_f3a"                       \
        "int index_bg_D3b"                       \
        "int index_bg_f3b"                       \
        "int index_bg_D3c"                       \
        "int index_bg_f3c"                       \
    )
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" "${filename}" | head -n 1)
    patch_class "${filename}" $((linenr + 1)) "insert" "below" "${new_lines[@]}"
    filename="${class_dir}/include/background.h"
    pattern='int +index_bg_f *;'
    new_lines=(                                                                       \
        "int index_bg_D2;   /**< second-order growth factor D2(a) */"                 \
        "int index_bg_f2;   /**< second-order growth rate f2(a) = [dlnD2]/[dlna] */"  \
        "int index_bg_D3a;  /**< third-order growth factor D3a(a) */"                 \
        "int index_bg_f3a;  /**< third-order growth rate f3a(a) = [dlnD3a]/[dlna] */" \
        "int index_bg_D3b;  /**< third-order growth factor D3b(a) */"                 \
        "int index_bg_f3b;  /**< third-order growth rate f3b(a) = [dlnD3b]/[dlna] */" \
        "int index_bg_D3c;  /**< third-order growth factor D3c(a) */"                 \
        "int index_bg_f3c;  /**< third-order growth rate f3c(a) = [dlnD3c]/[dlna] */" \
    )
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" "${filename}" | head -n 1)
    patch_class "${filename}" $((linenr + 1)) "insert" "below" "${new_lines[@]}"
    pattern='int +index_bi_D_prime *;'
    new_lines=(                                                                                     \
        "int index_bi_D2;        /**< {C} second-order growth factor D2(a) */"                      \
        "int index_bi_D2_prime;  /**< {C} D2 satisfies \\\\f$ D2''(\\\\tau) = -aHD2'(\\\\tau) +     \
3/2 a^2 \\\\rho_M (D2(\\\\tau) + D^2(\\\\tau)) \\\\f$ */"                                           \
        "int index_bi_D3a;        /**< {C} third-order growth factor D3a(a) */"                     \
        "int index_bi_D3a_prime;  /**< {C} D3a satisfies \\\\f$ D3a''(\\\\tau) = -aHD3a'(\\\\tau) + \
3/2 a^2 \\\\rho_M (D3a(\\\\tau) + 2D^3(\\\\tau)) \\\\f$ */"                                         \
        "int index_bi_D3b;        /**< {C} third-order growth factor D3b(a) */"                     \
        "int index_bi_D3b_prime;  /**< {C} D3b satisfies \\\\f$ D3b''(\\\\tau) = -aHD3b'(\\\\tau) + \
3/2 a^2 \\\\rho_M (D3b(\\\\tau) + 2D(\\\\tau)D2(\\\\tau) + 2D^3(\\\\tau)) \\\\f$ */"                \
        "int index_bi_D3c;        /**< {C} third-order growth factor D3c(a) */"                     \
        "int index_bi_D3c_prime;  /**< {C} D3c satisfies \\\\f$ D3c''(\\\\tau) = -aHD3c'(\\\\tau) + \
3/2 a^2 \\\\rho_M D^3(\\\\tau) \\\\f$ */"                                                           \
    )
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" "${filename}" | head -n 1)
    patch_class "${filename}" $((linenr + 1)) "insert" "below" "${new_lines[@]}"
    pattern='double *\\* *deg_ncdm'
    new_lines=(                                                                                  \
        "double * growthfac_contrib_ncdm;  /**< ncdm contribution factors for growth factors */" \
    )
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" "${filename}" | head -n 1)
    patch_class "${filename}" ${linenr} "insert" "below" "${new_lines[@]}"
    filename="${class_dir}/source/input.c"
    pattern='class_read_list_of_doubles_or_default *\\( *\"deg_ncdm\"'
    new_lines=(                                                                 \
    	"/* Read growth factor contribution of each ncdm species: */"           \
    	"class_read_list_of_doubles_or_default("                                \
    	"  \"growthfac_contrib_ncdm\",pba->growthfac_contrib_ncdm,0.0,N_ncdm);" \
    )
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" "${filename}" | head -n 1)
    patch_class "${filename}" $((linenr - 1)) "insert" "below" "${new_lines[@]}"
    pattern='pba *-> *deg_ncdm_default *='
    new_lines=(                               \
    	"pba->growthfac_contrib_ncdm = NULL;" \
    )
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" "${filename}" | head -n 1)
    patch_class "${filename}" ${linenr} "insert" "below" "${new_lines[@]}"
    filename="${class_dir}/source/background.c"
    pattern='free *\\( *pba *-> *deg_ncdm *\\)'
    new_lines=(                               \
    	"free(pba->growthfac_contrib_ncdm);" \
    )
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" "${filename}" | head -n 1)
    patch_class "${filename}" ${linenr} "insert" "below" "${new_lines[@]}"
    pattern='class_define_index *\\( *pba *-> *index_bg_f *,'
    new_lines=(                                                       \
        "/* Second-order growth factor and rate */"                   \
        "class_define_index(pba->index_bg_D2, _TRUE_, index_bg, 1);"  \
        "class_define_index(pba->index_bg_f2, _TRUE_, index_bg, 1);"  \
        "/* Third-order growth factors and rates */"                  \
        "class_define_index(pba->index_bg_D3a, _TRUE_, index_bg, 1);" \
        "class_define_index(pba->index_bg_f3a, _TRUE_, index_bg, 1);" \
        "class_define_index(pba->index_bg_D3b, _TRUE_, index_bg, 1);" \
        "class_define_index(pba->index_bg_f3b, _TRUE_, index_bg, 1);" \
        "class_define_index(pba->index_bg_D3c, _TRUE_, index_bg, 1);" \
        "class_define_index(pba->index_bg_f3c, _TRUE_, index_bg, 1);" \
    )
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" "${filename}" | head -n 1)
    patch_class "${filename}" $((linenr + 1)) "insert" "below" "${new_lines[@]}"
    pattern='class_define_index *\\( *pba *-> *index_bi_D_prime *,'
    new_lines=(                                                             \
        "/* -> Second-order equation for second-order growth factor */"     \
        "class_define_index(pba->index_bi_D2, _TRUE_, index_bi, 1);"        \
        "class_define_index(pba->index_bi_D2_prime, _TRUE_, index_bi, 1);"  \
        "/* -> Third-order equations for third-order growth factors */"     \
        "class_define_index(pba->index_bi_D3a, _TRUE_, index_bi, 1);"       \
        "class_define_index(pba->index_bi_D3a_prime, _TRUE_, index_bi, 1);" \
        "class_define_index(pba->index_bi_D3b, _TRUE_, index_bi, 1);"       \
        "class_define_index(pba->index_bi_D3b_prime, _TRUE_, index_bi, 1);" \
        "class_define_index(pba->index_bi_D3c, _TRUE_, index_bi, 1);"       \
        "class_define_index(pba->index_bi_D3c_prime, _TRUE_, index_bi, 1);" \
    )
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" "${filename}" | head -n 1)
    patch_class "${filename}" $((linenr + 1)) "insert" "below" "${new_lines[@]}"
    pattern='memcopy_result *= *memcpy *\\( *pba *-> *background_table'
    new_lines=(                                                                                              \
        "/* Normalise D2 by the squared factor used for D and construct f2 = D2_prime/(aHD) */"              \
        "pvecback[pba->index_bg_D2] = pData[i*pba->bi_size + pba->index_bi_D2]/"                             \
        "  pow(pData[(pba->bt_size - 1)*pba->bi_size + pba->index_bi_D], 2);"                                \
        "pvecback[pba->index_bg_f2] = pData[i*pba->bi_size + pba->index_bi_D2_prime]/"                       \
        "  (pData[i*pba->bi_size + pba->index_bi_D2]*pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H]);"  \
        "/* Normalise D3a by the cubed factor used for D and construct f3a = D3a_prime/(aHD) */"             \
        "pvecback[pba->index_bg_D3a] = pData[i*pba->bi_size + pba->index_bi_D3a]/"                           \
        "  pow(pData[(pba->bt_size - 1)*pba->bi_size + pba->index_bi_D], 3);"                                \
        "pvecback[pba->index_bg_f3a] = pData[i*pba->bi_size + pba->index_bi_D3a_prime]/"                     \
        "  (pData[i*pba->bi_size + pba->index_bi_D3a]*pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H]);" \
        "/* Normalise D3b by the cubed factor used for D and construct f3b = D3b_prime/(aHD) */"             \
        "pvecback[pba->index_bg_D3b] = pData[i*pba->bi_size + pba->index_bi_D3b]/"                           \
        "  pow(pData[(pba->bt_size - 1)*pba->bi_size + pba->index_bi_D], 3);"                                \
        "pvecback[pba->index_bg_f3b] = pData[i*pba->bi_size + pba->index_bi_D3b_prime]/"                     \
        "  (pData[i*pba->bi_size + pba->index_bi_D3b]*pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H]);" \
        "/* Normalise D3c by the cubed factor used for D and construct f3c = D3c_prime/(aHD) */"             \
        "pvecback[pba->index_bg_D3c] = pData[i*pba->bi_size + pba->index_bi_D3c]/"                           \
        "  pow(pData[(pba->bt_size - 1)*pba->bi_size + pba->index_bi_D], 3);"                                \
        "pvecback[pba->index_bg_f3c] = pData[i*pba->bi_size + pba->index_bi_D3c_prime]/"                     \
        "  (pData[i*pba->bi_size + pba->index_bi_D3c]*pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H]);" \
    )
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" "${filename}" | head -n 1)
    patch_class "${filename}" $((linenr - 1)) "insert" "above" "${new_lines[@]}"
    pattern='pvecback_integration *\\[ *pba *-> *index_bi_D_prime *\\] *='
    new_lines=(                                                                                       \
        "/* Use proper initial conditions for the growth factors */"                                  \
        "double Omega0_M = pba->Omega0_b;"                                                            \
        "if (pba->has_cdm == _TRUE_)"                                                                 \
        "  Omega0_M += pba->Omega0_cdm;"                                                              \
        "if (pba->has_dcdm == _TRUE_)"                                                                \
        "  Omega0_M += pba->Omega_ini_dcdm;  /* take dcdm into account */"                            \
        "double Omega0_R_eff = pow(a/pba->a_today, 4)*pow(pvecback[pba->index_bg_H]/pba->H0, 2);  \
/* take all relativistic species into account */"                                                     \
        "double eps = 3./2.*Omega0_M/Omega0_R_eff*a;"                                                 \
        "double aH = a*pvecback[pba->index_bg_H];"                                                    \
        "double C = 1.0;  /* arbitrary */"                                                            \
        "pvecback_integration[pba->index_bi_D]         = pow(C, 1)   *(1. + 1.*eps + 1./4.*eps*eps);" \
        "pvecback_integration[pba->index_bi_D_prime]   = pow(C, 1)*aH*(0. + 1.*eps + 1./2.*eps*eps);" \
        "pvecback_integration[pba->index_bi_D2]        = pow(C, 2)   *(0. + 1.*eps + 3./4.*eps*eps);" \
        "pvecback_integration[pba->index_bi_D2_prime]  = pow(C, 2)*aH*(0. + 1.*eps + 3./2.*eps*eps);" \
        "pvecback_integration[pba->index_bi_D3a]       = pow(C, 3)   *(0. + 2.*eps + 2./1.*eps*eps);" \
        "pvecback_integration[pba->index_bi_D3a_prime] = pow(C, 3)*aH*(0. + 2.*eps + 4./1.*eps*eps);" \
        "pvecback_integration[pba->index_bi_D3b]       = pow(C, 3)   *(0. + 2.*eps + 5./2.*eps*eps);" \
        "pvecback_integration[pba->index_bi_D3b_prime] = pow(C, 3)*aH*(0. + 2.*eps + 5./1.*eps*eps);" \
        "pvecback_integration[pba->index_bi_D3c]       = pow(C, 3)   *(0. + 1.*eps + 3./4.*eps*eps);" \
        "pvecback_integration[pba->index_bi_D3c_prime] = pow(C, 3)*aH*(0. + 1.*eps + 3./2.*eps*eps);" \
    )
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" "${filename}" | head -n 1)
    patch_class "${filename}" $((linenr + 1)) "insert" "below" "${new_lines[@]}"
    pattern='class_store_columntitle *\\( *titles *, *\"gr\\.fac\\. f\"'
    new_lines=(                                                     \
        "class_store_columntitle(titles, \"gr.fac. D2\", _TRUE_);"  \
        "class_store_columntitle(titles, \"gr.fac. f2\", _TRUE_);"  \
        "class_store_columntitle(titles, \"gr.fac. D3a\", _TRUE_);" \
        "class_store_columntitle(titles, \"gr.fac. f3a\", _TRUE_);" \
        "class_store_columntitle(titles, \"gr.fac. D3b\", _TRUE_);" \
        "class_store_columntitle(titles, \"gr.fac. f3b\", _TRUE_);" \
        "class_store_columntitle(titles, \"gr.fac. D3c\", _TRUE_);" \
        "class_store_columntitle(titles, \"gr.fac. f3c\", _TRUE_);" \
    )
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" "${filename}" | head -n 1)
    patch_class "${filename}" $((linenr + 1)) "insert" "below" "${new_lines[@]}"
    pattern='class_store_double *\\( *dataptr *, *pvecback *\\[ *pba *-> *index_bg_f *\\]'
    new_lines=(                                                                       \
        "class_store_double(dataptr, pvecback[pba->index_bg_D2], _TRUE_, storeidx);"  \
        "class_store_double(dataptr, pvecback[pba->index_bg_f2], _TRUE_, storeidx);"  \
        "class_store_double(dataptr, pvecback[pba->index_bg_D3a], _TRUE_, storeidx);" \
        "class_store_double(dataptr, pvecback[pba->index_bg_f3a], _TRUE_, storeidx);" \
        "class_store_double(dataptr, pvecback[pba->index_bg_D3b], _TRUE_, storeidx);" \
        "class_store_double(dataptr, pvecback[pba->index_bg_f3b], _TRUE_, storeidx);" \
        "class_store_double(dataptr, pvecback[pba->index_bg_D3c], _TRUE_, storeidx);" \
        "class_store_double(dataptr, pvecback[pba->index_bg_f3c], _TRUE_, storeidx);" \
    )
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" "${filename}" | head -n 1)
    patch_class "${filename}" $((linenr + 1)) "insert" "below" "${new_lines[@]}"
    pattern='rho_M *\\+='
    new_lines=(                                                                    \
        "/* Include dcdm in growth factors */"                                     \
        "if (pba->has_dcdm == _TRUE_)"                                             \
        "  rho_M += pvecback[pba->index_bg_rho_dcdm];"                             \
        "/**"                                                                      \
        " * Code for including the non-relativistic contribution from ncdm in"     \
        " * growth factors. For small (realistic for neutrinos) ncdm masses,"      \
        " * ncdm will only cluster on large, linear scales. Including ncdm"        \
        " * when computing the scale-independent growth factors will then add"     \
        " * a correction to all scales, which should really only be applied"       \
        " * to large scales, greater than the free-streaming scale,"               \
        " * see (96) in https://arxiv.org/abs/astro-ph/0603494"                    \
        " */"                                                                      \
        "double rho_ncdm, p_ncdm;"                                                 \
        "int n_ncdm;"                                                              \
        "if (pba->has_ncdm == _TRUE_) {"                                           \
        "  for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++) {"                     \
	    "    if (pba->growthfac_contrib_ncdm[n_ncdm] == 0.) {"                     \
	    "      continue;"                                                          \
	    "    }"                                                                    \
        "    class_call(background_ncdm_momenta("                                  \
        "      pba->q_ncdm_bg[n_ncdm],"                                            \
        "      pba->w_ncdm_bg[n_ncdm],"                                            \
        "      pba->q_size_ncdm_bg[n_ncdm],"                                       \
        "      pba->M_ncdm[n_ncdm],"                                               \
        "      pba->factor_ncdm[n_ncdm],"                                          \
        "      1./a - 1.,"                                                         \
        "      NULL,"                                                              \
        "      &rho_ncdm,"                                                         \
        "      &p_ncdm,"                                                           \
        "      NULL,"                                                              \
        "      NULL),"                                                             \
        "      pba->error_message,"                                                \
        "      pba->error_message);"                                               \
        "    rho_M += pba->growthfac_contrib_ncdm[n_ncdm]*(rho_ncdm - 3.*p_ncdm);" \
        "  }"                                                                      \
        "}"                                                                        \
    )
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" "${filename}" | head -n 1)
    patch_class "${filename}" $((linenr + 1)) "insert" "below" "${new_lines[@]}"
    pattern='dy *\\[ *pba *-> *index_bi_D_prime *\\] *='
    new_lines=(                                                                           \
        "/* Second-order growth factor */"                                                \
        "dy[pba->index_bi_D2] = y[pba->index_bi_D2_prime];"                               \
        "dy[pba->index_bi_D2_prime] = -a*H*y[pba->index_bi_D2_prime]"                     \
        "  + 1.5*a*a*rho_M*(y[pba->index_bi_D2] + pow(y[pba->index_bi_D], 2));"           \
        "/* Third-order growth factors */"                                                \
        "dy[pba->index_bi_D3a] = y[pba->index_bi_D3a_prime];"                             \
        "dy[pba->index_bi_D3a_prime] = -a*H*y[pba->index_bi_D3a_prime]"                   \
        "  + 1.5*a*a*rho_M*(y[pba->index_bi_D3a] + 2.*pow(y[pba->index_bi_D], 3));"       \
        "dy[pba->index_bi_D3b] = y[pba->index_bi_D3b_prime];"                             \
        "dy[pba->index_bi_D3b_prime] = -a*H*y[pba->index_bi_D3b_prime]"                   \
        "  + 1.5*a*a*rho_M*(y[pba->index_bi_D3b] "                                        \
        "    + 2.*y[pba->index_bi_D]*y[pba->index_bi_D2] + 2.*pow(y[pba->index_bi_D], 3)" \
        "  );"                                                                            \
        "dy[pba->index_bi_D3c] = y[pba->index_bi_D3c_prime];"                             \
        "dy[pba->index_bi_D3c_prime] = -a*H*y[pba->index_bi_D3c_prime]"                   \
        "  + 1.5*a*a*rho_M*pow(y[pba->index_bi_D], 3);"                                   \
    )
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" "${filename}" | head -n 1)
    patch_class "${filename}" $((linenr + 1)) "insert" "below" "${new_lines[@]}"
    # Include theta_tot in perturbation output
    linenr=$(new_class_perturbation_linenr "perturb_prepare_output")
    new_lines=(                                                               \
        "/* Include theta_tot in perturbation output */"                      \
        "class_store_columntitle(ppt->scalar_titles, \"theta_tot\", _TRUE_);" \
    )
    patch_class "${class_dir}/source/perturbations.c" ${linenr} "insert" "below" "${new_lines[@]}"
    linenr=$(new_class_perturbation_linenr "perturb_print_variables")
    new_lines=(                                                                                  \
        "/* Include theta_tot in perturbation output */"                                         \
        "double rho_plus_p_tot = -2./3.*pvecback[pba->index_bg_H_prime]/a + 2./3.*pba->K/(a*a);" \
        "double theta_tot = ppw->rho_plus_p_theta/rho_plus_p_tot;"                               \
        "class_store_double(dataptr, theta_tot, _TRUE_, storeidx);"                              \
    )
    patch_class "${class_dir}/source/perturbations.c" ${linenr} "insert" "below" "${new_lines[@]}"
    # Include h_prime in perturbation output
    linenr=$(new_class_perturbation_linenr "perturb_prepare_output")
    new_lines=(                                                                                \
        "/* Include h_prime in perturbation output */"                                         \
        "class_store_columntitle(ppt->scalar_titles, \"h_prime\", ppt->gauge == synchronous);" \
    )
    patch_class "${class_dir}/source/perturbations.c" ${linenr} "insert" "below" "${new_lines[@]}"
    linenr=$(new_class_perturbation_linenr "perturb_print_variables")
    new_lines=(                                                          \
        "/* Include h_prime in perturbation output */"                   \
        "class_store_double(dataptr, pvecmetric[ppw->index_mt_h_prime]," \
        "  ppt->gauge == synchronous, storeidx);"                        \
    )
    patch_class "${class_dir}/source/perturbations.c" ${linenr} "insert" "below" "${new_lines[@]}"
    # Include H_T_prime in perturbation output
    linenr=$(new_class_perturbation_linenr "perturb_prepare_output")
    new_lines=(                                                               \
        "/* Include H_T_prime (in N-body gauge) in perturbation output */"    \
        "class_store_columntitle(ppt->scalar_titles, \"H_T_prime\", _TRUE_);" \
    )
    patch_class "${class_dir}/source/perturbations.c" ${linenr} "insert" "below" "${new_lines[@]}"
    linenr=$(new_class_perturbation_linenr "perturb_print_variables")
    new_lines=(                                                                             \
        "/**"                                                                               \
        " * Include H_T_prime (in N-body gauge) in perturbation output."                    \
        " * Here we make use of rho_plus_p_tot defined earlier."                            \
        " */"                                                                               \
        "double p_tot_prime = 0.0;"                                                         \
        "/* Photons */"                                                                     \
        " p_tot_prime += -3.*a*H*(1. + 1./3.)*1./3.*pvecback[pba->index_bg_rho_g];"         \
        "/* Baryons have no pressure */"                                                    \
        "/* Ultra relativistic species */"                                                  \
        "if (pba->has_ur == _TRUE_)"                                                        \
        "  p_tot_prime += -3.*a*H*(1. + 1./3.)*1./3.*pvecback[pba->index_bg_rho_ur];"       \
        "/* Cold dark matter has no pressure */"                                            \
        "/* Non-cold dark matter */"                                                        \
        "if (pba->has_ncdm == _TRUE_) {"                                                    \
        "  for(n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++)"                                   \
        "    p_tot_prime += -a*H*(5.*pvecback[pba->index_bg_p_ncdm1+n_ncdm]"                \
        "    - pvecback[pba->index_bg_pseudo_p_ncdm1+n_ncdm]);"                             \
        "}"                                                                                 \
        "/* Decaying cold dark matter has no pressure */"                                   \
        "/* Decay radiation */"                                                             \
        "if (pba->has_dr == _TRUE_)"                                                        \
        "  p_tot_prime += -3.*a*H*(1. + 1./3.)*1./3.*pvecback[pba->index_bg_rho_dr]"        \
        "    + 1./3.*a*pba->Gamma_dcdm*pvecback[pba->index_bg_rho_dcdm];"                   \
        "/* Dark energy fluid */"                                                           \
        "if (pba->has_fld == _TRUE_) {"                                                     \
        "  p_tot_prime += a*H*pvecback[pba->index_bg_rho_fld]"                              \
        "    *(a*dw_over_da_fld - 3.*w_fld*(1. + w_fld));"                                  \
        "}"                                                                                 \
        "/* Scalar field */"                                                                \
        "if (pba->has_scf == _TRUE_) {"                                                     \
        "  p_tot_prime += -H/a*pvecback[pba->index_bg_phi_prime_scf]"                       \
        "    *pvecback[pba->index_bg_phi_prime_scf]"                                        \
        "    - 2./3.*pvecback[pba->index_bg_dV_scf]*pvecback[pba->index_bg_phi_prime_scf];" \
        "}"                                                                                 \
        "/* Lambda has constant pressure */"                                                \
        "double H_T_prime = 3.*a*H/rho_plus_p_tot*("                                        \
        "  - ppw->delta_p"                                                                  \
        "  + p_tot_prime*theta_tot/(k*k)"                                                   \
        "  + ppw->rho_plus_p_shear);"                                                       \
        "class_store_double(dataptr, H_T_prime, _TRUE_, storeidx);"                         \
    )
    patch_class "${class_dir}/source/perturbations.c" ${linenr} "insert" "below" "${new_lines[@]}"
    # Do not convert synchronous variables to Newtonian gauge
    pattern='converting *synchronous *variables *to *newtonian *ones'
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" \
        "${class_dir}/source/perturbations.c" | head -n 1)
    ((linenr += 1))
    new_lines=(                                              \
        "/* Do not convert to Newtonian gauge */"            \
        "if (0 == 1) {  /* (ppt->gauge == synchronous) { */" \
    )
    patch_class "${class_dir}/source/perturbations.c" ${linenr} "replace" "below" "${new_lines[@]}"
    # Add 'node', 'num_threads' and 'message' to the classy.Class
    # initialiser as optional arguments and store them in the
    # 'background' struct. The 'num_threads' integer will hold the
    # number of MPI processes on the local node, which signals the
    # number of OpenMP threads to use. All three variables will be used
    # to make CLASS print out status updates during perturbation
    # computations.
    linenr_1="$(grep -n 'struct background' "${class_dir}/include/background.h" | head -n 1)"
    linenr_1="${linenr_1%%:*}"
    i=0
    while :; do
        ((i += 1))
        linenr_2="$(grep -n '};' "${class_dir}/include/background.h" | head -n ${i} | tail -n 1)"
        linenr_2="${linenr_2%%:*}"
        if [ ${linenr_2} -gt ${linenr_1} ]; then
            break
        fi
    done
    new_lines=(                                                \
        "/**"                                                  \
        " * Used to set number of OpenMP threads and to print" \
        " * status updates during perturbation computations."  \
        " */"                                                  \
        "int node, num_threads;"                               \
        "char* message;"                                       \
    )
    patch_class "${class_dir}/include/background.h" ${linenr_2} "insert" "below" "${new_lines[@]}"
    linenr="$(grep -n "__cinit__" "${class_dir}/python/classy.pyx")"
    linenr="${linenr%%:*}"
    sed -i "${linenr}s/.*/    def __cinit__(self, \
        default=False, concept_class_call=False, node=0, num_threads=-1, message=''): \
        # Changed for use with CONCEPT/" "${class_dir}/python/classy.pyx"
    ((linenr += 1))
    new_lines=(                                                            \
        "import os"                                                        \
        "os.environ.pop('CONCEPT_CLASS_CALL', None)"                       \
        "if concept_class_call:"                                           \
        "    os.environ['CONCEPT_CLASS_CALL'] = '1'"                       \
        "self.ba.node = <int>node"                                         \
        "self.ba.num_threads = <int>num_threads"                           \
        "self.ba.message = <char*>malloc((len(message) + 1)*sizeof(char))" \
        "strcpy(self.ba.message, message.encode())"                        \
    )
    patch_class "${class_dir}/python/classy.pyx" ${linenr} "insert" "below" "${new_lines[@]}"
    linenr_1="$(grep -n 'cdef struct background:' "${class_dir}/python/cclassy.pxd")"
    linenr_1="${linenr_1%%:*}"
    i=0
    while :; do
        ((i += 1))
        linenr_2="$(grep -n 'cdef struct' "${class_dir}/python/cclassy.pxd" \
            | head -n ${i} | tail -n 1)"
        linenr_2="${linenr_2%%:*}"
        if [ ${linenr_2} -eq ${linenr_1} ]; then
            ((i += 1))
            linenr_2="$(grep -n 'cdef struct' "${class_dir}/python/cclassy.pxd" \
                | head -n ${i} | tail -n 1)"
            linenr_2="${linenr_2%%:*}"
            while :; do
                ((linenr_2 -= 1)) || :
                line="$(sed "${linenr_2}!d" "${class_dir}/python/cclassy.pxd")"
                if [ -n "${line}" ]; then
                    break
                fi
            done
            ((linenr_2 += 1))
            break
        fi
    done
    new_lines=(           \
        "int node"        \
        "int num_threads" \
        "char* message"   \
    )
    patch_class "${class_dir}/python/cclassy.pxd" ${linenr_2} "insert" "above" "${new_lines[@]}"
    pattern='pba *-> *shooting_failed *= *_FALSE_ *;'
    new_lines=(
        "char * CONCEPT_CLASS_CALL = getenv(\"CONCEPT_CLASS_CALL\");" \
        "if (CONCEPT_CLASS_CALL == NULL) {"                           \
        "  pba->node = 0;"                                            \
        "  pba->num_threads = -1;"                                    \
        "  pba->message = (char*)malloc(1*sizeof(char));"             \
        "  pba->message[0] = '\\\\0';"                                \
        "}"                                                           \
    )
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" "${class_dir}/source/input.c" | head -n 1)
    patch_class "${class_dir}/source/input.c" $((linenr + 1)) "insert" "below" "${new_lines[@]}"
    linenr="$(grep -n 'for (index_k = ppt->k_size\[index_md\]-1; index_k >=0; index_k--)' \
        "${class_dir}/source/perturbations.c")"
    linenr="${linenr%%:*}"
    ((linenr += 1))
    new_lines=(                                                     \
        "if ((abort == _FALSE_) && (pba->message[0] != '\\\\0')) {" \
        "  printf("                                                 \
        "    pba->message,"                                         \
        "    pba->node,"                                            \
        "    thread,"                                               \
        "    ppt->k[index_md][index_k],"                            \
        "    ppt->k_size[index_md] - 1 - index_k,"                  \
        "    ppt->k_size[index_md] - 1"                             \
        "  );"                                                      \
        "  fflush(stdout);"                                         \
        "}"                                                         \
    )
    patch_class "${class_dir}/source/perturbations.c" ${linenr} "insert" "below" "${new_lines[@]}"
    pattern='#pragma *omp *parallel'
    linenr=$(awk "\$0 ~ \"${pattern}\" {print NR}" "${class_dir}/source/perturbations.c" \
        | head -n 1)
    new_lines=(                                                                  \
        "if (pba->num_threads != -1) {"                                          \
        "  /**"                                                                  \
        "   * Explicitly set the number of OpenMP threads."                      \
        "   * Note that the value of OMP_NUM_THREADS is now completely ignored." \
        "   */"                                                                  \
        "  omp_set_num_threads(pba->num_threads);"                               \
        "}"                                                                      \
    )
    patch_class "${class_dir}/source/perturbations.c" ${linenr} "insert" "above" "${new_lines[@]}"
    # If CO𝘕CEPT is installed we now switch out the values of various
    # constants in the CLASS source code so that they match the values
    # used in CO𝘕CEPT.
    if [ "${concept_works}" == "True" ]; then
        printf "Patching physical constants in CLASS to match values used in ${esc_concept}\n"
        class_constants=(                                            \
            _Mpc_over_m_   "units.Mpc/units.m"                       \
            _Gyr_over_Mpc_ "light_speed*units.Gyr/units.Mpc"         \
            _c_            "light_speed*units.s/units.m"             \
            _G_            "G_Newton*units.kg*units.s**2/units.m**3" \
            _eV_           "units.eV/units.J"                        \
            _k_B_          "NotImplemented"                          \
            _h_P_          "2*π*ħ/(units.J*units.s)"                 \
        )
        for ((name_index = 0; name_index < ${#class_constants[@]}; name_index += 2)); do
            ((expr_index = name_index + 1))
            name=${class_constants[${name_index}]}
            expr="${class_constants[${expr_index}]}"
            if [ "${expr}" == "NotImplemented" ]; then
                # This CLASS constant has no equivalent in CO𝘕CEPT
                continue
            fi
            value="$(concept_print "${expr}")"
            redefine_class "${class_dir}/include/background.h" ${name} "${value}"
        done
    fi
    # Build CLASS, including the Python wrapper classy. We explicitly
    # specify Python 3 as the language level and replace -O4 with -O3
    # (which are equivalent on compilers that recognise -O4).
    # We try with and without the -ffast-math option, try different
    # compilers and various flags for correctly linking to OpenMP.
    # The classy Python wrapper is hard-coded to use gcc in setup.py to
    # test for existence of the mvec library, as well as to get the
    # compiler library directory.
    sed -i '1s/^/# cython: language_level=3\n/' "python/classy.pyx"
    sed -i 's/-O4/-O3/' "Makefile"
    cp "Makefile" "Makefile_ori"
    cp "python/setup.py" "python/setup.py_ori"
    class_install_func() {
        if [ "${hardcoded_gcc}" == "False" ] && [ -z "${CC}" ]; then
            return 1
        fi
        if [ -d "${blas_dir}" ]; then
            export LD_LIBRARY_PATH="${blas_dir}/lib:${LD_LIBRARY_PATH}"
        fi
        export PATH="${python_dir}/bin:${PATH}"
        cp "Makefile_ori" "Makefile"
        cp "python/setup.py_ori" "python/setup.py"
        make clean || :
        if [ "${fast_math}" == "False" ]; then
            sed -i '0,/-ffast-math/s/-ffast-math/#-ffast-math/' "Makefile"
        fi
        if [ -n "${CC}" ]; then
            sed -i "s/CC * =/CC = ${CC//\//\\/}  #/" "Makefile"
        fi
        if [ -n "${OMPFLAG}" ]; then
            sed -i "s/OMPFLAG *=/OMPFLAG = ${OMPFLAG}  #/" "Makefile"
        fi
        sed -i "s/extra_link_args *=/extra_link_args=\
'${gomp_in_extra_link_args} ${OMPFLAG}'.split(),  #/" "python/setup.py"
        if [ "${hardcoded_gcc}" == "False" ]; then
            # Remove hard-coding of gcc
            linenr_1=$(awk "\$0 ~ \"gcc\" {print NR}" "python/setup.py" \
                | head -n 1)
            linenr_2=$(awk "\$0 ~ \"mvec\" {print NR}" "python/setup.py" \
                | tail -n 1)
            sed -i "${linenr_1},${linenr_2}d" "python/setup.py"
            sed -i "${linenr_1}iliblist = ['class']" "python/setup.py"
            sed -i 's/GCCPATH//g' "python/setup.py"
            # Check for the mvec library
            rm -rf "tmp_mvec" || :
            mkdir -p "tmp_mvec"
            cd "tmp_mvec"
            echo "int main(void){ return 0; }" > test.c
            mvec_warning="$("${CC}" -lmvec test.c 2>&1 | grep mvec || :)"
            cd ..
            rm -rf "tmp_mvec" || :
            if [ -z "${mvec_warning}" ]; then
                # mvec found.
                # Add the mvec and m library.
                ((linenr = linenr_1 + 1))
                sed -i "${linenr}iliblist += ['mvec', 'm']" "python/setup.py"
            fi
        fi
        PYTHON="${python}" make ${make_jobs} 2>&1 || return 1
        "${python}" -c "import classy; classy.Class().compute()" || return 1
        # Test CLASS
        if [ "${do_tests}" == "True" ]; then
            if [ -d "${blas_dir}" ]; then
                export LD_LIBRARY_PATH="${blas_dir}/lib:${LD_LIBRARY_PATH}"
            fi
            # We do not use the built-in test_class.py, as this requires
            # tens of gigabytes of memory and uses the Nose and
            # parameterized Python packages. Instead we perform a simple
            # test of our own, which only does a background computation.
            if ! "${python}" -c "
import sys
from classy import Class
cosmo = Class()
cosmo.compute()
sys.exit(int(cosmo.get_background()['proper time [Gyr]'][-1]) != 13)
"; then
                test_success="False"
            fi
            if [ "${test_success}" == "False" ]; then
                echo "CLASS did not pass a simple background computation test" > "test_log"
            fi
        fi
        # Cleanup if installing in slim mode
        if [ "${slim}" == "True" ]; then
            (cd "${class_dir}" && make clean) || :
            rm -rf "${class_dir}/doc" || :
            rm -rf "${class_dir}/output/"* || :
            rm -f "${class_dir}/Makefile_ori" "${class_dir}/python/setup.py_ori" || :
        fi
    }
    install "CLASS" "no_init_install"               \
        gomp_in_extra_link_args "" "-gomp"      ";" \
        hardcoded_gcc "False" "True"            ";" \
        compiler "${compiler_possibilities[@]}" ";" \
        OMPFLAG "-fopenmp" "-openmp" "-qopenmp" ";" \
        fast_math "True" "False"                ";" \

fi