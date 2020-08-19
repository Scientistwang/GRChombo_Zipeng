/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "ADMFixedBGVars.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "simd.hpp"

class Potential
{
  public:
    struct params_t
    {
        double mass;
        double self_interaction;
    };

    // class params
    params_t m_params;

    // add alias for metric vars
    template <class data_t>
    using MetricVars = typename ADMFixedBGVars::template Vars<data_t>;
    //Zipeng edit
    // functions that specifies m-field distribution
        
    template <class data_t>
    double m_field(data_t x, double y, double z, double m_0) const
    {   
        //insert some meaningful m field functions here
        double m = m_0 * ( pow(x,-1.0) + pow(y,-1.0) + pow(z,-1.0) );
        return m;
    }   

    template <class data_t>
    double compute_partial_t_m (data_t x, double y, double z, double m_0) const
    {   
        //necesssary because Tensor only has 3 components?
        //need to write function that corresponds to the m_field above
        return 0;
    }   

    template <class data_t>
    void compute_partial_m (data_t x, double y, double z, double m_0, 
                            Tensor<1, data_t> &partial_m) const
    {   
        //need to write function that corresponds to the m_field above
        partial_m[0] = -1.0 * m_0 * pow(x,-2.0);
        partial_m[1] = -1.0 * m_0 * pow(y,-2.0);
        partial_m[3] = -1.0 * m_0 * pow(z,-2.0);
    }
    
  public:
    //! The constructor
    Potential(params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the proca field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &dVdA, data_t &dphidt,
		    	   Coordinates<data_t> coords, //Zipeng edit
                           const vars_t<data_t> &vars,
                           const vars_t<Tensor<1, data_t>> &d1,
                           const MetricVars<data_t> &metric_vars) const
    {
        // calculate full spatial christoffel symbols and gamma^ij
        using namespace TensorAlgebra;
        const auto gamma_UU = compute_inverse(metric_vars.gamma);
        const auto chris_phys =
            compute_christoffel(metric_vars.d1_gamma, gamma_UU);

        // for ease of reading
        double c4 = m_params.self_interaction;

	//zipeng-edit
	//coords dependent mass
	data_t coords_mass = m_field(coords.x, coords.y, coords.z, m_params.mass);

        // Here we are defining often used terms
        // DA[i][j] = D_i A_j
        Tensor<2, data_t> DA;
        FOR2(i, j)
        {
            DA[i][j] = d1.Avec[j][i];
            FOR1(k) { DA[i][j] += -chris_phys.ULL[k][i][j] * vars.Avec[k]; }
        }

        // DAScalar = D_i A^i
        data_t DA_scalar;
        DA_scalar = 0;
        FOR2(i, j) { DA_scalar += DA[i][j] * gamma_UU[i][j]; }

        // Xsquared = X^/mu X_/mu
        data_t Xsquared;
        Xsquared = -vars.phi * vars.phi;
        FOR2(i, j) { Xsquared += gamma_UU[i][j] * vars.Avec[j] * vars.Avec[i]; }

        // C = 1 + 4 c4 A^k A_k - 12 c4 phi^2
        data_t C = 1.0 - 12.0 * c4 * vars.phi * vars.phi;
        FOR2(i, j)
        {
            C += 4.0 * c4 * gamma_UU[i][j] * vars.Avec[j] * vars.Avec[i];
        }

        // dVdA = mu^2 ( 1 + 4 c4 (A^k A_k - phi^2))
        dVdA = pow(coords_mass, 2.0) * (1.0 + 4.0 * c4 * Xsquared);

        // dphidt - for now the whole thing is here since it depends mainly
        // on the form of the potential - except the advection term which is in
        // the ProcaField code
        dphidt = 0;
        FOR2(i, j)
        {
            dphidt += -gamma_UU[i][j] * vars.Avec[i] * metric_vars.d1_lapse[j];
        }
        // QUESTION: Should this be lapse * Z / C  or lapse * Z??
        dphidt += -metric_vars.lapse * vars.Z / C;
        FOR4(i, j, k, l)
        {
            dphidt += -8.0 * c4 * metric_vars.lapse / C * gamma_UU[i][k] *
                      gamma_UU[j][l] * vars.Avec[i] * vars.Avec[j] * DA[k][l];
        }
        dphidt += metric_vars.lapse / C * (1.0 + 4.0 * c4 * Xsquared) *
                  (metric_vars.K * vars.phi - DA_scalar);
        FOR1(i)
        {
            dphidt += 8.0 * c4 * vars.phi * metric_vars.lapse / C *
                      (vars.Evec[i] * vars.Avec[i]);
        }
        FOR4(i, j, k, l)
        {
            dphidt += 8.0 * c4 * vars.phi * metric_vars.lapse / C *
                      (-metric_vars.K_tensor[i][j] * vars.Avec[k] *
                       vars.Avec[l] * gamma_UU[i][k] * gamma_UU[j][l]);
        }
        FOR2(i, j)
        {
            dphidt += 8.0 * c4 * vars.phi * metric_vars.lapse / C *
                      (2.0 * vars.Avec[i] * d1.phi[j] * gamma_UU[i][j]);
        }

        //Zipeng Edit
        //questionable things: metric_vars.shift 

        Tensor<1,data_t> partial_m;
        compute_partial_m(coords.x, coords.y, coords.z, m_params.mass, partial_m);
        double partial_t_m = compute_partial_t_m (coords.x, coords.y, coords.z, m_params.mass);

        data_t DmA;
        // DmA = partial_mu m * A^mu
        DmA = 0;
        DmA += partial_t_m * vars.phi;
        FOR1(i)
        {
            DmA += partial_m[i] * vars.Avec[i];
        }

        data_t D_m_beta;
        // D_m_beta = partial_i m * beta^i
        D_m_beta = 0;
        FOR1(i)
        {
            D_m_beta += partial_m[i] * metric_vars.shift[i];
        }

	dphidt += -2.0 / m_field(coords.x, coords.y, coords.z, m_params.mass) 
		  * (metric_vars.lapse * DmA + partial_t_m * vars.phi - D_m_beta * vars.phi);
    }
};

#endif /* POTENTIAL_HPP_ */
