// 18O2 and 16O2 model for Gordon Holtgrieve and Daniel Schindler
// 2008-12-09 ZTA
// using ADMB v7.1.1, with MinGW g++ 3.4.4 at AFSC and MS VC++ 6.0 at UW
//
// updated January 2014 by ZTA
// respiration was split into a 'base' pool and a 'production' pool, and coefficients
// were updated in the calculations for the Schmidt number and O2 concentration at saturation
// using ADMB revision 1561 (compiled from source), with cygwin and g++ 4.8.2 on Windows 7 64-bit


DATA_SECTION

  init_number altitude;         // altitude
  init_number aspect;           // aspect
  init_int DST_adj;             // 0 - no daylight savings time, 1 - daylight savings time
  init_number latitude;         // latitude
  init_number longitude;        // longitude
  init_number slope;            // slope (of what?)
  init_number solar_const;      // solar constant
  init_number GST_adj;          // should be a whole number or a whole number + 0.5
  init_number transmiss;        // transmissivity
  init_number area;             // area (in m^2)
  init_number depth;            // depth (in m)
  init_number salinity;         // salinity

  init_number alphaP;           // ???
  init_number alphaGk;          // ???
  init_int S_meas_flag;         // 0 - do not use measured S (irradiance?), 1 - use measured S
  init_int Light_sat_flag;      // 0 - daylight is not saturated (linear increase), 1 - daylight is saturated (asymptotic increase), 2 - daylight is not saturated (hyperbolic tangent)

  init_number est_S_scale;      // scaling factor to scale estimated S to measured S

  init_int RK_step_size;        // step size for Runge-Kutta algorithm; should be between 2 and 10

  // NOTE:  the temperature data is the longest data set and most complete of the measured data

  init_int n_irr;               // number of irradiance measurements
  init_int n_tC;                // number of temperature (degC) measurements
  init_int n_O2conc;            // number of O2 concentration measurements
  init_int n_delta_18O_O2;      // number of delta 18O_O2 measurements

  init_vector irr_data(1,n_irr);
  init_vector irr_day_of_year(1,n_irr);
  init_vector irr_hour_of_day(1,n_irr);
  init_vector irr_time_interval(1,n_irr);

  init_vector tC_data(1,n_tC);
  init_vector tC_day_of_year(1,n_tC);
  init_vector tC_hour_of_day(1,n_tC);
  init_vector tC_time_interval(1,n_tC);

  init_vector O2conc_data(1,n_O2conc);
  init_vector O2conc_day_of_year(1,n_O2conc);
  init_vector O2conc_hour_of_day(1,n_O2conc);
  init_vector O2conc_time_interval(1,n_O2conc);

  init_vector delta_18O_O2_data(1,n_delta_18O_O2);
  init_vector delta_18O_O2_day_of_year(1,n_delta_18O_O2);
  init_vector delta_18O_O2_hour_of_day(1,n_delta_18O_O2);
  init_vector delta_18O_O2_time_interval(1,n_delta_18O_O2);

  init_int end_of_data_value_1;                 // indicator of end of values

  // map time intervals to tC_time_interval indices
  int ii;
  int jj;
  ivector irr_time_interval_map(1,n_irr);
  ivector O2conc_time_interval_map(1,n_O2conc);
  ivector delta_18O_O2_time_interval_map(1,n_delta_18O_O2);


  // flags for estimated parameters
!! ad_comm::change_datafile_name("BaMM.cfg");

  init_number trefC;                            // reference temperature in degrees Celsius, e.g., 15 deg C
  init_int base_respir_switch;                  // -1: off, 1 or 2: on
  init_int Eb_switch;                           // -1: off, 1 or 2: on
  init_vector Eb_params(1,4);                   // 1: min, 2: max, 3: mean for prior, 4: std dev for prior

  init_int Ep_switch;                           // -1: off, 1 or 2: on
  init_vector Ep_params(1,4);                   // 1: min, 2: max, 3: mean for prior, 4: std dev for prior
  init_int beta_prod_switch;                    // -1: off, 1 or 2: on
  init_vector beta_prod_params(1,4);            // 1: min, 2: max, 3: mean for prior, 4: std dev for prior
  init_int N_steps_for_C_pool;                  // number of preceding time steps over which the labile carbon pool is accumulated

  init_int isotope_switch;                      // -1: off, 1: on

  init_int k20_switch;                          // -1: off, 1 or 2: on
  init_vector k20_params(1,4);                  // 1: min, 2: max, 3: mean for prior, 4: std dev for prior
  init_int alphaR_switch;                       // -1: off, 1 or 2: on
  init_int alphaR_transform_switch;             // -1: off, 1: on
  init_vector alphaR_params(1,4);               // 1: min, 2: max, 3: mean for prior, 4: std dev for prior
  init_int delta_18O_H2O_switch;                // -1: off, 1 or 2: on
  init_vector delta_18O_H2O_params(1,4);        // 1: min, 2: max, 3: mean for prior, 4: std dev for prior
  init_int init_O2_conc_switch;                 // -1: off, 1 or 2: on
  init_vector init_O2_conc_params(1,4);         // 1: min, 2: max, 3: mean for prior, 4: std dev for prior
  init_int init_delta_18O_O2_switch;            // -1: off, 1 or 2: on
  init_vector init_delta_18O_O2_params(1,4);    // 1: min, 2: max, 3: mean for prior, 4: std dev for prior
  init_int sigma_O2_conc_switch;                // -1: off, 1 or 2: on
  init_vector sigma_O2_conc_params(1,4);        // 1: min, 2: max, 3: mean for prior, 4: std dev for prior
  init_int sigma_delta_18O_O2_switch;           // -1: off, 1 or 2: on
  init_vector sigma_delta_18O_O2_params(1,4);   // 1: min, 2: max, 3: mean for prior, 4: std dev for prior

  init_number delta_18O_O2_scale_factor;        // scale factor for delta_18O_O2 NLL component

  init_int end_of_data_value_2;                 // indicator of end of values

  int prod_respir_switch;

  int base_respir_phase;
  int Eb_phase;
  int Ep_phase;
  int beta_prod_phase;
  int k20_phase;
  int alphaR_phase;
  int delta_18O_H2O_phase;
  int init_O2_conc_phase;
  int init_delta_18O_O2_phase;
  int sigma_O2_conc_phase;
  int sigma_delta_18O_O2_phase;

  number Eb_min;
  number Eb_max;
  number Ep_min;
  number Ep_max;
  number beta_prod_min;
  number beta_prod_max;
  number k20_min;
  number k20_max;
  number alphaR_min;
  number alphaR_max;
  number delta_18O_H2O_min;
  number delta_18O_H2O_max;
  number init_O2_conc_min;
  number init_O2_conc_max;
  number init_delta_18O_O2_min;
  number init_delta_18O_O2_max;
  number sigma_O2_conc_min;
  number sigma_O2_conc_max;
  number sigma_delta_18O_O2_min;
  number sigma_delta_18O_O2_max;


  number R_air;
  number R_SMOW;


  int mcmc_iter                     // number of MCMC iterations


  number pi;
  number deg_to_rad;                // conversion factor from degrees to radians
  number rad_to_deg;                // conversion factor from radians to degrees
  number degC_to_degK;              // conversion factor from degrees Celsius to degrees Kelvin

  number tref;                      // reference temperature, e.g., 15 degrees Celsius, in degrees Kelvin
  number kk;                        // Boltzmann constant, in electronvolts/degrees Kelvin

  number atmos_corr;                // atmospheric pressure correction


  int Light_sat_Pmax_phase;         // phase depending upon Light_sat_flag
  int Light_sat_I_halfSat_phase;    // phase depending upon Light_sat_flag
  int Light_sat_alpha_P_I_phase;    // phase depending upon Light_sat_flag


  int debug_flag;


  !!#define NUM_DIFFEQ      2       // number of differential equations

  !!#define IDX_MIN         1       // first element of parameter array is the min value
  !!#define IDX_MAX         2       // second element of parameter array is the max value
  !!#define IDX_MEAN        3       // third element of parameter array is the mean value for the normal prior
  !!#define IDX_STD_DEV     4       // fourth element of parameter array is the std dev for the normal prior


  !! cout << endl;
  !! cout << "Data checks";
  !! if (end_of_data_value_1 == -999 && end_of_data_value_2 == -999)
  !! {
  !!     cout << ": okay" << endl;
  !! }
  !! else
  !! {
  !!     cout << ": mismatch" << endl;
  !!     cout << "Latitude\t" << latitude << endl;
  !!     cout << "Longitude\t" << longitude << endl;
  !!     cout << "alphaR switch\t" << alphaR_switch << endl;
  !!     cout << "delta_18O_O2_data\t" << delta_18O_O2_data << endl;
  !!     cout << "Config checks" << endl;
  !!     cout << "base respir switch\t" << base_respir_switch << endl;
  !!     cout << "isotope switch\t" << isotope_switch << endl;
  !!     cout << "sigma_delta_18O_O2 params\t" << sigma_delta_18O_O2_params << endl;
  !!     cout << endl;
  !!     exit(1);
  !! }


 LOCAL_CALCS

    // map irradiance_time_intervals to tC_time_intervals by index
    irr_time_interval_map = 0;
    jj = 0;
    for (ii = 1; ii <= n_irr; ii++)
    {
        do
        {
            jj++;
        }
        while (jj <= n_tC && irr_time_interval(ii) != tC_time_interval(jj));

        if (jj <= n_tC && irr_time_interval(ii) == tC_time_interval(jj))
        {
            irr_time_interval_map(ii) = jj;
        }
    }

    // map O2conc_time_interval to tC_time_intervals by index
    O2conc_time_interval_map = 0;
    jj = 0;
    for (ii = 1; ii <= n_O2conc; ii++)
    {
        do
        {
            jj++;
        }
        while (jj <= n_tC && O2conc_time_interval(ii) != tC_time_interval(jj));

        if (jj <= n_tC && O2conc_time_interval(ii) == tC_time_interval(jj))
        {
            O2conc_time_interval_map(ii) = jj;
        }
    }

    // map delta_18O_O2_time_interval to tC_time_intervals by index
    delta_18O_O2_time_interval_map = 0;
    jj = 0;
    for (ii = 1; ii <= n_delta_18O_O2; ii++)
    {
        do
        {
            jj++;
        }
        while (jj <= n_tC && delta_18O_O2_time_interval(ii) != tC_time_interval(jj));

        if (jj <= n_tC && delta_18O_O2_time_interval(ii) == tC_time_interval(jj))
        {
            delta_18O_O2_time_interval_map(ii) = jj;
        }
    }


    // R_air        = 0.002052924;
    R_air        = 0.002052322;
    R_SMOW       = 0.0020052;


    pi           = 4.0 * atan(1.0);
    deg_to_rad   = pi / 180.0;
    rad_to_deg   = 180.0 / pi;
    degC_to_degK = 273.15;

    tref         = trefC + degC_to_degK;    // reference temperature for base respiration calculation
    kk           = 8.6173324e-5;            // Boltzmann constant, in eV/degK; see http://physics.nist.gov/cgi-bin/cuu/Value?tkev

    atmos_corr   = pow(((288.15 - (0.0065 * altitude)) / 288.15), 5.256);    // correction for atmospheric pressure


    if (est_S_scale <= 0.0)
    {
        est_S_scale = 1.0;
    }


    mcmc_iter    = 0;


    if (Light_sat_flag == 1)
    {
        // if light saturation flag == 1, then estimate Pmax and I_halfSat and not alpha_P_I
        Light_sat_Pmax_phase      = 1;
        Light_sat_I_halfSat_phase = 1;
        Light_sat_alpha_P_I_phase = -1;
    }
    else if (Light_sat_flag == 2)
    {
        // if light saturation flag == 2, then estimate alpha_P_I and Pmax and not I_halfSat
        Light_sat_Pmax_phase      = 1;
        Light_sat_I_halfSat_phase = -1;
        Light_sat_alpha_P_I_phase = 1;
    }
    else
    {
        // if light saturation flag == 0 or other, then estimate alpha_P_I and not Pmax and I_halfSat
        Light_sat_Pmax_phase      = -1;
        Light_sat_I_halfSat_phase = -1;
        Light_sat_alpha_P_I_phase = 1;
    }


    base_respir_phase = 1;
    if (base_respir_switch == -1)
    {
        base_respir_phase = -1;
    }

    Eb_phase = 2;
    if (Eb_switch == -1)
    {
        Eb_phase = -2;
    }

    k20_phase = 1;
    if (k20_switch == -1)
    {
        k20_phase = -1;
    }

    alphaR_phase = 3;
    if (alphaR_switch == -1)
    {
        alphaR_phase = -3;
    }

    delta_18O_H2O_phase = 1;
    if (delta_18O_H2O_switch == -1 || isotope_switch == -1)
    {
        delta_18O_H2O_phase = -1;
    }

    init_O2_conc_phase = 1;
    if (init_O2_conc_switch == -1)
    {
        init_O2_conc_phase = -1;
    }

    init_delta_18O_O2_phase = 1;
    if (init_delta_18O_O2_switch == -1 || isotope_switch == -1)
    {
        init_delta_18O_O2_phase = -1;
    }

    sigma_O2_conc_phase = 5;
    if (sigma_O2_conc_switch == -1)
    {
        sigma_O2_conc_phase = -5;
    }

    sigma_delta_18O_O2_phase = 5;
    if (sigma_delta_18O_O2_switch == -1 || isotope_switch == -1)
    {
        sigma_delta_18O_O2_phase = -5;
    }

    Ep_phase = 3;
    if (Ep_switch == -1)
    {
        Ep_phase = -3;
    }

    beta_prod_phase = 4;
    if (beta_prod_switch == -1)
    {
        beta_prod_phase = -4;
    }

    // if any of the production-dependent respiration parameters are on, then include prod respir in calculations
    prod_respir_switch = max(Ep_switch,beta_prod_switch);

    if (N_steps_for_C_pool < 1)
    {
        N_steps_for_C_pool = 20;    // default number of time steps, to cover 3 to 5 hours
    }

    debug_flag   = 0;


 END_CALCS


INITIALIZATION_SECTION

  // using a PIN file instead


PARAMETER_SECTION

  init_number log_Pmax(Light_sat_Pmax_phase);
  init_number log_I_halfSat(Light_sat_I_halfSat_phase);
  init_number log_alpha_P_I(Light_sat_alpha_P_I_phase);


 LOCAL_CALCS

    Eb_min                  = Eb_params(IDX_MIN);
    Eb_max                  = Eb_params(IDX_MAX);

    Ep_min                  = Ep_params(IDX_MIN);
    Ep_max                  = Ep_params(IDX_MAX);

    beta_prod_min           = beta_prod_params(IDX_MIN);
    beta_prod_max           = beta_prod_params(IDX_MAX);

    k20_min                 = k20_params(IDX_MIN);
    k20_max                 = k20_params(IDX_MAX);

    alphaR_min              = alphaR_params(IDX_MIN);
    alphaR_max              = alphaR_params(IDX_MAX);

    delta_18O_H2O_min       = delta_18O_H2O_params(IDX_MIN);
    delta_18O_H2O_max       = delta_18O_H2O_params(IDX_MAX);

    init_O2_conc_min        = init_O2_conc_params(IDX_MIN);
    init_O2_conc_max        = init_O2_conc_params(IDX_MAX);

    init_delta_18O_O2_min   = init_delta_18O_O2_params(IDX_MIN);
    init_delta_18O_O2_max   = init_delta_18O_O2_params(IDX_MAX);

    sigma_O2_conc_min       = sigma_O2_conc_params(IDX_MIN);
    sigma_O2_conc_max       = sigma_O2_conc_params(IDX_MAX);

    sigma_delta_18O_O2_min  = sigma_delta_18O_O2_params(IDX_MIN);
    sigma_delta_18O_O2_max  = sigma_delta_18O_O2_params(IDX_MAX);

 END_CALCS


  init_number log_Rref(base_respir_phase);

  init_bounded_number log_Eb(Eb_min,Eb_max,Eb_phase);

  init_bounded_number log_Ep(Ep_min,Ep_max,Ep_phase);

  init_bounded_number log_beta_prod(beta_prod_min,beta_prod_max,beta_prod_phase);

  init_bounded_number log_k20(k20_min,k20_max,k20_phase);

  !! if (alphaR_transform_switch > 0)
  !! {
      init_bounded_number log_alphaR(alphaR_min,alphaR_max,alphaR_phase);
  !! }
  !! else
  !! {
      init_bounded_number alphaR(alphaR_min,alphaR_max,alphaR_phase);
  !! }

  init_bounded_number delta_18O_H2O(delta_18O_H2O_min,delta_18O_H2O_max,delta_18O_H2O_phase);

  init_bounded_number log_init_O2_conc(init_O2_conc_min,init_O2_conc_max,init_O2_conc_phase);
  init_bounded_number log_init_delta_18O_O2(init_delta_18O_O2_min,init_delta_18O_O2_max,init_delta_18O_O2_phase);

  init_bounded_number log_sigma_O2_conc(sigma_O2_conc_min,sigma_O2_conc_max,sigma_O2_conc_phase);
  init_bounded_number log_sigma_delta_18O_O2(sigma_delta_18O_O2_min,sigma_delta_18O_O2_max,sigma_delta_18O_O2_phase);

  vector est_irr(1,n_tC);
  vector est_O2conc(1,n_tC);
  vector est_delta_18O_O2(1,n_tC);

  vector O2_Sat(1,n_tC);     // a function of water temperature

  // replace this
  // sdreport_number Pmax;
  // sdreport_number I_halfSat;
  // sdreport_number alpha_P_I;
  // with this, and update est_val depending upon Light_sat_flag
  number Pmax;
  number I_halfSat;
  number alpha_P_I;
  sdreport_vector est_val(1,2);

  number Rref;
  number Eb;
  number Ep;
  number beta_prod;

  number k20;
  number alphaR_real;

  number init_O2_conc;
  number init_delta_18O_O2;

  vector P(1,n_tC);
  vector R(1,n_tC);
  vector Rb(1,n_tC);
  vector Rp(1,n_tC);
  vector G(1,n_tC);

  number sigma_O2_conc;
  number sigma_delta_18O_O2;

  number Sc20;                      // Schmidt number at 20 degC

  vector f(1,15);

  objective_function_value obj_fun;


PROCEDURE_SECTION

  calculate_initial_values();

  calculate_estimated_irradiance();

  calculate_O2();

  calculate_objective_function();

  if (mceval_phase())
  {
      output_mcmc_parameters();
  }


FUNCTION calculate_initial_values

    // Pmax, in units of mg O2 m-2 hr-1
    Pmax                = mfexp(log_Pmax);                      // Pmax is always positive

    // I @ 1/2 Sat., in units of uE s-1 m-2
    I_halfSat           = mfexp(log_I_halfSat);                 // I_halfSat is always positive

    // the ratio of Pmax to I_halfSat, when light saturation is off
    alpha_P_I           = mfexp(log_alpha_P_I);                 // alpha_P_I is always positive

    // Rref, in units of mg O2 m-2 hr-1
    Rref                = mfexp(log_Rref);                      // Rref is always positive

    // Eb, in units of eV
    Eb                  = mfexp(log_Eb);                        // Eb is always positive

    // Ep, in units of eV
    Ep                  = mfexp(log_Ep);                        // Ep is always positive

    // beta_prod, no units
    beta_prod           = -mfexp(log_beta_prod);                // beta_prod is always negative

    // k20, in units of m hr-1
    k20                 = mfexp(log_k20);                       // k20 is always positive

    // alphaR              = (1.0 / (1.0 + mfexp(-log_alphaR)));   // 0 < alphaR < 1
    if (alphaR_transform_switch > 0)
    {
        alphaR_real     = (((pi / 2.0) + atan(log_alphaR)) / pi); // 0 < alphaR < 1
    }
    else
    {
        alphaR_real     = alphaR;
    }

    // delta 18O-H2O, in units of ‰ vs SMOW
    // delta_18O_H2O       = -mfexp(log_delta_18O_H2O);            // delta_18O_H2O is always negative (true?)

    // initial O2 concentration, in units of mg liter-1
    init_O2_conc        = mfexp(log_init_O2_conc);              // init_O2_conc is always positive

    // initial delta 18O-O2, in units of ‰ vs SMOW
    init_delta_18O_O2   = mfexp(log_init_delta_18O_O2);         // init_delta_18O_O2 is always positive


    sigma_O2_conc       = mfexp(log_sigma_O2_conc);             // sigma_O2_conc is always positive
    sigma_delta_18O_O2  = mfexp(log_sigma_delta_18O_O2);        // sigma_delta_18O_O2 is always positive

    Sc20                = Schmidt_number(20.0);                 // Schmidt number for 20 degC

    if (Light_sat_flag == 1)
    {
        // if light saturation flag == 1, then estimate Pmax and I_halfSat and not alpha_P_I
        est_val(1) = Pmax;
        est_val(2) = I_halfSat;
    }
    else if (Light_sat_flag == 2)
    {
        // if light saturation flag == 2, then estimate alpha_P_I and Pmax and not I_halfSat
        est_val(1) = Pmax;
        est_val(2) = alpha_P_I;
    }
    else
    {
        // if light saturation flag == 0 or other, then estimate alpha_P_I and not Pmax and I_halfSat
        est_val(1) = alpha_P_I;
        est_val(2) = alpha_P_I;
    }

    for (int i = 1; i <= n_tC; i++)
    {
        O2_Sat(i) = O2_saturation(tC_data(i));
    }

    if (debug_flag) cout << "end of calculate_initial_values" << endl;


FUNCTION dvariable Schmidt_number(dvariable tC)

    dvariable Sc    = 0.0;
    dvariable m     = 0.0;
    dvariable sq_tC = square(tC);

    // the Schmidt number for oxygen; see O2_transfer PDF, equations 10 and 11

    // January 2014 - new coefficients from Raymond et al. 2012 L&O:Fl & Envir 2:41-53
    Sc = 1568.0 + (-86.04 * tC) + (2.142 * sq_tC) + (-0.0216 * tC * sq_tC);

    if (salinity > 0.0)
    {
        m   = 2.474e-3 + (3.286e-5 * tC);
        Sc *= (1.0 + (salinity * m));
    }

    if (debug_flag) cout << "end of Schmidt_number" << endl;

    return(Sc);


FUNCTION dvariable k_tC(dvariable tC)

    dvariable k = 0.0;

    k = k20 * sqrt(Sc20 / Schmidt_number(tC));

    if (debug_flag) cout << "end of k_tC" << endl;

    return(k);


FUNCTION dvariable O2_saturation(dvariable tC)

    dvariable O2_sat = 0.0;
    dvariable log_O2 = 0.0;
    dvariable tK     = tC + degC_to_degK;
    dvariable tS     = log((298.15 - tC) / tK);
    dvariable sq_tS  = square(tS);

    // January 2014 - updated per Garcia and Gordon 1992 L&O 37(6):1307-1312
    if (tK > 0.0)
    {
        // January 2014 - coefficients come from Benson and Krause data, per GWH
        log_O2 = 2.00907 + (3.22014 * tS) + (4.05010 * sq_tS) + (4.94457 * tS * sq_tS) + (-0.256847 * square(sq_tS)) + (3.88767 * tS * square(sq_tS)) +
                 (salinity * (-6.24523e-3 + (-7.37614e-3 * tS) + (-1.0341e-2 * sq_tS) + (-8.17083e-3 * tS * sq_tS))) + (-4.88682e-7 * square(salinity));
        O2_sat = atmos_corr * 1.42905 * mfexp(log_O2);
    }

    if (debug_flag) cout << "end of O2_saturation" << endl;

    return(O2_sat);


FUNCTION dvariable delta_to_ratio(dvariable delta_val)

    dvariable ratio_val = 0.0;

    // converts delta notation to isotopic ratios
    ratio_val = ((delta_val / 1000.0) + 1.0) * R_SMOW;

    if (debug_flag) cout << "end of delta_to_ratio" << endl;

    return(ratio_val);


FUNCTION dvariable ratio_to_delta(dvariable ratio_val)

    dvariable delta_val = 0.0;

    // converts isotopic ratios to delta notation
    delta_val = ((ratio_val / R_SMOW) - 1.0) * 1000.0;

    if (debug_flag) cout << "end of ratio_to_delta" << endl;

    return(delta_val);


FUNCTION dvariable S_at_surface(int day_of_year, dvariable hour_of_day)

    // calculate irradiance at the surface, in microEinsteins s-1 m-2

    dvariable S = 0.0;
    dvariable solar_hour = hour_of_day + (4.0 / 60.0 * longitude) - GST_adj - DST_adj;
    dvariable hour_angle = 15.0 * (solar_hour - 12.0);
    dvariable solar_declination, solar_altitude, solar_azimuth, solar_azimuth_corr;
    dvariable angle_incidence, optical_corr, S_outer, S_normal, S_direct, S_diffuse;

    solar_declination = -23.45 * cos(360.0 * (day_of_year + 10.0) / 365.0 * deg_to_rad);
    solar_altitude    = asin(sin(latitude * deg_to_rad) * sin(solar_declination * deg_to_rad) + cos(latitude * deg_to_rad) * cos(solar_declination * deg_to_rad) * cos(hour_angle * deg_to_rad));
    solar_azimuth     = sin(solar_declination * deg_to_rad) * cos(latitude * deg_to_rad) * sin(latitude * deg_to_rad) * cos(hour_angle * deg_to_rad) / cos(solar_altitude);

    if (solar_azimuth > 0.999)
    {
        solar_azimuth = 0.999 / deg_to_rad;
    }
    else
    {
        solar_azimuth = solar_azimuth / deg_to_rad;
    }

    if (hour_angle < 12.0)
    {
        solar_azimuth_corr = solar_azimuth;
    }
    else
    {
        solar_azimuth_corr = 360.0 - solar_azimuth;
    }

    angle_incidence = (sin(solar_altitude) * cos(slope * deg_to_rad)) + (cos(solar_altitude) * sin(slope * deg_to_rad) * cos((solar_azimuth_corr - aspect) * deg_to_rad));

    optical_corr    = atmos_corr * pow(transmiss,(sqrt(1229.0 + square(614.0 * sin(solar_altitude))) - (614.0 * sin(solar_altitude))));

    S_outer  = solar_const * (1.0 + (0.034 * cos(2.0 * pi * day_of_year / 365.0)));
    S_normal = S_outer * optical_corr;

    S_direct = S_normal * angle_incidence;
    if (S_direct < 0.0)
    {
        S_direct = 0.0;
    }

    S_diffuse = S_outer * (0.271 - (0.294 * optical_corr)) * sin(solar_altitude);
    if (S_diffuse < 0.0)
    {
        S_diffuse = 0.0;
    }

    S = 4.57 * 0.45 * est_S_scale * (S_direct + S_diffuse);

    if (debug_flag) cout << "end of S_at_surface" << endl;

    return(S);


FUNCTION calculate_estimated_irradiance

    // fill the est_irr array (dim 1:n_tC) using irr_data (dim 1:n_irr)
    // where n_irr <= n_tC

    int i;

    if (S_meas_flag == 0 || (S_meas_flag != 0 && n_irr != n_tC))
    {
        // calculate irradiance
        for (i = 1; i <= n_tC; i++)
        {
            est_irr(i) = S_at_surface(tC_day_of_year(i), tC_hour_of_day(i));
        }
    }
    else
    {
        // use measured irradiance
        est_irr = irr_data;

    }

    if (debug_flag) cout << "end of calculate_estimated_irradiance" << endl;


FUNCTION dvariable total_primary_production(dvar_vector local_irr_arr, dvariable local_Pmax, dvariable local_I_halfSat, dvar_vector local_tC_time_interval)

    dvariable Pval      = 0.0;
    dvariable p1        = 0.0;
    dvariable p2        = 0.0;
    dvariable time_step = 0.0;
    dvariable diel_dur  = 0.0;
    int i;

    int istart = local_irr_arr.indexmin();
    int iend   = local_irr_arr.indexmax();

    int itime_start = local_tC_time_interval.indexmin();
    int itime_end   = local_tC_time_interval.indexmax();

    if (istart != itime_start || iend != itime_end)
    {
        cout << "ERROR (in total_primary_production): irradiation array is not the same length as temperature time interval array" << endl;
    }

    dvar_vector sum_vector(istart,iend);    // for the whole time period

    sum_vector.initialize();

    p1 = local_I_halfSat / local_Pmax;

    if (Light_sat_flag == 1)
    {
        p2 = 1.0 / local_Pmax;
    }

    // find the element which at the end of the first 24-hour period
    int idx_end_of_first_day = -1;
    for (i = (istart+1); i <= iend && idx_end_of_first_day == -1; i++)
    {
        if (local_tC_time_interval(i) >= 24.0)
        {
            idx_end_of_first_day = i;
        }
    }

    if (idx_end_of_first_day == -1)
    {
        idx_end_of_first_day = iend;
    }

    // calculate total primary production in the ecosystem over each 24-hour time period
    for (i = (istart+1); i <= iend; i++)
    {
        time_step = local_tC_time_interval(i) - local_tC_time_interval(i-1);

        if (local_irr_arr(i) != -99999)
        {
            if (Light_sat_flag == 1)
            {
                // if using light saturation, then use Pmax and I_halfSat and asymptotic function and not alpha_P_I
                sum_vector(i) = (local_irr_arr(i) / (p1 + (p2 * local_irr_arr(i)))) * time_step;
            }
            else if (Light_sat_flag == 2)
            {
                // if not using light saturation, then use alpha_P_I and hyperbolic tangent function and not Pmax and I_halfSat
                sum_vector(i) = local_Pmax * tanh(alpha_P_I * local_irr_arr(i) / local_Pmax) * time_step;
            }
            else
            {
                // if not using light saturation, then use alpha_P_I and linear function and not Pmax and I_halfSat
                sum_vector(i) = (alpha_P_I * local_irr_arr(i)) * time_step;
            }
        }
    }

    // calculate the integrated value for each 24-hour period and then average the values
    Pval = calculate_integrated_moving_window_average(sum_vector, local_tC_time_interval, istart, iend, idx_end_of_first_day);

    if (debug_flag) cout << "end of total_primary_production" << endl;

    return(Pval);


FUNCTION dvariable total_respiration(dvar_vector local_tC_data, dvariable local_log_Rref, dvariable local_Eb, dvariable local_Ep, dvariable local_beta_prod, dvar_vector local_tC_time_interval, int respir_type_flag)

    // respir_type_flag ==  0 -> return integrated community respiration (R = Rb + Rp)
    // respir_type_flag ==  1 -> return integrated production-dependent respiration (Rp)
    // respir_type_flag == -1 -> return integrated base ecosystem respiration (Rb)

    dvariable Rval        = 0.0;
    dvariable time_step   = 0.0;
    dvariable diel_dur    = 0.0;
    dvariable base_respir = 0.0;
    dvariable prod_respir = 0.0;
    dvariable local_tK    = 0.0;
    dvariable tref_mult   = 0.0;
    int i, j, prev_idx;

    if (respir_type_flag < -1 || respir_type_flag > 1)
    {
        cout << "ERROR (in total_respiration): invalid respiration type flag " << respir_type_flag << endl;
        return(-1.0);
    }

    int istart = local_tC_data.indexmin();
    int iend   = local_tC_data.indexmax();

    int itime_start = local_tC_time_interval.indexmin();
    int itime_end   = local_tC_time_interval.indexmax();

    if (istart != itime_start || iend != itime_end)
    {
        cout << "ERROR (in total_respiration): temperature array is not the same length as temperature time interval array" << endl;
        return(-1.0);
    }

    dvar_vector sum_vector(istart,iend);    // for the whole time period

    sum_vector.initialize();

    // find the element which at the end of the first 24-hour period
    int idx_end_of_first_day = -1;
    for (i = (istart+1); i <= iend && idx_end_of_first_day == -1; i++)
    {
        if (local_tC_time_interval(i) >= 24.0)
        {
            idx_end_of_first_day = i;
        }
    }

    if (idx_end_of_first_day == -1)
    {
        idx_end_of_first_day = iend;
    }

    // calculate total community respiration in the ecosystem over each 24-hour time period
    for (i = (istart+1); i <= iend; i++)
    {
        time_step = local_tC_time_interval(i) - local_tC_time_interval(i-1);

        if (local_tC_data(i) != -99999)
        {
            // previous respiration equation
            // sum_vector(i) = local_R20 * time_step * pow(1.047,(local_tC_data(i) - 20.0));

            local_tK      = local_tC_data(i) + degC_to_degK;

            // tref_mult     = (1.0 / (kk * tref)) - (1.0 / (kk * local_tK));
            // REWRITE
            tref_mult     = (local_tK - tref) / (kk * tref * local_tK);

            // base ecosystem respiration equation from Yvon-Durocher et al. 2012 Nature
            if (respir_type_flag <= 0)
            {
                sum_vector(i) = mfexp((local_Eb * tref_mult) + local_log_Rref);
            }

            if (prod_respir_switch > 0 && i > N_steps_for_C_pool && respir_type_flag >= 0)
            {
                // production-dependent respiration equation
                for (j = 1; j <= N_steps_for_C_pool; j++)
                {
                    prev_idx = i - j;

                    // use legit index values only
                    if (prev_idx >= itime_start)
                    {
                        sum_vector(i) += (mfexp(local_Ep * tref_mult) * P(prev_idx) * mfexp(local_beta_prod * j));
                    }
                }
            }
        }

        sum_vector(i) *= time_step;
    }

    // calculate the integrated value for each 24-hour period and then average the values
    Rval = calculate_integrated_moving_window_average(sum_vector, local_tC_time_interval, istart, iend, idx_end_of_first_day);

    if (debug_flag) cout << "end of total_respiration" << endl;

    return(Rval);


FUNCTION dvariable total_mass_flux(dvar_vector local_tC_data, dvar_vector local_O2_Sat, dvar_vector local_est_O2conc, dvar_vector local_tC_time_interval)

    dvariable Gval = 0.0;
    dvariable time_step = 0.0;
    dvariable diel_dur  = 0.0;
    int i;

    int istart = local_tC_data.indexmin();
    int iend   = local_tC_data.indexmax();

    int itime_start = local_tC_time_interval.indexmin();
    int itime_end   = local_tC_time_interval.indexmax();

    if (istart != itime_start || iend != itime_end)
    {
        cout << "ERROR (in total_mass_flux): temperature array is not the same length as temperature time interval array" << endl;
    }

    dvar_vector sum_vector(istart,iend);    // for the whole time period

    sum_vector.initialize();

    // find the element which at the end of the first 24-hour period
    int idx_end_of_first_day = -1;
    for (i = (istart+1); i <= iend && idx_end_of_first_day == -1; i++)
    {
        if (local_tC_time_interval(i) >= 24.0)
        {
            idx_end_of_first_day = i;
        }
    }

    if (idx_end_of_first_day == -1)
    {
        idx_end_of_first_day = iend;
    }

    // calculate the magnutude (absolute value) of total mass flux of O2 by gas exchange
    // in the ecosystem over each 24-hour time period
    for (i = (istart+1); i <= iend; i++)
    {
        time_step = local_tC_time_interval(i) - local_tC_time_interval(i-1);

        if (local_tC_data(i) != -99999)
        {
            sum_vector(i) = fabs(k_tC(local_tC_data(i)) * (local_O2_Sat(i) - local_est_O2conc(i)) * 1000.0) * time_step;
        }
    }

    // calculate the integrated value for each 24-hour period and then average the values
    Gval = calculate_integrated_moving_window_average(sum_vector, local_tC_time_interval, istart, iend, idx_end_of_first_day);

    if (debug_flag) cout << "end of total_mass_flux" << endl;

    return(Gval);


FUNCTION dvariable calculate_integrated_moving_window_average(dvar_vector data_val, dvar_vector local_tC_time_interval, int istart, int iend, int idx_end_of_first_day)

    dvariable avg_val = 0.0;    // average value of each 24-hour window
    int i, j, num_windows;

    num_windows = iend - idx_end_of_first_day + 1;

    // sum the integrated values for each 24-hour period
    for (i = 0; i < num_windows; i++)
    {
        for (j = (istart+1+i); j <= (idx_end_of_first_day + i); j++)
        {
            // avg_val += (data_val(j) / (local_tC_time_interval(idx_end_of_first_day + i) - local_tC_time_interval(istart + i)));
            // 2009-03-01 - don't divide by the day length (in hours)
            avg_val += (data_val(j));
        }
    }

    // divide by number of 24-hour windows for the average
    avg_val /= double(num_windows);

    return(avg_val);


FUNCTION calculate_O2

    // this function assumes that irr_arr and tC_data are (1,n_tC)

    dvariable p1 = 0.0;
    dvariable p2 = 0.0;
    dvariable tK = 0.0;

    dvariable tref_mult = 0.0;

    dvariable K_O2, O2_Sat_mg, a, B, c, D, P_Iso, alphaGeq, R_H2O;

    // for solving the two differential equations for each time step
    dvar_vector x(1,NUM_DIFFEQ), y(1,NUM_DIFFEQ);

    int j, prev_idx;

    est_O2conc.initialize();
    est_delta_18O_O2.initialize();

    P.initialize();
    R.initialize();
    Rb.initialize();
    Rp.initialize();
    G.initialize();

    // calculate the est_O2conc and the est_delta_18O_O2
    p1 = I_halfSat / Pmax;

    if (Light_sat_flag == 1)
    {
        p2 = 1.0 / Pmax;
    }

    R_H2O = delta_to_ratio(delta_18O_H2O);

    // initialize the estimated values
    est_O2conc(1)       = init_O2_conc;
    est_delta_18O_O2(1) = init_delta_18O_O2;

    if (Light_sat_flag == 1)
    {
        // if using light saturation, then use Pmax and I_halfSat and asymptotic function and not alpha_P_I
        P(1) = (est_irr(1) / (p1 + (p2 * est_irr(1)))) * area;
    }
    else if (Light_sat_flag == 2)
    {
        // if not using light saturation, then use alpha_P_I and Pmax and hyperbolic tangent function and not I_halfSat
        P(1) = Pmax * tanh(alpha_P_I * est_irr(1) / Pmax) * area;
    }
    else
    {
        // if not using light saturation, then use alpha_P_I and linear function and not Pmax and I_halfSat
        P(1) = (alpha_P_I * est_irr(1)) * area;
    }

    // rate of respiration
    // previous respiration equation
    // R(1) = R20 * area * pow(1.047,(tC_data(1) - 20.0));

    tK = tC_data(1) + degC_to_degK;

    // tref_mult = (1.0 / (kk * tref)) - (1.0 / (kk * tK));
    // REWRITE
    tref_mult = (tK - tref) / (kk * tref * tK);

    // base ecosystem respiration equation
    Rb(1) = mfexp((Eb * tref_mult) + log_Rref) * area;

    if (prod_respir_switch > 0 && 1 > N_steps_for_C_pool)
    {
        // production-dependent respiration equation
        for (j = 1; j <= N_steps_for_C_pool; j++)
        {
            prev_idx = 1 - j;

            // use legit index values only
            if (prev_idx >= 1)
            {
                // NOPE, because this is the beginning of the time series, so there are no values for preceding time steps
                // Rp(1) += (mfexp(Ep * tref_mult) * P(prev_idx) * mfexp(beta_prod * j));
            }
        }
    }

    // total respiration = base respiration and prod respiration
    R(1) = Rb(1) + Rp(1);

    G(1) = (k_tC(tC_data(1)) * (O2_Sat(1) - est_O2conc(1)) * 1000.0) * area;

    for (int i = 2; i <= n_tC; i++)
    {
        x.initialize();
        y.initialize();

        // rate of photosynthesis
        if (Light_sat_flag == 1)
        {
            // if using light saturation, then use Pmax and I_halfSat and asymptotic function and not alpha_P_I
            P(i)  = (est_irr(i) / (p1 + (p2 * est_irr(i)))) * area;
        }
        else if (Light_sat_flag == 2)
        {
            // if not using light saturation, then use alpha_P_I and Pmax and hyperbolic tangent function and not I_halfSat
            P(i)  = Pmax * tanh(alpha_P_I * est_irr(i) / Pmax) * area;
        }
        else
        {
            // if not using light saturation, then use alpha_P_I and linear function and not Pmax and I_halfSat
            P(i)  = (alpha_P_I * est_irr(i)) * area;
        }

        // rate of respiration
        // previous respiration equation
        // R(i)      = R20 * area * pow(1.047,(tC_data(i) - 20.0));

        tK        = tC_data(i) + degC_to_degK;

        // tref_mult = (1.0 / (kk * tref)) - (1.0 / (kk * tK));
        // REWRITE
        tref_mult = (tK - tref) / (kk * tref * tK);

        // base ecosystem respiration equation
        Rb(i)     = mfexp((Eb * tref_mult) + log_Rref) * area;

        if (prod_respir_switch > 0 && i > N_steps_for_C_pool)
        {
            // production-dependent respiration equation
            for (j = 1; j <= N_steps_for_C_pool; j++)
            {
                prev_idx = i - j;

                // use legit index values only
                if (prev_idx >= 1)
                {
                    Rp(i) += (mfexp(Ep * tref_mult) * P(prev_idx) * mfexp(beta_prod * j));
                }
            }

            Rp(i) *= area;
        }

        // total respiration = base respiration and prod respiration
        R(i)      = Rb(i) + Rp(i);

        // atmos-water gas exchange coefficient
        K_O2      = k_tC(tC_data(i)) / depth;

        // convert to mg
        O2_Sat_mg = O2_Sat(i) * 1000.0 * depth * area;

        // parameters for solving system of two differential equations
        alphaGeq  = ((-0.72951 + (426.96 / (tC_data(i) + degC_to_degK))) / 1000.0) + 1.0;

        a         = P(i) - R(i) + (K_O2 * O2_Sat_mg);
        B         = alphaR_real * R(i);
        c         = K_O2 * alphaGk;
        D         = O2_Sat_mg * alphaGeq * R_air;

        P_Iso     = P(i) * alphaP * R_H2O;

        // for O2
        x(1)      =  est_O2conc(i - 1) * 1000.0 * depth * area;

        // for delta 18 O to 16 O
        x(2)      = x(1) * delta_to_ratio(est_delta_18O_O2(i - 1));

        // solve the two differential equations for this time step
        y = solve_diffeq_RungeKutta(tC_time_interval(i-1), tC_time_interval(i), RK_step_size, x, K_O2, P_Iso, a, B, c, D);

        // retrieve values
        est_O2conc(i)       = y(1) * 0.001 / (depth * area);
        est_delta_18O_O2(i) = ratio_to_delta(y(2) / y(1));

        // rate of total mass flux of O2 by gas exchange
        G(i)      = (k_tC(tC_data(i)) * (O2_Sat(i) - est_O2conc(i)) * 1000.0) * area;
    }

    if (debug_flag) cout << "end of calculate_O2" << endl;


FUNCTION dvariable discrete_diffeq(int diffeq_num, dvariable time_step, dvar_vector x, dvariable K_O2, dvariable P_Iso, dvariable a, dvariable B, dvariable c, dvariable D)

    // this function assumes that x is (1,NUM_DIFFEQ)

    dvariable diffeq_val = 0.0;

    // calculate value for discrete interval for system of differential equations
    if (diffeq_num == 1)
    {
        diffeq_val = a - (K_O2 * x(1));
    }
    else if (diffeq_num == 2 && x(1) > 0.0)
    {
        diffeq_val = P_Iso - (B * (x(2) / x(1))) + (c * (D - x(2)));
    }

    if (debug_flag) cout << "end of discrete_diffeq" << endl;

    return(diffeq_val);


FUNCTION dvar_vector solve_diffeq_RungeKutta(dvariable prev_time, dvariable curr_time, int step_size, dvar_vector x, dvariable K_O2, dvariable P_Iso, dvariable a, dvariable B, dvariable c, dvariable D)

    // this function assumes that x is (1,NUM_DIFFEQ)

    dvar_vector diffeq_arr(1,NUM_DIFFEQ);
    dvariable time_start, time_incr;
    dvar_vector yt(1,NUM_DIFFEQ), k1(1,NUM_DIFFEQ), k2(1,NUM_DIFFEQ), k3(1,NUM_DIFFEQ), k4(1,NUM_DIFFEQ);
    int i, j, k;

    diffeq_arr.initialize();
    diffeq_arr = x;

    time_incr = (curr_time - prev_time) / double(step_size);

    for (i = 0; i < step_size; i++)
    {
        time_start = prev_time + (double(i) * time_incr);

        for (j = 1; j <= NUM_DIFFEQ; j++)
        {
            k1(j) = time_incr * discrete_diffeq(j, time_start, diffeq_arr, K_O2, P_Iso, a, B, c, D);
        }

        for (j = 1; j <= NUM_DIFFEQ; j++)
        {
            yt(j) = diffeq_arr(j) + (0.5 * k1(j));
        }

        for (j = 1; j <= NUM_DIFFEQ; j++)
        {
            k2(j) = time_incr * discrete_diffeq(j, time_start + (time_incr * 0.5), yt, K_O2, P_Iso, a, B, c, D);
        }

        for (j = 1; j <= NUM_DIFFEQ; j++)
        {
            yt(j) = diffeq_arr(j) + (0.5 * k2(j));
        }

        for (j = 1; j <= NUM_DIFFEQ; j++)
        {
            k3(j) = time_incr * discrete_diffeq(j, time_start + (time_incr * 0.5), yt, K_O2, P_Iso, a, B, c, D);
        }

        for (j = 1; j <= NUM_DIFFEQ; j++)
        {
            yt(j) = diffeq_arr(j) + k3(j);
        }

        for (j = 1; j <= NUM_DIFFEQ; j++)
        {
            k4(j) = time_incr * discrete_diffeq(j, time_start + time_incr, yt, K_O2, P_Iso, a, B, c, D);
        }

        for (j = 1; j <= NUM_DIFFEQ; j++)
        {
            diffeq_arr(j) = diffeq_arr(j) + ((k1(j) + (2.0 * k2(j)) + (2.0 * k3(j)) + k4(j)) / 6.0);
        }
    }

    if (debug_flag) cout << "end of solve_diffeq_RungeKutta" << endl;

    return(diffeq_arr);


FUNCTION calculate_objective_function

    int i;

    f.initialize();

    // NLL components for O2 concentration
    for (i = 1; i <= n_O2conc; i++)
    {
        f(1) += square(O2conc_data(i) - est_O2conc(O2conc_time_interval_map(i)));
    }
    f(1) /= (2.0 * square(sigma_O2_conc));

    if (sigma_O2_conc_switch > 0)
    {
        f(2) = n_O2conc * log(sigma_O2_conc);
    }

    if (isotope_switch != -1)
    {
        // NLL components for delta_18O_O2
        for (i = 1; i <= n_delta_18O_O2; i++)
        {
            f(3) += square(delta_18O_O2_data(i) - est_delta_18O_O2(delta_18O_O2_time_interval_map(i)));
        }
        f(3) /= (2.0 * square(sigma_delta_18O_O2));

        if (sigma_delta_18O_O2_switch != -1)
        {
            f(4) = n_delta_18O_O2 * log(sigma_delta_18O_O2);
        }

        // scale factor for delta_18O_O2 NLL components
        f(3) *= delta_18O_O2_scale_factor;
        f(4) *= delta_18O_O2_scale_factor;
    }

    // normal priors on parameters
    if (Eb_switch == 2)
    {
        f(5) = square(Eb - Eb_params(IDX_MEAN)) / (2.0 * square(Eb_params(IDX_STD_DEV)));
    }

    if (Ep_switch == 2)
    {
        f(6) = square(Ep - Ep_params(IDX_MEAN)) / (2.0 * square(Ep_params(IDX_STD_DEV)));
    }

    if (k20_switch == 2)
    {
        f(7) = square(k20 - k20_params(IDX_MEAN)) / (2.0 * square(k20_params(IDX_STD_DEV)));
    }

    if (alphaR_switch == 2)
    {
        f(8) = square(alphaR_real - alphaR_params(IDX_MEAN)) / (2.0 * square(alphaR_params(IDX_STD_DEV)));
    }

    if (delta_18O_H2O_switch == 2 && isotope_switch != -1)
    {
        f(9) = square(delta_18O_H2O - delta_18O_H2O_params(IDX_MEAN)) / (2.0 * square(delta_18O_H2O_params(IDX_STD_DEV)));
    }

    if (init_O2_conc_switch == 2)
    {
        f(10) = square(init_O2_conc - init_O2_conc_params(IDX_MEAN)) / (2.0 * square(init_O2_conc_params(IDX_STD_DEV)));
    }

    if (init_delta_18O_O2_switch == 2 && isotope_switch != -1)
    {
        f(11) = square(init_delta_18O_O2 - init_delta_18O_O2_params(IDX_MEAN)) / (2.0 * square(init_delta_18O_O2_params(IDX_STD_DEV)));
    }

    if (sigma_O2_conc_switch == 2)
    {
        f(12) = square(sigma_O2_conc - sigma_O2_conc_params(IDX_MEAN)) / (2.0 * square(sigma_O2_conc_params(IDX_STD_DEV)));
    }

    if (sigma_delta_18O_O2_switch == 2 && isotope_switch != -1)
    {
        f(13) = square(sigma_delta_18O_O2 - sigma_delta_18O_O2_params(IDX_MEAN)) / (2.0 * square(sigma_delta_18O_O2_params(IDX_STD_DEV)));
    }

    if (beta_prod_switch == 2)
    {
        f(14) = square(beta_prod - beta_prod_params(IDX_MEAN)) / (2.0 * square(beta_prod_params(IDX_STD_DEV)));
    }

    obj_fun = sum(f);

    if (debug_flag) cout << "end of calculate_objective_function" << endl;


FUNCTION output_mcmc_parameters

    ofstream mcmcfile("BaMM_mcmc.dat", ios::out | ios::app);

    mcmc_iter++;

    if (mcmcfile.is_open())
    {
        mcmcfile << mcmc_iter << "\t" << obj_fun;

        if (Light_sat_flag == 1)
        {
            mcmcfile << "\tPmax\t" << Pmax << "\tI_halfSat\t" << I_halfSat;
        }
        else if (Light_sat_flag == 2)
        {
            mcmcfile << "\tPmax\t" << Pmax << "\talpha_P_I\t" << alpha_P_I;
        }
        else
        {
            mcmcfile << "\talpha_P_I\t" << alpha_P_I;
        }

        mcmcfile << "\tRref\t" << Rref << "\tEb\t" << Eb << "\tEp\t" << Ep << "\tbeta_prod\t" << beta_prod << "\tk20\t" << k20;

        if (alphaR_transform_switch > 0)
        {
            mcmcfile << "\tlog_alphaR\t" << log_alphaR;
        }

        mcmcfile << "\talphaR\t" << alphaR_real << "\tdelta_18O_H2O\t" << delta_18O_H2O << "\tinit_O2_conc\t" << init_O2_conc << "\tinit_delta_18O_O2\t" << init_delta_18O_O2 << "\tsigma_O2_conc\t" << sigma_O2_conc << "\tsigma_delta_18O_O2\t" << sigma_delta_18O_O2 << "\test_O2conc\t" << est_O2conc << "\test_delta_18O_O2\t" << est_delta_18O_O2 << "\tIntegrated PP\t" << total_primary_production(est_irr, Pmax, I_halfSat, tC_time_interval) << "\tIntegrated CR\t" << total_respiration(tC_data, log_Rref, Eb, Ep, beta_prod, tC_time_interval,0) << "\tIntegrated G\t" << total_mass_flux(tC_data, O2_Sat, est_O2conc, tC_time_interval) << "\tP\t" << P << "\tR\t" << R << "\tG\t" << G << "\tRb\t" << Rb << "\tRp\t" << Rp << endl;
    }

    mcmcfile.close();

    if (debug_flag) cout << "end of output_mcmc_parameters" << endl;


REPORT_SECTION

    report << "objective function" << endl;
    report << obj_fun << endl;
    report << "objective function components" << endl;
    report << f << endl;
    report << endl;

    report << "Estimated and derived values" << endl;

    if (Light_sat_flag == 1)
    {
        report << "Pmax (in mg O2 m-2 hr-1)\t" << Pmax << endl;
        report << "I_halfSat (in microE s-1 m-2)\t" << I_halfSat << endl;
    }
    else if (Light_sat_flag == 2)
    {
        report << "Pmax (in mg O2 m-2 hr-1)\t" << Pmax << endl;
        report << "alpha_P_I\t" << alpha_P_I << endl;
    }
    else
    {
        report << "alpha_P_I\t" << alpha_P_I << endl;
    }

    report << "Tref (in degrees Celsius)\t" << trefC << "\t(input, not estimated)" << endl;
    report << "Rref (in mg O2 m-2 hr-1)\t" << Rref << (base_respir_switch == -1 ? "\t(not estimated)" : "") << endl;
    report << "Eb (in eV)\t" << Eb << (Eb_switch == -1 ? "\t(not estimated)" : "") << endl;
    report << "Ep (in eV)\t" << Ep << (Ep_switch == -1 ? "\t(not estimated)" : "") << endl;
    report << "beta_prod\t" << beta_prod << (beta_prod_switch == -1 ? "\t(not estimated)" : "") << endl;
    report << "N (number of preceding time steps for carbon pool accumulation)\t" << N_steps_for_C_pool << "\t(not estimated)" << endl;
    report << "k20 (in m hr-1)\t" << k20 << (k20_switch == -1 ? "\t(not estimated)" : "") << endl;
    report << "alphaR\t" << alphaR_real << (alphaR_switch == -1 ? "\t(not estimated)" : "") << endl;
    report << "delta_18O_H2O (in ‰ vs SMOW)\t" << delta_18O_H2O << ((delta_18O_H2O_switch == -1 || isotope_switch == -1) ? "\t(not estimated)" : "") << endl;
    report << "init_O2_conc (in mg liter-1)\t" << init_O2_conc << (init_O2_conc_switch == -1 ? "\t(not estimated)" : "") << endl;
    report << "init_delta_18O_O2 (in ‰ vs SMOW)\t" << init_delta_18O_O2 << ((init_delta_18O_O2_switch == -1 || isotope_switch == -1) ? "\t(not estimated)" : "") << endl;
    report << "sigma_O2_conc\t" << sigma_O2_conc << (sigma_O2_conc_switch == -1 ? "\t(not estimated)" : "") << endl;
    report << "sigma_delta_18O_O2\t" << sigma_delta_18O_O2 << ((sigma_delta_18O_O2_switch == -1 || isotope_switch == -1) ? "\t(not estimated)" : "") << endl;
    report << endl;

    report << "Scale factor for estimated S\t" << est_S_scale << endl;
    report << endl;

    report << "Integrated primary production (P, in mg O2 m-2 d-1)\t" << total_primary_production(est_irr, Pmax, I_halfSat, tC_time_interval) << endl;
    report << "Integrated community respiration (R, in mg O2 m-2 d-1)\t" << total_respiration(tC_data, log_Rref, Eb, Ep, beta_prod, tC_time_interval,0) << endl;
    report << "Integrated magnitude of total mass flux of O2 by gas exchange (G, in mg O2 m-2 d-1)\t" << total_mass_flux(tC_data, O2_Sat, est_O2conc, tC_time_interval) << endl;
    report << "Integrated production-dependent respiration (Rp, in mg O2 m-2 d-1)\t" << total_respiration(tC_data, log_Rref, Eb, Ep, beta_prod, tC_time_interval,1) << endl;
    report << "Integrated base ecosystem respiration (Rb, in mg O2 m-2 d-1)\t" << total_respiration(tC_data, log_Rref, Eb, Ep, beta_prod, tC_time_interval,-1) << endl;
    report << endl;

    report << "Interval\tDay of year\tHour of day\tS (microE s-1 m-2)\tO2 conc (mg liter-1)\tdelta 18O_O2 (‰ vs SMOW)\tO2 sat\tPct Sat\tP\tR\tG\tRb\tRp" << endl;
    for (int i = 1; i <= n_tC; i++)
    {
        report << i << "\t" << tC_day_of_year(i) << "\t" << tC_hour_of_day(i) << "\t" << est_irr(i) << "\t" << est_O2conc(i) << "\t" << est_delta_18O_O2(i) << "\t" << O2_Sat(i) << "\t" << (100.0 * est_O2conc(i) / O2_Sat(i)) << "\t" << P(i) << "\t" << R(i) << "\t" << G(i) << "\t" << Rb(i) << "\t" << Rp(i) << endl;
    }
    report << endl;


TOP_OF_MAIN_SECTION

  gradient_structure::set_NUM_DEPENDENT_VARIABLES(2000);
  gradient_structure::set_MAX_NVAR_OFFSET(2000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(1000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(10000000);
  arrmblsize = 20000000; // use instead of gradient_structure::set_ARRAY_MEMBLOCK_SIZE


RUNTIME_SECTION

    convergence_criteria 1.e-4 1.e-4 1.e-4 1.e-4 1.e-5 1.e-5 1.e-7
    maximum_function_evaluations 5000 5000 5000 5000 5000 5000 5000



