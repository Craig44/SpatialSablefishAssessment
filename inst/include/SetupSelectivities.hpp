/*
*  SetupSelectivities.h
*
*  This script has functions that are used to map selectivity parameters and control variables
*  to develop an array (age x year) for each fishery and survey.
*
*/


// logistic ogive function parametersised with a_50 and a_to95
template <class Type>
vector<Type> logistic_ogive(vector<Type>& ages, Type& sel_50, Type& sel_95) {
  //std::cout << "logistic_ogive\n";
  int n_ages = ages.size();
  vector<Type> logis(n_ages);
  for (int age = 0;  age < n_ages; ++age) {
    logis[age] = Type(1.0) / (Type(1.0) + pow(Type(19.0), (sel_50 - ages[age]) / sel_95));
  }
  return logis;
}
// Alternative logistic ogive function parametersised with a_50 and delta
template <class Type>
vector<Type> logistic_alt_ogive(vector<Type>& ages, Type& sel_50, Type& delta) {
  //std::cout << "logistic_ogive\n";
  int n_ages = ages.size();
  vector<Type> logis(n_ages);
  for (int age = 0;  age < n_ages; ++age) {
    logis[age] = Type(1.0) / (Type(1.0) + exp(-delta * (ages[age] - sel_50)));
  }
  return logis;
}
// Double normal Punt et. al 1996 gamma parameterization
template <class Type>
vector<Type> double_normal_ogive(vector<Type>& ages, Type& sel_50, Type& delta) {
  //std::cout << "logistic_ogive\n";
  int n_ages = ages.size();
  vector<Type> d_norm(n_ages);
  for (int age = 0;  age < n_ages; ++age) {
    d_norm[age] = (pow(ages[age] / sel_50, sel_50 / (0.5*(sqrt(square(sel_50)+4*square(delta))-sel_50)))*exp((sel_50-ages[age])/(0.5*(sqrt(square(sel_50)+4*square(delta))-sel_50))));
  }
  return d_norm;
}

// Taken from the SS3 manual
template <class Type>
vector<Type> double_normal_SS3_ogive(vector<Type>& ages, Type& min_x_val_, Type& first_sel_val_, Type& p1trans_, Type& p2trans_, Type& p3trans_, Type& p4trans_, Type& p5trans_, Type& p6trans_) {
  //std::cout << "logistic_ogive\n";
  int n_ages = ages.size();
  vector<Type> d_norm(n_ages);
  // NOTE: NEED to get the model min and max age or min and max size, and the bin width (assumed to be 1 below)
  // min age is popn.state.col_min or get_int("min_age")
  // max age is popn.state.col_max or get_int("max_age")
  // length bins is popn.state.class_mins; min size is the value in the first bin and max size is the value in the last bin
  // from Adam Langley, sent on 2022-02-22
  for (int age_ndx = 0;  age_ndx < n_ages; ++age_ndx) {
    if (ages[age_ndx] > min_x_val_) {
      // see also https://github.com/nmfs-stock-synthesis/stock-synthesis/blob/main/SS_selex.tpl, lines 436 and 504
      Type x_p1trans   = ages[age_ndx] - p1trans_;
      Type x_p2trans   = ages[age_ndx] - p2trans_;
      Type asc         = exp(-1.0 * (x_p1trans * x_p1trans / p3trans_));
      Type asc_scaled  = (p5trans_ + (1.0 - p5trans_) * (asc - 0.0) / (1.0 - 0.0));
      Type desc        = exp(-1.0 * (x_p2trans * x_p2trans / p4trans_));
      Type stj         = exp(-1.0 * ((40.0 - p2trans_) * (40.0 - p2trans_) / p4trans_));
      Type desc_scaled = (1.0 + (p6trans_ - 1.0) * (desc - 1.0) /(stj - 1.0));
      Type join1       = 1.0 / (1.0 + exp(-1.0 * (20.0 * x_p1trans / (1.0 + fabs(x_p1trans)))));
      Type join2       = 1.0 / (1.0 + exp(-1.0 * (20.0 * x_p2trans / (1.0 + fabs(x_p2trans)))));
      d_norm[age_ndx] = asc_scaled * (1.0 - join1) + join1 * ((1.0 - join2) + (desc_scaled * join2));
    } else {
      d_norm[age_ndx] = exp(first_sel_val_);
    }
  }
  return d_norm;
}






/*
 *  BuildSelectivity is the function that will build a selectivity array that can have different time-blocks and a different ogive per block given parameters and control variables
 *  @param sel_params dim: parameters (untransformed) where n_time_blocks is the number of time-blocks n_time_blocks x n_params
 *  @param sel_type vector of length = n_time_blocks. THis will define how many columns are in sel_params
 *  0 = Specifies whether it is a logistic (= 0)
 *  1 = Specifies the gamma formuation which mimics the double normal with two parameters
 *  2 = power function not scaled to have a max = 1 like in the assessment
 *  3 = Alternative logistic parameterisation
 *  4 = exponential decay selectivity
 *  5 = double normal with 3 parameters
 *  6 = SS3 double normal with 8 parameters (usually only 3-4 are estimated)
 *  @param ages vector of ages.
 *  @param sel_array is a prepopulated array that will be fulled with the correct . Specifies which time-block selectivity to which year.
 *  @param normalise if true then it will scale the selectivity to have max value = 1
 */
template <class Type>
void BuildSelectivity(array<Type> sel_params, vector<int> sel_type, vector<Type> ages, array<Type>& sel_array, bool normalise) {
  int n_time_blocks = sel_params.dim(0);
  int n_ages = ages.size();

  for(int i = 0; i < n_time_blocks; ++i) {
    if(sel_type(i) == 0) {
      // logistic, expect two parameters
      sel_array.col(i) = logistic_alt_ogive(ages, sel_params(i, 0), sel_params(i, 1));
    } else if (sel_type(i) == 1) {
      // Double normal
      sel_array.col(i) = double_normal_ogive(ages, sel_params(i, 0), sel_params(i, 1));
    } else if (sel_type(i) == 2) {
      // Power function used for the GOA trawl selectivity
      for(int age_ndx = 0; age_ndx < n_ages; ++age_ndx)
        sel_array(age_ndx, i) = (1.0 / pow(ages[age_ndx],sel_params(i, 0)));
      // scale by max TODO: DELETE THIS after validated it shouldn't be here
      // A max call can cause AD issues. Or cause a split in the AD graph which I don't like at all.... consider using the expoential decay i.e., sel_type == 4 instead
    }  else if (sel_type(i) == 3) {
      // Alternative logistic parameterisation
      sel_array.col(i) = logistic_ogive(ages, sel_params(i, 0), sel_params(i, 1));
    } else if (sel_type(i) == 4) {
      // Exponential decay function
      for(int age_ndx = 0; age_ndx < n_ages; ++age_ndx)
        sel_array(age_ndx, i) = exp(-1.0 * ages[age_ndx] * sel_params(i, 0));
    } else if (sel_type(i) == 5) {
      // Double normal with 3 parameters
      // param order mu, sigma_r, sigma_l
      Type delta = 5.0;
      Type tmp = log(0.5);
      Type stmp = 0.0;
      for(int age_ndx = 0; age_ndx < n_ages; ++age_ndx) {
      stmp = 1.0 / (1.0 + exp(-delta * (ages[age_ndx] - sel_params(i, 0))));
      sel_array(age_ndx, i) = stmp * exp(tmp * square((ages[age_ndx] - sel_params(i, 0)) / sel_params(i, 1))) + (1.0 - stmp) *
        exp(tmp * square((ages[age_ndx] - sel_params(i, 0)) / sel_params(i, 2)));
      }
    }
    // scale by max TODO: DELETE THIS after validated it shouldn't be here
    // A max call can cause AD issues.
    if(normalise) {
      Type max_sel = max(vector<Type>(sel_array.col(i)));
      for(int age_ndx = 0; age_ndx < n_ages; ++age_ndx)
        sel_array(age_ndx, i) /= max_sel;
    }
  }
}



