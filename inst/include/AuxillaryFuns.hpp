
// Zerofun functions
template <class Type>
Type ZeroFun(Type x, Type delta) {
  if (x >= delta)
    return x;

  return delta / (2.0 - (x / delta));
}

/*
 * Updated the multinom fun because TMB's version cannot deal with p = 0
 * where as R's can
 */
template <class Type>
Type dmultinom_upd(vector<Type> x, vector<Type> p, int give_log=0)
{
  p /= p.sum(); // double check its nomalised
  Type logres = lgamma(x.sum() + Type(1));
  Type p_delta;
  Type delta = 0.00001;
  for(int i = 0; i < p.size(); ++i) {
    p_delta = ZeroFun(p(i), delta); // delta for the zerofun function to stop p = 0
    logres += x(i)*log(p_delta) - lgamma(x(i)+Type(1));
  }
  if(give_log)
    return logres;
  else
    return exp(logres);
}

/*
 *  an index folding method to help get the region and release year for tagged fish in the tagged partition
 *  @param region_ndx (starts at 0 goes to (n_regions - 1))
 *  @param release_year_ndx (0:n_years_to_retain_tagged_cohorts_for - 1). 0 is current year release, 1 is tag fished release last year, 2 is tag release 2 years before current
 *  @return an index to look up the tagged partition which covers both these indicies
 */
int get_tag_release_event_ndx(int region_ndx, int release_event_year_ndx, int n_regions) {
  return release_event_year_ndx * n_regions + region_ndx;
}
/*
 *  square of a variable
 */
template <class Type>
Type square(Type x){return x*x;}
/*
 *  square of a variable overloaded for vector input
 */
template <class Type>
vector<Type> square(vector<Type> x) {
  return x*x;
}

// Return the natural logarithm of one plus the specified value.
template <class Type>
Type log1p(Type x) {
  return log(1 + x);
}
//Calculates the log of 1 plus the exponential of the specified value without overflow.
template <class Type>
Type log1p_exp(Type a) {
  if (a > 0.0) {
    return a + log1p(exp(-a));
  }
  return log1p(exp(a));
}

// transform Y -Inf-Inf -> X bound lb - ub
template <class Type>
Type invlogit_general(Type& Y, Type& lb, Type& ub) {
  return(lb + (ub - lb) * (1 / (1 + exp(-Y))));
}
// transform Y -Inf-Inf -> X bound lb - ub
template <class Type>
vector<Type> invlogit_general(vector<Type>& Y, Type& lb, Type& ub) {
  vector<Type> X(Y.size());
  for(int i = 0; i < X.size(); ++i) {
    X(i) = lb + (ub - lb) * (1 / (1 + exp(-Y(i))));
  }
  return(X);
}

/* Parameter transform */
template <class Type>
Type f(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}


/*
 * Geometric mean
 */
template <class Type>
Type geo_mean(vector<Type>& x){
  return exp((log(x).sum())/x.size());
}

/*
 *
 */
template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x,eps,Type(0.01)*pow(x-eps,2),Type(0));
  return CppAD::CondExpGe(x,eps,x,eps/(Type(2)-x/eps));
}


template <class Type>
void get_covar(vector<Type> theta, matrix<Type>& covar, int n, int type) {
  // Source
  // https://github.com/glmmTMB/glmmTMB/blob/9cb23d77afc042be5bced3713c40529023d9ef7e/glmmTMB/src/glmmTMB.cpp#L361
  if (type == 1) {
    // case: diag_covstruct
    //std::cerr << "diag_covstruct" << std::endl;
    Type sd = exp(theta(0));
    for(int i = 0; i < n; ++i)
      covar(i,i) = sd * sd;
  } else if (type == 2) {
    // case: us_covstruct
    //std::cerr << "us_covstruct" << std::endl;
    vector<Type> sd = exp(theta.head(n));
    vector<Type> corr_transf = theta.tail(theta.size() - n);
    density::UNSTRUCTURED_CORR_t<Type> nldens(corr_transf);
    // for covariance
    //density::VECSCALE_t<density::UNSTRUCTURED_CORR_t<Type> > scnldens = density::VECSCALE(nldens, sd);
    for(int i=0; i<n; i++) {
      for(int j=0; j<n; j++) {
        covar(i,j) = nldens.cov()(i,j) * sd(i) * sd(j);
      }
    }
  } else if (type == 3){
    // case: cs_covstruct Compound Symmetry: Heterogenous.
    //std::cerr << "cs_covstruct" << std::endl;

    vector<Type> sd = exp(theta.head(n));
    Type corr_transf = theta(n);
    Type a = Type(1) / (n - Type(1)); // Correcting for violations of sphericity
    Type rho = invlogit(corr_transf) * (Type(1) + a) - a;
    for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
        covar(i,j) = (i==j ? Type(1) * sd(i) * sd(i) : rho * sd(i) * sd(j)) ;

  } else if (type == 4){
    // case: toep_covstruct Toeplitz: Heterogenous.
    //std::cerr << "toep_covstruct" << std::endl;
    vector<Type> sd = exp(theta.head(n));
    vector<Type> parms = theta.tail(n-1);              // Corr parms
    parms = parms / sqrt(Type(1.0) + parms * parms );  // Now in (-1,1)
    for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
        covar(i,j) = (i==j ? Type(1) * sd(i) * sd(i) : parms( (i > j ? i-j : j-i) - 1 ) * sd(i) * sd(j));
  } else if (type == 5){
    //std::cerr << "ar1_covstruct" << std::endl;
    // case: ar1_covstruct
    //  * NOTE: Valid parameter space is phi in [-1, 1]
    //  * NOTE: 'times' not used as we assume unit distance between consecutive time points.
    vector<Type> sd(n);
    sd.fill(exp(theta(0)));
    Type corr_transf = theta(1);
    Type phi = corr_transf / sqrt(1.0 + pow(corr_transf, 2));
    for(int i=0; i < n; i++){
      for(int j=0; j < n; j++){
        covar(i,j) = pow(phi, abs(i-j)) * sd(i) * sd(i);
      }
    }
  }

  return;
}



/*
 * centred log transform
 */
template <class Type>
vector<Type> crl(vector<Type>& x) {
  return log(x / geo_mean(x));
}
/*
 * rmultinomm - for simulate call
 */
template <class Type>
vector<Type> rmultinom(vector<Type> prob, Type N) {
  vector<Type> sim_X(prob.size());
  sim_X.setZero();
  // Now simulate using the uniform random variable
  Type rng_uniform;
  Type cumulative_expect;
  while(N > 0) {
    rng_uniform = runif(Type(0),Type(1));
    //std::cout << rng_uniform << " ";
    cumulative_expect = 0.0;
    for (unsigned i = 0; i < prob.size(); ++i) {
      cumulative_expect += prob[i];
      if (cumulative_expect >= rng_uniform) {
        sim_X[i] += 1.0;
        break;
      }
      //sim_X[prob.size() - 1] += 1.0;
    }
    N -= 1;
  }
  //std::cout << "\n";
  return(sim_X);
}

/*
 *  Simulate a single draw from a multinomial-dirichlet distribution
 */
template <class Type>
vector<Type> rdirichletmulti(vector<Type> fitted_props, Type& n_eff, Type& theta) {
  vector<Type> dirichlet_draw(fitted_props.size());
  for(int ndx = 0; ndx < fitted_props.size(); ndx++)
    dirichlet_draw(ndx) = rgamma(fitted_props(ndx) * theta * n_eff, (Type)1.0);// shape, rate = 1.0

  Type dirich_total = dirichlet_draw.sum();
  dirichlet_draw /= dirich_total;
  return(rmultinom(dirichlet_draw, n_eff));
}
/*
 * inverse centred log transform
 * Up to a constant so standardise
 */
template <class Type>
vector<Type> inv_crl(vector<Type>& y){
  return exp(y) / (exp(y)).sum();
}

// Beverton-Holt SR relationship function
template<class Type>
Type BevertonHolt(Type SSB, Type B0, Type h) {
  Type ssb_ratio = SSB / B0;
  Type part_2 = (1 - ((5*h - 1) / (4*h)) * ( 1 - ssb_ratio));
  return (ssb_ratio / part_2);
}

// Beverton-Holt SR relationship function without equilibrium assumptions
template<class Type>
Type BevertonHoltNoEquil(Type a, Type b, Type SSB) {
  Type Rt = (a + SSB) / (SSB * b);
  return Rt;
}

/*
 * Simplex functionality
 */
// simplex_jacobian - calculate the jacobian for estiamting the logistic simplex rather than the simplex
template <class Type>
Type simplex_jacobian(vector<Type> zk, vector<Type> X_orig) {
  std::cout << "Simplex Jacobian\n";
  Type Jacobian = 1.0;
  for (int k = 0; k < (zk.size() - 1); ++k) { // because I store zk in a matrix with xk in simplex_restore() the zk vector has one more null/empty element than it should have so we don't iterate over this.
    if (k == 0) {
      Jacobian *= (zk[k] * (1 - zk[k]));
    } else {
      Jacobian *= (zk[k] * (1 - zk[k])) * (1 - sum(vector<Type>(X_orig.segment(0,k))));
    }
  }
  return Jacobian;
}

// This class takes the transformed (logistic simplex) values, and their original unit vector and gives the
// recent update value.
// returns a matrix row 1 (C++ call: row(0)) is the zk - needed for the jacobian
// teh second row (C++ call: row(1)) is xk the untransformed values
template <class Type>
matrix<Type> Simplex_restore(vector<Type> logit_unit_vec, vector<Type> last_unit_vec) {
  std::cout << "Simplex_restore\n";
  int n_pars = last_unit_vec.size();
  matrix<Type> matrix_return(2,n_pars);
  Type total = 0.0;
  for (int k = 0; k < (logit_unit_vec.size()); ++k) {
    matrix_return(0,k) = invlogit(logit_unit_vec[k] + log(Type(1) / ((Type(logit_unit_vec.size()) + Type(1)) - Type(k + 1))));

    if (k == 0)
      matrix_return(1,k) = matrix_return(0,k);
    else
      matrix_return(1,k) = (Type(1) - sum(vector<Type>(last_unit_vec.segment(0,k)))) * matrix_return(0,k);
    total += matrix_return(1,k);
    //std::cerr << "k = " << k + 1 << " zk = " << matrix_return(0,k) << " xk = " << matrix_return(1,k) << std::endl; // Output warning
  }
  //std::cerr << vector<Type>(matrix_return.row(1)) << std::endl;
  matrix_return(1,n_pars - 1) = 1 - total;
  //std::cerr << vector<Type>(matrix_return.row(1)) << std::endl;

  return matrix_return;
}

