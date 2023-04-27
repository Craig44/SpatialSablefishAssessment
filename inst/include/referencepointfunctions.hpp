// Calculate spawner per recruit
/*
 * @param F F to test
 * @param sel vector of ages selectivity
 * @param natural_mortality natural mortality
 * @param waa vector of weights at age
 * @param paa vector of proportions at age
 * @param propZ_ssb scalar interpolating SSB
 * @param ages vector of ages
 * @return SPR
 */

template<class Type>
Type get_SPR(Type F, vector<Type>& sel, Type& natural_mortality, vector<Type>& waa, vector<Type>& paa, Type propZ_ssb, vector<Type>& ages) {
  Type Na = 1;
  vector<Type> Za = natural_mortality + sel * F;
  vector<Type> Sa = exp(-Za);
  Type SPR = Na * exp(-Za(0) * propZ_ssb) * waa(0) * paa(0);
  for(unsigned age_iter = 1; age_iter < (ages.size() * 4); ++age_iter) {
    if(age_iter >= ages.size() ) {
      SPR += Na * exp(-Za(ages.size() - 1) * propZ_ssb) * waa(ages.size() - 1) * paa(ages.size() - 1);
      Na *= Sa(ages.size() - 1);
    } else {
      SPR += Na * exp(-Za(age_iter) * propZ_ssb) * waa(age_iter) * paa(age_iter);
      Na *= Sa(age_iter);
    }
  }
  return SPR;
}

// Calculate Yeild per recruit
/*
 * @param F F to test
 * @param sel vector of ages selectivity
 * @param natural_mortality natural mortality
 * @param waa vector of weights at age
 * @param ages vector of ages
 * @return SPR
 */
template<class Type>
Type get_YPR(Type F, vector<Type>& sel, Type& natural_mortality, vector<Type>& waa, vector<Type>& ages) {
  Type Na = 1;
  Type YPR = (sel(0) * F) / (natural_mortality + sel(0) * F) * (Type(1) - exp(-(natural_mortality + sel(0) * F))) * Na * waa(0);
  for(unsigned age_iter = 1; age_iter < ages.size() * 4; ++age_iter) {
    if(age_iter >= ages.size() ) {
      Na *= exp(-(natural_mortality + sel(ages.size() - 1) * F));
      YPR += (sel(ages.size() - 1) * F) / (natural_mortality + sel(ages.size() - 1) * F) * (Type(1) - exp(-(natural_mortality + sel(ages.size() - 1) * F))) * Na * waa(ages.size() - 1);
    } else {
      Na *= exp(-(natural_mortality + sel(age_iter) * F));
      YPR += (sel(age_iter) * F) / (natural_mortality + sel(age_iter) * F) * (Type(1) - exp(-(natural_mortality + sel(age_iter) * F))) * Na * waa(age_iter);

    }
  }
  return YPR;
}
// Calculate derivative for Yeild per recruit
/*
 * @param F F to test
 * @param sel vector of ages selectivity
 * @param natural_mortality natural mortality
 * @param waa vector of weights at age
 * @param ages vector of ages
 * @return derivative YPR
 */
template<class Type>
Type get_dYPR(Type F, vector<Type>& sel, Type& natural_mortality, vector<Type>& waa, vector<Type>& ages) {
  Type ln_F = log(F);
  Type h = 0.001;
  Type v = -get_YPR(exp(ln_F + 2.0 * h), sel, natural_mortality, waa, ages) + 8.0 * get_YPR(exp(ln_F + h), sel, natural_mortality, waa, ages) - 8.0 * get_YPR(exp(ln_F - h), sel, natural_mortality, waa, ages) + get_YPR(exp(ln_F - 2.0 * h), sel, natural_mortality, waa, ages);
  Type g = v / (12.0 * h);
  return(g / F);
}
// Calculate Equilbrium SSB for a given F
/*
 * @param F F to test
 * @param sel vector of ages selectivity
 * @param natural_mortality natural mortality
 * @param waa vector of weights at age
 * @param ages vector of ages
 * @return SPR
 */
template<class Type>
Type get_SSBe(Type& F, vector<Type>& N_init, vector<Type>& sel, Type& natural_mortality, vector<Type>& waa, vector<Type>& ages, vector<Type>& paa, Type propZ_ssb, int n_runs, int maxAgePlusGroup) {

  vector<Type> Fa = F * sel;
  vector<Type> Za = Fa + natural_mortality;
  vector<Type> Sa = exp(-Za);
  array<Type> N_equilibrium(ages.size(), n_runs);
  N_equilibrium.col(0) = N_init;

  // Run annual cycle
  for(int year_ndx = 1; year_ndx < n_runs; ++year_ndx) {
    // Initial recruitment
    N_equilibrium(0, year_ndx) = N_init(0);
    // Ageing + Z
    for(int age_ndx = 1; age_ndx < ages.size() ; ++age_ndx)
      N_equilibrium(age_ndx, year_ndx) = N_equilibrium(age_ndx - 1, year_ndx - 1) * Sa(age_ndx - 1);
    if(maxAgePlusGroup == 1) {
      N_equilibrium(ages.size() - 1, year_ndx) = N_equilibrium(ages.size()  - 2, year_ndx - 1) * Sa(ages.size()  - 2) +
        N_equilibrium(ages.size()  - 1, year_ndx - 1) * Sa(ages.size()  - 1);
    }
  }

  // Calculate SSBs an interpolation bewtween the year, starting with previous years Paritition
  Type ssbe = 0;
  for(int age_ndx = 0; age_ndx < ages.size(); ++age_ndx)
    ssbe += N_equilibrium(age_ndx, n_runs - 1) * exp(-Za(age_ndx) * propZ_ssb) * paa(age_ndx) * waa(age_ndx);

  return(ssbe);
}
