// Initially generated by TMBtools: modify with care!!

#define TMB_LIB_INIT R_init_SpatialSablefishAssessment_TMBExports

#include <TMB.hpp>

#include <AuxillaryFuns.hpp>

#include <SetupSelectivities.hpp>

#include "TagIntegratedValidate.hpp"

#include "TagIntegrated.hpp"

#include "TagIntegratedAgeBasedMovement.hpp"

#include "CurrentAssessment.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "TagIntegrated") {
    return TagIntegrated(this);
  } else if(model == "TagIntegratedValidate") {
    return TagIntegratedValidate(this);
  } else if(model == "TagIntegratedAgeBasedMovement") {
    return TagIntegratedAgeBasedMovement(this);
  } else if(model == "Assessment") {
    return CurrentAssessment(this);
  } else {
    Rf_error("Unknown model.");
  }
  return 0;
}
