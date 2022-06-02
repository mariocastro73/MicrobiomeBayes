#ifndef sdeModel_h
#define sdeModel_h 1

// Stochastic Logistic model

// class definition
class sdeModel {
  public:
    static const int nParams = 3; // number of model parameters
    static const int nDims = 1; // number of sde dimensions
    static const bool diagDiff = true; // whether diffusion function is diagonal
    // static const bool sdDiff = true; // whether diffusion is on sd or var scale
    static const bool sdDiff = false; // whether diffusion is on sd or var scale
    void sdeDr(double *dr, double *x, double *theta); // drift function
    void sdeDf(double *df, double *x, double *theta); // diffusion function
    bool isValidParams(double *theta); // parameter validator
    bool isValidData(double *x, double *theta); // data validator
};

// drift function
inline void sdeModel::sdeDr(double *dr, double *x, double *theta) {
  dr[0] = theta[0]*x[0]*(1-x[0]/theta[1]); // r*x(1-x/K) 
  return;
}

inline void sdeModel::sdeDf(double *df, double *x, double *theta) {
  df[0] = theta[2]*x[0]; // matrix element (1,1)
  return;
}

// parameter validator
inline bool sdeModel::isValidParams(double *theta) {
  bool val = theta[0] > 0.0;
  val = val && theta[1] > 0.0;
  val = val && theta[2] > 0.0;
  return val;
}

// data validator
inline bool sdeModel::isValidData(double *x, double *theta) {
  return (x[0] >= 0.0);
}

#endif