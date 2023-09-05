#ifndef sdeModel_h
#define sdeModel_h 1

// Lotka-Volterra Predator-Prey model

// class definition
class sdeModel {
  public:
    static const int nParams = 13; // number of model parameters
    static const int nDims = 3; // number of sde dimensions
    // static const bool diagDiff = false; // whether diffusion function is diagonal
    static const bool diagDiff = true; // whether diffusion function is diagonal
    static const bool sdDiff = true; // whether diffusion is on sd or var scale
    // static const bool sdDiff = false; // whether diffusion is on sd or var scale
    void sdeDr(double *dr, double *x, double *theta); // drift function
    void sdeDf(double *df, double *x, double *theta); // diffusion function
    bool isValidParams(double *theta); // parameter validator
    bool isValidData(double *x, double *theta); // data validator
};

// drift function
inline void sdeModel::sdeDr(double *dr, double *x, double *theta) {
  dr[0] = x[0]*(theta[0]+theta[3]*x[0]+theta[4]*x[1]+theta[5]*x[2]);
  dr[1] = x[1]*(theta[1]+theta[6]*x[0]+theta[7]*x[1]+theta[8]*x[2]);
  dr[2] = x[2]*(theta[2]+theta[9]*x[0]+theta[10]*x[1]+theta[11]*x[2]);
    
  return;
}


inline void sdeModel::sdeDf(double *df, double *x, double *theta) {
  df[0] = theta[12];
  df[1] = theta[12];
  df[2] = theta[12];
  return;
}

// parameter validator
inline bool sdeModel::isValidParams(double *theta) {
  bool val = theta[12] > 0.0;
  return val;
}

// data validator
inline bool sdeModel::isValidData(double *x, double *theta) {
  return (x[0] >= 0.0) && (x[1] >= 0.0) && (x[2] >= 0.0);
}

#endif
