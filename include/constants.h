#ifndef CONSTANTS_H
#define CONSTANTS_H

/* Gravitational constant [m^3/kg/s^2] */
#define CONST_G         6.67e-11

/* Ceres mass [kg] */
#define CONST_M_CERES   9.39e20

/* Ceres radius [m] */
#define CONST_R_CERES   476.2e3

/* Boltzmann constant [J/K] */
#define CONST_K_B       1.38e-23

/* Avogadro number [1/mol] */
#define CONST_MOLE      6.02e23

/* H2O molar mass [kg/mol] */
#define CONST_MMASS_H2O 18e-3

/* H2O molecular mass [kg] */
#define CONST_MMOL_H2O  (CONST_MMASS_H2O / CONST_MOLE)

/* Escape velocity at surface [m/s] */
#define CONST_VE_CERES  (sqrt(2.0 * CONST_G * CONST_M_CERES / CONST_R_CERES))

/* Rotation period [s] (9.07417 hours) */
#define CONST_P_CERES   (9.07417 * 3600.0)

/* Photodissociation lifetime [s] at 2.8 AU */
#define CONST_LIFETIME  (1.0 / 1.02e-5 * (2.8 * 2.8))

/* Collision cross-section [m^2] */
#define CONST_CROSS_SECTION 1e-19

/* Stefan-Boltzmann constant [W/m^2/K^4] */
#define CONST_SIGMA     5.6704e-8

/* Solar constant at 1 AU [W/m^2] */
#define CONST_SOLAR_1AU 1368.0

/* Pi */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#endif /* CONSTANTS_H */
