
abs(Z_1-Z_2)<0.001 && Z_1>0 && Z_2>0

LOGMASS-log10(Kroupa_MILES_stellar_mass)

log10(AGE*1e9)-log10(Kroupa_MILES_age_lightW)


LOGMASS-log10(Kroupa_MILES_stellar_mass)

abs(pow(10, LOGMASS)-Kroupa_MILES_stellar_mass) /sqrt(pow((Kroupa_MILES_stellar_mass_up - Kroupa_MILES_stellar_mass_low)/2,2)+ pow((pow(10,MAXLOGMASS)-pow(10,MINLOGMASS))/2.,2))


(Kroupa_MILES_stellar_mass_up - Kroupa_MILES_stellar_mass_low)/2 / (pow(10,MAXLOGMASS)-pow(10,MINLOGMASS))/2.

(pow(10, LOGMASS)-Kroupa_MILES_stellar_mass) /sqrt(pow((Kroupa_MILES_stellar_mass_up - Kroupa_MILES_stellar_mass_low)/2,2)+ pow((pow(10,MAXLOGMASS)-pow(10,MINLOGMASS))/2.,2)) 

(pow(10, LOGMASS)-Kroupa_STELIB_stellar_mass) /sqrt(pow((Kroupa_STELIB_stellar_mass_up - Kroupa_STELIB_stellar_mass_low)/2,2)+ pow((pow(10,MAXLOGMASS)-pow(10,MINLOGMASS))/2.,2)) 

_3 && abs(E_BV-Kroupa_MILES_spm_EBV)<0.02
_3 && abs(E_BV-Kroupa_STELIB_spm_EBV)<0.02

exp(-x*x/2.)*330000

(Kroupa_STELIB_stellar_mass_up - Kroupa_MILES_stellar_mass_low)/Kroupa_MILES_stellar_mass


