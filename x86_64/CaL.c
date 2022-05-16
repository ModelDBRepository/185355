/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__lca
#define _nrn_initial _nrn_initial__lca
#define nrn_cur _nrn_cur__lca
#define _nrn_current _nrn_current__lca
#define nrn_jacob _nrn_jacob__lca
#define nrn_state _nrn_state__lca
#define _net_receive _net_receive__lca 
#define _f_trates _f_trates__lca 
#define rates rates__lca 
#define states states__lca 
#define trates trates__lca 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define glcabar _p[0]
#define glcabar_columnindex 0
#define glca _p[1]
#define glca_columnindex 1
#define ilca _p[2]
#define ilca_columnindex 2
#define einf _p[3]
#define einf_columnindex 3
#define etau _p[4]
#define etau_columnindex 4
#define e _p[5]
#define e_columnindex 5
#define De _p[6]
#define De_columnindex 6
#define elca _p[7]
#define elca_columnindex 7
#define eexp _p[8]
#define eexp_columnindex 8
#define _g _p[9]
#define _g_columnindex 9
#define _ion_elca	*_ppvar[0]._pval
#define _ion_ilca	*_ppvar[1]._pval
#define _ion_dilcadv	*_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_rates(void);
 static void _hoc_states(void);
 static void _hoc_trates(void);
 static void _hoc_vtrap(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_lca", _hoc_setdata,
 "rates_lca", _hoc_rates,
 "states_lca", _hoc_states,
 "trates_lca", _hoc_trates,
 "vtrap_lca", _hoc_vtrap,
 0, 0
};
#define vtrap vtrap_lca
 extern double vtrap( double , double );
 /* declare global and static user variables */
#define usetable usetable_lca
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_lca", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "glcabar_lca", "mho/cm2",
 "glca_lca", "mho/cm2",
 "ilca_lca", "mA/cm2",
 "etau_lca", "ms",
 0,0
};
 static double delta_t = 1;
 static double e0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "usetable_lca", &usetable_lca,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"lca",
 "glcabar_lca",
 0,
 "glca_lca",
 "ilca_lca",
 "einf_lca",
 "etau_lca",
 0,
 "e_lca",
 0,
 0};
 static Symbol* _lca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 10, _prop);
 	/*initialize range parameters*/
 	glcabar = 0;
 	_prop->param = _p;
 	_prop->param_size = 10;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_lca_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* elca */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ilca */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dilcadv */
 
}
 static void _initlists();
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _CaL_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("lca", 2.0);
 	_lca_sym = hoc_lookup("lca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 10, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "lca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "lca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "lca_ion");
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 lca /home/hines/tmp/hack/185355/CaL.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96520.0;
 static double R = 8.3134;
 static double _zq10 ;
 static double *_t_einf;
 static double *_t_eexp;
 static double *_t_etau;
static int _reset;
static char *modelname = "CaL.mod L-type Cav channel";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_trates(double);
static int rates(double);
static int states();
static int trates(double);
 static void _n_trates(double);
 
static int  states (  ) {
   trates ( _threadargscomma_ v ) ;
   e = e + eexp * ( einf - e ) ;
   
/*VERBATIM*/
        return 0;
  return 0; }
 
static void _hoc_states(void) {
  double _r;
   _r = 1.;
 states (  );
 hoc_retpushx(_r);
}
 
static int  rates (  double _lv ) {
   double _lalpha , _lbeta , _lsum ;
 _zq10 = pow( 3.0 , ( ( celsius - 6.3 ) / 10.0 ) ) ;
   _lalpha = - 15.69 * vtrap ( _threadargscomma_ _lv - 81.5 , - 10.0 ) ;
   _lbeta = 0.29 * exp ( - _lv / 10.86 ) ;
   _lsum = _lalpha + _lbeta ;
   etau = 1.0 / _lsum ;
   einf = _lalpha / _lsum ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   _r = 1.;
 rates (  *getarg(1) );
 hoc_retpushx(_r);
}
 static double _mfac_trates, _tmin_trates;
 static void _check_trates();
 static void _check_trates() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_dt;
  static double _sav_celsius;
  if (!usetable) {return;}
  if (_sav_dt != dt) { _maktable = 1;}
  if (_sav_celsius != celsius) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_trates =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_trates)/200.; _mfac_trates = 1./_dx;
   for (_i=0, _x=_tmin_trates; _i < 201; _x += _dx, _i++) {
    _f_trates(_x);
    _t_einf[_i] = einf;
    _t_eexp[_i] = eexp;
    _t_etau[_i] = etau;
   }
   _sav_dt = dt;
   _sav_celsius = celsius;
  }
 }

 static int trates(double _lv){ _check_trates();
 _n_trates(_lv);
 return 0;
 }

 static void _n_trates(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_trates(_lv); return; 
}
 _xi = _mfac_trates * (_lv - _tmin_trates);
 if (isnan(_xi)) {
  einf = _xi;
  eexp = _xi;
  etau = _xi;
  return;
 }
 if (_xi <= 0.) {
 einf = _t_einf[0];
 eexp = _t_eexp[0];
 etau = _t_etau[0];
 return; }
 if (_xi >= 200.) {
 einf = _t_einf[200];
 eexp = _t_eexp[200];
 etau = _t_etau[200];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 einf = _t_einf[_i] + _theta*(_t_einf[_i+1] - _t_einf[_i]);
 eexp = _t_eexp[_i] + _theta*(_t_eexp[_i+1] - _t_eexp[_i]);
 etau = _t_etau[_i] + _theta*(_t_etau[_i+1] - _t_etau[_i]);
 }

 
static int  _f_trates (  double _lv ) {
   double _ltinc ;
 rates ( _threadargscomma_ _lv ) ;
   _ltinc = - dt * _zq10 ;
   eexp = 1.0 - exp ( _ltinc / etau ) ;
    return 0; }
 
static void _hoc_trates(void) {
  double _r;
    _r = 1.;
 trates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap (  double _lx , double _ly ) {
   double _lvtrap;
 if ( fabs ( _lx / _ly ) < 1e-6 ) {
     _lvtrap = _ly * ( 1.0 - _lx / _ly / 2.0 ) ;
     }
   else {
     _lvtrap = _lx / ( exp ( _lx / _ly ) - 1.0 ) ;
     }
   
return _lvtrap;
 }
 
static void _hoc_vtrap(void) {
  double _r;
   _r =  vtrap (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ hoc_execerror("lca", "cannot be used with CVODE"); return 0;}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_lca_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_lca_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_lca_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  e = e0;
 {
   trates ( _threadargscomma_ v ) ;
   e = einf ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  elca = _ion_elca;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   glca = glcabar * e * e ;
   ilca = glca * ( v - elca ) ;
   }
 _current += ilca;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  elca = _ion_elca;
 _g = _nrn_current(_v + .001);
 	{ double _dilca;
  _dilca = ilca;
 _rhs = _nrn_current(_v);
  _ion_dilcadv += (_dilca - ilca)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ilca += ilca ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  elca = _ion_elca;
 { error =  states();
 if(error){fprintf(stderr,"at line 60 in file CaL.mod:\n	SOLVE states\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
   _t_einf = makevector(201*sizeof(double));
   _t_eexp = makevector(201*sizeof(double));
   _t_etau = makevector(201*sizeof(double));
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/home/hines/tmp/hack/185355/CaL.mod";
static const char* nmodl_file_text = 
  "TITLE CaL.mod L-type Cav channel\n"
  "COMMENT\n"
  "\n"
  "Mod file by A. Hanuschkin <AH, 2011> for:\n"
  "Yim MY, Hanuschkin A, Wolfart J (2015) Hippocampus 25:297-308.\n"
  "http://onlinelibrary.wiley.com/doi/10.1002/hipo.22373/abstract\n"
  "\n"
  "Mod File history:\n"
  "- fitted H-H parameter N-Ca from Jaffe DB, Ross WN, Lisman JE,  Lasser-Ross N, Miyakawa H, Johnston D (1994) Journal of Neurophysiology, Vol. 71 no. 3, 1065-1077\n"
  "- Ca ion & L/T/N-Ca channels model of  Aradi I, Holmes WR (1999) J Comput Neurosci 6:215-35\n"
  "- checked and adapted by Hanuschkin in 2011\n"
  "- Note that eCa is calculated during simulation by ccanl.mod. ecat, ecal values set in Santhakumar are not used in our model scripts.\n"
  " \n"
  "ENDCOMMENT\n"
  " \n"
  "UNITS {\n"
  "        (mA) =		(milliamp)\n"
  "        (mV) =		(millivolt)\n"
  "        (uF) = 		(microfarad)\n"
  "	(molar) = 	(1/liter)\n"
  "	(nA) = 		(nanoamp)\n"
  "	(mM) = 		(millimolar)\n"
  "	(um) = 		(micron)\n"
  "	FARADAY = 96520 (coul)\n"
  "	R = 8.3134	(joule/degC)\n"
  "}\n"
  " \n"
  "NEURON { \n"
  "SUFFIX lca\n"
  "USEION lca READ elca WRITE ilca VALENCE 2 \n"
  "RANGE  glca\n"
  "RANGE glcabar\n"
  "RANGE einf, etau, ilca\n"
  "}\n"
  " \n"
  "INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}\n"
  " \n"
  "PARAMETER {\n"
  "        v (mV) \n"
  "        celsius = 6.3 (degC)\n"
  "        dt (ms) \n"
  "	glcabar (mho/cm2)\n"
  "}\n"
  " \n"
  "STATE {\n"
  "	e\n"
  "}\n"
  " \n"
  "ASSIGNED {\n"
  "	glca (mho/cm2)\n"
  "	ilca (mA/cm2)\n"
  "	elca (mV)\n"
  "\n"
  "	einf \n"
  "	etau (ms) \n"
  "	eexp      \n"
  "} \n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states\n"
  "        glca = glcabar*e*e\n"
  "	ilca = glca*(v-elca)\n"
  "}\n"
  " \n"
  "UNITSOFF\n"
  " \n"
  "INITIAL {\n"
  "	trates(v)\n"
  "	e = einf\n"
  "}\n"
  "\n"
  "PROCEDURE states() {	:Computes state variables e \n"
  "        trates(v)	:      at the current v and dt.\n"
  "	e = e + eexp*(einf-e)\n"
  "        VERBATIM\n"
  "        return 0;\n"
  "        ENDVERBATIM\n"
  "}\n"
  " \n"
  "LOCAL q10\n"
  "\n"
  "PROCEDURE rates(v) {  :Computes rate and other constants at current v.\n"
  "                      :Call once from HOC to initialize inf at resting v.\n"
  "        LOCAL  alpha, beta, sum\n"
  "        q10 = 3^((celsius - 6.3)/10) : q10=1 for 6.3 celcius\n"
  "                :\"e\" LCa activation system\n"
  "        alpha = -15.69*vtrap(v-81.5,-10)	\n"
  "	beta = 0.29*exp(-v/10.86)	\n"
  "	sum = alpha+beta        \n"
  "	etau = 1/sum      einf = alpha/sum\n"
  "                :no LCa inactivation system\n"
  "}\n"
  " \n"
  "PROCEDURE trates(v) {  :Computes rate and other constants at current v.\n"
  "                      :Call once from HOC to initialize inf at resting v.\n"
  "	LOCAL tinc\n"
  "        TABLE  einf, eexp, etau\n"
  "	DEPEND dt, celsius FROM -100 TO 100 WITH 200\n"
  "                           \n"
  "	rates(v)	: not consistently executed from here if usetable_hh == 1\n"
  "		: so don't expect the tau values to be tracking along with\n"
  "		: the inf values in hoc\n"
  "\n"
  "	tinc = -dt * q10\n"
  "	eexp = 1 - exp(tinc/etau)\n"
  "}\n"
  " \n"
  "FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.\n"
  "        if (fabs(x/y) < 1e-6) {\n"
  "                vtrap = y*(1 - x/y/2)\n"
  "        }else{  \n"
  "                vtrap = x/(exp(x/y) - 1)\n"
  "        }\n"
  "}\n"
  " \n"
  "UNITSON\n"
  "\n"
  ;
#endif
