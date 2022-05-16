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
 
#define nrn_init _nrn_init__tca
#define _nrn_initial _nrn_initial__tca
#define nrn_cur _nrn_cur__tca
#define _nrn_current _nrn_current__tca
#define nrn_jacob _nrn_jacob__tca
#define nrn_state _nrn_state__tca
#define _net_receive _net_receive__tca 
#define _f_trates _f_trates__tca 
#define rates rates__tca 
#define states states__tca 
#define trates trates__tca 
 
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
#define gcatbar _p[0]
#define gcatbar_columnindex 0
#define gtca _p[1]
#define gtca_columnindex 1
#define itca _p[2]
#define itca_columnindex 2
#define ainf _p[3]
#define ainf_columnindex 3
#define binf _p[4]
#define binf_columnindex 4
#define atau _p[5]
#define atau_columnindex 5
#define btau _p[6]
#define btau_columnindex 6
#define a _p[7]
#define a_columnindex 7
#define b _p[8]
#define b_columnindex 8
#define Da _p[9]
#define Da_columnindex 9
#define Db _p[10]
#define Db_columnindex 10
#define etca _p[11]
#define etca_columnindex 11
#define aexp _p[12]
#define aexp_columnindex 12
#define bexp _p[13]
#define bexp_columnindex 13
#define _g _p[14]
#define _g_columnindex 14
#define _ion_etca	*_ppvar[0]._pval
#define _ion_itca	*_ppvar[1]._pval
#define _ion_ditcadv	*_ppvar[2]._pval
 
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
 "setdata_tca", _hoc_setdata,
 "rates_tca", _hoc_rates,
 "states_tca", _hoc_states,
 "trates_tca", _hoc_trates,
 "vtrap_tca", _hoc_vtrap,
 0, 0
};
#define vtrap vtrap_tca
 extern double vtrap( double , double );
 /* declare global and static user variables */
#define usetable usetable_tca
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_tca", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gcatbar_tca", "mho/cm2",
 "gtca_tca", "mho/cm2",
 "itca_tca", "mA/cm2",
 "atau_tca", "ms",
 "btau_tca", "ms",
 0,0
};
 static double a0 = 0;
 static double b0 = 0;
 static double delta_t = 1;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "usetable_tca", &usetable_tca,
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
"tca",
 "gcatbar_tca",
 0,
 "gtca_tca",
 "itca_tca",
 "ainf_tca",
 "binf_tca",
 "atau_tca",
 "btau_tca",
 0,
 "a_tca",
 "b_tca",
 0,
 0};
 static Symbol* _tca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 15, _prop);
 	/*initialize range parameters*/
 	gcatbar = 0;
 	_prop->param = _p;
 	_prop->param_size = 15;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_tca_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* etca */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* itca */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_ditcadv */
 
}
 static void _initlists();
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _CaT_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("tca", 2.0);
 	_tca_sym = hoc_lookup("tca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 15, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "tca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "tca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "tca_ion");
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 tca /home/hines/tmp/hack/185355/CaT.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96520.0;
 static double R = 8.3134;
 static double _zq10 ;
 static double *_t_ainf;
 static double *_t_aexp;
 static double *_t_binf;
 static double *_t_bexp;
 static double *_t_atau;
 static double *_t_btau;
static int _reset;
static char *modelname = "CaT.mod T-type Cav channel";

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
   a = a + aexp * ( ainf - a ) ;
   b = b + bexp * ( binf - b ) ;
   
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
   _lalpha = - 0.2 * vtrap ( _threadargscomma_ _lv - 19.26 , - 10.0 ) ;
   _lbeta = 0.009 * exp ( - _lv / 22.03 ) ;
   _lsum = _lalpha + _lbeta ;
   atau = 1.0 / _lsum ;
   ainf = _lalpha / _lsum ;
   _lalpha = 1e-6 * exp ( - _lv / 16.26 ) ;
   _lbeta = 1.0 / ( exp ( ( 29.79 - _lv ) / 10.0 ) + 1.0 ) ;
   _lsum = _lalpha + _lbeta ;
   btau = 1.0 / _lsum ;
   binf = _lalpha / _lsum ;
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
    _t_ainf[_i] = ainf;
    _t_aexp[_i] = aexp;
    _t_binf[_i] = binf;
    _t_bexp[_i] = bexp;
    _t_atau[_i] = atau;
    _t_btau[_i] = btau;
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
  ainf = _xi;
  aexp = _xi;
  binf = _xi;
  bexp = _xi;
  atau = _xi;
  btau = _xi;
  return;
 }
 if (_xi <= 0.) {
 ainf = _t_ainf[0];
 aexp = _t_aexp[0];
 binf = _t_binf[0];
 bexp = _t_bexp[0];
 atau = _t_atau[0];
 btau = _t_btau[0];
 return; }
 if (_xi >= 200.) {
 ainf = _t_ainf[200];
 aexp = _t_aexp[200];
 binf = _t_binf[200];
 bexp = _t_bexp[200];
 atau = _t_atau[200];
 btau = _t_btau[200];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 ainf = _t_ainf[_i] + _theta*(_t_ainf[_i+1] - _t_ainf[_i]);
 aexp = _t_aexp[_i] + _theta*(_t_aexp[_i+1] - _t_aexp[_i]);
 binf = _t_binf[_i] + _theta*(_t_binf[_i+1] - _t_binf[_i]);
 bexp = _t_bexp[_i] + _theta*(_t_bexp[_i+1] - _t_bexp[_i]);
 atau = _t_atau[_i] + _theta*(_t_atau[_i+1] - _t_atau[_i]);
 btau = _t_btau[_i] + _theta*(_t_btau[_i+1] - _t_btau[_i]);
 }

 
static int  _f_trates (  double _lv ) {
   double _ltinc ;
 rates ( _threadargscomma_ _lv ) ;
   _ltinc = - dt * _zq10 ;
   aexp = 1.0 - exp ( _ltinc / atau ) ;
   bexp = 1.0 - exp ( _ltinc / btau ) ;
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
 
static int _ode_count(int _type){ hoc_execerror("tca", "cannot be used with CVODE"); return 0;}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_tca_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_tca_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_tca_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  a = a0;
  b = b0;
 {
   trates ( _threadargscomma_ v ) ;
   a = ainf ;
   b = binf ;
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
  etca = _ion_etca;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   gtca = gcatbar * a * a * b ;
   itca = gtca * ( v - etca ) ;
   }
 _current += itca;

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
  etca = _ion_etca;
 _g = _nrn_current(_v + .001);
 	{ double _ditca;
  _ditca = itca;
 _rhs = _nrn_current(_v);
  _ion_ditcadv += (_ditca - itca)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_itca += itca ;
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
  etca = _ion_etca;
 { error =  states();
 if(error){fprintf(stderr,"at line 59 in file CaT.mod:\n	SOLVE states\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
   _t_ainf = makevector(201*sizeof(double));
   _t_aexp = makevector(201*sizeof(double));
   _t_binf = makevector(201*sizeof(double));
   _t_bexp = makevector(201*sizeof(double));
   _t_atau = makevector(201*sizeof(double));
   _t_btau = makevector(201*sizeof(double));
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/home/hines/tmp/hack/185355/CaT.mod";
static const char* nmodl_file_text = 
  "TITLE CaT.mod T-type Cav channel\n"
  "COMMENT\n"
  "\n"
  "Mod File by A. Hanuschkin <AH, 2011> for:\n"
  "Yim MY, Hanuschkin A, Wolfart J (2015) Hippocampus 25:297-308.\n"
  "http://onlinelibrary.wiley.com/doi/10.1002/hipo.22373/abstract\n"
  "\n"
  "Mod File history:\n"
  "- fitted H-H parameter N-Ca from Jaffe DB, Ross WN, Lisman JE,  Lasser-Ross N, Miyakawa H, Johnston D (1994) Journal of Neurophysiology, Vol. 71 no. 3, 1065-1077\n"
  "- Ca ion & L/T/N-Ca channels model of  Aradi I, Holmes WR (1999) J Comput Neurosci 6:215-35\n"
  "- Note that eCa is calculated during simulation by ccanl.mod. ecat, ecal values set in Santhakumar are not used in our model scripts.\n"
  "\n"
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
  "SUFFIX tca\n"
  "USEION tca READ etca WRITE itca VALENCE 2 \n"
  "RANGE gtca\n"
  "RANGE gcatbar\n"
  "RANGE ainf, atau, binf, btau, itca\n"
  "}\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}\n"
  " \n"
  "PARAMETER {\n"
  "        v (mV) \n"
  "        celsius = 6.3 (degC)\n"
  "        dt (ms) \n"
  "	gcatbar (mho/cm2)\n"
  "}\n"
  " \n"
  "STATE {\n"
  "	a b\n"
  "}\n"
  " \n"
  "ASSIGNED {\n"
  "        gtca (mho/cm2)\n"
  "	itca (mA/cm2)\n"
  "	etca (mV)\n"
  "\n"
  "	ainf binf\n"
  "	atau (ms) btau (ms) \n"
  "	aexp bexp      \n"
  "} \n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states\n"
  "        gtca = gcatbar*a*a*b\n"
  "	itca = gtca*(v-etca)\n"
  "}\n"
  " \n"
  "UNITSOFF\n"
  " \n"
  "INITIAL {\n"
  "	trates(v)\n"
  "	a = ainf\n"
  "	b = binf\n"
  "}\n"
  "\n"
  "PROCEDURE states() {	:Computes state variables a and b \n"
  "        trates(v)	:      at the current v and dt.\n"
  "	a = a + aexp*(ainf-a) : i.e. a_{t+1} = a_t*exp(-dt/atau)+ainf*(1-exp(-dt/atau)); da/dt = 1/atau*(ainf-a)\n"
  "	b = b + bexp*(binf-b)\n"
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
  "        q10 = 3^((celsius - 6.3)/10) : q10 = 1 for 6.3 celsius\n"
  "                :\"a\" TCa activation system\n"
  "        alpha = -0.2*vtrap(v-19.26,-10)		\n"
  "	beta = 0.009*exp(-v/22.03)		\n"
  "	sum = alpha+beta        \n"
  "	atau = 1/sum      ainf = alpha/sum\n"
  "                :\"b\" TCa inactivation system\n"
  "	alpha = 1e-6*exp(-v/16.26)		\n"
  "	beta = 1/(exp((29.79-v)/10)+1)		\n"
  "	sum = alpha+beta        \n"
  "	btau = 1/sum      binf = alpha/sum\n"
  "}\n"
  " \n"
  "PROCEDURE trates(v) {  :Computes rate and other constants at current v.\n"
  "                      :Call once from HOC to initialize inf at resting v.\n"
  "	LOCAL tinc\n"
  "        TABLE  ainf, aexp, binf, bexp, atau, btau\n"
  "	DEPEND dt, celsius FROM -100 TO 100 WITH 200\n"
  "                           \n"
  "	rates(v)	: not consistently executed from here if usetable_hh == 1\n"
  "		: so don't expect the tau values to be tracking along with\n"
  "		: the inf values in hoc\n"
  "\n"
  "	       tinc = -dt * q10\n"
  "	aexp = 1 - exp(tinc/atau)\n"
  "	bexp = 1 - exp(tinc/btau)\n"
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
