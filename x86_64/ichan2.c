/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
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
 
#define nrn_init _nrn_init__ichan2
#define _nrn_initial _nrn_initial__ichan2
#define nrn_cur _nrn_cur__ichan2
#define _nrn_current _nrn_current__ichan2
#define nrn_jacob _nrn_jacob__ichan2
#define nrn_state _nrn_state__ichan2
#define _net_receive _net_receive__ichan2 
#define rates rates__ichan2 
#define states states__ichan2 
#define trates trates__ichan2 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gnatbar _p[0]
#define gnatbar_columnindex 0
#define gkfbar _p[1]
#define gkfbar_columnindex 1
#define gksbar _p[2]
#define gksbar_columnindex 2
#define gl _p[3]
#define gl_columnindex 3
#define el _p[4]
#define el_columnindex 4
#define ggabaa _p[5]
#define ggabaa_columnindex 5
#define egabaa _p[6]
#define egabaa_columnindex 6
#define ina _p[7]
#define ina_columnindex 7
#define ik _p[8]
#define ik_columnindex 8
#define il _p[9]
#define il_columnindex 9
#define igabaa _p[10]
#define igabaa_columnindex 10
#define m _p[11]
#define m_columnindex 11
#define h _p[12]
#define h_columnindex 12
#define nf _p[13]
#define nf_columnindex 13
#define ns _p[14]
#define ns_columnindex 14
#define ena _p[15]
#define ena_columnindex 15
#define ek _p[16]
#define ek_columnindex 16
#define Dm _p[17]
#define Dm_columnindex 17
#define Dh _p[18]
#define Dh_columnindex 18
#define Dnf _p[19]
#define Dnf_columnindex 19
#define Dns _p[20]
#define Dns_columnindex 20
#define gna _p[21]
#define gna_columnindex 21
#define gkf _p[22]
#define gkf_columnindex 22
#define gks _p[23]
#define gks_columnindex 23
#define minf _p[24]
#define minf_columnindex 24
#define hinf _p[25]
#define hinf_columnindex 25
#define nfinf _p[26]
#define nfinf_columnindex 26
#define nsinf _p[27]
#define nsinf_columnindex 27
#define mtau _p[28]
#define mtau_columnindex 28
#define htau _p[29]
#define htau_columnindex 29
#define nftau _p[30]
#define nftau_columnindex 30
#define nstau _p[31]
#define nstau_columnindex 31
#define mexp _p[32]
#define mexp_columnindex 32
#define hexp _p[33]
#define hexp_columnindex 33
#define nfexp _p[34]
#define nfexp_columnindex 34
#define nsexp _p[35]
#define nsexp_columnindex 35
#define q10 _p[36]
#define q10_columnindex 36
#define v _p[37]
#define v_columnindex 37
#define _g _p[38]
#define _g_columnindex 38
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
#define _ion_ek	*_ppvar[3]._pval
#define _ion_ik	*_ppvar[4]._pval
#define _ion_dikdv	*_ppvar[5]._pval
 
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
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
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
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_ichan2", _hoc_setdata,
 "rates_ichan2", _hoc_rates,
 "states_ichan2", _hoc_states,
 "trates_ichan2", _hoc_trates,
 "vtrap_ichan2", _hoc_vtrap,
 0, 0
};
#define vtrap vtrap_ichan2
 extern double vtrap( _threadargsprotocomma_ double , double );
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gnatbar_ichan2", "mho/cm2",
 "gkfbar_ichan2", "mho/cm2",
 "gksbar_ichan2", "mho/cm2",
 "gl_ichan2", "mho/cm2",
 "el_ichan2", "mV",
 "ggabaa_ichan2", "mho/cm2",
 "egabaa_ichan2", "mV",
 "ina_ichan2", "mA/cm2",
 "ik_ichan2", "mA/cm2",
 "il_ichan2", "mA/cm2",
 "igabaa_ichan2", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 static double ns0 = 0;
 static double nf0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
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
"ichan2",
 "gnatbar_ichan2",
 "gkfbar_ichan2",
 "gksbar_ichan2",
 "gl_ichan2",
 "el_ichan2",
 "ggabaa_ichan2",
 "egabaa_ichan2",
 0,
 "ina_ichan2",
 "ik_ichan2",
 "il_ichan2",
 "igabaa_ichan2",
 0,
 "m_ichan2",
 "h_ichan2",
 "nf_ichan2",
 "ns_ichan2",
 0,
 0};
 static Symbol* _na_sym;
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 39, _prop);
 	/*initialize range parameters*/
 	gnatbar = 0;
 	gkfbar = 0;
 	gksbar = 0;
 	gl = 0;
 	el = 0;
 	ggabaa = 0;
 	egabaa = 0;
 	_prop->param = _p;
 	_prop->param_size = 39;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 6, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[3]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[4]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[5]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static void _initlists();
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _ichan2_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("na", 1.0);
 	ion_reg("k", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 39, 6);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "k_ion");
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 ichan2 /home/hines/tmp/hack/185355/ichan2.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96520.0;
 static double R = 8.3134;
static int _reset;
static char *modelname = "ichan2.mod combination of Nav and Kv channels for Hodgin-Huxely-type action potential mechanism";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_threadargsprotocomma_ double);
static int states(_threadargsproto_);
static int trates(_threadargsprotocomma_ double);
 
static int  states ( _threadargsproto_ ) {
   trates ( _threadargscomma_ v ) ;
   m = m + mexp * ( minf - m ) ;
   h = h + hexp * ( hinf - h ) ;
   nf = nf + nfexp * ( nfinf - nf ) ;
   ns = ns + nsexp * ( nsinf - ns ) ;
    return 0; }
 
static void _hoc_states(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 states ( _p, _ppvar, _thread, _nt );
 hoc_retpushx(_r);
}
 
static int  rates ( _threadargsprotocomma_ double _lv ) {
   double _lalpha , _lbeta , _lsum ;
 q10 = pow( 3.0 , ( ( celsius - 6.3 ) / 10.0 ) ) ;
   _lalpha = - 0.3 * vtrap ( _threadargscomma_ ( _lv + 60.0 - 17.0 ) , - 5.0 ) ;
   _lbeta = 0.3 * vtrap ( _threadargscomma_ ( _lv + 60.0 - 45.0 ) , 5.0 ) ;
   _lsum = _lalpha + _lbeta ;
   mtau = 1.0 / _lsum ;
   minf = _lalpha / _lsum ;
   _lalpha = 0.23 / exp ( ( _lv + 60.0 + 5.0 ) / 20.0 ) ;
   _lbeta = 3.33 / ( 1.0 + exp ( ( _lv + 60.0 - 47.5 ) / - 10.0 ) ) ;
   _lsum = _lalpha + _lbeta ;
   htau = 1.0 / _lsum ;
   hinf = _lalpha / _lsum ;
   _lalpha = - 0.028 * vtrap ( _threadargscomma_ ( _lv + 65.0 - 35.0 ) , - 6.0 ) ;
   _lbeta = 0.1056 / exp ( ( _lv + 65.0 - 10.0 ) / 40.0 ) ;
   _lsum = _lalpha + _lbeta ;
   nstau = 1.0 / _lsum ;
   nsinf = _lalpha / _lsum ;
   _lalpha = - 0.07 * vtrap ( _threadargscomma_ ( _lv + 65.0 - 47.0 ) , - 6.0 ) ;
   _lbeta = 0.264 / exp ( ( _lv + 65.0 - 22.0 ) / 40.0 ) ;
   _lsum = _lalpha + _lbeta ;
   nftau = 1.0 / _lsum ;
   nfinf = _lalpha / _lsum ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 rates ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int  trates ( _threadargsprotocomma_ double _lv ) {
   double _ltinc ;
 rates ( _threadargscomma_ _lv ) ;
   _ltinc = - dt * q10 ;
   mexp = 1.0 - exp ( _ltinc / mtau ) ;
   hexp = 1.0 - exp ( _ltinc / htau ) ;
   nfexp = 1.0 - exp ( _ltinc / nftau ) ;
   nsexp = 1.0 - exp ( _ltinc / nstau ) ;
    return 0; }
 
static void _hoc_trates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 trates ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap ( _threadargsprotocomma_ double _lx , double _ly ) {
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
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  vtrap ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ hoc_execerror("ichan2", "cannot be used with CVODE"); return 0;}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_k_sym, _ppvar, 3, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 4, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 5, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  h = h0;
  m = m0;
  ns = ns0;
  nf = nf0;
 {
   trates ( _threadargscomma_ v ) ;
   m = minf ;
   h = hinf ;
   nf = nfinf ;
   ns = nsinf ;
   }
 
}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
  ena = _ion_ena;
  ek = _ion_ek;
 initmodel(_p, _ppvar, _thread, _nt);
  }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   gna = gnatbar * m * m * m * h ;
   gkf = gkfbar * nf * nf * nf * nf ;
   gks = gksbar * ns * ns * ns * ns ;
   ina = gna * ( v - ena ) ;
   ik = gkf * ( v - ek ) + gks * ( v - ek ) ;
   il = gl * ( v - el ) ;
   igabaa = ggabaa * ( v - egabaa ) ;
   }
 _current += ina;
 _current += ik;
 _current += il;
 _current += igabaa;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
  ena = _ion_ena;
  ek = _ion_ek;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dik;
 double _dina;
  _dina = ina;
  _dik = ik;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dinadv += (_dina - ina)/.001 ;
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 
}
 
}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
  ena = _ion_ena;
  ek = _ion_ek;
 {  { states(_p, _ppvar, _thread, _nt); }
  }  }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/hines/tmp/hack/185355/ichan2.mod";
static const char* nmodl_file_text = 
  "TITLE ichan2.mod combination of Nav and Kv channels for Hodgin-Huxely-type action potential mechanism\n"
  " \n"
  "COMMENT\n"
  "\n"
  "Original Mod File:\n"
  "Original name 'ichan2.mod'\n"
  "Santhakumar V, Aradi I, Soltesz I (2005) J Neurophysiol 93:437-53 \n"
  "https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=51781&file=%2fdentategyrusnet2005%2fichan2.mod\n"
  "Morgan RJ, Soltesz I (2008) Proc Natl Acad Sci U S A 105:6179-84\n"
  "Morgan RJ, Santhakumar V, Soltesz I (2007) Prog Brain Res 163:639-58\n"
  "Cutsuridis V, Cobb S, Graham BP (2009) Hippocampus 20(3):423-46 \n"
  "\n"
  "Current version by A. Hanuschkin <AH, 2011> for:\n"
  "Yim MY, Hanuschkin A, Wolfart J (2015) Hippocampus 25:297-308.\n"
  "http://onlinelibrary.wiley.com/doi/10.1002/hipo.22373/abstract\n"
  "\n"
  "Changes in current versus original version:\n"
  " - added a tonic (leak) GABAA conductance to be modified during epilepsy (see Young CC, Stegen M, Bernard R, Muller M, Bischofberger J, Veh RW, Haas CA, Wolfart J (2009) J Physiol 587:4213-4233)\n"
  " - checked, simplified (reduced to single k Ion, Ekf=Eks=Ek) & commented by A. Hanuschkin 2011,2012\n"
  "\n"
  "Mod File history:\n"
  "I_Na: (equivalent to Santhakumar et al 2005)\n"
  "* modified parameters (shifted voltage dependence by 68mV) \n"
  "  Aradi I and Soltesz I (2002) J Physiol. 538(Pt 1):227-51.\n"
  "* dynamics (NOTE in this paper (at least) sign error in alpha_m, beta_m and alpha_n! (1-exp()) -> should be (exp()-1) <ah>)\n"
  "Aradi and Soltesz (2002)\n"
  "* modified from \n"
  "  Yuen GL, Durand D, (1991) Neuroscience 41(2-3):411-23.\n"
  "* Aradi and Soltesz (2002) parameters =  \n"
  "  Aradi I, Holmes WR (1999) J Comput Neurosci 6:215-35\n"
  "parameters = \n"
  "  Yuen and Durand (1991) parameters + V shift of 16mV\n"
  "\n"
  "I_Kf: (equivalent to Santhakumar et al 2005)\n"
  "* modified parameters Aradi and Soltesz (2002)/Aradi & Holmes (1999) (shifted voltage dependence by 65mV) \n"
  "* NOTE typo in formular beta_n in Aradi and Soltesz (2002): beta_n = 0.264/exp((v-22)/4) -> should be beta_n = 0.264/exp((v-22)/40) <ah>\n"
  "* modified parameters from Yuen and Durand (1991)\n"
  "\n"
  "I_Ks: (equivalent to Santhakumar et al 2005)\n"
  "* modified parameters Aradi & Holmes (1999) (shifted voltage dependence by 65mV) \n"
  "\n"
  "I_leak: (equivalent to Santhakumar et al 2005)\n"
  "\n"
  "I_GABAA: (tonic GABAA leak (see above), added in Yim et al (2015))\n"
  "* replicated from I_leak\n"
  " \n"
  "A. Hanuschkin(c) 2011,2012\n"
  "\n"
  "ENDCOMMENT\n"
  " \n"
  "UNITS {\n"
  "        (mA) =(milliamp)\n"
  "        (mV) =(millivolt)\n"
  "        (uF) = (microfarad)\n"
  "	(molar) = (1/liter)\n"
  "	(nA) = (nanoamp)\n"
  "	(mM) = (millimolar)\n"
  "	(um) = (micron)\n"
  "	FARADAY = 96520 (coul)\n"
  "	R = 8.3134	(joule/degC)\n"
  "}\n"
  " \n"
  "NEURON { \n"
  ": \"suffix marks the mechanism to be distributed and whose variables & parameters are identified in hoc by a particular suffix\" The Neuron Book Chap 9.5\n"
  "SUFFIX ichan2\n"
  "\n"
  ": ION usage block\n"
  "USEION na READ ena WRITE ina VALENCE 1			: Na current\n"
  "USEION k READ ek WRITE ik  				: K current\n"
  "NONSPECIFIC_CURRENT il, igabaa 				: leak current\n"
  "\n"
  ": range variable definition block,\n"
  ": i.e. variables that might change with space along a compartment / could be declared global in this case\n"
  "RANGE gnatbar, gkfbar, gksbar				: gbar values for Na, K(slow/fast)\n"
  "RANGE gl, el, ina, ik, il, ggabaa, igabaa, egabaa	: gbar and reversal poti for leak current	\n"
  "}\n"
  "\n"
  ": The INDEPENDENT statement was omitted; INDEPENDENT statements are irrelevant to NEURON. http://www.neuron.yale.edu/phpbb/viewtopic.php?f=16&t=2351 \n"
  ": INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}\n"
  " \n"
  ": Variables whose values are normally specified by the user are parameters, and are declared in a PARAMETER block.\n"
  ": Variables in the parameter section will have global scope\n"
  "PARAMETER {						\n"
  "        :v (mV) \n"
  "        :celsius = 6.3 (degC)\n"
  "        :dt (ms) \n"
  "\n"
  "        gnatbar (mho/cm2)   				: Na (gbar and reversal poti)\n"
  "        ena  	(mV)	\n"
  "		\n"
  "	gkfbar 	(mho/cm2)				: K  (gbar(slow/fast), reversal is ek)\n"
  "	gksbar = 0 (mho/cm2)	                        : init to 0 (not included in BC, HIPP and MC) <ah>\n"
  "        ek     	(mV)                      \n"
  "\n"
  "	gl 	(mho/cm2)    				: leak (gbar and reversal poti)\n"
  " 	el 	(mV)\n"
  "\n"
  "	ggabaa 	(mho/cm2)    				: GABAA (gbar and reversal poti)\n"
  " 	egabaa 	(mV)\n"
  "}\n"
  "\n"
  "\n"
  ": \"If a model involves differential equations [..] their dependent variables or unknowns are to be listed in the STATE block\" The Neuron Book Chap 9.5\n"
  "STATE {\n"
  "	m h nf ns\n"
  "}\n"
  "\n"
  ": The ASSIGNED block is used to declare two kinds of variables\n"
  ": 1) those given values outside the mod file (variables potentially available to every mechanism (e.g. v, celsius,t..)\n"
  ": 2) left hand side of assignment statements (unknowns in set of equations, dependent variables in differential euqtions ...)\n"
  "ASSIGNED {		\n"
  ": 1)\n"
  "        v (mV) \n"
  "        celsius (degC)\n"
  "        dt (ms) \n"
  "	\n"
  ": 2) \n"
  "        gna (mho/cm2) 					: Na\n"
  "        ina (mA/cm2)\n"
  "	\n"
  "        gkf (mho/cm2)					: K\n"
  "        gks (mho/cm2)\n"
  "	ik (mA/cm2)\n"
  "\n"
  "	il (mA/cm2)					: leak \n"
  "\n"
  "	igabaa (mA/cm2)					: GABAA \n"
  "\n"
  "	minf hinf nfinf nsinf				: left hand side of differential equations\n"
  " 	mtau (ms) htau (ms) nftau (ms) nstau (ms)	: and other assignment variables\n"
  "	mexp hexp nfexp nsexp\n"
  "	q10\n"
  "} \n"
  "\n"
  ": This block is evaluated every time step. \n"
  "BREAKPOINT {\n"
  "	SOLVE states					: here the state variables are updated \n"
  "        gna = gnatbar*m*m*m*h  			: calculated g at timepoint t\n"
  "        gkf = gkfbar*nf*nf*nf*nf\n"
  "        gks = gksbar*ns*ns*ns*ns\n"
  "\n"
  "        ina = gna*(v - ena)				: calculated currents flowing\n"
  "       	ik = gkf*(v-ek) + gks*(v-ek)\n"
  "	il = gl*(v-el)\n"
  "	igabaa = ggabaa*(v-egabaa)\n"
  "}\n"
  " \n"
  ": UNITSOFF\n"
  " \n"
  ": Called from Neuron during initializing the model\n"
  "INITIAL {\n"
  "	trates(v)\n"
  "	\n"
  "	m = minf\n"
  "	h = hinf\n"
  "        nf = nfinf\n"
  "        ns = nsinf\n"
  "}\n"
  "\n"
  ": discreticed versions of the differential equations, hence a PROCEDURE and not DERIVATIVE block\n"
  "PROCEDURE states() {	: Computes state variables m, h, and n \n"
  "        trates(v)	: at the current v and dt.\n"
  "        m = m + mexp*(minf-m)\n"
  "        h = h + hexp*(hinf-h)\n"
  "        nf = nf + nfexp*(nfinf-nf)\n"
  "        ns = ns + nsexp*(nsinf-ns)\n"
  "}\n"
  "\n"
  ": moved this to assign block <ah> \n"
  ": LOCAL q10\n"
  "\n"
  ":Computes rate and other constants at current v.\n"
  "PROCEDURE rates(v) {  \n"
  "        LOCAL  alpha, beta, sum\n"
  "        q10 = 3^((celsius - 6.3)/10)\n"
  "                :\"m\" sodium activation system - act and inact cross at -40	: shifted by 68mV compared to in Aradi 1999/2002\n"
  "	alpha = -0.3*vtrap((v+60-17),-5)		: in Aradi 1999: alpha = -0.3*vtrap((v-25),-5); in Aradi 2002: alpha = 0.3*vtrap((v-25),-5) <ah>\n"
  "	beta = 0.3*vtrap((v+60-45),5)			: in Aradi 1999: beta = 0.3*vtrap((v-53),5);  in Aradi 2002:  beta = -0.3*vtrap((v-53),5) <ah>\n"
  "	sum = alpha+beta        \n"
  "	mtau = 1/sum      minf = alpha/sum\n"
  "                :\"h\" sodium inactivation system		: shifted by 68mV compared to in Aradi 1999/2002\n"
  "	alpha = 0.23/exp((v+60+5)/20)			: in Aradi 1999/2002:  alpha = 0.23/exp((v-3)/20) <ah>\n"
  "	beta = 3.33/(1+exp((v+60-47.5)/-10))		: in Aradi 1999/2002:  beta = 3.33/(1+exp((v-55.5)/-10)) <ah>\n"
  "	sum = alpha+beta\n"
  "	htau = 1/sum \n"
  "        hinf = alpha/sum \n"
  "\n"
  "\n"
  "             :\"ns\" sKDR activation system		: shifted by 65mV compared to Aradi 1999 <ah>\n"
  "        alpha = -0.028*vtrap((v+65-35),-6)		: in Aradi 1999: alpha = -0.028*vtrap((v-35),-6) \n"
  "	beta = 0.1056/exp((v+65-10)/40)			: in Aradi 1999: beta = 0.1056/exp((v-10)/40)   \n"
  "	sum = alpha+beta        			\n"
  "	nstau = 1/sum      nsinf = alpha/sum		\n"
  "            :\"nf\" fKDR activation system		: shifted by 65mV compared to Aradi 1999/2002 <ah>\n"
  "        alpha = -0.07*vtrap((v+65-47),-6)		: in Aradi 1999: alpha = -0.07*vtrap((v-47),-6); in Aradi 2002: alpha = 0.07*vtrap((v-47),-6) <ah>\n"
  "	beta = 0.264/exp((v+65-22)/40)			: in Aradi 1999/2002: beta = 0.264/exp((v-22)/40)  // probably typo in Aradi & Soltez 2002 there: beta = 0.264/exp((v-22)/4) <ah>\n"
  "	sum = alpha+beta        \n"
  "	nftau = 1/sum      nfinf = alpha/sum\n"
  "}\n"
  "\n"
  ": Computes rate and other constants at current v. \n"
  "PROCEDURE trates(v) {  \n"
  "	LOCAL tinc\n"
  "        : TABLE minf, mexp, hinf, hexp, nfinf, nfexp, nsinf, nsexp, mtau, htau, nftau, nstau   : <ah>\n"
  "	: DEPEND dt, celsius FROM -100 TO 100 WITH 200					       : <ah>\n"
  "                           \n"
  "	rates(v)	: not consistently executed from here if usetable_hh == 1\n"
  "			: so don't expect the tau values to be tracking along with\n"
  "			: the inf values in hoc\n"
  "\n"
  "        tinc = -dt * q10\n"
  "        mexp = 1 - exp(tinc/mtau)\n"
  "        hexp = 1 - exp(tinc/htau)\n"
  "	nfexp = 1 - exp(tinc/nftau)\n"
  "	nsexp = 1 - exp(tinc/nstau)\n"
  "}\n"
  "\n"
  ":Traps for 0 in denominator of rate eqns.\n"
  "FUNCTION vtrap(x,y) {  \n"
  "        if (fabs(x/y) < 1e-6) {\n"
  "                vtrap = y*(1 - x/y/2)\n"
  "        }else{  \n"
  "                vtrap = x/(exp(x/y) - 1)\n"
  "        }\n"
  "}\n"
  " \n"
  ":UNITSON\n"
  "\n"
  ;
#endif
