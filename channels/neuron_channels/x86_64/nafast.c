/* Created by Language version: 6.2.0 */
/* NOT VECTORIZED */
#include <stdio.h>
#include <math.h>
#include "scoplib.h"
#undef PI
 
#include "md1redef.h"
#include "section.h"
#include "md2redef.h"

#if METHOD3
extern int _method3;
#endif

#undef exp
#define exp hoc_Exp
extern double hoc_Exp();
 
#define _threadargscomma_ /**/
#define _threadargs_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gnabar _p[0]
#define ina _p[1]
#define m _p[2]
#define h _p[3]
#define ena _p[4]
#define Dm _p[5]
#define Dh _p[6]
#define _g _p[7]
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 /* declaration of user functions */
 static int _hoc_alp();
 static int _hoc_bet();
 static int _hoc_expM1();
 static int _hoc_rates();
 static int _mechtype;
extern int nrn_get_mechtype();
 static _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range();
 _prop = hoc_getdata_range(_mechtype);
 _p = _prop->param; _ppvar = _prop->dparam;
 ret(1.);
}
 /* connect user functions to hoc names */
 static IntFunc hoc_intfunc[] = {
 "setdata_Na_mit_usb", _hoc_setdata,
 "alp_Na_mit_usb", _hoc_alp,
 "bet_Na_mit_usb", _hoc_bet,
 "expM1_Na_mit_usb", _hoc_expM1,
 "rates_Na_mit_usb", _hoc_rates,
 0, 0
};
#define alp alp_Na_mit_usb
#define bet bet_Na_mit_usb
#define expM1 expM1_Na_mit_usb
 extern double alp();
 extern double bet();
 extern double expM1();
 /* declare global and static user variables */
#define htau htau_Na_mit_usb
 double htau = 0;
#define hinf hinf_Na_mit_usb
 double hinf = 0;
#define mtau mtau_Na_mit_usb
 double mtau = 0;
#define minf minf_Na_mit_usb
 double minf = 0;
#define usetable usetable_Na_mit_usb
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "gnabar_Na_mit_usb", 0, 1e+09,
 "usetable_Na_mit_usb", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "mtau_Na_mit_usb", "ms",
 "htau_Na_mit_usb", "ms",
 "gnabar_Na_mit_usb", "mho/cm2",
 "ina_Na_mit_usb", "mA/cm2",
 0,0
};
 static double delta_t = 1;
 static double h0 = 0;
 static double m0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "minf_Na_mit_usb", &minf_Na_mit_usb,
 "hinf_Na_mit_usb", &hinf_Na_mit_usb,
 "mtau_Na_mit_usb", &mtau_Na_mit_usb,
 "htau_Na_mit_usb", &htau_Na_mit_usb,
 "usetable_Na_mit_usb", &usetable_Na_mit_usb,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(), nrn_init(), nrn_state();
 static void nrn_cur(), nrn_jacob();
 
static int _ode_count(), _ode_map(), _ode_spec(), _ode_matsol();
 
#define _cvode_ieq _ppvar[3]._i
 /* connect range variables in _p that hoc is supposed to know about */
 static char *_mechanism[] = {
 "6.2.0",
"Na_mit_usb",
 "gnabar_Na_mit_usb",
 0,
 "ina_Na_mit_usb",
 0,
 "m_Na_mit_usb",
 "h_Na_mit_usb",
 0,
 0};
 static Symbol* _na_sym;
 
static void nrn_alloc(_prop)
	Prop *_prop;
{
	Prop *prop_ion, *need_memb();
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 8, _prop);
 	/*initialize range parameters*/
 	gnabar = 0.12;
 	_prop->param = _p;
 	_prop->param_size = 8;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 
}
 static _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 _nafast_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("na", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
  hoc_register_dparam_size(_mechtype, 4);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 Na_mit_usb /home/aditya/aditya_OB_model/channels/neuron_channels/x86_64/nafast.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_minf;
 static double *_t_hinf;
 static double *_t_mtau;
 static double *_t_htau;
static int _reset;
static char *modelname = "HH fast sodium channel";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static _modl_cleanup(){ _match_recurse=1;}
static _f_rates();
static rates();
 
static int _ode_spec1(), _ode_matsol1();
 static _n_rates();
 static int _slist1[2], _dlist1[2];
 static int states();
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / mtau ;
   Dh = ( hinf - h ) / htau ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rates ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / htau )) ;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / mtau)))*(- ( ( ( minf ) ) / mtau ) / ( ( ( ( - 1.0) ) ) / mtau ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / htau)))*(- ( ( ( hinf ) ) / htau ) / ( ( ( ( - 1.0) ) ) / htau ) - h) ;
   }
  return 0;
}
 
double alp (  _lv , _li )  
	double _lv , _li ;
 {
   double _lalp;
 if ( _li  == 0.0 ) {
     _lalp = 0.32 * expM1 ( _threadargscomma_ - ( _lv * 1.0 + 42.0 ) , 4.0 ) ;
     }
   else if ( _li  == 1.0 ) {
     _lalp = 0.128 / ( exp ( ( _lv * 1.0 + 38.0 ) / 18.0 ) ) ;
     }
   
return _lalp;
 }
 
static int _hoc_alp() {
  double _r;
   _r =  alp (  *getarg(1) , *getarg(2) ) ;
 ret(_r);
}
 
double bet (  _lv , _li )  
	double _lv , _li ;
 {
   double _lbet;
 if ( _li  == 0.0 ) {
     _lbet = 0.28 * expM1 ( _threadargscomma_ _lv * 1.0 + 15.0 , 5.0 ) ;
     }
   else if ( _li  == 1.0 ) {
     _lbet = 4.0 / ( exp ( - ( _lv * 1.0 + 15.0 ) / 5.0 ) + 1.0 ) ;
     }
   
return _lbet;
 }
 
static int _hoc_bet() {
  double _r;
   _r =  bet (  *getarg(1) , *getarg(2) ) ;
 ret(_r);
}
 
double expM1 (  _lx , _ly )  
	double _lx , _ly ;
 {
   double _lexpM1;
 if ( fabs ( _lx / _ly ) < 1e-6 ) {
     _lexpM1 = _ly * ( 1.0 - _lx / _ly / 2.0 ) ;
     }
   else {
     _lexpM1 = _lx / ( exp ( _lx / _ly ) - 1.0 ) ;
     }
   
return _lexpM1;
 }
 
static int _hoc_expM1() {
  double _r;
   _r =  expM1 (  *getarg(1) , *getarg(2) ) ;
 ret(_r);
}
 static double _mfac_rates, _tmin_rates;
 static _check_rates();
 static _check_rates() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  if (!usetable) {return;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_rates =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_rates)/200.; _mfac_rates = 1./_dx;
   for (_i=0, _x=_tmin_rates; _i < 201; _x += _dx, _i++) {
    _f_rates(_x);
    _t_minf[_i] = minf;
    _t_hinf[_i] = hinf;
    _t_mtau[_i] = mtau;
    _t_htau[_i] = htau;
   }
  }
 }

 static rates(double _lv){ _check_rates();
 _n_rates(_lv);
 return;
 }

 static _n_rates(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_rates(_lv); return; 
}
 _xi = _mfac_rates * (_lv - _tmin_rates);
 _i = (int) _xi;
 if (_xi <= 0.) {
 minf = _t_minf[0];
 hinf = _t_hinf[0];
 mtau = _t_mtau[0];
 htau = _t_htau[0];
 return; }
 if (_i >= 200) {
 minf = _t_minf[200];
 hinf = _t_hinf[200];
 mtau = _t_mtau[200];
 htau = _t_htau[200];
 return; }
 _theta = _xi - (double)_i;
 minf = _t_minf[_i] + _theta*(_t_minf[_i+1] - _t_minf[_i]);
 hinf = _t_hinf[_i] + _theta*(_t_hinf[_i+1] - _t_hinf[_i]);
 mtau = _t_mtau[_i] + _theta*(_t_mtau[_i+1] - _t_mtau[_i]);
 htau = _t_htau[_i] + _theta*(_t_htau[_i+1] - _t_htau[_i]);
 }

 
static int  _f_rates (  _lv )  
	double _lv ;
 {
   double _la , _lb ;
 _la = alp ( _threadargscomma_ _lv , 0.0 ) ;
   _lb = bet ( _threadargscomma_ _lv , 0.0 ) ;
   mtau = 1.0 / ( _la + _lb ) ;
   minf = _la / ( _la + _lb ) ;
   _la = alp ( _threadargscomma_ _lv , 1.0 ) ;
   _lb = bet ( _threadargscomma_ _lv , 1.0 ) ;
   htau = 1.0 / ( _la + _lb ) ;
   hinf = _la / ( _la + _lb ) ;
    return 0; }
 
static int _hoc_rates() {
  double _r;
    _r = 1.;
 rates (  *getarg(1) ) ;
 ret(_r);
}
 
static int _ode_count(_type) int _type;{ return 2;}
 
static int _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
     _ode_spec1 ();
  }}
 
static int _ode_map(_ieq, _pv, _pvdot, _pp, _ppd, _atol, _type) int _ieq, _type; double** _pv, **_pvdot, *_pp, *_atol; Datum* _ppd; { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static int _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
 _ode_matsol1 ();
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  m = m0;
 {
   rates ( _threadargscomma_ v ) ;
   m = minf ;
   h = hinf ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
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
  ena = _ion_ena;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   ina = gnabar * m * m * m * h * ( v - ena ) ;
   }
 _current += ina;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
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
  ena = _ion_ena;
 _g = _nrn_current(_v + .001);
 	{ double _dina;
  _dina = ina;
 _rhs = _nrn_current(_v);
  _ion_dinadv += (_dina - ina)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
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

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
 double _break, _save;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
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
 _break = t + .5*dt; _save = t;
 v=_v;
{
  ena = _ion_ena;
 { {
 for (; t < _break; t += dt) {
 error =  states();
 if(error){fprintf(stderr,"at line 44 in file nafast.mod:\n	SOLVE states METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 
}}
 t = _save;
 } }}

}

static terminal(){}

static _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
   _t_minf = makevector(201*sizeof(double));
   _t_hinf = makevector(201*sizeof(double));
   _t_mtau = makevector(201*sizeof(double));
   _t_htau = makevector(201*sizeof(double));
_first = 0;
}
