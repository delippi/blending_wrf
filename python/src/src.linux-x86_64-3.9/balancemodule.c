/* File: balancemodule.c
 * This file is auto-generated with f2py (version:1.23.5).
 * f2py is a Fortran to Python Interface Generator (FPIG), Second Edition,
 * written by Pearu Peterson <pearu@cens.ioc.ee>.
 * Generation date: Fri Apr 14 16:56:03 2023
 * Do not edit this file directly unless you know what you are doing!!!
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef PY_SSIZE_T_CLEAN
#define PY_SSIZE_T_CLEAN
#endif /* PY_SSIZE_T_CLEAN */

/* Unconditionally included */
#include <Python.h>
#include <numpy/npy_os.h>

/*********************** See f2py2e/cfuncs.py: includes ***********************/
#include <stdarg.h>
#include "fortranobject.h"
#include <math.h>

/**************** See f2py2e/rules.py: mod_rules['modulebody'] ****************/
static PyObject *balance_error;
static PyObject *balance_module;

/*********************** See f2py2e/cfuncs.py: typedefs ***********************/
/*need_typedefs*/

/****************** See f2py2e/cfuncs.py: typedefs_generated ******************/
/*need_typedefs_generated*/

/********************** See f2py2e/cfuncs.py: cppmacros **********************/
#if defined(PREPEND_FORTRAN)
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F
#else
#define F_FUNC(f,F) _##f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F##_
#else
#define F_FUNC(f,F) _##f##_
#endif
#endif
#else
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif
#endif
#if defined(UNDERSCORE_G77)
#define F_FUNC_US(f,F) F_FUNC(f##_,F##_)
#else
#define F_FUNC_US(f,F) F_FUNC(f,F)
#endif

#define rank(var) var ## _Rank
#define shape(var,dim) var ## _Dims[dim]
#define old_rank(var) (PyArray_NDIM((PyArrayObject *)(capi_ ## var ## _tmp)))
#define old_shape(var,dim) PyArray_DIM(((PyArrayObject *)(capi_ ## var ## _tmp)),dim)
#define fshape(var,dim) shape(var,rank(var)-dim-1)
#define len(var) shape(var,0)
#define flen(var) fshape(var,0)
#define old_size(var) PyArray_SIZE((PyArrayObject *)(capi_ ## var ## _tmp))
/* #define index(i) capi_i ## i */
#define slen(var) capi_ ## var ## _len
#define size(var, ...) f2py_size((PyArrayObject *)(capi_ ## var ## _tmp), ## __VA_ARGS__, -1)

#define CHECKSCALAR(check,tcheck,name,show,var)\
    if (!(check)) {\
        char errstring[256];\
        sprintf(errstring, "%s: "show, "("tcheck") failed for "name, var);\
        PyErr_SetString(balance_error,errstring);\
        /*goto capi_fail;*/\
    } else 
#ifdef DEBUGCFUNCS
#define CFUNCSMESS(mess) fprintf(stderr,"debug-capi:"mess);
#define CFUNCSMESSPY(mess,obj) CFUNCSMESS(mess) \
    PyObject_Print((PyObject *)obj,stderr,Py_PRINT_RAW);\
    fprintf(stderr,"\n");
#else
#define CFUNCSMESS(mess)
#define CFUNCSMESSPY(mess,obj)
#endif

#ifndef max
#define max(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) ((a < b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a < b) ? (a) : (b))
#endif


/************************ See f2py2e/cfuncs.py: cfuncs ************************/
static int f2py_size(PyArrayObject* var, ...)
{
  npy_int sz = 0;
  npy_int dim;
  npy_int rank;
  va_list argp;
  va_start(argp, var);
  dim = va_arg(argp, npy_int);
  if (dim==-1)
    {
      sz = PyArray_SIZE(var);
    }
  else
    {
      rank = PyArray_NDIM(var);
      if (dim>=1 && dim<=rank)
        sz = PyArray_DIM(var, dim-1);
      else
        fprintf(stderr, "f2py_size: 2nd argument value=%d fails to satisfy 1<=value<=%d. Result will be 0.\n", dim, rank);
    }
  va_end(argp);
  return sz;
}

static int
double_from_pyobj(double* v, PyObject *obj, const char *errmess)
{
    PyObject* tmp = NULL;
    if (PyFloat_Check(obj)) {
        *v = PyFloat_AsDouble(obj);
        return !(*v == -1.0 && PyErr_Occurred());
    }

    tmp = PyNumber_Float(obj);
    if (tmp) {
        *v = PyFloat_AsDouble(tmp);
        Py_DECREF(tmp);
        return !(*v == -1.0 && PyErr_Occurred());
    }

    if (PyComplex_Check(obj)) {
        PyErr_Clear();
        tmp = PyObject_GetAttrString(obj,"real");
    }
    else if (PyBytes_Check(obj) || PyUnicode_Check(obj)) {
        /*pass*/;
    }
    else if (PySequence_Check(obj)) {
        PyErr_Clear();
        tmp = PySequence_GetItem(obj, 0);
    }

    if (tmp) {
        if (double_from_pyobj(v,tmp,errmess)) {Py_DECREF(tmp); return 1;}
        Py_DECREF(tmp);
    }
    {
        PyObject* err = PyErr_Occurred();
        if (err==NULL) err = balance_error;
        PyErr_SetString(err,errmess);
    }
    return 0;
}

static int
int_from_pyobj(int* v, PyObject *obj, const char *errmess)
{
    PyObject* tmp = NULL;

    if (PyLong_Check(obj)) {
        *v = Npy__PyLong_AsInt(obj);
        return !(*v == -1 && PyErr_Occurred());
    }

    tmp = PyNumber_Long(obj);
    if (tmp) {
        *v = Npy__PyLong_AsInt(tmp);
        Py_DECREF(tmp);
        return !(*v == -1 && PyErr_Occurred());
    }

    if (PyComplex_Check(obj)) {
        PyErr_Clear();
        tmp = PyObject_GetAttrString(obj,"real");
    }
    else if (PyBytes_Check(obj) || PyUnicode_Check(obj)) {
        /*pass*/;
    }
    else if (PySequence_Check(obj)) {
        PyErr_Clear();
        tmp = PySequence_GetItem(obj, 0);
    }

    if (tmp) {
        if (int_from_pyobj(v, tmp, errmess)) {
            Py_DECREF(tmp);
            return 1;
        }
        Py_DECREF(tmp);
    }

    {
        PyObject* err = PyErr_Occurred();
        if (err == NULL) {
            err = balance_error;
        }
        PyErr_SetString(err, errmess);
    }
    return 0;
}


/********************* See f2py2e/cfuncs.py: userincludes *********************/
/*need_userincludes*/

/********************* See f2py2e/capi_rules.py: usercode *********************/


/* See f2py2e/rules.py */
extern void F_FUNC(balance,BALANCE)(int*,int*,int*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);
/*eof externroutines*/

/******************** See f2py2e/capi_rules.py: usercode1 ********************/


/******************* See f2py2e/cb_rules.py: buildcallback *******************/
/*need_callbacks*/

/*********************** See f2py2e/rules.py: buildapi ***********************/

/********************************** balance **********************************/
static char doc_f2py_rout_balance_balance[] = "\
ph,mu_2d = balance(qw,qvapor,ph1,phb,t,tw,mu1,psw,mub_2d,hgt_2d,dnw_1d,znu_1d,rdnw_1d,rdn_1d,ptop,ph,mu_2d,[ml,nl,kk1,overwrite_ph,overwrite_mu_2d])\n\nWrapper for ``balance``.\
\n\nParameters\n----------\n"
"qw : input rank-3 array('d') with bounds (ml,nl,-1 + kk1)\n"
"qvapor : input rank-3 array('d') with bounds (ml,nl,-1 + kk1)\n"
"ph1 : input rank-3 array('d') with bounds (ml,nl,kk1)\n"
"phb : input rank-3 array('d') with bounds (ml,nl,kk1)\n"
"t : input rank-3 array('d') with bounds (ml,nl,-1 + kk1)\n"
"tw : input rank-3 array('d') with bounds (ml,nl,-1 + kk1)\n"
"mu1 : input rank-2 array('d') with bounds (ml,nl)\n"
"psw : input rank-2 array('d') with bounds (ml,nl)\n"
"mub_2d : input rank-2 array('d') with bounds (ml,nl)\n"
"hgt_2d : input rank-2 array('d') with bounds (ml,nl)\n"
"dnw_1d : input rank-1 array('d') with bounds (-1 + kk1)\n"
"znu_1d : input rank-1 array('d') with bounds (-1 + kk1)\n"
"rdnw_1d : input rank-1 array('d') with bounds (-1 + kk1)\n"
"rdn_1d : input rank-1 array('d') with bounds (-1 + kk1)\n"
"ptop : input float\n"
"ph : input rank-3 array('d') with bounds (ml,nl,kk1)\n"
"mu_2d : input rank-2 array('d') with bounds (ml,nl)\n"
"\nOther Parameters\n----------------\n"
"ml : input int, optional\n    Default: shape(qw, 0)\n"
"nl : input int, optional\n    Default: shape(qw, 1)\n"
"kk1 : input int, optional\n    Default: 1 + shape(qw, 2)\n"
"overwrite_ph : input int, optional\n    Default: 0\n"
"overwrite_mu_2d : input int, optional\n    Default: 0\n"
"\nReturns\n-------\n"
"ph : rank-3 array('d') with bounds (ml,nl,kk1)\n"
"mu_2d : rank-2 array('d') with bounds (ml,nl)";
/* extern void F_FUNC(balance,BALANCE)(int*,int*,int*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*); */
static PyObject *f2py_rout_balance_balance(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(int*,int*,int*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*)) {
    PyObject * volatile capi_buildvalue = NULL;
    volatile int f2py_success = 1;
/*decl*/

    int ml = 0;
    PyObject *ml_capi = Py_None;
    int nl = 0;
    PyObject *nl_capi = Py_None;
    int kk1 = 0;
    PyObject *kk1_capi = Py_None;
    double *qw = NULL;
    npy_intp qw_Dims[3] = {-1, -1, -1};
    const int qw_Rank = 3;
    PyArrayObject *capi_qw_tmp = NULL;
    int capi_qw_intent = 0;
    PyObject *qw_capi = Py_None;
    double *qvapor = NULL;
    npy_intp qvapor_Dims[3] = {-1, -1, -1};
    const int qvapor_Rank = 3;
    PyArrayObject *capi_qvapor_tmp = NULL;
    int capi_qvapor_intent = 0;
    PyObject *qvapor_capi = Py_None;
    double *ph1 = NULL;
    npy_intp ph1_Dims[3] = {-1, -1, -1};
    const int ph1_Rank = 3;
    PyArrayObject *capi_ph1_tmp = NULL;
    int capi_ph1_intent = 0;
    PyObject *ph1_capi = Py_None;
    double *phb = NULL;
    npy_intp phb_Dims[3] = {-1, -1, -1};
    const int phb_Rank = 3;
    PyArrayObject *capi_phb_tmp = NULL;
    int capi_phb_intent = 0;
    PyObject *phb_capi = Py_None;
    double *t = NULL;
    npy_intp t_Dims[3] = {-1, -1, -1};
    const int t_Rank = 3;
    PyArrayObject *capi_t_tmp = NULL;
    int capi_t_intent = 0;
    PyObject *t_capi = Py_None;
    double *tw = NULL;
    npy_intp tw_Dims[3] = {-1, -1, -1};
    const int tw_Rank = 3;
    PyArrayObject *capi_tw_tmp = NULL;
    int capi_tw_intent = 0;
    PyObject *tw_capi = Py_None;
    double *mu1 = NULL;
    npy_intp mu1_Dims[2] = {-1, -1};
    const int mu1_Rank = 2;
    PyArrayObject *capi_mu1_tmp = NULL;
    int capi_mu1_intent = 0;
    PyObject *mu1_capi = Py_None;
    double *psw = NULL;
    npy_intp psw_Dims[2] = {-1, -1};
    const int psw_Rank = 2;
    PyArrayObject *capi_psw_tmp = NULL;
    int capi_psw_intent = 0;
    PyObject *psw_capi = Py_None;
    double *mub_2d = NULL;
    npy_intp mub_2d_Dims[2] = {-1, -1};
    const int mub_2d_Rank = 2;
    PyArrayObject *capi_mub_2d_tmp = NULL;
    int capi_mub_2d_intent = 0;
    PyObject *mub_2d_capi = Py_None;
    double *hgt_2d = NULL;
    npy_intp hgt_2d_Dims[2] = {-1, -1};
    const int hgt_2d_Rank = 2;
    PyArrayObject *capi_hgt_2d_tmp = NULL;
    int capi_hgt_2d_intent = 0;
    PyObject *hgt_2d_capi = Py_None;
    double *dnw_1d = NULL;
    npy_intp dnw_1d_Dims[1] = {-1};
    const int dnw_1d_Rank = 1;
    PyArrayObject *capi_dnw_1d_tmp = NULL;
    int capi_dnw_1d_intent = 0;
    PyObject *dnw_1d_capi = Py_None;
    double *znu_1d = NULL;
    npy_intp znu_1d_Dims[1] = {-1};
    const int znu_1d_Rank = 1;
    PyArrayObject *capi_znu_1d_tmp = NULL;
    int capi_znu_1d_intent = 0;
    PyObject *znu_1d_capi = Py_None;
    double *rdnw_1d = NULL;
    npy_intp rdnw_1d_Dims[1] = {-1};
    const int rdnw_1d_Rank = 1;
    PyArrayObject *capi_rdnw_1d_tmp = NULL;
    int capi_rdnw_1d_intent = 0;
    PyObject *rdnw_1d_capi = Py_None;
    double *rdn_1d = NULL;
    npy_intp rdn_1d_Dims[1] = {-1};
    const int rdn_1d_Rank = 1;
    PyArrayObject *capi_rdn_1d_tmp = NULL;
    int capi_rdn_1d_intent = 0;
    PyObject *rdn_1d_capi = Py_None;
    double ptop = 0;
    PyObject *ptop_capi = Py_None;
    double *ph = NULL;
    npy_intp ph_Dims[3] = {-1, -1, -1};
    const int ph_Rank = 3;
    PyArrayObject *capi_ph_tmp = NULL;
    int capi_ph_intent = 0;
    int capi_overwrite_ph = 0;
    PyObject *ph_capi = Py_None;
    double *mu_2d = NULL;
    npy_intp mu_2d_Dims[2] = {-1, -1};
    const int mu_2d_Rank = 2;
    PyArrayObject *capi_mu_2d_tmp = NULL;
    int capi_mu_2d_intent = 0;
    int capi_overwrite_mu_2d = 0;
    PyObject *mu_2d_capi = Py_None;
    static char *capi_kwlist[] = {"qw","qvapor","ph1","phb","t","tw","mu1","psw","mub_2d","hgt_2d","dnw_1d","znu_1d","rdnw_1d","rdn_1d","ptop","ph","mu_2d","ml","nl","kk1","overwrite_ph","overwrite_mu_2d",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
    if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
        "OOOOOOOOOOOOOOOOO|OOOii:balance.balance",\
        capi_kwlist,&qw_capi,&qvapor_capi,&ph1_capi,&phb_capi,&t_capi,&tw_capi,&mu1_capi,&psw_capi,&mub_2d_capi,&hgt_2d_capi,&dnw_1d_capi,&znu_1d_capi,&rdnw_1d_capi,&rdn_1d_capi,&ptop_capi,&ph_capi,&mu_2d_capi,&ml_capi,&nl_capi,&kk1_capi,&capi_overwrite_ph,&capi_overwrite_mu_2d))
        return NULL;
/*frompyobj*/
    /* Processing variable ptop */
        f2py_success = double_from_pyobj(&ptop,ptop_capi,"balance.balance() 15th argument (ptop) can't be converted to double");
    if (f2py_success) {
    /* Processing variable qw */
    ;
    capi_qw_intent |= F2PY_INTENT_IN;
    capi_qw_tmp = array_from_pyobj(NPY_DOUBLE,qw_Dims,qw_Rank,capi_qw_intent,qw_capi);
    if (capi_qw_tmp == NULL) {
        PyObject *exc, *val, *tb;
        PyErr_Fetch(&exc, &val, &tb);
        PyErr_SetString(exc ? exc : balance_error,"failed in converting 1st argument `qw' of balance.balance to C/Fortran array" );
        npy_PyErr_ChainExceptionsCause(exc, val, tb);
    } else {
        qw = (double *)(PyArray_DATA(capi_qw_tmp));

    /* Processing variable ml */
    if (ml_capi == Py_None) ml = shape(qw, 0); else
        f2py_success = int_from_pyobj(&ml,ml_capi,"balance.balance() 1st keyword (ml) can't be converted to int");
    if (f2py_success) {
    CHECKSCALAR(shape(qw, 0) == ml,"shape(qw, 0) == ml","1st keyword ml","balance:ml=%d",ml) {
    /* Processing variable nl */
    if (nl_capi == Py_None) nl = shape(qw, 1); else
        f2py_success = int_from_pyobj(&nl,nl_capi,"balance.balance() 2nd keyword (nl) can't be converted to int");
    if (f2py_success) {
    CHECKSCALAR(shape(qw, 1) == nl,"shape(qw, 1) == nl","2nd keyword nl","balance:nl=%d",nl) {
    /* Processing variable kk1 */
    if (kk1_capi == Py_None) kk1 = 1 + shape(qw, 2); else
        f2py_success = int_from_pyobj(&kk1,kk1_capi,"balance.balance() 3rd keyword (kk1) can't be converted to int");
    if (f2py_success) {
    CHECKSCALAR(shape(qw, 2) == -1 + kk1,"shape(qw, 2) == -1 + kk1","3rd keyword kk1","balance:kk1=%d",kk1) {
    /* Processing variable qvapor */
    qvapor_Dims[0]=ml,qvapor_Dims[1]=nl,qvapor_Dims[2]=-1 + kk1;
    capi_qvapor_intent |= F2PY_INTENT_IN;
    capi_qvapor_tmp = array_from_pyobj(NPY_DOUBLE,qvapor_Dims,qvapor_Rank,capi_qvapor_intent,qvapor_capi);
    if (capi_qvapor_tmp == NULL) {
        PyObject *exc, *val, *tb;
        PyErr_Fetch(&exc, &val, &tb);
        PyErr_SetString(exc ? exc : balance_error,"failed in converting 2nd argument `qvapor' of balance.balance to C/Fortran array" );
        npy_PyErr_ChainExceptionsCause(exc, val, tb);
    } else {
        qvapor = (double *)(PyArray_DATA(capi_qvapor_tmp));

    /* Processing variable t */
    t_Dims[0]=ml,t_Dims[1]=nl,t_Dims[2]=-1 + kk1;
    capi_t_intent |= F2PY_INTENT_IN;
    capi_t_tmp = array_from_pyobj(NPY_DOUBLE,t_Dims,t_Rank,capi_t_intent,t_capi);
    if (capi_t_tmp == NULL) {
        PyObject *exc, *val, *tb;
        PyErr_Fetch(&exc, &val, &tb);
        PyErr_SetString(exc ? exc : balance_error,"failed in converting 5th argument `t' of balance.balance to C/Fortran array" );
        npy_PyErr_ChainExceptionsCause(exc, val, tb);
    } else {
        t = (double *)(PyArray_DATA(capi_t_tmp));

    /* Processing variable tw */
    tw_Dims[0]=ml,tw_Dims[1]=nl,tw_Dims[2]=-1 + kk1;
    capi_tw_intent |= F2PY_INTENT_IN;
    capi_tw_tmp = array_from_pyobj(NPY_DOUBLE,tw_Dims,tw_Rank,capi_tw_intent,tw_capi);
    if (capi_tw_tmp == NULL) {
        PyObject *exc, *val, *tb;
        PyErr_Fetch(&exc, &val, &tb);
        PyErr_SetString(exc ? exc : balance_error,"failed in converting 6th argument `tw' of balance.balance to C/Fortran array" );
        npy_PyErr_ChainExceptionsCause(exc, val, tb);
    } else {
        tw = (double *)(PyArray_DATA(capi_tw_tmp));

    /* Processing variable ph1 */
    ph1_Dims[0]=ml,ph1_Dims[1]=nl,ph1_Dims[2]=kk1;
    capi_ph1_intent |= F2PY_INTENT_IN;
    capi_ph1_tmp = array_from_pyobj(NPY_DOUBLE,ph1_Dims,ph1_Rank,capi_ph1_intent,ph1_capi);
    if (capi_ph1_tmp == NULL) {
        PyObject *exc, *val, *tb;
        PyErr_Fetch(&exc, &val, &tb);
        PyErr_SetString(exc ? exc : balance_error,"failed in converting 3rd argument `ph1' of balance.balance to C/Fortran array" );
        npy_PyErr_ChainExceptionsCause(exc, val, tb);
    } else {
        ph1 = (double *)(PyArray_DATA(capi_ph1_tmp));

    /* Processing variable phb */
    phb_Dims[0]=ml,phb_Dims[1]=nl,phb_Dims[2]=kk1;
    capi_phb_intent |= F2PY_INTENT_IN;
    capi_phb_tmp = array_from_pyobj(NPY_DOUBLE,phb_Dims,phb_Rank,capi_phb_intent,phb_capi);
    if (capi_phb_tmp == NULL) {
        PyObject *exc, *val, *tb;
        PyErr_Fetch(&exc, &val, &tb);
        PyErr_SetString(exc ? exc : balance_error,"failed in converting 4th argument `phb' of balance.balance to C/Fortran array" );
        npy_PyErr_ChainExceptionsCause(exc, val, tb);
    } else {
        phb = (double *)(PyArray_DATA(capi_phb_tmp));

    /* Processing variable mu1 */
    mu1_Dims[0]=ml,mu1_Dims[1]=nl;
    capi_mu1_intent |= F2PY_INTENT_IN;
    capi_mu1_tmp = array_from_pyobj(NPY_DOUBLE,mu1_Dims,mu1_Rank,capi_mu1_intent,mu1_capi);
    if (capi_mu1_tmp == NULL) {
        PyObject *exc, *val, *tb;
        PyErr_Fetch(&exc, &val, &tb);
        PyErr_SetString(exc ? exc : balance_error,"failed in converting 7th argument `mu1' of balance.balance to C/Fortran array" );
        npy_PyErr_ChainExceptionsCause(exc, val, tb);
    } else {
        mu1 = (double *)(PyArray_DATA(capi_mu1_tmp));

    /* Processing variable psw */
    psw_Dims[0]=ml,psw_Dims[1]=nl;
    capi_psw_intent |= F2PY_INTENT_IN;
    capi_psw_tmp = array_from_pyobj(NPY_DOUBLE,psw_Dims,psw_Rank,capi_psw_intent,psw_capi);
    if (capi_psw_tmp == NULL) {
        PyObject *exc, *val, *tb;
        PyErr_Fetch(&exc, &val, &tb);
        PyErr_SetString(exc ? exc : balance_error,"failed in converting 8th argument `psw' of balance.balance to C/Fortran array" );
        npy_PyErr_ChainExceptionsCause(exc, val, tb);
    } else {
        psw = (double *)(PyArray_DATA(capi_psw_tmp));

    /* Processing variable mub_2d */
    mub_2d_Dims[0]=ml,mub_2d_Dims[1]=nl;
    capi_mub_2d_intent |= F2PY_INTENT_IN;
    capi_mub_2d_tmp = array_from_pyobj(NPY_DOUBLE,mub_2d_Dims,mub_2d_Rank,capi_mub_2d_intent,mub_2d_capi);
    if (capi_mub_2d_tmp == NULL) {
        PyObject *exc, *val, *tb;
        PyErr_Fetch(&exc, &val, &tb);
        PyErr_SetString(exc ? exc : balance_error,"failed in converting 9th argument `mub_2d' of balance.balance to C/Fortran array" );
        npy_PyErr_ChainExceptionsCause(exc, val, tb);
    } else {
        mub_2d = (double *)(PyArray_DATA(capi_mub_2d_tmp));

    /* Processing variable hgt_2d */
    hgt_2d_Dims[0]=ml,hgt_2d_Dims[1]=nl;
    capi_hgt_2d_intent |= F2PY_INTENT_IN;
    capi_hgt_2d_tmp = array_from_pyobj(NPY_DOUBLE,hgt_2d_Dims,hgt_2d_Rank,capi_hgt_2d_intent,hgt_2d_capi);
    if (capi_hgt_2d_tmp == NULL) {
        PyObject *exc, *val, *tb;
        PyErr_Fetch(&exc, &val, &tb);
        PyErr_SetString(exc ? exc : balance_error,"failed in converting 10th argument `hgt_2d' of balance.balance to C/Fortran array" );
        npy_PyErr_ChainExceptionsCause(exc, val, tb);
    } else {
        hgt_2d = (double *)(PyArray_DATA(capi_hgt_2d_tmp));

    /* Processing variable dnw_1d */
    dnw_1d_Dims[0]=-1 + kk1;
    capi_dnw_1d_intent |= F2PY_INTENT_IN;
    capi_dnw_1d_tmp = array_from_pyobj(NPY_DOUBLE,dnw_1d_Dims,dnw_1d_Rank,capi_dnw_1d_intent,dnw_1d_capi);
    if (capi_dnw_1d_tmp == NULL) {
        PyObject *exc, *val, *tb;
        PyErr_Fetch(&exc, &val, &tb);
        PyErr_SetString(exc ? exc : balance_error,"failed in converting 11st argument `dnw_1d' of balance.balance to C/Fortran array" );
        npy_PyErr_ChainExceptionsCause(exc, val, tb);
    } else {
        dnw_1d = (double *)(PyArray_DATA(capi_dnw_1d_tmp));

    /* Processing variable znu_1d */
    znu_1d_Dims[0]=-1 + kk1;
    capi_znu_1d_intent |= F2PY_INTENT_IN;
    capi_znu_1d_tmp = array_from_pyobj(NPY_DOUBLE,znu_1d_Dims,znu_1d_Rank,capi_znu_1d_intent,znu_1d_capi);
    if (capi_znu_1d_tmp == NULL) {
        PyObject *exc, *val, *tb;
        PyErr_Fetch(&exc, &val, &tb);
        PyErr_SetString(exc ? exc : balance_error,"failed in converting 12nd argument `znu_1d' of balance.balance to C/Fortran array" );
        npy_PyErr_ChainExceptionsCause(exc, val, tb);
    } else {
        znu_1d = (double *)(PyArray_DATA(capi_znu_1d_tmp));

    /* Processing variable rdnw_1d */
    rdnw_1d_Dims[0]=-1 + kk1;
    capi_rdnw_1d_intent |= F2PY_INTENT_IN;
    capi_rdnw_1d_tmp = array_from_pyobj(NPY_DOUBLE,rdnw_1d_Dims,rdnw_1d_Rank,capi_rdnw_1d_intent,rdnw_1d_capi);
    if (capi_rdnw_1d_tmp == NULL) {
        PyObject *exc, *val, *tb;
        PyErr_Fetch(&exc, &val, &tb);
        PyErr_SetString(exc ? exc : balance_error,"failed in converting 13rd argument `rdnw_1d' of balance.balance to C/Fortran array" );
        npy_PyErr_ChainExceptionsCause(exc, val, tb);
    } else {
        rdnw_1d = (double *)(PyArray_DATA(capi_rdnw_1d_tmp));

    /* Processing variable rdn_1d */
    rdn_1d_Dims[0]=-1 + kk1;
    capi_rdn_1d_intent |= F2PY_INTENT_IN;
    capi_rdn_1d_tmp = array_from_pyobj(NPY_DOUBLE,rdn_1d_Dims,rdn_1d_Rank,capi_rdn_1d_intent,rdn_1d_capi);
    if (capi_rdn_1d_tmp == NULL) {
        PyObject *exc, *val, *tb;
        PyErr_Fetch(&exc, &val, &tb);
        PyErr_SetString(exc ? exc : balance_error,"failed in converting 14th argument `rdn_1d' of balance.balance to C/Fortran array" );
        npy_PyErr_ChainExceptionsCause(exc, val, tb);
    } else {
        rdn_1d = (double *)(PyArray_DATA(capi_rdn_1d_tmp));

    /* Processing variable ph */
    capi_ph_intent |= (capi_overwrite_ph?0:F2PY_INTENT_COPY);
    ph_Dims[0]=ml,ph_Dims[1]=nl,ph_Dims[2]=kk1;
    capi_ph_intent |= F2PY_INTENT_IN|F2PY_INTENT_OUT;
    capi_ph_tmp = array_from_pyobj(NPY_DOUBLE,ph_Dims,ph_Rank,capi_ph_intent,ph_capi);
    if (capi_ph_tmp == NULL) {
        PyObject *exc, *val, *tb;
        PyErr_Fetch(&exc, &val, &tb);
        PyErr_SetString(exc ? exc : balance_error,"failed in converting 16th argument `ph' of balance.balance to C/Fortran array" );
        npy_PyErr_ChainExceptionsCause(exc, val, tb);
    } else {
        ph = (double *)(PyArray_DATA(capi_ph_tmp));

    /* Processing variable mu_2d */
    capi_mu_2d_intent |= (capi_overwrite_mu_2d?0:F2PY_INTENT_COPY);
    mu_2d_Dims[0]=ml,mu_2d_Dims[1]=nl;
    capi_mu_2d_intent |= F2PY_INTENT_IN|F2PY_INTENT_OUT;
    capi_mu_2d_tmp = array_from_pyobj(NPY_DOUBLE,mu_2d_Dims,mu_2d_Rank,capi_mu_2d_intent,mu_2d_capi);
    if (capi_mu_2d_tmp == NULL) {
        PyObject *exc, *val, *tb;
        PyErr_Fetch(&exc, &val, &tb);
        PyErr_SetString(exc ? exc : balance_error,"failed in converting 17th argument `mu_2d' of balance.balance to C/Fortran array" );
        npy_PyErr_ChainExceptionsCause(exc, val, tb);
    } else {
        mu_2d = (double *)(PyArray_DATA(capi_mu_2d_tmp));

/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
                (*f2py_func)(&ml,&nl,&kk1,qw,qvapor,ph1,phb,t,tw,mu1,psw,mub_2d,hgt_2d,dnw_1d,znu_1d,rdnw_1d,rdn_1d,&ptop,ph,mu_2d);
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
        if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
        CFUNCSMESS("Building return value.\n");
        capi_buildvalue = Py_BuildValue("NN",capi_ph_tmp,capi_mu_2d_tmp);
/*closepyobjfrom*/
/*end of closepyobjfrom*/
        } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
    }  /*if (capi_mu_2d_tmp == NULL) ... else of mu_2d*/
    /* End of cleaning variable mu_2d */
    }  /*if (capi_ph_tmp == NULL) ... else of ph*/
    /* End of cleaning variable ph */
    if((PyObject *)capi_rdn_1d_tmp!=rdn_1d_capi) {
        Py_XDECREF(capi_rdn_1d_tmp); }
    }  /*if (capi_rdn_1d_tmp == NULL) ... else of rdn_1d*/
    /* End of cleaning variable rdn_1d */
    if((PyObject *)capi_rdnw_1d_tmp!=rdnw_1d_capi) {
        Py_XDECREF(capi_rdnw_1d_tmp); }
    }  /*if (capi_rdnw_1d_tmp == NULL) ... else of rdnw_1d*/
    /* End of cleaning variable rdnw_1d */
    if((PyObject *)capi_znu_1d_tmp!=znu_1d_capi) {
        Py_XDECREF(capi_znu_1d_tmp); }
    }  /*if (capi_znu_1d_tmp == NULL) ... else of znu_1d*/
    /* End of cleaning variable znu_1d */
    if((PyObject *)capi_dnw_1d_tmp!=dnw_1d_capi) {
        Py_XDECREF(capi_dnw_1d_tmp); }
    }  /*if (capi_dnw_1d_tmp == NULL) ... else of dnw_1d*/
    /* End of cleaning variable dnw_1d */
    if((PyObject *)capi_hgt_2d_tmp!=hgt_2d_capi) {
        Py_XDECREF(capi_hgt_2d_tmp); }
    }  /*if (capi_hgt_2d_tmp == NULL) ... else of hgt_2d*/
    /* End of cleaning variable hgt_2d */
    if((PyObject *)capi_mub_2d_tmp!=mub_2d_capi) {
        Py_XDECREF(capi_mub_2d_tmp); }
    }  /*if (capi_mub_2d_tmp == NULL) ... else of mub_2d*/
    /* End of cleaning variable mub_2d */
    if((PyObject *)capi_psw_tmp!=psw_capi) {
        Py_XDECREF(capi_psw_tmp); }
    }  /*if (capi_psw_tmp == NULL) ... else of psw*/
    /* End of cleaning variable psw */
    if((PyObject *)capi_mu1_tmp!=mu1_capi) {
        Py_XDECREF(capi_mu1_tmp); }
    }  /*if (capi_mu1_tmp == NULL) ... else of mu1*/
    /* End of cleaning variable mu1 */
    if((PyObject *)capi_phb_tmp!=phb_capi) {
        Py_XDECREF(capi_phb_tmp); }
    }  /*if (capi_phb_tmp == NULL) ... else of phb*/
    /* End of cleaning variable phb */
    if((PyObject *)capi_ph1_tmp!=ph1_capi) {
        Py_XDECREF(capi_ph1_tmp); }
    }  /*if (capi_ph1_tmp == NULL) ... else of ph1*/
    /* End of cleaning variable ph1 */
    if((PyObject *)capi_tw_tmp!=tw_capi) {
        Py_XDECREF(capi_tw_tmp); }
    }  /*if (capi_tw_tmp == NULL) ... else of tw*/
    /* End of cleaning variable tw */
    if((PyObject *)capi_t_tmp!=t_capi) {
        Py_XDECREF(capi_t_tmp); }
    }  /*if (capi_t_tmp == NULL) ... else of t*/
    /* End of cleaning variable t */
    if((PyObject *)capi_qvapor_tmp!=qvapor_capi) {
        Py_XDECREF(capi_qvapor_tmp); }
    }  /*if (capi_qvapor_tmp == NULL) ... else of qvapor*/
    /* End of cleaning variable qvapor */
    } /*CHECKSCALAR(shape(qw, 2) == -1 + kk1)*/
    } /*if (f2py_success) of kk1*/
    /* End of cleaning variable kk1 */
    } /*CHECKSCALAR(shape(qw, 1) == nl)*/
    } /*if (f2py_success) of nl*/
    /* End of cleaning variable nl */
    } /*CHECKSCALAR(shape(qw, 0) == ml)*/
    } /*if (f2py_success) of ml*/
    /* End of cleaning variable ml */
    if((PyObject *)capi_qw_tmp!=qw_capi) {
        Py_XDECREF(capi_qw_tmp); }
    }  /*if (capi_qw_tmp == NULL) ... else of qw*/
    /* End of cleaning variable qw */
    } /*if (f2py_success) of ptop*/
    /* End of cleaning variable ptop */
/*end of cleanupfrompyobj*/
    if (capi_buildvalue == NULL) {
/*routdebugfailure*/
    } else {
/*routdebugleave*/
    }
    CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
    return capi_buildvalue;
}
/******************************* end of balance *******************************/
/*eof body*/

/******************* See f2py2e/f90mod_rules.py: buildhooks *******************/
/*need_f90modhooks*/

/************** See f2py2e/rules.py: module_rules['modulebody'] **************/

/******************* See f2py2e/common_rules.py: buildhooks *******************/

/*need_commonhooks*/

/**************************** See f2py2e/rules.py ****************************/

static FortranDataDef f2py_routine_defs[] = {
    {"balance",-1,{{-1}},0,(char *)F_FUNC(balance,BALANCE),(f2py_init_func)f2py_rout_balance_balance,doc_f2py_rout_balance_balance},

/*eof routine_defs*/
    {NULL}
};

static PyMethodDef f2py_module_methods[] = {

    {NULL,NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "balance",
    NULL,
    -1,
    f2py_module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit_balance(void) {
    int i;
    PyObject *m,*d, *s, *tmp;
    m = balance_module = PyModule_Create(&moduledef);
    Py_SET_TYPE(&PyFortran_Type, &PyType_Type);
    import_array();
    if (PyErr_Occurred())
        {PyErr_SetString(PyExc_ImportError, "can't initialize module balance (failed to import numpy)"); return m;}
    d = PyModule_GetDict(m);
    s = PyUnicode_FromString("1.23.5");
    PyDict_SetItemString(d, "__version__", s);
    Py_DECREF(s);
    s = PyUnicode_FromString(
        "This module 'balance' is auto-generated with f2py (version:1.23.5).\nFunctions:\n"
"    ph,mu_2d = balance(qw,qvapor,ph1,phb,t,tw,mu1,psw,mub_2d,hgt_2d,dnw_1d,znu_1d,rdnw_1d,rdn_1d,ptop,ph,mu_2d,ml=shape(qw, 0),nl=shape(qw, 1),kk1=1 + shape(qw, 2),overwrite_ph=0,overwrite_mu_2d=0)\n"
".");
    PyDict_SetItemString(d, "__doc__", s);
    Py_DECREF(s);
    s = PyUnicode_FromString("1.23.5");
    PyDict_SetItemString(d, "__f2py_numpy_version__", s);
    Py_DECREF(s);
    balance_error = PyErr_NewException ("balance.error", NULL, NULL);
    /*
     * Store the error object inside the dict, so that it could get deallocated.
     * (in practice, this is a module, so it likely will not and cannot.)
     */
    PyDict_SetItemString(d, "_balance_error", balance_error);
    Py_DECREF(balance_error);
    for(i=0;f2py_routine_defs[i].name!=NULL;i++) {
        tmp = PyFortranObject_NewAsAttr(&f2py_routine_defs[i]);
        PyDict_SetItemString(d, f2py_routine_defs[i].name, tmp);
        Py_DECREF(tmp);
    }

/*eof initf2pywraphooks*/
/*eof initf90modhooks*/

/*eof initcommonhooks*/


#ifdef F2PY_REPORT_ATEXIT
    if (! PyErr_Occurred())
        on_exit(f2py_report_on_exit,(void*)"balance");
#endif
    return m;
}
#ifdef __cplusplus
}
#endif
