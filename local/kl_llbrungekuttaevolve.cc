/* FILE: kl_llbrungekuttaevolve.cc                 -*-Mode: c++-*-
 *
 * For description see the header file.
 *
 */

#include <float.h>
#include <string>

#include "nb.h"
#include "director.h"
#include "kl_timedriver.h"     // KL(m)
#include "simstate.h"
#include "kl_llbrungekuttaevolve.h" // KL(m)
#include "kl_llb_util.h"            // KL(m)
#include "oxswarn.h"                // KL(m) for messages
#include "key.h"
#include "energy.h"             // Needed to make MSVC++ 5 happy
#include "kl_llb_util.h"

#if OOMMF_THREADS
# include "oxsthread.h"
#endif // OOMMF_THREADS

OC_USE_STRING;

// Oxs_Ext registration support
OXS_EXT_REGISTER(Klm_LLB_RKEvolve);

/* End includes */

// Revision information, set via CVS keyword substitution
static const Oxs_WarningMessageRevisionInfo revision_info
  (__FILE__,
   "$Revision: 1.00 $",
   "$Date: 2010/04/01 14:08:00 $",
   "$Author: lebecki $",
   "Kristof M. Lebecki (lebecki(ot)uni-konstanz.de)");

/////////////////////////////////////////////////////////////////////////
///////// DEBUGGING DEBUGGING DEBUGGING DEBUGGING DEBUGGING /////////////
// KL(m)
// These functions are defined in rungekuttaevovle.cc. Just re-use them.
// Other option would be to copy their definition here and preceed with
// "static" keyword.
extern void SpinDiff(const char*,const Oxs_SimState&);
extern void VecDiff(const char*,const Oxs_MeshValue<ThreeVector>&);
extern void RealDiff(const char*,const Oxs_MeshValue<OC_REAL8m>&);
///////// DEBUGGING DEBUGGING DEBUGGING DEBUGGING DEBUGGING /////////////
/////////////////////////////////////////////////////////////////////////

// KL(m)
// Additional debugging
#define KL_DEBUG 1
// For small dM/dt we see unusuall (wrong?) behavior
#define SMALL_DM_DT 1.5e11

// KL(m)
#define FALSE  0
#define TRUE   1
// This is "usual" magnetization saturation. There are some places
/// in the "original" code, where I had to change some settings from
/// m to M, from dm_dt to dM_dt, and so on.
/// Like: allowed_error, start_dm, min_timestep ("Timestep lower bound").
/// Following I define "common" value used for this transition,
/// assuming 1Tesla=800,000A/m to be "common".
/// Units: A/m
#define COMMON_MS 800000.0

#if OOMMF_THREADS

////////////////////////////////////////////////////////////////////////
///////// THREAD A  THREAD A  THREAD A  THREAD A  THREAD A /////////////
////////////////////////////////////////////////////////////////////////
class _Klm_LLB_RKEvolve_RKFBase54_ThreadA : public Oxs_ThreadRunObj {
public:
  static Oxs_Mutex job_control;
  static OC_INT4m offset;

  // Imports
  const Oxs_Mesh* mesh_ptr;
  const Oxs_MeshValue<OC_REAL8m>* Ms_ptr;
  const Oxs_MeshValue<OC_REAL8m>* gamma_ptr;
  const Oxs_MeshValue<OC_REAL8m>* alpha_ptr;
  const Oxs_MeshValue<ThreeVector>* mxH_ptr;
  const Oxs_MeshValue<ThreeVector>* spin_ptr;

  // Exports
  Oxs_MeshValue<ThreeVector>* dm_dt_ptr;
  /// Note: mxH and dm_dt may be same Oxs_MeshValue objects.
  OC_REAL8m dE_dt_sum;
!!! change m->M
  OC_REAL8m max_dm_dt;

  // More imports
  OC_INT4m do_precess;

  // Job control (imports)
  OC_INT4m vecsize;
  OC_INT4m block_size;

  _Klm_LLB_RKEvolve_RKFBase54_ThreadA()
    : mesh_ptr(0), Ms_ptr(0), gamma_ptr(0), alpha_ptr(0), mxH_ptr(0),
!!! change m->M
!!! do we need H_ptr? What for?
      spin_ptr(0), dm_dt_ptr(0), dE_dt_sum(0.0), max_dm_dt(0.0),
      vecsize(0), block_size(0) {}

  void Cmd(int threadnumber, void* data);
};

Oxs_Mutex _Klm_LLB_RKEvolve_RKFBase54_ThreadA::job_control;
OC_INT4m     _Klm_LLB_RKEvolve_RKFBase54_ThreadA::offset(0);

void _Klm_LLB_RKEvolve_RKFBase54_ThreadA::Cmd(int /* threadnumber */,
                                         void* /* data */)
{ // Imports
  const Oxs_Mesh& mesh = *mesh_ptr;
  const Oxs_MeshValue<OC_REAL8m>& Ms    = *Ms_ptr;
  const Oxs_MeshValue<OC_REAL8m>& gamma = *gamma_ptr;
  const Oxs_MeshValue<OC_REAL8m>& alpha = *alpha_ptr;
  const Oxs_MeshValue<ThreeVector>& mxH = *mxH_ptr;
  const Oxs_MeshValue<ThreeVector>& spin = *spin_ptr;

  // Exports
  Oxs_MeshValue<ThreeVector>& dm_dt = *dm_dt_ptr;
  /// Note: mxH and dm_dt may be same Oxs_MeshValue objects.
  dE_dt_sum = 0.0;
  max_dm_dt = 0.0;

  // Support variables
  OC_REAL8m max_dm_dt_sq = 0.0;

  while(1) {
    job_control.Lock();
    OC_INT4m istart = offset;
    OC_INT4m istop = ( offset += block_size );
    job_control.Unlock();

    if(istart>=vecsize) break;
    if(istop>vecsize) istop=vecsize;

    OC_INT4m j;
    const OC_INT4m jstart = istart;
    const OC_INT4m jstop  = istop;

    for(j=jstart;j<jstop;++j) {
      if(Ms[j]==0) {
        dm_dt[j].Set(0.0,0.0,0.0);
      } else {
        OC_REAL8m coef1 = gamma[j];
        OC_REAL8m coef2 = -1*alpha[j]*coef1;
        if(!do_precess) coef1 = 0.0;
        
        ThreeVector scratch1 = mxH[j];
        ThreeVector scratch2 = mxH[j];
        // Note: mxH may be same as dm_dt

        scratch1 *= coef1;   // coef1 == 0 if do_precess if false
        OC_REAL8m mxH_sq = scratch2.MagSq();
        OC_REAL8m dm_dt_sq = mxH_sq*(coef1*coef1+coef2*coef2);
        if(dm_dt_sq>max_dm_dt_sq) {
          max_dm_dt_sq = dm_dt_sq;
        }
!!! tu tez jest rownanie LLB!
        dE_dt_sum += mxH_sq * Ms[j] * mesh.Volume(j) * coef2;
        scratch2 ^= spin[j]; // ((mxH)xm)
        scratch2 *= coef2;  // -alpha.gamma((mxH)xm) = alpha.gamma(mx(mxH))
        dm_dt[j] = scratch1 + scratch2;
      }
    }
  }
  max_dm_dt = sqrt(max_dm_dt_sq);  // Note: This result may slightly
  /// differ from max_i |dm_dt[i]| by rounding errors.
}

void Klm_LLB_RKEvolve::Initialize_Threaded_Calculate_dm_dt
(const Oxs_SimState& state,
 const Oxs_MeshValue<ThreeVector>& mxH,
 Oxs_MeshValue<ThreeVector>& dm_dt,
 vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadA>& thread_data)
{
  const Oxs_Mesh* mesh = state.mesh;
  const OC_INT4m vecsize = mesh->Size();

  // Fill out alpha and gamma meshvalue arrays, as necessary.
  if(mesh_id != mesh->Id() || !gamma.CheckMesh(mesh)
     || !alpha.CheckMesh(mesh)) {
    UpdateMeshArrays(mesh);
  }


  const OC_INT4m MaxThreadCount = Oc_GetMaxThreadCount();
  thread_data.resize(MaxThreadCount);

  _Klm_LLB_RKEvolve_RKFBase54_ThreadA::job_control.Lock();
  _Klm_LLB_RKEvolve_RKFBase54_ThreadA::offset = 0;
  _Klm_LLB_RKEvolve_RKFBase54_ThreadA::job_control.Unlock();

  OC_INT4m ibs = (vecsize + MaxThreadCount - 1)/MaxThreadCount;
  ibs = (ibs + 7)/8;
  if(ibs%8) { ibs += 8 - (ibs%8); }
  if(ibs> static_cast<OC_INT4m>((vecsize*3+3)/4)) ibs = vecsize;

  for(OC_INT4m ithread=0;ithread<MaxThreadCount;++ithread) {
    thread_data[ithread].mesh_ptr = mesh;
    thread_data[ithread].Ms_ptr = state.Ms;
    thread_data[ithread].gamma_ptr = &gamma;
    thread_data[ithread].alpha_ptr = &alpha;
    thread_data[ithread].mxH_ptr = &mxH;
    thread_data[ithread].spin_ptr = &state.spin;
    thread_data[ithread].dm_dt_ptr = &dm_dt;
    thread_data[ithread].do_precess = do_precess;
    thread_data[ithread].vecsize = vecsize;
    thread_data[ithread].block_size = ibs;
  }
}

void Klm_LLB_RKEvolve::Finalize_Threaded_Calculate_dm_dt
(const vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadA>& thread_data, // Import
 OC_REAL8m pE_pt,          // Import
 OC_REAL8m& max_dm_dt,     // Export
 OC_REAL8m& dE_dt,         // Export
 OC_REAL8m& min_timestep_export) // Export
{
  max_dm_dt = 0.0;
  dE_dt = 0.0;
  const OC_INT4m thread_count = thread_data.size();
  for(OC_INT4m ithread=0;ithread<thread_count;++ithread) {
    if(thread_data[ithread].max_dm_dt>max_dm_dt) {
      max_dm_dt = thread_data[ithread].max_dm_dt;
    }
    dE_dt += thread_data[ithread].dE_dt_sum;
  }

  dE_dt = -1 * MU0 * dE_dt + pE_pt;
  /// The first term is (partial E/partial M)*dM/dt, the
  /// second term is (partial E/partial t)*dt/dt.  Note that,
  /// provided Ms_[i]>=0, that by constructions dE_dt_sum above
  /// is always non-negative, so dE_dt_ can only be made positive
  /// by positive pE_pt_.

  // Get bound on smallest stepsize that would actually
  // change spin new_max_dm_dt_index:
  min_timestep_export = DBL_MAX/64.;
  if(max_dm_dt>1 || OC_REAL8_EPSILON<min_timestep_export*max_dm_dt) {
    min_timestep_export = OC_REAL8_EPSILON/max_dm_dt;
    // A timestep of size min_timestep will be hopelessly lost
    // in roundoff error.  So increase a bit, based on an empirical
    // fudge factor.  This fudge factor can be tested by running a
    // problem with start_dM = 0.  If the evolver can't climb its
    // way out of the stepsize=0 hole, then this fudge factor is too
    // small.  So far, the most challenging examples have been
    // problems with small cells with nearly aligned spins, e.g., in
    // a remanent state with an external field is applied at t=0.
    // Damping ratio doesn't seem to have much effect, either way.
    min_timestep_export *= 64;
  } else {
    // Degenerate case: max_dm_dt_ must be exactly or very nearly
    // zero.  Punt.
    min_timestep_export = 1.0;
  }
}

////////////////////////////////////////////////////////////////////////
///////// THREAD B  THREAD B  THREAD B  THREAD B  THREAD B /////////////
////////////////////////////////////////////////////////////////////////
class _Klm_LLB_RKEvolve_RKFBase54_ThreadB : public Oxs_ThreadRunObj {
public:
  static Oxs_Mutex job_control;
  static OC_INT4m offset;

  // Imports
  const Oxs_Mesh* mesh_ptr;
  const Oxs_MeshValue<OC_REAL8m>* Ms_ptr;
  const Oxs_MeshValue<OC_REAL8m>* gamma_ptr;
  const Oxs_MeshValue<OC_REAL8m>* alpha_ptr;
  const Oxs_MeshValue<ThreeVector>* mxH_ptr;
  const Oxs_MeshValue<ThreeVector>* base_spin_ptr;
  const vector<OC_REAL8m>* b_ptr;       // State advancement
  const vector<const Oxs_MeshValue<ThreeVector>*>* A_ptr; // State advancement

  // Exports
  Oxs_MeshValue<ThreeVector>* dm_dt_ptr;
  /// Note: mxH and dm_dt may be same Oxs_MeshValue objects.
  Oxs_MeshValue<ThreeVector>* work_spin_ptr; // Import and export
  Oxs_MeshValue<ThreeVector>* kn_ptr; // State advancement; kn = b*A
  OC_REAL8m dE_dt_sum;
  OC_REAL8m max_dm_dt;
  OC_REAL8m norm_error;  // State advancement
  // Also sets spins in work_state to base_state.spin + mstep*(b*A)

  // More imports
  OC_REAL8m mstep; // State advancement
  OC_REAL8m b_dm_dt; // State advancement

  OC_INT4m do_precess;

  // Job control (imports)
  OC_INT4m vecsize;
  OC_INT4m block_size;

  _Klm_LLB_RKEvolve_RKFBase54_ThreadB()
    : mesh_ptr(0), Ms_ptr(0), gamma_ptr(0), alpha_ptr(0), mxH_ptr(0),
      base_spin_ptr(0), b_ptr(0), A_ptr(0),
      dm_dt_ptr(0), work_spin_ptr(0), kn_ptr(0),
      dE_dt_sum(0.0), max_dm_dt(0.0),
      norm_error(0.0), mstep(0.0), b_dm_dt(0.0),
      do_precess(1), vecsize(0), block_size(0) {}

  void Cmd(int threadnumber, void* data);
};

Oxs_Mutex _Klm_LLB_RKEvolve_RKFBase54_ThreadB::job_control;
OC_INT4m     _Klm_LLB_RKEvolve_RKFBase54_ThreadB::offset(0);

void _Klm_LLB_RKEvolve_RKFBase54_ThreadB::Cmd(int /* threadnumber */,
                                         void* /* data */)
{
  // Imports
  const Oxs_Mesh& mesh = *mesh_ptr;
  const Oxs_MeshValue<OC_REAL8m>& Ms    = *Ms_ptr;
  const Oxs_MeshValue<OC_REAL8m>& gamma = *gamma_ptr;
  const Oxs_MeshValue<OC_REAL8m>& alpha = *alpha_ptr;
  const Oxs_MeshValue<ThreeVector>& mxH = *mxH_ptr;
  const Oxs_MeshValue<ThreeVector>& base_spin = *base_spin_ptr;
  const vector<OC_REAL8m>& b = *b_ptr;
  const vector<const Oxs_MeshValue<ThreeVector>*>& A = *A_ptr;

  const OC_INDEX term_count = b.size();

  // Exports
  Oxs_MeshValue<ThreeVector>& dm_dt = *dm_dt_ptr;
  Oxs_MeshValue<ThreeVector>& work_spin = *work_spin_ptr;
  /// Note: mxH and dm_dt may be same Oxs_MeshValue objects.
  Oxs_MeshValue<ThreeVector>& kn = *kn_ptr;

  dE_dt_sum = 0.0;
  max_dm_dt = 0.0;
  norm_error = 0.0;

  // Support variables
  OC_REAL8m max_dm_dt_sq = 0.0;

  OC_REAL8m min_normsq0  = DBL_MAX;
  OC_REAL8m max_normsq0  = 0.0;

  OC_REAL8m min_normsq1  = DBL_MAX;
  OC_REAL8m max_normsq1  = 0.0;

  while(1) {
    job_control.Lock();
    OC_INT4m istart = offset;
    OC_INT4m istop = ( offset += block_size );
    job_control.Unlock();

    if(istart>=vecsize) break;
    if(istop>vecsize) istop=vecsize;

    OC_INT4m j;
    const OC_INT4m jstart = istart;
    const OC_INT4m jstop  = istop;

    for(j=jstart; j < jstart + (jstop-jstart)%2 ; ++j) {
      ThreeVector dm_dt_result(0,0,0);
      if(Ms[j]!=0.0) {
        OC_REAL8m coef1 = gamma[j];
        OC_REAL8m coef2 = -1*alpha[j]*coef1;
        if(!do_precess) coef1 = 0.0;
        
        ThreeVector scratch1 = mxH[j];
        ThreeVector scratch2 = mxH[j];
        // Note: mxH may be same as dm_dt

        scratch1 *= coef1;   // coef1 == 0 if do_precess if false
        OC_REAL8m mxH_sq = scratch2.MagSq();
        OC_REAL8m dm_dt_sq = mxH_sq*(coef1*coef1+coef2*coef2);
        if(dm_dt_sq>max_dm_dt_sq) {
          max_dm_dt_sq = dm_dt_sq;
        }
!!! tu tez jest rownanie LLB!
        dE_dt_sum += mxH_sq * Ms[j] * mesh.Volume(j) * coef2;
        scratch2 ^= work_spin[j]; // ((mxH)xm)
        scratch2 *= coef2;  // -alpha.gamma((mxH)xm) = alpha.gamma(mx(mxH))
        dm_dt_result = scratch1 + scratch2;
      }
      dm_dt[j] = dm_dt_result;

      // State advance
      ThreeVector tempspin = dm_dt_result;
      tempspin *= b_dm_dt;
      for(OC_INDEX k=0;k<term_count;++k) {
        tempspin += b[k] * (*(A[k]))[j];
      }
      kn[j] = tempspin;  // Note: kn may be one of the A[k].
      tempspin *= mstep;
      tempspin += base_spin[j];
      OC_REAL8m magsq = tempspin.MakeUnit();
      if(magsq<min_normsq0) min_normsq0=magsq;
      if(magsq>max_normsq0) max_normsq0=magsq;
      work_spin[j] = tempspin;
    }

    for(;j<jstop;j+=2) {
      ThreeVector dm_dt_result0(0,0,0);
      if(Ms[j]!=0.0) {
        OC_REAL8m coef1 = gamma[j];
        OC_REAL8m coef2 = -1*alpha[j]*coef1;
        if(!do_precess) coef1 = 0.0;
        
        ThreeVector scratch1 = mxH[j];
        ThreeVector scratch2 = mxH[j];
        // Note: mxH may be same as dm_dt

        scratch1 *= coef1;   // coef1 == 0 if do_precess if false
        OC_REAL8m mxH_sq = scratch2.MagSq();
        OC_REAL8m dm_dt_sq = mxH_sq*(coef1*coef1+coef2*coef2);
        if(dm_dt_sq>max_dm_dt_sq) {
          max_dm_dt_sq = dm_dt_sq;
        }
!!! tu tez jest rownanie LLB!
        dE_dt_sum += mxH_sq * Ms[j] * mesh.Volume(j) * coef2;
        scratch2 ^= work_spin[j]; // ((mxH)xm)
        scratch2 *= coef2;  // -alpha.gamma((mxH)xm) = alpha.gamma(mx(mxH))
        dm_dt_result0 = scratch1 + scratch2;
      }

      ThreeVector dm_dt_result1(0,0,0);
      if(Ms[j+1]!=0.0) {
        OC_REAL8m coef1 = gamma[j+1];
        OC_REAL8m coef2 = -1*alpha[j+1]*coef1;
        if(!do_precess) coef1 = 0.0;
        
        ThreeVector scratch1 = mxH[j+1];
        ThreeVector scratch2 = mxH[j+1];
        // Note: mxH may be same as dm_dt

        scratch1 *= coef1;   // coef1 == 0 if do_precess if false
        OC_REAL8m mxH_sq = scratch2.MagSq();
        OC_REAL8m dm_dt_sq = mxH_sq*(coef1*coef1+coef2*coef2);
        if(dm_dt_sq>max_dm_dt_sq) {
          max_dm_dt_sq = dm_dt_sq;
        }
!!! tu tez jest rownanie LLB!
        dE_dt_sum += mxH_sq * Ms[j+1] * mesh.Volume(j+1) * coef2;
        scratch2 ^= work_spin[j+1]; // ((mxH)xm)
        scratch2 *= coef2;  // -alpha.gamma((mxH)xm) = alpha.gamma(mx(mxH))
        dm_dt_result1 = scratch1 + scratch2;
      }

      // State advance
      // Note: The code below is designed to keep rounding errors the
      //       same between j and j+1, especially in the case where
      //       extra precision is carried in intermediate floating point
      //       operations.
      ThreeVector tempsum0(0,0,0);
      for(OC_INDEX k=0;k<term_count;++k) {
        const OC_REAL8m bk = b[k];
        tempsum0.Accum(bk,(*(A[k]))[j]);
      }
      tempsum0.Accum(b_dm_dt,dm_dt_result0);
      ThreeVector wspin0 = base_spin[j];
      wspin0.Accum(mstep,tempsum0);

      kn[j]   = tempsum0;
      dm_dt[j]   = dm_dt_result0;
      OC_REAL8m magsq0 = wspin0.MakeUnit();
      work_spin[j]   = wspin0;

      if(magsq0<min_normsq0) min_normsq0=magsq0;
      if(magsq0>max_normsq0) max_normsq0=magsq0;

      ThreeVector tempsum1(0,0,0);
      for(OC_INDEX k=0;k<term_count;++k) {
        const OC_REAL8m bk = b[k];
        tempsum1.Accum(bk,(*(A[k]))[j+1]);
      }
      tempsum1.Accum(b_dm_dt,dm_dt_result1);
      ThreeVector wspin1 = base_spin[j+1];
      wspin1.Accum(mstep,tempsum1);

      kn[j+1] = tempsum1;
      dm_dt[j+1] = dm_dt_result1;
      OC_REAL8m magsq1 = wspin1.MakeUnit();
      work_spin[j+1] = wspin1;

      if(magsq1<min_normsq1) min_normsq1=magsq1; 

      if(magsq1>max_normsq1) max_normsq1=magsq1;
    }
  }

  // Thread-level exports
  if(min_normsq0 <= max_normsq0) {
    // At least one "j" was processed
    if(min_normsq1<min_normsq0) min_normsq0=min_normsq1;
    if(max_normsq1>max_normsq0) max_normsq0=max_normsq1;
    max_dm_dt = sqrt(max_dm_dt_sq);  // Note: This result may slightly
    /// differ from max_j |dm_dt[j]| by rounding errors.
    norm_error = OC_MAX(sqrt(max_normsq0)-1.0,1.0 - sqrt(min_normsq0));
  }

VecDiff("FINI THREAD B spin",work_spin);
VecDiff("             dm_dt",dm_dt);
VecDiff("                kn",kn);
}

void Klm_LLB_RKEvolve::Initialize_Threaded_Calculate_dm_dt
(const Oxs_SimState& base_state, // Import
 Oxs_SimState& work_state,       // Import and export
 const Oxs_MeshValue<ThreeVector>& mxH,
 Oxs_MeshValue<ThreeVector>& dm_dt,
 vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadB>& thread_data,
 OC_REAL8m mstep,
 OC_REAL8m b_dm_dt,
 const vector<OC_REAL8m>& b,
 vector<const Oxs_MeshValue<ThreeVector>*>& A,
 Oxs_MeshValue<ThreeVector>& kn) // export = b*A
{
  const Oxs_Mesh* mesh = base_state.mesh;
  const OC_INT4m vecsize = mesh->Size();

  // Fill out alpha and gamma meshvalue arrays, as necessary.
  if(mesh_id != mesh->Id() || !gamma.CheckMesh(mesh)
     || !alpha.CheckMesh(mesh)) {
    UpdateMeshArrays(mesh);
  }

  const OC_INT4m MaxThreadCount = Oc_GetMaxThreadCount();
  thread_data.resize(MaxThreadCount);

  _Klm_LLB_RKEvolve_RKFBase54_ThreadB::job_control.Lock();
  _Klm_LLB_RKEvolve_RKFBase54_ThreadB::offset = 0;
  _Klm_LLB_RKEvolve_RKFBase54_ThreadB::job_control.Unlock();

  OC_INT4m ibs = (vecsize + MaxThreadCount - 1)/MaxThreadCount;
  ibs = (ibs + 7)/8;
  if(ibs%8) { ibs += 8 - (ibs%8); }
  if(ibs> static_cast<OC_INT4m>((vecsize*3+3)/4)) ibs = vecsize;

SpinDiff("INIT THREAD B",base_state);
RealDiff(" gamma",gamma); /**/
RealDiff(" alpha",alpha); /**/
VecDiff(" mxH",mxH); /**/

  for(OC_INT4m ithread=0;ithread<MaxThreadCount;++ithread) {
    thread_data[ithread].mesh_ptr = mesh;
    thread_data[ithread].Ms_ptr = base_state.Ms;
    thread_data[ithread].gamma_ptr = &gamma;
    thread_data[ithread].alpha_ptr = &alpha;
    thread_data[ithread].mxH_ptr = &mxH;
    thread_data[ithread].base_spin_ptr = &base_state.spin;

    thread_data[ithread].b_ptr = &b;
    thread_data[ithread].A_ptr = &A;

    thread_data[ithread].dm_dt_ptr = &dm_dt;
    thread_data[ithread].work_spin_ptr = &work_state.spin;

    thread_data[ithread].kn_ptr = &kn;

    thread_data[ithread].mstep = mstep;
    thread_data[ithread].b_dm_dt = b_dm_dt;

    thread_data[ithread].do_precess = do_precess;
    thread_data[ithread].vecsize = vecsize;
    thread_data[ithread].block_size = ibs;
  }

  // Adjust state code
  work_state.ClearDerivedData();
  const OC_INDEX term_count = b.size();
  if(term_count != A.size()) {
    OXS_THROW(Oxs_BadParameter,"Failure in"
              " Klm_LLB_RKEvolve::Initialize_Threaded_Calculate_dm_dt:"
              " imports b and A are different lengths.");

  }
  if(term_count == 0) {
    OXS_THROW(Oxs_BadParameter,"Failure in"
              " Klm_LLB_RKEvolve::Initialize_Threaded_Calculate_dm_dt:"
              " imports b and A are empty.");

  }

  if(static_cast<OC_INDEX>(vecsize) != A[0]->Size()) {
    OXS_THROW(Oxs_BadParameter,"Failure in"
              " Klm_LLB_RKEvolve::Initialize_Threaded_Calculate_dm_dt:"
              " imports A[] and mesh are different lengths.");

  }

  kn.AdjustSize(mesh);
}

void Klm_LLB_RKEvolve::Finalize_Threaded_Calculate_dm_dt
(const Oxs_SimState& base_state, // Import
 const vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadB>& thread_data, // Import
 OC_REAL8m pE_pt,             // Import
 OC_REAL8m hstep,             // Import
 Oxs_SimState& work_state, // Export
 OC_REAL8m& max_dm_dt,        // Export
 OC_REAL8m& dE_dt,            // Export
 OC_REAL8m& min_timestep_export,    // Export
 OC_REAL8m& norm_error)       // Export
{
  max_dm_dt = 0.0;
  dE_dt = 0.0;
  norm_error = 0.0;
  const OC_INT4m thread_count = thread_data.size();
  for(OC_INT4m ithread=0;ithread<thread_count;++ithread) {
    if(thread_data[ithread].max_dm_dt>max_dm_dt) {
      max_dm_dt = thread_data[ithread].max_dm_dt;
    }
    dE_dt += thread_data[ithread].dE_dt_sum;
    if(thread_data[ithread].norm_error>norm_error) {
      norm_error = thread_data[ithread].norm_error;
    }
  }

  dE_dt = -1 * MU0 * dE_dt + pE_pt;
  // KL(m) This statement will most probably not be valid any more
  /// The first term is (partial E/partial M)*dM/dt, the
  /// second term is (partial E/partial t)*dt/dt.  Note that,
  /// provided Ms_[i]>=0, that by constructions dE_dt_sum above
  /// is always non-negative, so dE_dt_ can only be made positive
  /// by positive pE_pt_.

  // Get bound on smallest stepsize that would actually
  // change spin new_max_dm_dt_index:
  min_timestep_export = DBL_MAX/64.;
  if(max_dm_dt>1 || OC_REAL8_EPSILON<min_timestep_export*max_dm_dt) {
    min_timestep_export = OC_REAL8_EPSILON/max_dm_dt;
    // A timestep of size min_timestep will be hopelessly lost
    // in roundoff error.  So increase a bit, based on an empirical
    // fudge factor.  This fudge factor can be tested by running a
    // problem with start_dm = 0.  If the evolver can't climb its
    // way out of the stepsize=0 hole, then this fudge factor is too
    // small.  So far, the most challenging examples have been
    // problems with small cells with nearly aligned spins, e.g., in
    // a remanent state with an external field is applied at t=0.
    // Damping ratio doesn't seem to have much effect, either way.
    min_timestep_export *= 64;
  } else {
    // Degenerate case: max_dm_dt_ must be exactly or very nearly
    // zero.  Punt.
    min_timestep_export = 1.0;
  }

  // Finalize state adjustments
  work_state.last_timestep=hstep;

  // Adjust time and iteration fields in new_state
  if(base_state.stage_number != work_state.stage_number) {
    // New stage
    work_state.stage_start_time = base_state.stage_start_time
                                + base_state.stage_elapsed_time;
    work_state.stage_elapsed_time = work_state.last_timestep;
  } else {
    work_state.stage_start_time = base_state.stage_start_time;
    work_state.stage_elapsed_time = base_state.stage_elapsed_time
                                  + work_state.last_timestep;
  }
  // Don't touch iteration counts. (?!)  The problem is that one call
  // to Klm_LLB_RKEvolve::Step() takes 2 half-steps, and it is the
  // result from these half-steps that are used as the export state.
  // If we increment the iteration count each time through here, then
  // the iteration count goes up by 2 for each call to Step().  So
  // instead, we leave iteration count at whatever value was filled
  // in by the Klm_LLB_RKEvolve::NegotiateTimeStep() method.
}

////////////////////////////////////////////////////////////////////////
///////// THREAD C  THREAD C  THREAD C  THREAD C  THREAD C /////////////
////////////////////////////////////////////////////////////////////////
class _Klm_LLB_RKEvolve_RKFBase54_ThreadC : public Oxs_ThreadRunObj {
public:
  static Oxs_Mutex job_control;
  static OC_INT4m offset;

  // Imports
  const Oxs_Mesh* mesh_ptr;
  const Oxs_MeshValue<OC_REAL8m>* Ms_ptr;
  const Oxs_MeshValue<OC_REAL8m>* gamma_ptr;
  const Oxs_MeshValue<OC_REAL8m>* alpha_ptr;
  const Oxs_MeshValue<ThreeVector>* mxH_ptr;
  const Oxs_MeshValue<ThreeVector>* base_spin_ptr;
  const vector<OC_REAL8m>* b1_ptr;      // State advancement
  const vector<OC_REAL8m>* b2_ptr;      // State advancement
  const vector<const Oxs_MeshValue<ThreeVector>*>* A_ptr; // State advancement

  // Exports
  Oxs_MeshValue<ThreeVector>* dm_dt_ptr;
  /// Note: mxH and dm_dt may be same Oxs_MeshValue objects.
  Oxs_MeshValue<ThreeVector>* work_spin_ptr; // Import and export
  Oxs_MeshValue<ThreeVector>* kn_ptr; // State advancement; kn = b2*A
  OC_REAL8m dE_dt_sum;
  OC_REAL8m max_dm_dt;
  OC_REAL8m norm_error;  // State advancement
  // Also sets spins in work_state to base_state.spin + mstep*(b1*A)

  // More imports
  OC_REAL8m mstep; // State advancement
  OC_REAL8m b1_dm_dt;  // State advancement
  OC_REAL8m b2_dm_dt; // State advancement

  OC_INT4m do_precess;

  // Job control (imports)
  OC_INT4m vecsize;
  OC_INT4m block_size; 

  _Klm_LLB_RKEvolve_RKFBase54_ThreadC()
    : mesh_ptr(0), Ms_ptr(0), gamma_ptr(0), alpha_ptr(0), mxH_ptr(0),
      base_spin_ptr(0), b1_ptr(0), b2_ptr(0), A_ptr(0),
      dm_dt_ptr(0), work_spin_ptr(0), kn_ptr(0),
      dE_dt_sum(0.0), max_dm_dt(0.0),
      norm_error(0.0), mstep(0.0), b1_dm_dt(0.0), b2_dm_dt(0.0),
      do_precess(1), vecsize(0), block_size(0) {}

  void Cmd(int threadnumber, void* data);
};

Oxs_Mutex _Klm_LLB_RKEvolve_RKFBase54_ThreadC::job_control;
OC_INT4m     _Klm_LLB_RKEvolve_RKFBase54_ThreadC::offset(0);

void _Klm_LLB_RKEvolve_RKFBase54_ThreadC::Cmd(int /* threadnumber */,
                                         void* /* data */)
{
  // Imports
  const Oxs_Mesh& mesh = *mesh_ptr;
  const Oxs_MeshValue<OC_REAL8m>& Ms    = *Ms_ptr;
  const Oxs_MeshValue<OC_REAL8m>& gamma = *gamma_ptr;
  const Oxs_MeshValue<OC_REAL8m>& alpha = *alpha_ptr;
  const Oxs_MeshValue<ThreeVector>& mxH = *mxH_ptr;
  const Oxs_MeshValue<ThreeVector>& base_spin = *base_spin_ptr;
  const vector<OC_REAL8m>& b1 = *b1_ptr;
  const vector<OC_REAL8m>& b2 = *b2_ptr;
  const vector<const Oxs_MeshValue<ThreeVector>*>& A = *A_ptr;

  const OC_INDEX term_count = b1.size();

  // Exports
  Oxs_MeshValue<ThreeVector>& dm_dt = *dm_dt_ptr;
  Oxs_MeshValue<ThreeVector>& work_spin = *work_spin_ptr;
  /// Note: mxH and dm_dt may be same Oxs_MeshValue objects.
  Oxs_MeshValue<ThreeVector>& kn = *kn_ptr;

  dE_dt_sum = 0.0;
  max_dm_dt = 0.0;
  norm_error = 0.0;

  // Support variables
  OC_REAL8m max_dm_dt_sq = 0.0;

  OC_REAL8m min_normsq  = DBL_MAX;
  OC_REAL8m max_normsq  = 0.0;

  while(1) {
    job_control.Lock();
    OC_INT4m istart = offset;
    OC_INT4m istop = ( offset += block_size );
    job_control.Unlock();

    if(istart>=vecsize) break;
    if(istop>vecsize) istop=vecsize;

    OC_INT4m j;
    const OC_INT4m jstart = istart;
    const OC_INT4m jstop  = istop;

    for(j=jstart; j<jstop ; ++j) {
      ThreeVector dm_dt_result(0,0,0);
      if(Ms[j]!=0.0) {
        OC_REAL8m coef1 = gamma[j];
        OC_REAL8m coef2 = -1*alpha[j]*coef1;
        if(!do_precess) coef1 = 0.0;
        
        ThreeVector scratch1 = mxH[j];
        ThreeVector scratch2 = mxH[j];
        // Note: mxH may be same as dm_dt

        scratch1 *= coef1;   // coef1 == 0 if do_precess if false
        OC_REAL8m mxH_sq = scratch2.MagSq();
        OC_REAL8m dm_dt_sq = mxH_sq*(coef1*coef1+coef2*coef2);
        if(dm_dt_sq>max_dm_dt_sq) {
          max_dm_dt_sq = dm_dt_sq;
        }
!!! tu tez jest rownanie LLB!
        dE_dt_sum += mxH_sq * Ms[j] * mesh.Volume(j) * coef2;
        scratch2 ^= work_spin[j]; // ((mxH)xm)
        scratch2 *= coef2;  // -alpha.gamma((mxH)xm) = alpha.gamma(mx(mxH))
        dm_dt_result = scratch1 + scratch2;
      }

      // State advance
      ThreeVector tempspin1(0,0,0);
      ThreeVector tempspin2(0,0,0);
      for(OC_INDEX k=0;k<term_count;++k) {
        const ThreeVector& tempvec = (*(A[k]))[j]; // note
        tempspin2.Accum(b2[k],tempvec);
        tempspin1.Accum(b1[k],tempvec);
      }
      tempspin1.Accum(b1_dm_dt,dm_dt_result);
      tempspin1 *= mstep;
      tempspin1 += base_spin[j];
      OC_REAL8m magsq = tempspin1.MakeUnit();
      work_spin[j] = tempspin1;

      dm_dt[j] = dm_dt_result;
      tempspin2.Accum(b2_dm_dt,dm_dt_result);
      kn[j] = tempspin2;

      if(magsq<min_normsq) min_normsq=magsq;
      if(magsq>max_normsq) max_normsq=magsq;
    }
  }

  // Thread-level exports
  if(min_normsq <= max_normsq) {
    // At least one "j" was processed
    max_dm_dt = sqrt(max_dm_dt_sq);  // Note: This result may slightly
    /// differ from max_j |dm_dt[j]| by rounding errors.
    norm_error = OC_MAX(sqrt(max_normsq)-1.0,1.0-sqrt(min_normsq));
  }
VecDiff("FINI THREAD C spin",work_spin);
VecDiff("             dm_dt",dm_dt);
VecDiff("                kn",kn);
}

void Klm_LLB_RKEvolve::Initialize_Threaded_Calculate_dm_dt
(const Oxs_SimState& base_state, // Import
 Oxs_SimState& work_state,       // Import and export
 const Oxs_MeshValue<ThreeVector>& mxH,
 Oxs_MeshValue<ThreeVector>& dm_dt,
 vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadC>& thread_data,
 OC_REAL8m mstep,
 OC_REAL8m b1_dm_dt,
 OC_REAL8m b2_dm_dt,
 const vector<OC_REAL8m>& b1,
 const vector<OC_REAL8m>& b2,
 vector<const Oxs_MeshValue<ThreeVector>*>& A,
 Oxs_MeshValue<ThreeVector>& kn) // Export kn = b2*A
  /// Also, work_state spin is set to base_state + mstep*(b1*A);
{
  const Oxs_Mesh* mesh = base_state.mesh;
  const OC_INT4m vecsize = mesh->Size();

  // Fill out alpha and gamma meshvalue arrays, as necessary.
  if(mesh_id != mesh->Id() || !gamma.CheckMesh(mesh)
     || !alpha.CheckMesh(mesh)) {
    UpdateMeshArrays(mesh);
  }

  const OC_INT4m MaxThreadCount = Oc_GetMaxThreadCount();
  thread_data.resize(MaxThreadCount);

  _Klm_LLB_RKEvolve_RKFBase54_ThreadC::job_control.Lock();
  _Klm_LLB_RKEvolve_RKFBase54_ThreadC::offset = 0;
  _Klm_LLB_RKEvolve_RKFBase54_ThreadC::job_control.Unlock();

  OC_INT4m ibs = (vecsize + MaxThreadCount - 1)/MaxThreadCount;
  ibs = (ibs + 7)/8;
  if(ibs%8) { ibs += 8 - (ibs%8); }
  if(ibs> static_cast<OC_INT4m>((vecsize*3+3)/4)) ibs = vecsize;

  for(OC_INT4m ithread=0;ithread<MaxThreadCount;++ithread) {
    thread_data[ithread].mesh_ptr = mesh;
    thread_data[ithread].Ms_ptr = base_state.Ms;
    thread_data[ithread].gamma_ptr = &gamma;
    thread_data[ithread].alpha_ptr = &alpha;
    thread_data[ithread].mxH_ptr = &mxH;
    thread_data[ithread].base_spin_ptr = &base_state.spin;

    thread_data[ithread].b1_ptr = &b1;
    thread_data[ithread].b2_ptr = &b2;
    thread_data[ithread].A_ptr = &A;

    thread_data[ithread].dm_dt_ptr = &dm_dt;
    thread_data[ithread].work_spin_ptr = &work_state.spin;

    thread_data[ithread].kn_ptr = &kn;

    thread_data[ithread].mstep = mstep;
    thread_data[ithread].b1_dm_dt = b1_dm_dt;
    thread_data[ithread].b2_dm_dt = b2_dm_dt;

    thread_data[ithread].do_precess = do_precess;
    thread_data[ithread].vecsize = vecsize;
    thread_data[ithread].block_size = ibs;
  }

  // Adjust state code
  work_state.ClearDerivedData();
  const OC_INDEX term_count = b1.size();
  if(term_count != b2.size()) {
    OXS_THROW(Oxs_BadParameter,"Failure in"
              " Klm_LLB_RKEvolve::Initialize_Threaded_Calculate_dm_dt:"
              " imports b1 and b2 are different lengths.");

  }
  if(term_count != A.size()) {
    OXS_THROW(Oxs_BadParameter,"Failure in"
              " Klm_LLB_RKEvolve::Initialize_Threaded_Calculate_dm_dt:"
              " imports b1 and A are different lengths.");

  }
  if(term_count == 0) {
    OXS_THROW(Oxs_BadParameter,"Failure in"
              " Klm_LLB_RKEvolve::Initialize_Threaded_Calculate_dm_dt:"
              " imports b and A are empty.");

  }

  if(static_cast<OC_INDEX>(vecsize) != A[0]->Size()) {
    OXS_THROW(Oxs_BadParameter,"Failure in"
              " Klm_LLB_RKEvolve::Initialize_Threaded_Calculate_dm_dt:"
              " imports A[] and mesh are different lengths.");

  }

  kn.AdjustSize(mesh);
}

void Klm_LLB_RKEvolve::Finalize_Threaded_Calculate_dm_dt
(const Oxs_SimState& base_state, // Import
 const vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadC>& thread_data, // Import
 OC_REAL8m pE_pt,             // Import
 OC_REAL8m hstep,             // Import
 Oxs_SimState& work_state, // Export
 OC_REAL8m& max_dm_dt,        // Export
 OC_REAL8m& dE_dt,            // Export
 OC_REAL8m& min_timestep_export,    // Export
 OC_REAL8m& norm_error)       // Export
{
  max_dm_dt = 0.0;
  dE_dt = 0.0;
  norm_error = 0.0;
  const OC_INT4m thread_count = thread_data.size();
  for(OC_INT4m ithread=0;ithread<thread_count;++ithread) {
    if(thread_data[ithread].max_dm_dt>max_dm_dt) {
      max_dm_dt = thread_data[ithread].max_dm_dt;
    }
    dE_dt += thread_data[ithread].dE_dt_sum;
    if(thread_data[ithread].norm_error>norm_error) {
      norm_error = thread_data[ithread].norm_error;
    }
  }

  dE_dt = -1 * MU0 * dE_dt + pE_pt;
  // KL(m) (probably) not valid any more
  /// The first term is (partial E/partial M)*dM/dt, the
  /// second term is (partial E/partial t)*dt/dt.  Note that,
  /// provided Ms_[i]>=0, that by constructions dE_dt_sum above
  /// is always non-negative, so dE_dt_ can only be made positive
  /// by positive pE_pt_.

  // Get bound on smallest stepsize that would actually
  // change spin new_max_dm_dt_index:
  min_timestep_export = DBL_MAX/64.;
  if(max_dm_dt>1 || OC_REAL8_EPSILON<min_timestep_export*max_dm_dt) {
    min_timestep_export = OC_REAL8_EPSILON/max_dm_dt;
    // A timestep of size min_timestep will be hopelessly lost
    // in roundoff error.  So increase a bit, based on an empirical
    // fudge factor.  This fudge factor can be tested by running a
    // problem with start_dm = 0.  If the evolver can't climb its
    // way out of the stepsize=0 hole, then this fudge factor is too
    // small.  So far, the most challenging examples have been
    // problems with small cells with nearly aligned spins, e.g., in
    // a remanent state with an external field is applied at t=0.
    // Damping ratio doesn't seem to have much effect, either way.
    min_timestep_export *= 64;
  } else {
    // Degenerate case: max_dm_dt_ must be exactly or very nearly
    // zero.  Punt.
    min_timestep_export = 1.0;
  }

  // Finalize state adjustments
  work_state.last_timestep=hstep;

  // Adjust time and iteration fields in new_state
  if(base_state.stage_number != work_state.stage_number) {
    // New stage
    work_state.stage_start_time = base_state.stage_start_time
                                + base_state.stage_elapsed_time;
    work_state.stage_elapsed_time = work_state.last_timestep;
  } else {
    work_state.stage_start_time = base_state.stage_start_time;
    work_state.stage_elapsed_time = base_state.stage_elapsed_time
                                  + work_state.last_timestep;
  }
  // Don't touch iteration counts. (?!)  The problem is that one call
  // to Klm_LLB_RKEvolve::Step() takes 2 half-steps, and it is the
  // result from these half-steps that are used as the export state.
  // If we increment the iteration count each time through here, then
  // the iteration count goes up by 2 for each call to Step().  So
  // instead, we leave iteration count at whatever value was filled
  // in by the Klm_LLB_RKEvolve::NegotiateTimeStep() method.
}

////////////////////////////////////////////////////////////////////////
///////// THREAD D  THREAD D  THREAD D  THREAD D  THREAD D /////////////
////////////////////////////////////////////////////////////////////////
class _Klm_LLB_RKEvolve_RKFBase54_ThreadD : public Oxs_ThreadRunObj {
public:
  static Oxs_Mutex job_control;
  static OC_INT4m offset;

  // Imports
  const Oxs_Mesh* mesh_ptr;
  const Oxs_MeshValue<OC_REAL8m>* Ms_ptr;
  const Oxs_MeshValue<OC_REAL8m>* gamma_ptr;
  const Oxs_MeshValue<OC_REAL8m>* alpha_ptr;
  const Oxs_MeshValue<ThreeVector>* mxH_ptr;
  const Oxs_MeshValue<ThreeVector>* work_spin_ptr;
  const Oxs_MeshValue<ThreeVector>* dD13456_ptr;

  // Exports
  Oxs_MeshValue<ThreeVector>* dm_dt_ptr;
  /// Note: mxH and dm_dt may be same Oxs_MeshValue objects.
  OC_REAL8m dE_dt_sum;
  OC_REAL8m max_dm_dt;
  OC_REAL8m max_dD_magsq;

  // More imports
  OC_REAL8m dc7;

  OC_INT4m do_precess;

  // Job control (imports)
  OC_INT4m vecsize;
  OC_INT4m block_size;

  _Klm_LLB_RKEvolve_RKFBase54_ThreadD()
    : mesh_ptr(0),Ms_ptr(0), gamma_ptr(0), alpha_ptr(0), mxH_ptr(0),
      work_spin_ptr(0), dD13456_ptr(0), dm_dt_ptr(0),
      dE_dt_sum(0.0), max_dm_dt(0.0),max_dD_magsq(0.0),
      dc7(0.0),
      do_precess(1), vecsize(0), block_size(0) {}

  void Cmd(int threadnumber, void* data);
};

Oxs_Mutex _Klm_LLB_RKEvolve_RKFBase54_ThreadD::job_control;
OC_INT4m     _Klm_LLB_RKEvolve_RKFBase54_ThreadD::offset(0);

void _Klm_LLB_RKEvolve_RKFBase54_ThreadD::Cmd(int /* threadnumber */,
                                         void* /* data */)
{ // This routine computes max_dD_magsq.  This involves computing
  // dm_dt at each mesh point, but that result is not stored.

  // Imports
  const Oxs_Mesh& mesh = *mesh_ptr;
  const Oxs_MeshValue<OC_REAL8m>& Ms    = *Ms_ptr;
  const Oxs_MeshValue<OC_REAL8m>& gamma = *gamma_ptr;
  const Oxs_MeshValue<OC_REAL8m>& alpha = *alpha_ptr;
  const Oxs_MeshValue<ThreeVector>& mxH = *mxH_ptr;
  const Oxs_MeshValue<ThreeVector>& work_spin = *work_spin_ptr;
  const Oxs_MeshValue<ThreeVector>& dD13456 = *dD13456_ptr;

  // Exports
  Oxs_MeshValue<ThreeVector>& dm_dt = *dm_dt_ptr;
  dE_dt_sum = 0.0;
  max_dm_dt = 0.0;
  max_dD_magsq = 0.0;

  // Support variables
  OC_REAL8m max_dm_dt_sq = 0.0;

  while(1) {
    job_control.Lock();
    OC_INT4m istart = offset;
    OC_INT4m istop = ( offset += block_size );
    job_control.Unlock();

    if(istart>=vecsize) break;
    if(istop>vecsize) istop=vecsize;

    OC_INT4m j;
    const OC_INT4m jstart = istart;
    const OC_INT4m jstop  = istop;

    for(j=jstart;j<jstop;++j) {
      ThreeVector dm_dt_result(0,0,0);
      if(Ms[j]!=0.0) {
        OC_REAL8m coef1 = gamma[j];
        OC_REAL8m coef2 = -1*alpha[j]*coef1;
        if(!do_precess) coef1 = 0.0;
        
        ThreeVector scratch1 = mxH[j];
        ThreeVector scratch2 = mxH[j];
        // Note: mxH may be same as dm_dt

        scratch1 *= coef1;   // coef1 == 0 if do_precess if false
        OC_REAL8m mxH_sq = scratch2.MagSq();
        OC_REAL8m dm_dt_sq = mxH_sq*(coef1*coef1+coef2*coef2);
        if(dm_dt_sq>max_dm_dt_sq) {
          max_dm_dt_sq = dm_dt_sq;
        }
!!! tu tez jest rownanie LLB!
        dE_dt_sum += mxH_sq * Ms[j] * mesh.Volume(j) * coef2;
        scratch2 ^= work_spin[j]; // ((mxH)xm)
        scratch2 *= coef2;  // -alpha.gamma((mxH)xm) = alpha.gamma(mx(mxH))
        dm_dt_result = scratch1 + scratch2;
      }
      ThreeVector dD = dm_dt_result;
      dD *= dc7;
      dD += dD13456[j];
      dm_dt[j] = dm_dt_result;
      OC_REAL8m dD_magsq = dD.MagSq();
      if(dD_magsq > max_dD_magsq) {
        max_dD_magsq = dD_magsq;
      }
    }
  }

  // Thread-level exports
  max_dm_dt = sqrt(max_dm_dt_sq);  // Note: This result may slightly
  /// differ from max_j |dm_dt[j]| by rounding errors.

}

void Klm_LLB_RKEvolve::Initialize_Threaded_Calculate_dm_dt
(const Oxs_SimState& work_state,       // Import
 const Oxs_MeshValue<ThreeVector>& mxH,
 Oxs_MeshValue<ThreeVector>& dm_dt,
 vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadD>& thread_data,
 OC_REAL8m dc7,
 const Oxs_MeshValue<ThreeVector>& dD13456)
{
  const Oxs_Mesh* mesh = work_state.mesh;
  const OC_INT4m vecsize = mesh->Size();

  // Fill out alpha and gamma meshvalue arrays, as necessary.
  if(mesh_id != mesh->Id() || !gamma.CheckMesh(mesh)
     || !alpha.CheckMesh(mesh)) {
    UpdateMeshArrays(mesh);
  }

  const OC_INT4m MaxThreadCount = Oc_GetMaxThreadCount();
  thread_data.resize(MaxThreadCount);

  _Klm_LLB_RKEvolve_RKFBase54_ThreadD::job_control.Lock();
  _Klm_LLB_RKEvolve_RKFBase54_ThreadD::offset = 0;
  _Klm_LLB_RKEvolve_RKFBase54_ThreadD::job_control.Unlock();

  OC_INT4m ibs = (vecsize + MaxThreadCount - 1)/MaxThreadCount;
  ibs = (ibs + 7)/8;
  if(ibs%8) { ibs += 8 - (ibs%8); }
  if(ibs> static_cast<OC_INT4m>((vecsize*3+3)/4)) ibs = vecsize;

  for(OC_INT4m ithread=0;ithread<MaxThreadCount;++ithread) {
    thread_data[ithread].mesh_ptr = mesh;
    thread_data[ithread].Ms_ptr = work_state.Ms;
    thread_data[ithread].gamma_ptr = &gamma;
    thread_data[ithread].alpha_ptr = &alpha;
    thread_data[ithread].mxH_ptr = &mxH;
    thread_data[ithread].work_spin_ptr = &work_state.spin;
    thread_data[ithread].dm_dt_ptr = &dm_dt;

    thread_data[ithread].dc7 = dc7;
    thread_data[ithread].dD13456_ptr = &dD13456;

    thread_data[ithread].do_precess = do_precess;
    thread_data[ithread].vecsize = vecsize;
    thread_data[ithread].block_size = ibs;
  }
}

void Klm_LLB_RKEvolve::Finalize_Threaded_Calculate_dm_dt
(const vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadD>& thread_data, // Import
 OC_REAL8m pE_pt,          // Import
 OC_REAL8m& max_dm_dt,      // Export
 OC_REAL8m& dE_dt,          // Export
 OC_REAL8m& min_timestep_export,  // Export
 OC_REAL8m& max_dD_magsq)   // Export
{
  max_dm_dt = 0.0;
  dE_dt = 0.0;
  max_dD_magsq = 0.0;
  const OC_INT4m thread_count = thread_data.size();
  for(OC_INT4m ithread=0;ithread<thread_count;++ithread) {
    if(thread_data[ithread].max_dm_dt>max_dm_dt) {
      max_dm_dt = thread_data[ithread].max_dm_dt;
    }
    dE_dt += thread_data[ithread].dE_dt_sum;
    if(thread_data[ithread].max_dD_magsq>max_dD_magsq) {
      max_dD_magsq = thread_data[ithread].max_dD_magsq;
    }
  }


  dE_dt = -1 * MU0 * dE_dt + pE_pt;
  // KL(m) (probably) not any more
  /// The first term is (partial E/partial M)*dM/dt, the
  /// second term is (partial E/partial t)*dt/dt.  Note that,
  /// provided Ms_[i]>=0, that by constructions dE_dt_sum above
  /// is always non-negative, so dE_dt_ can only be made positive
  /// by positive pE_pt_.

  // Get bound on smallest stepsize that would actually
  // change spin new_max_dm_dt_index:
  min_timestep_export = DBL_MAX/64.;
  if(max_dm_dt>1 || OC_REAL8_EPSILON<min_timestep_export*max_dm_dt) {
    min_timestep_export = OC_REAL8_EPSILON/max_dm_dt;
    // A timestep of size min_timestep will be hopelessly lost
    // in roundoff error.  So increase a bit, based on an empirical
    // fudge factor.  This fudge factor can be tested by running a
    // problem with start_dm = 0.  If the evolver can't climb its
    // way out of the stepsize=0 hole, then this fudge factor is too
    // small.  So far, the most challenging examples have been
    // problems with small cells with nearly aligned spins, e.g., in
    // a remanent state with an external field is applied at t=0.
    // Damping ratio doesn't seem to have much effect, either way.
    min_timestep_export *= 64;
  } else {
    // Degenerate case: max_dm_dt_ must be exactly or very nearly
    // zero.  Punt.
    min_timestep_export = 1.0;
  }
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#endif // OOMMF_THREADS

// Constructor
Klm_LLB_RKEvolve::Klm_LLB_RKEvolve(
  const char* name,     // Child instance id
  Oxs_Director* newdtr, // App director
  const char* argstr)   // MIF input block parameters
  : Klm_TimeEvolverVarMs(name,newdtr,argstr),
    mesh_id(0),
    max_step_decrease(0.03125), max_step_increase_limit(4.0),
    max_step_increase_adj_ratio(1.9),
    reject_goal(0.05), reject_ratio(0.05),
    energy_state_id(0),next_timestep(0.),
    rkstep_ptr(NULL)
{
  // Process arguments
  // KL(m)
  if( !HasInitValue("relative_temperature") )
    throw Oxs_Ext::Error(this,"\nRequired argument missing:"
                              " relative_temperature");
  relative_temperature = GetRealInitValue("relative_temperature",0.);
  if(relative_temperature<0 || relative_temperature>1) {
    char buf[4096];
    Oc_Snprintf(buf,sizeof(buf),
                "\nInvalid parameter value:"
                " Specified relative_temperature is %g"
                " (should be >=0. and <=1.)",
                static_cast<double>(relative_temperature));
    throw Oxs_ExtError(this,buf);
  }
  //
  min_timestep=GetRealInitValue("min_timestep",0.);
  max_timestep=GetRealInitValue("max_timestep",1e-10);
  if(max_timestep<=0.0) {
    char buf[4096];
    Oc_Snprintf(buf,sizeof(buf),
                "Invalid parameter value:"
                " Specified max time step is %g (should be >0.)",
                static_cast<double>(max_timestep));
    throw Oxs_ExtError(this,buf);
  }

  allowed_error_rate = GetRealInitValue("error_rate",1.0*(PI*1e9/180.)*COMMON_MS);
  // KL(m)
  // Changed units, now: A/m*s
  // Default value in radians remains the same as before: 1.0*PI*1e9/180.
  //if(allowed_error_rate>0.0) {
  //  allowed_error_rate *= PI*1e9/180.; // Convert from deg/ns to rad/s
  //}
  
  allowed_absolute_step_error
    = GetRealInitValue("absolute_step_error",0.2*(PI/180.)*COMMON_MS);
  // KL(m)
  // Changed units, now: A/m ( In some sense (A/m)*radian ).
  // Default value in radians remains the same as before: 0.2*PI/180.
  //if(allowed_absolute_step_error>0.0) {
  //  allowed_absolute_step_error *= PI/180.; // Convert from deg to rad
  //}
  allowed_relative_step_error
    = GetRealInitValue("relative_step_error",0.01);

  expected_energy_precision =
    GetRealInitValue("energy_precision",1e-10);

  // KL(m)
  expected_m_max_precision =
    GetRealInitValue("m_max_precision",1e-6);
  expected_m_min =
    GetRealInitValue("m_min",1e-100);

  reject_goal = GetRealInitValue("reject_goal",0.05);
  if(reject_goal<0.) {
    throw Oxs_ExtError(this,"Invalid initialization detected:"
         " \"reject_goal\" value must be non-negative.");
  }

  min_step_headroom = GetRealInitValue("min_step_headroom",0.33);
  if(min_step_headroom<0.) {
    throw Oxs_ExtError(this,"Invalid initialization detected:"
         " \"min_step_headroom\" value must be bigger than 0.");
  }

  max_step_headroom = GetRealInitValue("max_step_headroom",0.95);
  if(max_step_headroom<0.) {
    throw Oxs_ExtError(this,"Invalid initialization detected:"
         " \"max_step_headroom\" value must be bigger than 0.");
  }

  if(min_step_headroom>max_step_headroom) {
    throw Oxs_ExtError(this,"Invalid initialization detected:"
         " \"min_step_headroom\" value must not be larger than"
         " \"max_step_headroom\".");
  }


  if(HasInitValue("alpha")) {
    OXS_GET_INIT_EXT_OBJECT("alpha",Oxs_ScalarField,alpha_init);
  } else {
    alpha_init.SetAsOwner(dynamic_cast<Oxs_ScalarField *>
                          (MakeNew("Oxs_UniformScalarField",director,
                                   "value 0.5")));
  }

  // User may specify either gamma_G (Gilbert) or
  // gamma_LL (Landau-Lifshitz).
  gamma_style = GS_INVALID;
  if(HasInitValue("gamma_G") && HasInitValue("gamma_LL")) {
    throw Oxs_ExtError(this,"Invalid Specify block; "
                         "both gamma_G and gamma_LL specified.");
  } else if(HasInitValue("gamma_G")) {
    gamma_style = GS_G;
    OXS_GET_INIT_EXT_OBJECT("gamma_G",Oxs_ScalarField,gamma_init);
  } else if(HasInitValue("gamma_LL")) {
    gamma_style = GS_LL;
    OXS_GET_INIT_EXT_OBJECT("gamma_LL",Oxs_ScalarField,gamma_init);
  } else {
    gamma_style = GS_G;
    gamma_init.SetAsOwner(dynamic_cast<Oxs_ScalarField *>
                          (MakeNew("Oxs_UniformScalarField",director,
                                   "value -2.211e5")));
  }
  allow_signed_gamma = GetIntInitValue("allow_signed_gamma",0);
  do_precess = GetIntInitValue("do_precess",1);

  // KL(m)
  start_dM = GetRealInitValue("start_dM",0.01*(PI/180.)*COMMON_MS);
  // Units changed, now: (A/m). In some sense it is radian*A/m
  // Default value remains the same in units of radians: 0.01*PI/180
  //start_dM *= PI/180.; // Convert from deg*A/m to rad*A/m

  start_dt = GetRealInitValue("start_dt",max_timestep/8.);

  if(start_dM<0. && start_dt<0.) {
    throw Oxs_ExtError(this,"Invalid initialization detected:"
                       " at least one of \"start_dM\" and \"start_dt\""
                       "  must be nonnegative.");
  }

  stage_init_step_control = SISC_AUTO;  // Safety
  String stage_start = GetStringInitValue("stage_start","auto");
  Oxs_ToLower(stage_start);
  if(stage_start.compare("start_conditions")==0) {
    stage_init_step_control = SISC_START_COND;
  } else if(stage_start.compare("continuous")==0) {
    stage_init_step_control = SISC_CONTINUOUS;
  } else if(stage_start.compare("auto")==0
            || stage_start.compare("automatic")==0) {
    stage_init_step_control = SISC_AUTO;
  } else {
    throw Oxs_ExtError(this,"Invalid initialization detected:"
                         " \"stage_start\" value must be one of"
                         " start_conditions, continuous, or auto.");
  }

  String method = GetStringInitValue("method","rkf54");
  Oxs_ToLower(method); // Do case insensitive match
  if(method.compare("rkf54m")==0) {
    rkstep_ptr = &Klm_LLB_RKEvolve::TakeRungeKuttaFehlbergStep54M;
  } else if(method.compare("rkf54s")==0) {
    rkstep_ptr = &Klm_LLB_RKEvolve::TakeRungeKuttaFehlbergStep54S;
  } else if(method.compare("rkf54")==0) {
    rkstep_ptr = &Klm_LLB_RKEvolve::TakeRungeKuttaFehlbergStep54;
  } else {
    throw Oxs_ExtError(this,"Invalid initialization detected:"
                         " \"method\" value must be one of"
                         " rkf54, rkf54m, or rkf54s.");
  }

  // Setup outputs
  max_dM_dt_output.Setup(this,InstanceName(),"Max dM/dt","A/m*s",0,
     &Klm_LLB_RKEvolve::UpdateDerivedOutputs);
  dE_dt_output.Setup(this,InstanceName(),"dE/dt","J/s",0,
     &Klm_LLB_RKEvolve::UpdateDerivedOutputs);
  delta_E_output.Setup(this,InstanceName(),"Delta E","J",0,
     &Klm_LLB_RKEvolve::UpdateDerivedOutputs);
  // KL(m)
/*  
  min_m_output.Setup(this,InstanceName(),"Min m","",0,
     &Klm_LLB_RKEvolve::UpdateDerivedOutputs);
  max_m_output.Setup(this,InstanceName(),"Max m","",0,
     &Klm_LLB_RKEvolve::UpdateDerivedOutputs);
*/
  //
  dM_dt_output.Setup(this,InstanceName(),"dM/dt","A/m*s",1,
     &Klm_LLB_RKEvolve::UpdateDerivedOutputs);
  mxH_output.Setup(this,InstanceName(),"mxH","A/m",1,
     &Klm_LLB_RKEvolve::UpdateDerivedOutputs);

  max_dM_dt_output.Register(director,-5);
  dE_dt_output.Register(director,-5);
  delta_E_output.Register(director,-5);
  // KL(m)
//  min_m_output.Register(director,-5);
//  max_m_output.Register(director,-5);
  //
  dM_dt_output.Register(director,-5);
  mxH_output.Register(director,-5);

  // dM_dt and mxH output caches are used for intermediate storage,
  // so enable caching.
  dM_dt_output.CacheRequestIncrement(1);
  mxH_output.CacheRequestIncrement(1);

  // LLB consistency
/*
  if(!Activated_LLB_Mode(TRUE)) {
    char buf[4096];
    Oc_Snprintf(buf,sizeof(buf),
                "Evolver (LLB-type) is not consistent"
                " with the rest of the MIF-specification,"
                " like some of the energy terms.");
    throw Oxs_ExtError(this,buf);
  }
*/

  VerifyAllInitArgsUsed();

  // Reserve space for temp_state; see Step() method below
  director->ReserveSimulationStateRequest(1);

#if REPORT_TIME_RKDEVEL
  timer.resize(10);
  timer_counts.resize(10);
#endif
}

OC_BOOL Klm_LLB_RKEvolve::Init()
{
  Klm_TimeEvolverVarMs::Init();

  // LLB_Term consistency
  if(!Get_Use_LLB_Term()) {
    String msg=String("\nYou use LLB driver+evolver. "
      " In such a case you must also use Klm_LLB_Term.");
    throw Oxs_ExtError(this,msg.c_str());
  }
  
#if REPORT_TIME_RKDEVEL
  Oc_TimeVal cpu,wall;
  for(unsigned int ti=0;ti<timer.size();++ti) {
    timer[ti].GetTimes(cpu,wall);
    if(double(wall)>0.0) {
      fprintf(stderr,"               timer %2u ...   %7.2f cpu /%7.2f wall,"
              " (%s/%s)\n",
              ti,double(cpu),double(wall),InstanceName(),
              timer_counts[ti].name.c_str());
      if(timer_counts[ti].pass_count>0) {
        fprintf(stderr,"                \\---> passes = %d,"
                " bytes=%.2f MB, %.2f GB/sec\n",
                timer_counts[ti].pass_count,
                double(timer_counts[ti].bytes)/double(1024*1024),
                double(timer_counts[ti].bytes)
                /(double(1024*1024*1024)*double(wall)));
      }
    }
    timer[ti].Reset(); 
    timer_counts[ti].Reset();
 }
#endif // REPORT_TIME_RKDEVEL

  // Free scratch space allocated by previous problem (if any)
  vtmpA.Release();   vtmpB.Release();
  vtmpC.Release();   vtmpD.Release();

  // KL(m)
  // Free scratch space allocated by previous problem (if any)
  kltmpA.Release();   kltmpB.Release();
  kltmpC.Release();   kltmpD.Release();
  kltmpX.Release();

  // Free memory used by LLG gamma and alpha
  mesh_id = 0;
  alpha.Release();
  gamma.Release();

  max_step_increase = max_step_increase_limit;

  // Setup step_headroom and reject_ratio
  step_headroom = max_step_headroom;
  reject_ratio = reject_goal;

  energy_state_id=0;   // Mark as invalid state
  next_timestep=0.;    // Dummy value
  
  return 1;
}

// Destructor
Klm_LLB_RKEvolve::~Klm_LLB_RKEvolve()
{
#if REPORT_TIME_RKDEVEL
  Oc_TimeVal cpu,wall;
  for(unsigned int ti=0;ti<timer.size();++ti) {
    timer[ti].GetTimes(cpu,wall);
    if(double(wall)>0.0) {
      fprintf(stderr,"               timer %2u ...   %7.2f cpu /%7.2f wall,"
              " (%s/%s)\n",
              ti,double(cpu),double(wall),InstanceName(),
              timer_counts[ti].name.c_str());
      if(timer_counts[ti].pass_count>0) {
        fprintf(stderr,"                \\---> passes = %d,"
                " bytes=%.2f MB, %.2f GB/sec\n",
                timer_counts[ti].pass_count,
                double(timer_counts[ti].bytes)/double(1024*1024),
                double(timer_counts[ti].bytes)
                /(double(1024*1024*1024)*double(wall)));
      }
    }
 }
#endif // REPORT_TIME_RKDEVEL
  
  // LLB_Term consistency
  // We clear the usage flag.
  Clear_Use_LLB_Term();
}

void Klm_LLB_RKEvolve::UpdateMeshArrays(const Oxs_Mesh* mesh)
{
  mesh_id = 0; // Mark update in progress
  const OC_INDEX size = mesh->Size();
  OC_INDEX i;

  alpha_init->FillMeshValue(mesh,alpha);
  gamma_init->FillMeshValue(mesh,gamma);

  if(gamma_style == GS_G) { // Convert to LL form
    for(i=0;i<size;++i) {
      OC_REAL8m cell_alpha = alpha[i];
      gamma[i] /= (1+cell_alpha*cell_alpha);
    }
  }
  if(!allow_signed_gamma) {
    for(i=0;i<size;++i) gamma[i] = -1*fabs(gamma[i]);
  }

  mesh_id = mesh->Id();
}

//void Klm_LLB_RKEvolve::Calculate_dm_dt
void Klm_LLB_RKEvolve::Calculate_dM_dt // See KLnote1
(const Oxs_SimState& state_,
 const Oxs_MeshValue<ThreeVector>& mxH_,
 const Oxs_MeshValue<ThreeVector>& H_,  // KL(m). Additionally needed by LLB
 OC_REAL8m pE_pt_,
 Oxs_MeshValue<ThreeVector>& dM_dt_,    // KL(m) m->M, see KLnote1
 OC_REAL8m& max_dM_dt_,OC_REAL8m& dE_dt_,OC_REAL8m& min_timestep_export)
{ // Imports: state, mxH_, H_, pE_pt_
  // Exports: dM_dt_, max_dM_dt_, dE_dt_
  // NOTE: dM_dt_ is allowed, and in fact is encouraged,
  //   to be the same as mxH_.  In this case, mxH_ is
  //   overwritten by dM_dt on return.
  const Oxs_Mesh* mesh_ = state_.mesh;
  const Oxs_MeshValue<OC_REAL8m>& Ms_ = *(state_.Ms);
  const Oxs_MeshValue<ThreeVector>& spin_ = state_.spin;
  const OC_INDEX size = mesh_->Size(); // Assume import data are compatible
  // KL(m)
  const Oxs_MeshValue<OC_REAL8m>& Ms_T0_ = *(state_.GetPtr_Ms_T0());
#if KL_DEBUG
  if( state_.GetPtr_Ms_T0()==NULL ) {
    String msg="PROGRAMMING ERROR:"
      " state_.GetPtr_Ms_T0()==NULL in Klm_LLB_RKEvolve::Calculate_dM_dt\n";
    throw Oxs_ExtError(this,msg);
  }
#endif // KL_DEBUG

  ThreeVector scratch;

  // Fill out alpha and gamma meshvalue arrays, as necessary.
  if(mesh_id != mesh_->Id() || !gamma.CheckMesh(mesh_)
     || !alpha.CheckMesh(mesh_)) {
    UpdateMeshArrays(mesh_);
  }

  dM_dt_.AdjustSize(mesh_);  // For case &dm_dt_ != &mxH_

  OC_INDEX i;
  OC_REAL8m dE_dt_sum=0.0;
  OC_REAL8m max_dM_dt_sq = 0.0;
  const OC_REAL8m reltemp13 = 1.-relative_temperature/3.; // KL(m)
  const OC_REAL8m reltemp23 = relative_temperature*2./3.; // KL(m)
  for(i=0;i<size;i++) {
    if(Ms_[i]==0.0) {
      dM_dt_[i].Set(0,0,0);
    } else {
#if KL_DEBUG
      const Oxs_ThreeVector mxH_as_in_input = mxH_[i];
#endif
      // KL(m) old LLG and new LLB formulas:
      // LLG  dm/dt  =    PRECESSION -|gamma|alpha            .mx(mxH)
      // LLB  dM/dt  =|M|.PRECESSION -|gamma|alpha(1-t/3)MS_T0.mx(mxH)
      //                             +|gamma|alpha(2t/3) MS_T0.m.(m.H)
      // dE_dt_sum += V.H.dM/dt, where V = cell_volume
      // LLB  V.H.dM/dt = V.{ |gamma|alpha(1-t/3)MS(mxH)^2 + m.H(d|M|/dt) }
      // (as dE_dt_sum will be later multiplied by -MU0)
      //
      // Following comes (almost) original LLG code
      OC_REAL8m coef1 = gamma[i]; // Usually: -fabs(gamma)
      const OC_REAL8m coef2 = -1*alpha[i]*coef1; // |gamma|alpha
      if(!do_precess) coef1 = 0.0;
      ThreeVector scratch1 = mxH_[i];
      ThreeVector scratch2 = mxH_[i];
      // Note: mxH may be same as dm_dt (dM_dt?)
      scratch1 *= coef1;    // -|gamma|mxH
      scratch2 ^= spin_[i]; // ((mxH)xm)
      scratch2 *= coef2;    // -alpha.|gamma|(mx(mxH))
      // End of original code
      const OC_REAL8m mxH_MagSq = mxH_[i].MagSq();
      scratch2 *= reltemp13*Ms_T0_[i];          // KL(m)
      dM_dt_[i] = scratch1*Ms_[i] + scratch2;   // dM/dt, LLG-only part
      const OC_REAL8m mH = spin_[i]*H_[i];         // m.H
      ThreeVector scratch3 = spin_[i];
      scratch3 *= coef2*reltemp23*Ms_T0_[i]*mH; // dM/dt, LLB-only part
              /// alpha|gamma|.(2t/3).MS_T0.m.(m.H)
      dM_dt_[i] += scratch3;                    // dM/dt, whole LLB
      // Tracking maximum change
      const OC_REAL8m dM_dt_sq = dM_dt_[i].MagSq();
      if(dM_dt_sq>max_dM_dt_sq) {
        max_dM_dt_sq = dM_dt_sq;
      }
      // Energy change
      dE_dt_sum += mesh_->Volume(i)*( coef2*reltemp13*Ms_T0_[i]*mxH_MagSq
                                     +H_[i]*scratch3 );
#if KL_DEBUG
      // Consistency checks. They fail to work for small dM/dt
      if(sqrt(dM_dt_[i].MagSq())>SMALL_DM_DT) {
        // ### Other, possible method of energy-difference computation:
        const OC_REAL8m tmp = mesh_->Volume(i)*H_[i]*dM_dt_[i];
        OC_REAL8m rel_err;
        // rel_err will be normalized by tmp.
        // In case tmp is very "small" forget about the whole test.
        // "Small" means: compared to fabs(dE_dt_sum)/(i+1)
        //   (I did not found anything better :(. This cane be zero!)
        if( 0.0==tmp || fabs(tmp)<SMALL_VAL*fabs(dE_dt_sum)/(i+1)) 
          rel_err= 0.0;
        else rel_err= 
          (mesh_->Volume(i)*( coef2*reltemp13*Ms_T0_[i]*mxH_MagSq
                                         +H_[i]*scratch3 ) - tmp)/tmp;
        // KL(m) Maybe: 8 significant digits
        if( fabs(rel_err)>1e-8 ) {
          char buf[1024];
          Oc_Snprintf(buf,sizeof(buf), "Mismatch in dE/dt value (wrong formula?)\n"
            "Usually it would mean: wrong formula. But I have seen it for T=0\n"
	    "simulations. Is it maye the case?"
            "Stage %u, iteration %u, cell index %u, elapsed time %g, dM/dt %g\n"
            "dE/dt as coded:  %.20g\ndE/dt reference: %.20g\n",
            state_.stage_number, state_.iteration_count, i,
            state_.stage_start_time+state_.stage_elapsed_time, sqrt(dM_dt_[i].MagSq()),
            mesh_->Volume(i)*( coef2*reltemp13*Ms_T0_[i]*mxH_MagSq
                              +H_[i]*scratch3 ), 
            tmp);
          static Oxs_WarningMessage foo(3); // how often repeat it
          foo.Send(revision_info,OC_STRINGIFY(__LINE__), buf);
        }
        // mxH computation check
        const ThreeVector scratch9 = spin_[i]^H_[i]; // mxH_manually_computed
        OC_REAL8m rel_err2;
        // rel_err2 will be normalized by mxH_as_in_input.MagSq()
        // In case, where mxH is very "small" we forget about the whole test.
        // "Small" means: compared to m.MagSq()*H.MagSq(), where m.MagSq()==1.
        if( mxH_as_in_input.MagSq() < H_[i].MagSq()*SMALL_VAL*SMALL_VAL ) 
          rel_err2= 0.0;
        else rel_err2= (scratch9-mxH_as_in_input).MagSq()/mxH_as_in_input.MagSq();
        
        // KL(m) Seems to be: 5 significant digits (???), this is a square
        if( rel_err2>1e-10 ) {
          char buf[1024];
          Oc_Snprintf(buf,sizeof(buf), "Mismatch in mxH value (wrong H?)\n"
            "Stage %u, iteration %u, cell index %u, elapsed time %g, dM/dt %g\n"
            "m   %g %g %g\n"
            "H   %g %g %g\n"
            "mxH got from cache (is):\n"
            " %.9g %.9g %.9g\n"
            "mxH manually computed (should be):\n"
            " %.9g %.9g %.9g\n"
            "rel_err2: %g\n",
            state_.stage_number, state_.iteration_count, i,
            state_.stage_start_time+state_.stage_elapsed_time, sqrt(dM_dt_[i].MagSq()),
            spin_[i].x,spin_[i].y,spin_[i].z,
            H_[i].x,H_[i].y,H_[i].z,
            mxH_as_in_input.x,mxH_as_in_input.y,mxH_as_in_input.z,
            scratch9.x,scratch9.y,scratch9.z,
            rel_err2);
          static Oxs_WarningMessage foo(3); // how often repeat it
          foo.Send(revision_info,OC_STRINGIFY(__LINE__), buf);
        }
      }
#endif // KL_DEBUG
    }
  }
  max_dM_dt_ = sqrt(max_dM_dt_sq);

  dE_dt_ = -1 * MU0 * dE_dt_sum + pE_pt_;
  /// The first term is (partial E/partial M)*dM/dt, the
  /// second term is (partial E/partial t)*dt/dt.  Note that,
  /// provided Ms_[i]>=0, that by constructions dE_dt_sum above
  /// is always non-negative, so dE_dt_ can only be made positive
  /// by positive pE_pt_.
  /// KL(m) last sentence is not valid any more

  // Get bound on smallest stepsize that would actually
  // change spin new_max_dm_dt_index:
  min_timestep_export = DBL_MAX/64.;
  if((max_dM_dt_/COMMON_MS)>1 
   || OC_REAL8_EPSILON<min_timestep_export*(max_dM_dt_/COMMON_MS)) {
    min_timestep_export = OC_REAL8_EPSILON/(max_dM_dt_/COMMON_MS);
    // KL(m)
    // Normalized with COMMON_MS. We thus ev. increase mint_timestep
    // This should not be a problem, since MJD writes (28-4-2010):
    //   I don't think it really matters too much how you handle this.  
    //   It is okay to err on the lenient (i.e., bigger value for 
    //   timestep_lower_bound) because if you are up against
    //   this your stepsize is usually on the order of something 
    //   like 1e-20 s, so the configuration is essentially stagnant.  
    //   In real world systems thermal effects would kick it out 
    //   in some random direction, so pick some method and don't 
    //   worry too much about it.
    //
    // A timestep of size min_timestep will be hopelessly lost
    // in roundoff error.  So increase a bit, based on an empirical
    // fudge factor.  This fudge factor can be tested by running a
    // problem with start_dM = 0.  If the evolver can't climb its
    // way out of the stepsize=0 hole, then this fudge factor is too
    // small.  So far, the most challenging examples have been
    // problems with small cells with nearly aligned spins, e.g., in
    // a remanent state with an external field is applied at t=0.
    // Damping ratio doesn't seem to have much effect, either way.
    min_timestep_export *= 64;
  } else {
    // Degenerate case: max_dm_dt_ must be exactly or very nearly
    // zero.  Punt.
    min_timestep_export = 1.0;
  }
  return;
}

void Klm_LLB_RKEvolve::CheckCache(const Oxs_SimState& cstate)
{
  // Pull cached values out from cstate.
  // If cstate.Id() == energy_state_id, then cstate has been run
  // through either this method or UpdateDerivedOutputs.  Either
  // way, all derived state data should be stored in cstate,
  // except currently the "energy" mesh value array, which is
  // stored independently inside *this.  Eventually that should
  // probably be moved in some fashion into cstate too.
  if(energy_state_id != cstate.Id()) {
    // cached data out-of-date
    UpdateDerivedOutputs(cstate);
  }
  OC_BOOL cache_good = 1;
  OC_REAL8m max_dM_dt,dE_dt,delta_E,pE_pt;
//  OC_REAL8m min_m, max_m; // KL(m)
  OC_REAL8m timestep_lower_bound;  // Smallest timestep that can actually
  /// change spin with max_dM_dt (due to OC_REAL8_EPSILON restrictions).
  /// The next timestep is based on the error from the last step.  If
  /// there is no last step (either because this is the first step,
  /// or because the last state handled by this routine is different
  /// from the incoming current_state), then timestep is calculated
  /// so that max_dM_dt * timestep = start_dM  ( KL(m) changed ).

  cache_good &= cstate.GetDerivedData("Max dM/dt",max_dM_dt);
  cache_good &= cstate.GetDerivedData("dE/dt",dE_dt);
  cache_good &= cstate.GetDerivedData("Delta E",delta_E);
  // KL(m)
//  cache_good &= cstate.GetDerivedData("Min m",min_m);
//  cache_good &= cstate.GetDerivedData("Max m",max_m);
  //
  cache_good &= cstate.GetDerivedData("pE/pt",pE_pt);
  cache_good &= cstate.GetDerivedData("Timestep lower bound",
                                      timestep_lower_bound);
  cache_good &= (energy_state_id == cstate.Id());
  cache_good &= (dM_dt_output.cache.state_id == cstate.Id());

  if(!cache_good) {
    throw Oxs_ExtError(this,
       "Klm_LLB_RKEvolve::CheckCache: Invalid data cache.");
  }
}

#if !OOMMF_THREADS
void
Klm_LLB_RKEvolve::AdjustState
// KL(m) Here tables: spin, Ms, Ms_inverse are copied/filled.
// Here also are handled cases related to weird max_m, min_m.
(OC_REAL8m hstep,
 OC_REAL8m mstep,
 const Oxs_SimState& old_state,
// KL(m) See KLnote1
// const Oxs_MeshValue<ThreeVector>& dm_dt,
 const Oxs_MeshValue<ThreeVector>& dM_dt,
 Oxs_SimState& new_state,
 OC_REAL8m& norm_error) const
{
  new_state.ClearDerivedData();

  const Oxs_MeshValue<ThreeVector>& old_spin = old_state.spin;
  Oxs_MeshValue<ThreeVector>& new_spin = new_state.spin;

  // KL(m)
  const Oxs_MeshValue<OC_REAL8m>& old_Ms = (*(old_state.Ms.GetPtr()));
  // We nead write-access to new_Ms
  Oxs_MeshValue<OC_REAL8m>& new_Ms =
    *(const_cast<Oxs_MeshValue<OC_REAL8m>*>(new_state.Ms.GetPtr()));
  Oxs_MeshValue<OC_REAL8m>& new_Ms_inverse =
    *(const_cast<Oxs_MeshValue<OC_REAL8m>*>(new_state.Ms_inverse.GetPtr()));
  //

  if(!dM_dt.CheckMesh(old_state.mesh)) {
    throw Oxs_ExtError(this,
                         "Klm_LLB_RKEvolve::AdjustState:"
                         " Import spin and dM_dt are different sizes.");
  }
  new_spin.AdjustSize(old_state.mesh);
  new_Ms.AdjustSize(old_state.mesh); // KL(m)
  new_Ms_inverse.AdjustSize(old_state.mesh); // KL(m)
  const OC_INDEX size = old_state.mesh->Size();

  ThreeVector tempM; // KL(m) See KLnote1
  OC_INDEX i;
  for(i=0;i<size;++i) {
    tempM = old_spin[i] * old_Ms[i];
    tempM.Accum(mstep,dM_dt[i]);
    const OC_REAL8m newmag = sqrt(tempM.MakeUnit());
    new_spin[i] = tempM; // |tempM|=1
    new_Ms[i] = newmag;
    if(newmag!=0.0)
      new_Ms_inverse[i] = 1./newmag;
    else
      new_Ms_inverse[i] = 0.0;
  }
  norm_error = 0.0; // See KLnote1

  // Adjust time and iteration fields in new_state
  new_state.last_timestep=hstep;
  if(old_state.stage_number != new_state.stage_number) {
    // New stage
    new_state.stage_start_time = old_state.stage_start_time
                                + old_state.stage_elapsed_time;
    new_state.stage_elapsed_time = new_state.last_timestep;
  } else {
    new_state.stage_start_time = old_state.stage_start_time;
    new_state.stage_elapsed_time = old_state.stage_elapsed_time
                                  + new_state.last_timestep;
  }

  // Don't touch iteration counts. (?!)  The problem is that one call
  // to Klm_LLB_RKEvolve::Step() takes 2 half-steps, and it is the
  // result from these half-steps that are used as the export state.
  // If we increment the iteration count each time through here, then
  // the iteration count goes up by 2 for each call to Step().  So
  // instead, we leave iteration count at whatever value was filled
  // in by the Klm_LLB_RKEvolve::NegotiateTimeStep() method.
}

#else // !OOMMF_THREADS

class _Klm_LLB_RKEvolve_RKFBase54_AdjustState : public Oxs_ThreadRunObj {
public:
  static Oxs_Mutex job_control;
  static OC_INT4m offset;

  const Oxs_MeshValue<ThreeVector>* dm_dt_ptr;    // Import
  const Oxs_MeshValue<ThreeVector>* old_spin_ptr; // Import
  Oxs_MeshValue<ThreeVector>* new_spin_ptr;       // Export
  OC_REAL8m mstep;                               // Import
  OC_REAL8m norm_error;                          // Export
  
  // Job control (imports)
  OC_INT4m vecsize;
  OC_INT4m block_size;

  _Klm_LLB_RKEvolve_RKFBase54_AdjustState()
    : dm_dt_ptr(0), old_spin_ptr(0), new_spin_ptr(0),
      mstep(0.0), norm_error(0.0),
      vecsize(0), block_size(0) {}

  void Cmd(int threadnumber, void* data);
};

Oxs_Mutex _Klm_LLB_RKEvolve_RKFBase54_AdjustState::job_control;
OC_INT4m     _Klm_LLB_RKEvolve_RKFBase54_AdjustState::offset(0);

void _Klm_LLB_RKEvolve_RKFBase54_AdjustState::Cmd(int /* threadnumber */,
                                              void* /* data */)
{ // Imports
? change m->M ?
  const Oxs_MeshValue<ThreeVector>& dm_dt = *dm_dt_ptr;
  const Oxs_MeshValue<ThreeVector>& old_spin = *old_spin_ptr;

  // Exports
  Oxs_MeshValue<ThreeVector>& new_spin = *new_spin_ptr;

  // Support variables
  OC_REAL8m min_normsq = DBL_MAX;
  OC_REAL8m max_normsq = 0.0;

  OC_REAL8m min_normsq0 = DBL_MAX;
  OC_REAL8m max_normsq0 = 0.0;

  OC_REAL8m min_normsq1 = DBL_MAX;
  OC_REAL8m max_normsq1 = 0.0;

  while(1) {
    job_control.Lock();
    OC_INT4m istart = offset;
    OC_INT4m istop = ( offset += block_size );
    job_control.Unlock();

    if(istart>=vecsize) break;
    if(istop>vecsize) istop=vecsize;

    OC_INT4m j;
    const OC_INT4m jstart = istart;
    const OC_INT4m jstop  = istop;

    for(j=jstart; j < jstart + (jstop-jstart)%2 ; ++j) {
      ThreeVector tempspin = old_spin[j];
      tempspin.Accum(mstep,dm_dt[j]);
      OC_REAL8m magsq = tempspin.MakeUnit();
      if(magsq<min_normsq) min_normsq=magsq;
      if(magsq>max_normsq) max_normsq=magsq;
      new_spin[j] = tempspin;
    }

    for(;j<jstop;j+=2) {
      ThreeVector tempspin0 = old_spin[j];
      tempspin0.Accum(mstep,dm_dt[j]);
      OC_REAL8m magsq0 = tempspin0.MakeUnit();

      ThreeVector tempspin1 = old_spin[j+1];
      tempspin1.Accum(mstep,dm_dt[j+1]);
      OC_REAL8m magsq1 = tempspin1.MakeUnit();

      new_spin[j]   = tempspin0;
      new_spin[j+1] = tempspin1;

      if(magsq0<min_normsq0) min_normsq0=magsq0;
      if(magsq0>max_normsq0) max_normsq0=magsq0;

      if(magsq1<min_normsq1) min_normsq1=magsq1;
      if(magsq1>max_normsq1) max_normsq1=magsq1;
    }
  }

  if(max_normsq<max_normsq0) max_normsq = max_normsq0;
  if(max_normsq<max_normsq1) max_normsq = max_normsq1;

  if(min_normsq>min_normsq0) min_normsq = min_normsq0;
  if(min_normsq>min_normsq1) min_normsq = min_normsq1;
  
  norm_error = OC_MAX(sqrt(max_normsq)-1.0,1.0 - sqrt(min_normsq));
}

void
Klm_LLB_RKEvolve::AdjustState  // Threaded version
(OC_REAL8m hstep,
 OC_REAL8m mstep,
 const Oxs_SimState& old_state,
 const Oxs_MeshValue<ThreeVector>& dm_dt,
 Oxs_SimState& new_state,
 OC_REAL8m& norm_error) const
{
  new_state.ClearDerivedData();

  const Oxs_MeshValue<ThreeVector>& old_spin = old_state.spin;
  Oxs_MeshValue<ThreeVector>& new_spin = new_state.spin;

!!! You have to complete this coding as in non-threaded version
!!! Or maybe in _Klm_LLB_RKEvolve_RKFBase54_AdjustState ?
  const Oxs_MeshValue<OC_REAL8m>& old_Ms = (*(old_state.Ms.GetPtr()));
  const Oxs_MeshValue<OC_REAL8m>& old_Ms_inverse =
    (*(old_state.Ms_inverse.GetPtr()));
  // We nead write-access to new_Ms
  Oxs_MeshValue<OC_REAL8m>& new_Ms =
    *(const_cast<Oxs_MeshValue<OC_REAL8m>*>(new_state.Ms.GetPtr()));
  Oxs_MeshValue<OC_REAL8m>& new_Ms_inverse =
    *(const_cast<Oxs_MeshValue<OC_REAL8m>*>(new_state.Ms_inverse.GetPtr()));
  //

  if(!dm_dt.CheckMesh(old_state.mesh)) {
    throw Oxs_ExtError(this,
                         "Klm_LLB_RKEvolve::AdjustState:"
                         " Import spin and dm_dt are different sizes.");
  }
  new_spin.AdjustSize(old_state.mesh);
  const OC_INDEX vecsize = old_state.mesh->Size();

  vector<_Klm_LLB_RKEvolve_RKFBase54_AdjustState> thread_data;
  const OC_INT4m MaxThreadCount = Oc_GetMaxThreadCount();
  thread_data.resize(MaxThreadCount);

  _Klm_LLB_RKEvolve_RKFBase54_AdjustState::job_control.Lock();
  _Klm_LLB_RKEvolve_RKFBase54_AdjustState::offset = 0;
  _Klm_LLB_RKEvolve_RKFBase54_AdjustState::job_control.Unlock();

  OC_INT4m ibs = (vecsize + MaxThreadCount - 1)/MaxThreadCount;
  ibs = (ibs + 7)/8;
  if(ibs%8) { ibs += 8 - (ibs%8); }
  if(ibs> static_cast<OC_INT4m>((vecsize*3+3)/4)) ibs = vecsize;

  Oxs_ThreadTree threadtree;
  for(OC_INT4m ithread=0;ithread<MaxThreadCount;++ithread) {
    thread_data[ithread].dm_dt_ptr = &dm_dt;
    thread_data[ithread].old_spin_ptr = &old_spin;
    thread_data[ithread].new_spin_ptr = &new_spin;
    thread_data[ithread].mstep = mstep;
    thread_data[ithread].norm_error = 0.0;
    thread_data[ithread].vecsize= vecsize;
    thread_data[ithread].block_size= ibs;
    if(ithread!=0) threadtree.Launch(thread_data[ithread],0);
  }
  threadtree.LaunchRoot(thread_data[0],0);

  norm_error = -1.0;
  for(OC_INT4m ithread=0;ithread<MaxThreadCount;++ithread) {
    if(thread_data[ithread].norm_error > norm_error) {
      norm_error = thread_data[ithread].norm_error;
    }
  }

  // Adjust time and iteration fields in new_state
  new_state.last_timestep=hstep;
  if(old_state.stage_number != new_state.stage_number) {
    // New stage
    new_state.stage_start_time = old_state.stage_start_time
                                + old_state.stage_elapsed_time;
    new_state.stage_elapsed_time = new_state.last_timestep;
  } else {
    new_state.stage_start_time = old_state.stage_start_time;
    new_state.stage_elapsed_time = old_state.stage_elapsed_time
                                  + new_state.last_timestep;
  }

  // Don't touch iteration counts. (?!)  The problem is that one call
  // to Klm_LLB_RKEvolve::Step() takes 2 half-steps, and it is the
  // result from these half-steps that are used as the export state.
  // If we increment the iteration count each time through here, then
  // the iteration count goes up by 2 for each call to Step().  So
  // instead, we leave iteration count at whatever value was filled
  // in by the Klm_LLB_RKEvolve::NegotiateTimeStep() method.
}
#endif // OOMMF_THREADS

void Klm_LLB_RKEvolve::UpdateTimeFields
(const Oxs_SimState& cstate,
 Oxs_SimState& nstate,
 OC_REAL8m stepsize) const
{
  nstate.last_timestep=stepsize;
  if(cstate.stage_number != nstate.stage_number) {
    // New stage
    nstate.stage_start_time = cstate.stage_start_time
                              + cstate.stage_elapsed_time;
    nstate.stage_elapsed_time = stepsize;
  } else {
    nstate.stage_start_time = cstate.stage_start_time;
    nstate.stage_elapsed_time = cstate.stage_elapsed_time
                                + stepsize;
  }
}

void Klm_LLB_RKEvolve::NegotiateTimeStep
(const Klm_TimeDriver* driver,
 const Oxs_SimState&  cstate,
 Oxs_SimState& nstate,
 OC_REAL8m stepsize,
 OC_BOOL use_start_cond,
 OC_BOOL& force_step,
 OC_BOOL& driver_set_step) const
{ // This routine negotiates with driver over the proper step size.
  // As a side-effect, also initializes the nstate data structure,
  // where nstate is the "next state".

  // Pull needed cached values out from cstate.
  OC_REAL8m max_dM_dt;
  if(!cstate.GetDerivedData("Max dM/dt",max_dM_dt)) {
    throw Oxs_ExtError(this,
       "Klm_LLB_RKEvolve::NegotiateTimeStep: max_dM_dt not cached.");
  }
  OC_REAL8m timestep_lower_bound=0.;  // Smallest timestep that can actually
  /// change spin with max_dM_dt (due to OC_REAL8_EPSILON restrictions).
  if(!cstate.GetDerivedData("Timestep lower bound",
                            timestep_lower_bound)) {
    throw Oxs_ExtError(this,
       "Klm_LLB_RKEvolve::NegotiateTimeStep: "
       " timestep_lower_bound not cached.");
  }
  if(timestep_lower_bound<=0.0) {
    throw Oxs_ExtError(this,
       "Klm_LLB_RKEvolve::NegotiateTimeStep: "
       " cached timestep_lower_bound value not positive.");
  }

      // KL(m)! checks, to-do
      // - |M| is "very small"?
      // basing on it some desition?


  // The next timestep is based on the error from the last step.  If
  // there is no last step (either because this is the first step,
  // or because the last state handled by this routine is different
  // from the incoming current_state), then timestep is calculated
  // so that max_dM_dt * timestep = start_dM  ( changed by KL(m) ),
  // or timestep = start_dt.
  if(use_start_cond || stepsize<=0.0) {
    if(start_dM>=0.0) {
      if(start_dM < sqrt(DBL_MAX/4) * max_dM_dt) {
        stepsize = step_headroom * start_dM / max_dM_dt;
      } else {
        stepsize = sqrt(DBL_MAX/4);
      }
    }
    if(start_dt>=0.0) {
      if(start_dM<0.0 || stepsize>start_dt) {
        stepsize = start_dt;
      }
    }
  }

  // Insure step is not outside requested step bounds
  if(!use_start_cond && stepsize<timestep_lower_bound) {
    // Check for this before max_timestep, so max_timestep can
    // override.  On the one hand, if the timestep is too small to
    // move any spins, then taking the step just advances the
    // simulation time without changing the state; instead what would
    // be small accumulations are just lost.  On the other hand, if
    // the applied field is changing with time, then perhaps all that
    // is needed to get the magnetization moving is to advance the
    // simulation time.  In general, it is hard to tell which is
    // which, so we assume that the user knews what he was doing when
    // he set the max_timestep value (or accepted the default), and it
    // is up to him to notice if the simulation time is advancing
    // without any changes to the magnetization pattern.
    stepsize = timestep_lower_bound;
  }
  if(stepsize>max_timestep) stepsize = max_timestep;
  if(stepsize<min_timestep) stepsize = min_timestep;

  // Negotiate with driver over size of next step
  driver->FillStateMemberData(cstate,nstate);
  UpdateTimeFields(cstate,nstate,stepsize);

  // Update iteration count
  nstate.iteration_count = cstate.iteration_count + 1;
  nstate.stage_iteration_count = cstate.stage_iteration_count + 1;

  // Additional timestep control
  driver->FillStateSupplemental(nstate);

  // Check for forced step.
  // Note: The driver->FillStateSupplemental call may increase the
  //       timestep slightly to match a stage time termination
  //       criteria.  We should tweak the timestep size check
  //       to recognize that changes smaller than a few epsilon
  //       of the stage time are below our effective timescale
  //       resolution.
  force_step = 0;
  OC_REAL8m timestepcheck = nstate.last_timestep
                         - 4*OC_REAL8_EPSILON*nstate.stage_elapsed_time;
  if(timestepcheck<=min_timestep || timestepcheck<=timestep_lower_bound) {
    // Either driver wants to force stepsize, or else step size can't
    // be reduced any further.
    force_step=1;
  }

  // Is driver responsible for stepsize?
  if(nstate.last_timestep == stepsize) driver_set_step = 0;
  else                                 driver_set_step = 1;
}


OC_BOOL
Klm_LLB_RKEvolve::CheckError
(OC_REAL8m global_error_order,
 OC_REAL8m error,     // units: A/m
 OC_REAL8m stepsize,
 OC_REAL8m reference_stepsize,
 OC_REAL8m max_dM_dt, // units: A/m*s
 OC_REAL8m& new_stepsize)
{ // Returns 1 if step is good, 0 if error is too large.
  // Export new_stepsize is set to suggested stepsize
  // for next step.
  //
  // new_stepsize is initially filled with a relative stepsize
  // adjustment ratio, e.g., 1.0 means no change in stepsize.
  // At the end of this routine this ratio is multiplied by
  // stepsize to get the actual absolute stepsize.
  //
  // The import stepsize is the size (in seconds) of the actual
  // step.  The new_stepsize is computed from this based on the
  // various error estimates.  An upper bound is placed on the
  // size of new_stepsize relative to the imports stepsize and
  // reference_stepsize.  reference_stepsize has no effect if
  // it is smaller than max_step_increase*stepsize.  It is
  // usually used only in the case where the stepsize was
  // artificially reduced by, for example, a stage stopping
  // criterion.
  //
  // NOTE: This routine assumes the local error order is
  //     global_error_order + 1.

  OC_BOOL good_step = 1;
  OC_BOOL error_checked=0;

  if(allowed_relative_step_error>=0. || allowed_error_rate>=0.) {
    // Determine tighter rate bound.
    OC_REAL8m rate_error = 0.0;

    // KL(m)
    /// allowed_absolute_step_error in A/m
    ///       allowed_error_rate    in A/m
    ///       rate_error            in A/m
    ///       error                 in A/m

    if(allowed_relative_step_error<0.) {
      rate_error = allowed_error_rate;
    } else if(allowed_error_rate<0.) {
      rate_error = allowed_relative_step_error * max_dM_dt;
    } else {
      rate_error = allowed_relative_step_error * max_dM_dt;
      if(rate_error>allowed_error_rate) {
        rate_error = allowed_error_rate;
      }
    }
    rate_error *= stepsize;

    // Rate check
    if(error>rate_error) {
      good_step = 0;
      new_stepsize = pow(rate_error/error,1.0/global_error_order);
    } else {
      OC_REAL8m ratio = 0.125*DBL_MAX;
      // KL(m)? original line can be left here, I suppose.
      /// I was mainly concerned about the test error>=1
      if(error>=1 || rate_error<ratio*error) {
        OC_REAL8m test_ratio = rate_error/error;
        if(test_ratio<ratio) ratio = test_ratio;
      }
      new_stepsize = pow(ratio,1.0/global_error_order);
    }
    error_checked = 1;
  }

  // Absolute error check
  if(allowed_absolute_step_error>=0.0) {
    OC_REAL8m test_stepsize = 0.0;
    OC_REAL8m local_error_order = global_error_order + 1.0;
    if(error>allowed_absolute_step_error) {
      good_step = 0;
      test_stepsize = pow(allowed_absolute_step_error/error,
                          1.0/local_error_order);
    } else {
      OC_REAL8m ratio = 0.125*DBL_MAX;
      // KL(m)? original line can be left here, I suppose.
      /// I was mainly concerned about the test error>=1
      if(error>=1 || allowed_absolute_step_error<ratio*error) {
        OC_REAL8m test_ratio = allowed_absolute_step_error/error;
        if(test_ratio<ratio) ratio = test_ratio;
      }
      test_stepsize = pow(ratio,1.0/local_error_order);
    }
    if(!error_checked || test_stepsize<new_stepsize) {
      new_stepsize = test_stepsize;
    }
    error_checked = 1;
  }

  if(error_checked) {
    new_stepsize *= step_headroom;
    if(new_stepsize<max_step_decrease) {
      new_stepsize = max_step_decrease*stepsize;
    } else {
      new_stepsize *= stepsize;
      OC_REAL8m step_bound = stepsize * max_step_increase;
      const OC_REAL8m refrat = 0.85;  // Ad hoc value
      if(stepsize<reference_stepsize*refrat) {
        step_bound = OC_MIN(step_bound,reference_stepsize);
      } else if(stepsize<reference_stepsize) {
        OC_REAL8m ref_bound = reference_stepsize + (max_step_increase-1)
          *(stepsize-reference_stepsize*refrat)/(1-refrat);
        /// If stepsize = reference_stepsize*refrat,
        ///     then ref_bound = reference_stepsize
        /// If stepsize = reference_stepsize,
        ///     then ref_bound = step_bound
        /// Otherwise, linear interpolation on stepsize.
        step_bound = OC_MIN(step_bound,ref_bound);
      }
      if(new_stepsize>step_bound) new_stepsize = step_bound;
    }
  } else {
    new_stepsize = stepsize;
  }

  return good_step;
}

void Klm_LLB_RKEvolve::RungeKuttaFehlbergBase54
(RKF_SubType method,
 OC_REAL8m stepsize,
 Oxs_ConstKey<Oxs_SimState> current_state_key,
 // KL(m)!, see note below. ! - just as a mark here, nothing spec to do
 const Oxs_MeshValue<ThreeVector>& current_dM_dt, // See KLnote1
 Oxs_Key<Oxs_SimState>& next_state_key,
 OC_REAL8m& error_estimate, // KL(m) Note: this is in A/m units now. See KLnote1
 OC_REAL8m& global_error_order,
 OC_REAL8m& norm_error, // KL(m) Equal to zero.
 OC_BOOL& new_energy_and_dmdt_computed) // means:
 ///   new_energy_and_dMdt_computed   :)
{ // Runge-Kutta-Fehlberg routine with combined 4th and 5th
  // order Runge-Kutta steps.  The difference between the
  // two results (4th vs. 5th) is used to estimate the error.
  // The largest detected error detected cellsize is returned
  // in export error_estimate.  The units on this are radians.
  // The export global_error_order is set to 4 by this routine.
  // (The local error order is one better, i.e., 5.)  The norm_error
  // export is set to the cellwise maximum deviation from unit norm
  // across all the spins in the final state, before renormalization.

#if REPORT_TIME_RKDEVEL
timer[1].Start(); /**/
timer[2].Start(); /**/
#endif // REPORT_TIME_RKDEVEL

  // The following coefficients appear in
  //
  //   J. R. Dormand and P. J. Prince, ``A family of embedded
  //   Runge-Kutta formulae,'' J. Comp. Appl. Math., 6, 19--26
  //   (1980).
  //
  // They are also listed in J. Stoer and R. Bulirsch's book,
  // ``Introduction to Numerical Analysis,'' Springer, 2nd edition,
  // Sec. 7.2.5, p 454, but there are a number of errors in the S&B
  // account; the reader is recommended to refer directly to the D&P
  // paper.  A follow-up paper,
  //
  //   J. R. Dormand and P. J. Prince, ``A reconsideration of some 
  //   embedded Runge-Kutta formulae,'' J. Comp. Appl. Math., 15,
  //   203--211 (1986)
  //
  // provides formulae with improved stability and higher order.
  // See also
  //
  //   J. H. Williamson, ``Low-Storage Runge-Kutta Schemes,''
  //   J. Comp. Phys., 35, 48--56 (1980).
  //
  // FORMULAE for RK5(4)7FM:
  //
  //     dm_dt1 = dm_dt(t1,m1)
  //     dm_dt2 = dm_dt(t1+ (1/5)*h, m1+h*k1);
  //     dm_dt3 = dm_dt(t1+(3/10)*h, m1+h*k2);
  //     dm_dt4 = dm_dt(t1+ (4/5)*h, m1+h*k3);
  //     dm_dt5 = dm_dt(t1+ (8/9)*h, m1+h*k4);
  //     dm_dt6 = dm_dt(t1+     1*h, m1+h*k5);
  //     dm_dt7 = dm_dt(t1+     1*h, m1+h*k6);
  //  where
  //     k1 = dm_dt1*1/5
  //     k2 = dm_dt1*3/40       + dm_dt2*9/40
  //     k3 = dm_dt1*44/45      - dm_dt2*56/15      + dm_dt3*32/9
  //     k4 = dm_dt1*19372/6561 - dm_dt2*25360/2187 + dm_dt3*64448/6561
  //               - dm_dt4*212/729
  //     k5 = dm_dt1*9017/3168  - dm_dt2*355/33     + dm_dt3*46732/5247
  //               + dm_dt4*49/176   - dm_dt5*5103/18656
  //     k6 = dm_dt1*35/384     + 0                 + dm_dt3*500/1113
  //               + dm_dt4*125/192  - dm_dt5*2187/6784  + dm_dt6*11/84
  // Then
  //     Da = dm_dt1*35/384     + 0 + dm_dt3*500/1113   + dm_dt4*125/192
  //              - dm_dt5*2187/6784      + dm_dt6*11/84
  //     Db = dm_dt1*5179/57600 + 0 + dm_dt3*7571/16695 + dm_dt4*393/640
  //              - dm_dt5*92097/339200   + dm_dt6*187/2100   + dm_dt7*1/40
  // and
  //     m2a = m1 + h*Da
  //     m2b = m1 + h*Db.
  //
  // where m2a is the 5th order estimate, which is the candidate
  // for the next state, and m2b is the 4th order estimate used
  // only for error estimation/stepsize control.  Note that the
  // 4th order estimate uses more dm/dt evaluations (7) than the
  // 5th order method (6).  This is intentional; the coefficients
  // are selected to minimize error (see D&P paper cited above).
  // The extra calculation costs are minimal, because the 7th dm_dt
  // evaluation is at the candidate next state, so it is re-used
  // on the next step unless the step rejected.
  //
  // The error estimate is
  // 
  //     E = |m2b-m2a| = h*|Db-Da| = C*h^6
  //
  // where C is a constant that can be estimated in terms of E.
  // Note that we don't need to know C in order to adjust the
  // stepsize, because stepsize adjustment involves only the
  // ratio of the old stepsize to the new, so C drops out.

  // Scratch space usage:
  //  The import next_state is used for intermediate state
  // storage for all dm_dt computations.  The final computation
  // is for dm_dt7 at m1+h*k6 = m1+h*Da, which is the candidate
  // next state.  (Da=k6; see FSAL note below.)

  // Four temporary arrays, A-D, are used:
  //
  // Step \  Temp
  // Index \ Array:  A         B         C         D
  // ------+---------------------------------------------
  //   1   |      dm_dt2       -         -         -
  //   2   |      dm_dt2      k2         -         -
  //   3   |      dm_dt2     dm_dt3      -         -
  //   4   |      dm_dt2     dm_dt3     k3         -
  //   5   |      dm_dt2     dm_dt3    dm_dt4      -
  //   6   |      dm_dt2     dm_dt3    dm_dt4     k4
  //   7   |      dm_dt2     dm_dt3    dm_dt4    dm_dt5
  //   8   |        k5       dm_dt3    dm_dt4    dm_dt5
  //   9   |      dm_dt6     dm_dt3    dm_dt4    dm_dt5
  //  10   |      k6(3,6)    dD(3,6)   dm_dt4    dm_dt5
  //  11   |      dm_dt7     dD(3,6)   dm_dt4    dm_dt5
  //  12   |      dm_dt7       dD      dm_dt4    dm_dt5
  //
  // Here dD is Db-Da.  We don't need to compute Db directly.
  // The parenthesized numbers, e.g., k6(3,6) indicates
  // a partially formed value.  The total value k6 depends
  // upon dm_dt1, dm_dt3, dm_dt4, dm_dt5, and dm_dt6.  But
  // if we form k6 directly at step 11 in array A, then we
  // lose dm_dt6 which is needed to compute dD.  Rather than
  // use an additional array, we accumulate partial results
  // into arrays A and B for k6 and dD as indicated.
  //   Note that Da = k6.  Further, note that dm_dt7 is
  // dm_dt for the next state candidate.  (This is the 'F',
  // short for 'FSAL' ("first same as last"?) in the method
  // name, RK5(4)7FM.  The 'M' means minimized error norm,
  // 7 is the number of stages, and 5(4) is the main/subsidiary
  // integration formula order.  See the D&P 1986 paper for
  // details and additional references.)

  // Coefficient arrays, a, b, dc, defined by:
  //
  //   dm_dtN = dm_dt(t1+aN*h,m1+h*kN)
  //       kN = \sum_{M=1}^{M=N} dm_dtM*bNM
  //  Db - Da = \sum dm_dtM*dcM
  //

  // KL(m). Note KLnote1
  // Due to idea |Magnetization(time)|<>const I had to introduce one
  // change. Instead of dm_dt=dm/dt there will be dM_dt=dM/dt.
  // Thus, the whole Runge-Kutta machinery will work now with M,
  // instead of m.
  // This affects many (see below) variables / procedures (methods):
  //  (this list is not complete!)
  //   *old*              ->  *new*
  //   current_dm_dt          current_dM_dt
  //   AdjustState(...dM_dt...)
  //   mstep                  <untouched, it has nothing to do with m>
  //   error_norm             0.0
  //   start_dm               start_dM

  OC_REAL8m a1,a2,a3,a4; // a5 and a6 are 1.0
  OC_REAL8m b11,b21,b22,b31,b32,b33,b41,b42,b43,b44,
    b51,b52,b53,b54,b55,b61,b63,b64,b65,b66; // b62 is 0.0
  OC_REAL8m dc1,dc3,dc4,dc5,dc6,dc7;  // c[k] = b6k, and c^[2]=c[2]=0.0,
  /// where c are the coeffs for Da, c^ for Db, and dcM = c^[M]-c[M].

  switch(method) {
  case RK547FC:
    /////////////////////////////////////////////////////////////////
    // Coefficients for Dormand & Prince RK5(4)7FC
    a1 = OC_REAL8m(1.)/OC_REAL8m(5.);
    a2 = OC_REAL8m(3.)/OC_REAL8m(10.);
    a3 = OC_REAL8m(6.)/OC_REAL8m(13.);
    a4 = OC_REAL8m(2.)/OC_REAL8m(3.);
    // a5 and a6 are 1.0

    b11 =      OC_REAL8m(1.)/OC_REAL8m(5.);
  
    b21 =      OC_REAL8m(3.)/OC_REAL8m(40.);
    b22 =      OC_REAL8m(9.)/OC_REAL8m(40.);
  
    b31 =    OC_REAL8m(264.)/OC_REAL8m(2197.);
    b32 =    OC_REAL8m(-90.)/OC_REAL8m(2197.);
    b33 =    OC_REAL8m(840.)/OC_REAL8m(2197.);
  
    b41 =    OC_REAL8m(932.)/OC_REAL8m(3645.);
    b42 =    OC_REAL8m(-14.)/OC_REAL8m(27.);
    b43 =   OC_REAL8m(3256.)/OC_REAL8m(5103.);
    b44 =   OC_REAL8m(7436.)/OC_REAL8m(25515.);
  
    b51 =   OC_REAL8m(-367.)/OC_REAL8m(513.);
    b52 =     OC_REAL8m(30.)/OC_REAL8m(19.);
    b53 =   OC_REAL8m(9940.)/OC_REAL8m(5643.);
    b54 = OC_REAL8m(-29575.)/OC_REAL8m(8208.);
    b55 =   OC_REAL8m(6615.)/OC_REAL8m(3344.);
  
    b61 =     OC_REAL8m(35.)/OC_REAL8m(432.);
    b63 =   OC_REAL8m(8500.)/OC_REAL8m(14553.);
    b64 = OC_REAL8m(-28561.)/OC_REAL8m(84672.);
    b65 =    OC_REAL8m(405.)/OC_REAL8m(704.);
    b66 =     OC_REAL8m(19.)/OC_REAL8m(196.);
    // b62 is 0.0

    // Coefs for error calculation (c^[k] - c[k]).
    // Note that c[k] = b6k, and c^[2]=c[2]=0.0
    dc1 =     OC_REAL8m(11.)/OC_REAL8m(108.)    - b61;
    dc3 =   OC_REAL8m(6250.)/OC_REAL8m(14553.)  - b63;
    dc4 =  OC_REAL8m(-2197.)/OC_REAL8m(21168.)  - b64;
    dc5 =     OC_REAL8m(81.)/OC_REAL8m(176.)    - b65;
    dc6 =    OC_REAL8m(171.)/OC_REAL8m(1960.)   - b66;
    dc7 =      OC_REAL8m(1.)/OC_REAL8m(40.);
    break;

  case RK547FM:
    /////////////////////////////////////////////////////////////////
    // Coefficients for Dormand & Prince RK5(4)7FM
    a1 = OC_REAL8m(1.)/OC_REAL8m(5.);
    a2 = OC_REAL8m(3.)/OC_REAL8m(10.);
    a3 = OC_REAL8m(4.)/OC_REAL8m(5.);
    a4 = OC_REAL8m(8.)/OC_REAL8m(9.);
    // a5 and a6 are 1.0

    b11 =      OC_REAL8m(1.)/OC_REAL8m(5.);
  
    b21 =      OC_REAL8m(3.)/OC_REAL8m(40.);
    b22 =      OC_REAL8m(9.)/OC_REAL8m(40.);
  
    b31 =     OC_REAL8m(44.)/OC_REAL8m(45.);
    b32 =    OC_REAL8m(-56.)/OC_REAL8m(15.);
    b33 =     OC_REAL8m(32.)/OC_REAL8m(9.);
  
    b41 =  OC_REAL8m(19372.)/OC_REAL8m(6561.);
    b42 = OC_REAL8m(-25360.)/OC_REAL8m(2187.);
    b43 =  OC_REAL8m(64448.)/OC_REAL8m(6561.);
    b44 =   OC_REAL8m(-212.)/OC_REAL8m(729.);
  
    b51 =   OC_REAL8m(9017.)/OC_REAL8m(3168.);
    b52 =   OC_REAL8m(-355.)/OC_REAL8m(33.);
    b53 =  OC_REAL8m(46732.)/OC_REAL8m(5247.);
    b54 =     OC_REAL8m(49.)/OC_REAL8m(176.);
    b55 =  OC_REAL8m(-5103.)/OC_REAL8m(18656.);
  
    b61 =     OC_REAL8m(35.)/OC_REAL8m(384.);
    b63 =    OC_REAL8m(500.)/OC_REAL8m(1113.);
    b64 =    OC_REAL8m(125.)/OC_REAL8m(192.);
    b65 =  OC_REAL8m(-2187.)/OC_REAL8m(6784.);
    b66 =     OC_REAL8m(11.)/OC_REAL8m(84.);
    // b62 is 0.0

    // Coefs for error calculation (c^[k] - c[k]).
    // Note that c[k] = b6k, and c^[2]=c[2]=0.0
    dc1 =   OC_REAL8m(5179.)/OC_REAL8m(57600.)  - b61;
    dc3 =   OC_REAL8m(7571.)/OC_REAL8m(16695.)  - b63;
    dc4 =    OC_REAL8m(393.)/OC_REAL8m(640.)    - b64;
    dc5 = OC_REAL8m(-92097.)/OC_REAL8m(339200.) - b65;
    dc6 =    OC_REAL8m(187.)/OC_REAL8m(2100.)   - b66;
    dc7 =      OC_REAL8m(1.)/OC_REAL8m(40.);
    break;
  case RK547FS:
    /////////////////////////////////////////////////////////////////
    // Coefficients for Dormand & Prince RK5(4)7FS
    a1 = OC_REAL8m(2.)/OC_REAL8m(9.);
    a2 = OC_REAL8m(1.)/OC_REAL8m(3.);
    a3 = OC_REAL8m(5.)/OC_REAL8m(9.);
    a4 = OC_REAL8m(2.)/OC_REAL8m(3.);
    // a5 and a6 are 1.0

    b11 =      OC_REAL8m(2.)/OC_REAL8m(9.);
  
    b21 =      OC_REAL8m(1.)/OC_REAL8m(12.);
    b22 =      OC_REAL8m(1.)/OC_REAL8m(4.);
  
    b31 =     OC_REAL8m(55.)/OC_REAL8m(324.);
    b32 =    OC_REAL8m(-25.)/OC_REAL8m(108.);
    b33 =     OC_REAL8m(50.)/OC_REAL8m(81.);
  
    b41 =     OC_REAL8m(83.)/OC_REAL8m(330.);
    b42 =    OC_REAL8m(-13.)/OC_REAL8m(22.);
    b43 =     OC_REAL8m(61.)/OC_REAL8m(66.);
    b44 =      OC_REAL8m(9.)/OC_REAL8m(110.);
  
    b51 =    OC_REAL8m(-19.)/OC_REAL8m(28.);
    b52 =      OC_REAL8m(9.)/OC_REAL8m(4.);
    b53 =      OC_REAL8m(1.)/OC_REAL8m(7.);
    b54 =    OC_REAL8m(-27.)/OC_REAL8m(7.);
    b55 =     OC_REAL8m(22.)/OC_REAL8m(7.);
  
    b61 =     OC_REAL8m(19.)/OC_REAL8m(200.);
    b63 =      OC_REAL8m(3.)/OC_REAL8m(5.);
    b64 =   OC_REAL8m(-243.)/OC_REAL8m(400.);
    b65 =     OC_REAL8m(33.)/OC_REAL8m(40.);
    b66 =      OC_REAL8m(7.)/OC_REAL8m(80.);
    // b62 is 0.0

    // Coefs for error calculation (c^[k] - c[k]).
    // Note that c[k] = b6k, and c^[2]=c[2]=0.0
    dc1 =    OC_REAL8m(431.)/OC_REAL8m(5000.)  - b61;
    dc3 =    OC_REAL8m(333.)/OC_REAL8m(500.)   - b63;
    dc4 =  OC_REAL8m(-7857.)/OC_REAL8m(10000.) - b64;
    dc5 =    OC_REAL8m(957.)/OC_REAL8m(1000.)  - b65;
    dc6 =    OC_REAL8m(193.)/OC_REAL8m(2000.)  - b66;
    dc7 =     OC_REAL8m(-1.)/OC_REAL8m(50.);
    break;
  default:
    throw Oxs_ExtError(this,
                 "Klm_LLB_RKEvolve::RungeKuttaFehlbergBase54:"
                 " Programming error; Invalid sub-type.");
  }

#ifndef NDEBUG
  // COEFFICIENT CHECKS ////////////////////////////////////////
  // Try to catch some simple typing errors.  Oc_Nop calls below force
  // evaluation in order shown; otherwise, some compiler optimizations
  // reorder sums into form with less accuracy.
#define EPS (8*OC_REAL8_EPSILON)  // 6*OC_REAL8_EPSILON should be enough,
  /// but include a little slack compilers with bad numeric taste.
  if(fabs(b11-a1)>EPS ||
     fabs(b21+b22-a2)>EPS ||
     fabs(b31+b32+b33-a3)>EPS ||
     fabs(b41+b42+b43+b44-a4)>EPS ||
     fabs(Oc_Nop(b51+b52+b53) -1.0 + Oc_Nop(b54+b55))>EPS ||
     fabs(b61+b63+b64+b65+b66-1.0)>EPS) {
    char buf[512];
    Oc_Snprintf(buf,sizeof(buf),
                "Coefficient check failed:\n"
                "1: %g\n2: %g\n3: %g\n4: %g\n5: %g\n6: %g",
                static_cast<double>(b11-a1),
                static_cast<double>(b21+b22-a2),
                static_cast<double>(b31+b32+b33-a3),
                static_cast<double>(b41+b42+b43+b44-a4),
                static_cast<double>(b51+b52+b53+b54+b55-1.0),
                static_cast<double>(b61+b63+b64+b65+b66-1.0));
    throw Oxs_ExtError(this,buf);
  }
#endif // NDEBUG

  const Oxs_SimState* cstate = &(current_state_key.GetReadReference());
SpinDiff("\nPt A",*cstate); /**/

  OC_REAL8m pE_pt,max_dM_dt,dE_dt,timestep_lower_bound;
  OC_REAL8m dummy_error;

#if REPORT_TIME_RKDEVEL
timer[2].Stop(); /**/
++timer_counts[2].pass_count;
 timer_counts[2].name = "RKFB54 setup";
#endif // REPORT_TIME_RKDEVEL

#if !OOMMF_THREADS
  // Step 1
  AdjustState(stepsize*a1,stepsize*b11,*cstate,current_dM_dt,
              next_state_key.GetWriteReference(),dummy_error);
  GetEnergyDensity(next_state_key.GetReadReference(),temp_energy,
                   &vtmpA,&kltmpA,pE_pt);
  Calculate_dM_dt(next_state_key.GetReadReference(),vtmpA,kltmpA,pE_pt,
                  vtmpA,max_dM_dt,dE_dt,timestep_lower_bound);

  // Step 2
  vtmpB = current_dM_dt;
  vtmpB *= b21;
  vtmpB.Accumulate(b22,vtmpA);
  AdjustState(stepsize*a2,stepsize,*cstate,vtmpB,
              next_state_key.GetWriteReference(),dummy_error);

  // Step 3
  GetEnergyDensity(next_state_key.GetReadReference(),temp_energy,
                   &vtmpB,&kltmpB,pE_pt);
  Calculate_dM_dt(next_state_key.GetReadReference(),vtmpB,kltmpB,pE_pt,
                  vtmpB,max_dM_dt,dE_dt,timestep_lower_bound);

  // Step 4
  vtmpC = current_dM_dt;
  vtmpC *= b31;
  vtmpC.Accumulate(b32,vtmpA);
  vtmpC.Accumulate(b33,vtmpB);
  AdjustState(stepsize*a3,stepsize,*cstate,vtmpC,
              next_state_key.GetWriteReference(),dummy_error);

  // Step 5
  GetEnergyDensity(next_state_key.GetReadReference(),temp_energy,
                   &vtmpC,&kltmpC,pE_pt);
  Calculate_dM_dt(next_state_key.GetReadReference(),vtmpC,kltmpC,pE_pt,
                  vtmpC,max_dM_dt,dE_dt,timestep_lower_bound);

  // Step 6
  vtmpD = current_dM_dt;
  vtmpD *= b41;
  vtmpD.Accumulate(b42,vtmpA);
  vtmpD.Accumulate(b43,vtmpB);
  vtmpD.Accumulate(b44,vtmpC);
  AdjustState(stepsize*a4,stepsize,*cstate,vtmpD,
              next_state_key.GetWriteReference(),dummy_error);

  // Step 7
  GetEnergyDensity(next_state_key.GetReadReference(),temp_energy,
                   &vtmpD,&kltmpD,pE_pt);
  Calculate_dM_dt(next_state_key.GetReadReference(),vtmpD,kltmpD,pE_pt,
                  vtmpD,max_dM_dt,dE_dt,timestep_lower_bound);
  // Array holdings: A=dm_dt2   B=dm_dt3   C=dm_dt4   D=dm_dt5

  // Step 8
  vtmpA *= b52;
  vtmpA.Accumulate(b51,current_dM_dt);
  vtmpA.Accumulate(b53,vtmpB);
  vtmpA.Accumulate(b54,vtmpC);
  vtmpA.Accumulate(b55,vtmpD);
  AdjustState(stepsize,stepsize,*cstate,vtmpA,
              next_state_key.GetWriteReference(),dummy_error); // a5=1.0

  // Step 9
  GetEnergyDensity(next_state_key.GetReadReference(),temp_energy,
                   &vtmpA,&kltmpA,pE_pt);
  Calculate_dM_dt(next_state_key.GetReadReference(),vtmpA,kltmpA,pE_pt,
                  vtmpA,max_dM_dt,dE_dt,timestep_lower_bound);
  // Array holdings: A=dm_dt6   B=dm_dt3   C=dm_dt4   D=dm_dt5

  // Step 10
  OC_INDEX i;
  const OC_INDEX size = cstate->mesh->Size();
  for(i=0;i<size;i++) {
    ThreeVector dM_dt3 = vtmpB[i];
    ThreeVector dM_dt6 = vtmpA[i];
    vtmpA[i] = b63 * dM_dt3  +  b66 * dM_dt6;  // k6(3,6)
    vtmpB[i] = dc3 * dM_dt3  +  dc6 * dM_dt6;  // dD(3,6)
  }
  // Array holdings: A=k6(3,6)   B=dD(3,6)   C=dm_dt4   D=dm_dt5

  // Step 11
  vtmpA.Accumulate(b61,current_dM_dt);   // Note: b62=0.0
  vtmpA.Accumulate(b64,vtmpC);
  vtmpA.Accumulate(b65,vtmpD);
  AdjustState(stepsize,stepsize,*cstate,vtmpA,
              next_state_key.GetWriteReference(),norm_error); // a6=1.0
  const Oxs_SimState& endstate
    = next_state_key.GetReadReference(); // Candidate next state
  OC_REAL8m total_E;
  GetEnergyDensity(endstate,temp_energy,&mxH_output.cache.value,
                   &kltmpX, // KL(m)
                   pE_pt,total_E);
  mxH_output.cache.state_id=endstate.Id();
  Calculate_dM_dt(endstate,mxH_output.cache.value,kltmpX,pE_pt,
                  vtmpA,max_dM_dt,dE_dt,timestep_lower_bound);

  if(!endstate.AddDerivedData("Timestep lower bound",
                                timestep_lower_bound) ||
     !endstate.AddDerivedData("Max dM/dt",max_dM_dt) ||
     !endstate.AddDerivedData("pE/pt",pE_pt) ||
     !endstate.AddDerivedData("Total E",total_E) ||
     !endstate.AddDerivedData("dE/dt",dE_dt)) {
    throw Oxs_ExtError(this,
                 "Klm_LLB_RKEvolve::RungeKuttaFehlbergBase54:"
                 " Programming error; data cache already set.");
  }
  // Array holdings: A=dm_dt7   B=dD(3,6)   C=dm_dt4   D=dm_dt5

  // Step 12
  OC_REAL8m max_dD_sq=0.0;
  vtmpB.Accumulate(dc1,current_dM_dt);
  vtmpB.Accumulate(dc4,vtmpC);
  vtmpB.Accumulate(dc5,vtmpD);
  vtmpB.Accumulate(dc7,vtmpA);
  // Array holdings: A=dm_dt7   B=dD   C=dm_dt4   D=dm_dt5

  // next_state holds candidate next state, normalized and
  // with proper time field settings; see Step 11.  Note that
  // Step 11 also set norm_error.
  // KL(m) note: norm_error==0 now

  // Error estimate is max|m2a-m2b| = h*max|dD|
  // KL(m) note: this is related to "M" now, not to "m"
  for(i=0;i<size;i++) {
    OC_REAL8m magsq = vtmpB[i].MagSq();
    if(magsq>max_dD_sq) max_dD_sq = magsq;
  }

#else // OOMMF_THREADS
#if REPORT_TIME_RKDEVEL
timer[3].Start(); /**/
#endif // REPORT_TIME_RKDEVEL
  // Steps 1 and 2
  AdjustState(stepsize*a1,stepsize*b11,*cstate,current_dm_dt,
              next_state_key.GetWriteReference(),dummy_error);
#if REPORT_TIME_RKDEVEL
timer[3].Stop(); /**/
++timer_counts[3].pass_count;
 timer_counts[3].bytes += (cstate->mesh->Size())*(9*sizeof(OC_REAL8m));
 timer_counts[3].name = "RKF54 step 1";
#endif // REPORT_TIME_RKDEVEL

SpinDiff("Pt B",next_state_key.GetReadReference()); /**/
  GetEnergyDensity(next_state_key.GetReadReference(),temp_energy,
                   &vtmpA,NULL,pE_pt);
VecDiff("Pt B1 vtmpA",vtmpA); /**/
RealDiff("Pt B1 temp enery",temp_energy); /**/

#if REPORT_TIME_RKDEVEL
timer[4].Start(); /**/
#endif // REPORT_TIME_RKDEVEL
  {
    vector<OC_REAL8m> vb; vector<const Oxs_MeshValue<ThreeVector>*> vs;
    vb.push_back(b21); vs.push_back(&current_dm_dt);
    vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadB> thread_data;
    Oxs_SimState& work_state = next_state_key.GetWriteReference();
    Initialize_Threaded_Calculate_dm_dt(*cstate,work_state,
                                        vtmpA,vtmpA,thread_data,
                                        stepsize /* mstep */,
                                        b22 /* b_dm_dt */,
                                        vb,vs,vtmpB);

    Oxs_ThreadTree threadtree;
    const OC_INT4m thread_count = thread_data.size();
    for(OC_INT4m ithread=1;ithread<thread_count;++ithread) {
      threadtree.Launch(thread_data[ithread],0);
    }
    threadtree.LaunchRoot(thread_data[0],0);
    Finalize_Threaded_Calculate_dm_dt(*cstate,thread_data,pE_pt,
                                      stepsize*a2 /* hstep */,
                                      work_state,
                                      !!! m->M
                                      max_dm_dt,dE_dt,
                                      timestep_lower_bound,dummy_error);
  }
#if REPORT_TIME_RKDEVEL
timer[4].Stop(); /**/
++timer_counts[4].pass_count;
 timer_counts[4].bytes += (cstate->mesh->Size())*(3*(7+1)*sizeof(OC_REAL8m));
 timer_counts[4].name = "RKF54 step 2";
#endif // REPORT_TIME_RKDEVEL

  // Steps 3 and 4
SpinDiff("Pt C",next_state_key.GetReadReference()); /**/
  GetEnergyDensity(next_state_key.GetReadReference(),temp_energy,
                   &vtmpB,NULL,pE_pt);
VecDiff("Pt C1 vtmpB",vtmpB); /**/
#if REPORT_TIME_RKDEVEL
timer[5].Start(); /**/
#endif // REPORT_TIME_RKDEVEL
  {
    vector<OC_REAL8m> vb; vector<const Oxs_MeshValue<ThreeVector>*> vs;
    vb.push_back(b31); vs.push_back(&current_dm_dt);
    vb.push_back(b32); vs.push_back(&vtmpA);
    vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadB> thread_data;
    Oxs_SimState& work_state = next_state_key.GetWriteReference();
    Initialize_Threaded_Calculate_dm_dt(*cstate,work_state,
                                        vtmpB,vtmpB,thread_data,
                                        stepsize /* mstep */,
                                        b33 /* b_dm_dt */,
                                        vb,vs,vtmpC);

    Oxs_ThreadTree threadtree;
    const OC_INT4m thread_count = thread_data.size();
    for(OC_INT4m ithread=1;ithread<thread_count;++ithread) {
      threadtree.Launch(thread_data[ithread],0);
    }
    threadtree.LaunchRoot(thread_data[0],0);
    Finalize_Threaded_Calculate_dm_dt(*cstate,thread_data,pE_pt,
                                      stepsize*a3 /* hstep */,
                                      work_state,
                                      !!! m->M
                                      max_dm_dt,dE_dt,
                                      timestep_lower_bound,dummy_error);
  }
#if REPORT_TIME_RKDEVEL
timer[5].Stop(); /**/
++timer_counts[5].pass_count;
 timer_counts[5].bytes += (cstate->mesh->Size())*(3*(7+2)*sizeof(OC_REAL8m));
 timer_counts[5].name = "RKF54 step 4";
#endif // REPORT_TIME_RKDEVEL

  // Steps 5 and 6
SpinDiff("Pt D",next_state_key.GetReadReference()); /**/
  GetEnergyDensity(next_state_key.GetReadReference(),temp_energy,
                   &vtmpC,NULL,pE_pt);
#if REPORT_TIME_RKDEVEL
timer[6].Start(); /**/
#endif // REPORT_TIME_RKDEVEL
  {
    vector<OC_REAL8m> vb; vector<const Oxs_MeshValue<ThreeVector>*> vs;
    vb.push_back(b41); vs.push_back(&current_dm_dt);
    vb.push_back(b42); vs.push_back(&vtmpA);
    vb.push_back(b43); vs.push_back(&vtmpB);
    vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadB> thread_data;
    Oxs_SimState& work_state = next_state_key.GetWriteReference();
    Initialize_Threaded_Calculate_dm_dt(*cstate,work_state,
                                        vtmpC,vtmpC,thread_data,
                                        stepsize /* mstep */,
                                        b44 /* b_dm_dt */,
                                        vb,vs,vtmpD);

    Oxs_ThreadTree threadtree;
    const OC_INT4m thread_count = thread_data.size();
    for(OC_INT4m ithread=1;ithread<thread_count;++ithread) {
      threadtree.Launch(thread_data[ithread],0);
    }
    threadtree.LaunchRoot(thread_data[0],0);
    Finalize_Threaded_Calculate_dm_dt(*cstate,thread_data,pE_pt,
                                      stepsize*a4 /* hstep */,
                                      work_state,
                                      !!! m->M
                                      max_dm_dt,dE_dt,
                                      timestep_lower_bound,dummy_error);
  }
#if REPORT_TIME_RKDEVEL
timer[6].Stop(); /**/
++timer_counts[6].pass_count;
 timer_counts[6].bytes += (cstate->mesh->Size())*(3*(7+3)*sizeof(OC_REAL8m));
 timer_counts[6].name = "RKF54 step 6";
#endif // REPORT_TIME_RKDEVEL

  // Steps 7 and 8
SpinDiff("Pt E",next_state_key.GetReadReference()); /**/
  GetEnergyDensity(next_state_key.GetReadReference(),temp_energy,
                   &vtmpD,NULL,pE_pt);
#if REPORT_TIME_RKDEVEL
timer[7].Start(); /**/
#endif // REPORT_TIME_RKDEVEL
  {
    vector<OC_REAL8m> vb; vector<const Oxs_MeshValue<ThreeVector>*> vs;
    vb.push_back(b51); vs.push_back(&current_dm_dt);
    vb.push_back(b52); vs.push_back(&vtmpA);
    vb.push_back(b53); vs.push_back(&vtmpB);
    vb.push_back(b54); vs.push_back(&vtmpC);
    vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadB> thread_data;
    Oxs_SimState& work_state = next_state_key.GetWriteReference();
    Initialize_Threaded_Calculate_dm_dt(*cstate,work_state,
                                        vtmpD,vtmpD,thread_data,
                                        stepsize /* mstep */,
                                        b55 /* b_dm_dt */,
                                        vb,vs,vtmpA);

    Oxs_ThreadTree threadtree;
    const OC_INT4m thread_count = thread_data.size();
    for(OC_INT4m ithread=1;ithread<thread_count;++ithread) {
      threadtree.Launch(thread_data[ithread],0);
    }
    threadtree.LaunchRoot(thread_data[0],0);
    Finalize_Threaded_Calculate_dm_dt(*cstate,thread_data,pE_pt,
                                      stepsize,   /* hstep; a5==1 */
                                      work_state,
                                      !!! m->M
                                      max_dm_dt,dE_dt,
                                      timestep_lower_bound,dummy_error);
  }
#if REPORT_TIME_RKDEVEL
timer[7].Stop(); /**/
++timer_counts[7].pass_count;
 timer_counts[7].bytes += (cstate->mesh->Size())*(3*(7+4)*sizeof(OC_REAL8m));
 timer_counts[7].name = "RKF54 step 8";
#endif // REPORT_TIME_RKDEVEL

SpinDiff("Pt F",next_state_key.GetReadReference()); /**/
  GetEnergyDensity(next_state_key.GetReadReference(),temp_energy,
                   &vtmpA,NULL,pE_pt);
VecDiff("Pt F1 vtmpA",vtmpA); /**/
#if REPORT_TIME_RKDEVEL
timer[8].Start(); /**/
#endif // REPORT_TIME_RKDEVEL
  {
    vector<OC_REAL8m> vb1;
    vector<OC_REAL8m> vb2;
    vector<const Oxs_MeshValue<ThreeVector>*> vs;
    vb1.push_back(b61); vb2.push_back(dc1); vs.push_back(&current_dm_dt);
    vb1.push_back(b63); vb2.push_back(dc3); vs.push_back(&vtmpB);
    vb1.push_back(b64); vb2.push_back(dc4); vs.push_back(&vtmpC);
    vb1.push_back(b65); vb2.push_back(dc5); vs.push_back(&vtmpD);

    vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadC> thread_data;
    Oxs_SimState& work_state = next_state_key.GetWriteReference();
    Initialize_Threaded_Calculate_dm_dt(*cstate,work_state,
                                        vtmpA,vtmpA,thread_data,
                                        stepsize /* mstep */,
                                        b66 /* b1_dm_dt */,
                                        dc6 /* b2_dm_dt */,
                                        vb1,vb2,vs,vtmpB);

    Oxs_ThreadTree threadtree;
    const OC_INT4m thread_count = thread_data.size();
    for(OC_INT4m ithread=1;ithread<thread_count;++ithread) {
      threadtree.Launch(thread_data[ithread],0);
    }
    threadtree.LaunchRoot(thread_data[0],0);
    Finalize_Threaded_Calculate_dm_dt(*cstate,thread_data,pE_pt,
                                      stepsize,   /* hstep; a6==1 */
                                      work_state,
                                      !!! m->M
                                      max_dm_dt,dE_dt,
                                      timestep_lower_bound,norm_error);
    // Array holdings: A=dm_dt6   B=dD(1,3,4,5,6)   C=dm_dt4   D=dm_dt5
  }
  const Oxs_SimState& endstate
    = next_state_key.GetReadReference(); // Candidate next state
#if REPORT_TIME_RKDEVEL
timer[8].Stop(); /**/
++timer_counts[8].pass_count;
 timer_counts[8].bytes += (cstate->mesh->Size())*(3*(7+4)*sizeof(OC_REAL8m));
 timer_counts[8].name = "RKF54 step 9";
#endif // REPORT_TIME_RKDEVEL
  OC_REAL8m total_E;
SpinDiff("Pt G",endstate); /**/
  GetEnergyDensity(endstate,temp_energy,&mxH_output.cache.value,
                   NULL,pE_pt,total_E);
  mxH_output.cache.state_id=endstate.Id();
  OC_REAL8m max_dD_sq=0.0;

#if REPORT_TIME_RKDEVEL
timer[9].Start(); /**/
#endif // REPORT_TIME_RKDEVEL
  {
    vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadD> thread_data;
    Initialize_Threaded_Calculate_dm_dt(endstate,mxH_output.cache.value,vtmpA,
                                        thread_data,dc7,vtmpB);
    Oxs_ThreadTree threadtree;
    const OC_INT4m thread_count = thread_data.size();
    for(OC_INT4m ithread=1;ithread<thread_count;++ithread) {
      threadtree.Launch(thread_data[ithread],0);
    }
    threadtree.LaunchRoot(thread_data[0],0);
    Finalize_Threaded_Calculate_dm_dt(thread_data,pE_pt,
                                      !!! m->M
                                      max_dm_dt,dE_dt,
                                      timestep_lower_bound,max_dD_sq);
  }
  // Array holdings: A=dm_dt7   B=dD(1,3,4,5,6)   C=dm_dt4   D=dm_dt5
#if REPORT_TIME_RKDEVEL
timer[9].Stop(); /**/
++timer_counts[9].pass_count;
 timer_counts[9].bytes += (cstate->mesh->Size())*(3*5*sizeof(OC_REAL8m));
 timer_counts[9].name = "RKF54 step 10";
#endif // REPORT_TIME_RKDEVEL

  if(!endstate.AddDerivedData("Timestep lower bound",
                                timestep_lower_bound) ||
                                      !!! m->M
     !endstate.AddDerivedData("Max dm/dt",max_dm_dt) ||
     !endstate.AddDerivedData("pE/pt",pE_pt) ||
     !endstate.AddDerivedData("Total E",total_E) ||
     !endstate.AddDerivedData("dE/dt",dE_dt)) {
    throw Oxs_ExtError(this,
                 "Klm_LLB_RKEvolve::RungeKuttaFehlbergBase54:"
                 " Programming error; data cache already set.");
  }

#endif // OOMMF_THREADS

  error_estimate = stepsize * sqrt(max_dD_sq);
  global_error_order = 5.0;

  new_energy_and_dmdt_computed = 1;

#if REPORT_TIME_RKDEVEL
timer[1].Stop(); /**/
++timer_counts[1].pass_count;
 timer_counts[1].name = "RKFB54 total";
#endif // REPORT_TIME_RKDEVEL

}

void Klm_LLB_RKEvolve::TakeRungeKuttaFehlbergStep54
(OC_REAL8m stepsize,
 Oxs_ConstKey<Oxs_SimState> current_state_key,
 const Oxs_MeshValue<ThreeVector>& current_dM_dt,
 Oxs_Key<Oxs_SimState>& next_state_key,
 OC_REAL8m& error_estimate,OC_REAL8m& global_error_order,
 OC_REAL8m& norm_error,
 OC_BOOL& new_energy_and_dmdt_computed)
{
  RungeKuttaFehlbergBase54(RK547FC,stepsize,
     current_state_key,current_dM_dt,next_state_key,
     error_estimate,global_error_order,norm_error,
     new_energy_and_dmdt_computed);
}

void Klm_LLB_RKEvolve::TakeRungeKuttaFehlbergStep54M
(OC_REAL8m stepsize,
 Oxs_ConstKey<Oxs_SimState> current_state_key,
 const Oxs_MeshValue<ThreeVector>& current_dM_dt,
 Oxs_Key<Oxs_SimState>& next_state_key,
 OC_REAL8m& error_estimate,OC_REAL8m& global_error_order,
 OC_REAL8m& norm_error,
 OC_BOOL& new_energy_and_dmdt_computed)
{
  RungeKuttaFehlbergBase54(RK547FM,stepsize,
     current_state_key,current_dM_dt,next_state_key,
     error_estimate,global_error_order,norm_error,
     new_energy_and_dmdt_computed);
}

void Klm_LLB_RKEvolve::TakeRungeKuttaFehlbergStep54S
(OC_REAL8m stepsize,
 Oxs_ConstKey<Oxs_SimState> current_state_key,
 const Oxs_MeshValue<ThreeVector>& current_dM_dt,
 Oxs_Key<Oxs_SimState>& next_state_key,
 OC_REAL8m& error_estimate,OC_REAL8m& global_error_order,
 OC_REAL8m& norm_error,
 OC_BOOL& new_energy_and_dmdt_computed)
{
  RungeKuttaFehlbergBase54(RK547FS,stepsize,
     current_state_key,current_dM_dt,next_state_key,
     error_estimate,global_error_order,norm_error,
     new_energy_and_dmdt_computed);
}

OC_REAL8m Klm_LLB_RKEvolve::MaxDiff_mul
                    (const Oxs_MeshValue<ThreeVector>& vecA,
                     const Oxs_MeshValue<OC_REAL8m>&      vecA_mul,
		     const Oxs_MeshValue<ThreeVector>& vecB,
                     const Oxs_MeshValue<OC_REAL8m>&      vecB_mul)
{
  OC_INDEX size = vecA.Size();
  if(vecB.Size()!=size || vecB_mul.Size()!=size || vecA_mul.Size()!=size) {
    throw Oxs_ExtError(this,
                 "Klm_LLB_RKEvolve::MaxDiff_mul:"
                 " Import MeshValues incompatible (different lengths).");
  }
  OC_REAL8m max_magsq = 0.0;
  for(OC_INDEX i=0;i<size;i++) {
    ThreeVector vtemp = vecB_mul[i]*vecB[i] - vecA_mul[i]*vecA[i];
    OC_REAL8m magsq = vtemp.MagSq();
    if(magsq>max_magsq) max_magsq = magsq;
  }
  return sqrt(max_magsq);
}

void Klm_LLB_RKEvolve::AdjustStepHeadroom(OC_INT4m step_reject)
{ // step_reject should be 0 or 1, reflecting whether the current
  // step was rejected or not.  This routine updates reject_ratio
  // and adjusts step_headroom appropriately.

  // First adjust reject_ratio, weighing mostly the last
  // thirty or so results.
  reject_ratio = (31*reject_ratio + step_reject)/32.;

  // Adjust step_headroom
  if(reject_ratio>reject_goal && step_reject>0) {
    // Reject ratio too high and getting worse
    step_headroom *= 0.925;
  }
  if(reject_ratio<reject_goal && step_reject<1) {
    // Reject ratio too small and getting smaller
    step_headroom *= 1.075;
  }

  if(step_headroom>max_step_headroom) step_headroom=max_step_headroom;
  if(step_headroom<min_step_headroom) step_headroom=min_step_headroom;
}

////////////////////////////////////////////////////////////////////////
/// Klm_LLB_RKEvolve::ComputeEnergyChange  /////////////////////////
///    Threaded and non-threaded versions      /////////////////////////
////////////////////////////////////////////////////////////////////////
#if !OOMMF_THREADS

void Klm_LLB_RKEvolve::ComputeEnergyChange
(const Oxs_Mesh* mesh,
 const Oxs_MeshValue<OC_REAL8m>& current_energy,
 const Oxs_MeshValue<OC_REAL8m>& candidate_energy,
 OC_REAL8m& dE,OC_REAL8m& var_dE,OC_REAL8m& total_E)
{ // Computes cellwise difference between energies, and variance.
  // Export total_E is "current" energy (used for stepsize control).
  Nb_Xpfloat dE_xp      = 0.0;
  Nb_Xpfloat var_dE_xp  = 0.0;
  Nb_Xpfloat total_E_xp = 0.0;
  const OC_INDEX size = mesh->Size();
  for(OC_INDEX i=0;i<size;++i) {
    OC_REAL8m vol = mesh->Volume(i);
    OC_REAL8m e = vol*current_energy[i];
    OC_REAL8m new_e = vol*candidate_energy[i];
    total_E_xp += e;
    dE_xp += new_e - e;
    var_dE_xp += new_e*new_e + e*e;
  }
  total_E = total_E_xp.GetValue();
  dE      = dE_xp.GetValue();
  var_dE  = var_dE_xp.GetValue();
}

#else

class _Klm_LLBRKEvolve_RKFBase54_ComputeEnergyChange : public Oxs_ThreadRunObj {
public:
  static Oxs_Mutex job_control;
  static OC_INT4m offset;

  const Oxs_Mesh* mesh;                               // Import
  const Oxs_MeshValue<OC_REAL8m>* current_energy_ptr;    // Import
  const Oxs_MeshValue<OC_REAL8m>* candidate_energy_ptr;  // Import
  OC_REAL8m dE;                                          // Export
  OC_REAL8m var_dE;                                      // Export
  OC_REAL8m total_E;                                     // Export

  // Job control (imports)
  OC_INT4m vecsize;
  OC_INT4m block_size;

  _Klm_LLBRKEvolve_RKFBase54_ComputeEnergyChange()
    : mesh(0),current_energy_ptr(0),candidate_energy_ptr(0),
      dE(0.0), var_dE(0.0), total_E(0.0),
      vecsize(0), block_size(0) {}

  void Cmd(int threadnumber, void* data);
};

Oxs_Mutex _Klm_LLBRKEvolve_RKFBase54_ComputeEnergyChange::job_control;
OC_INT4m     _Klm_LLBRKEvolve_RKFBase54_ComputeEnergyChange::offset(0);

void _Klm_LLBRKEvolve_RKFBase54_ComputeEnergyChange::Cmd(int /* threadnumber */,
                                                      void* /* data */)
{ // Imports
  const Oxs_MeshValue<OC_REAL8m>& current_energy   = *current_energy_ptr;
  const Oxs_MeshValue<OC_REAL8m>& candidate_energy = *candidate_energy_ptr;

  // Exports
  dE = var_dE = total_E = 0.0;

  Nb_Xpfloat dE0=0.0, var_dE0=0.0, total_E0=0.0;
  Nb_Xpfloat dE1=0.0, var_dE1=0.0, total_E1=0.0;
  Nb_Xpfloat dE2=0.0, var_dE2=0.0, total_E2=0.0;
  Nb_Xpfloat dE3=0.0, var_dE3=0.0, total_E3=0.0;

  while(1) {
    job_control.Lock();
    OC_INT4m istart = offset;
    OC_INT4m istop = ( offset += block_size );
    job_control.Unlock();

    if(istart>=vecsize) break;
    if(istop>vecsize) istop=vecsize;

    OC_INT4m j;
    const OC_INT4m jstart = istart;
    const OC_INT4m jstop  = istop;

    for(j=jstart; j < jstart + (jstop-jstart)%4 ; ++j) {
      OC_REAL8m vol = mesh->Volume(j);
      OC_REAL8m e = vol*current_energy[j];
      OC_REAL8m new_e = vol*candidate_energy[j];
      total_E0 += e;
      dE0 += new_e - e;
      var_dE0 += new_e*new_e + e*e;
    }

    for(;j<jstop;j+=4) {
      OC_REAL8m vol0 = mesh->Volume(j);
      OC_REAL8m vol1 = mesh->Volume(j+1);

      OC_REAL8m e0 = vol0*current_energy[j];
      OC_REAL8m e1 = vol1*current_energy[j+1];

      OC_REAL8m new_e0 = vol0*candidate_energy[j];
      OC_REAL8m new_e1 = vol1*candidate_energy[j+1];

      total_E0 += e0;
      dE0 += new_e0 - e0;
      var_dE0 += new_e0*new_e0 + e0*e0;

      total_E1 += e1;
      dE1 += new_e1 - e1;
      var_dE1 += new_e1*new_e1 + e1*e1;

      OC_REAL8m vol2 = mesh->Volume(j+2);
      OC_REAL8m vol3 = mesh->Volume(j+3);

      OC_REAL8m e2 = vol2*current_energy[j+2];
      OC_REAL8m e3 = vol3*current_energy[j+3];
      OC_REAL8m new_e2 = vol2*candidate_energy[j+2];
      OC_REAL8m new_e3 = vol3*candidate_energy[j+3];

      total_E2 += e2;
      dE2 += new_e2 - e2;
      var_dE2 += new_e2*new_e2 + e2*e2;

      total_E3 += e3;
      dE3 += new_e3 - e3;
      var_dE3 += new_e3*new_e3 + e3*e3;
    }
  }

  dE = dE0 + dE1 + dE2 + dE3;
  var_dE = var_dE0 + var_dE1 + var_dE2 + var_dE3;
  total_E = total_E0 + total_E1 + total_E2 + total_E3;
}

void
Klm_LLB_RKEvolve::ComputeEnergyChange // Threaded version
(const Oxs_Mesh* mesh,
 const Oxs_MeshValue<OC_REAL8m>& current_energy,
 const Oxs_MeshValue<OC_REAL8m>& candidate_energy,
 OC_REAL8m& dE,OC_REAL8m& var_dE,OC_REAL8m& total_E)
{ // Computes cellwise difference between energies, and variance.
  // Export total_E is "current" energy (used for stepsize control).

  const OC_INDEX vecsize = mesh->Size();

  vector<_Klm_LLBRKEvolve_RKFBase54_ComputeEnergyChange> thread_data;
  const OC_INT4m MaxThreadCount = Oc_GetMaxThreadCount();
  thread_data.resize(MaxThreadCount);

  _Klm_LLBRKEvolve_RKFBase54_ComputeEnergyChange::job_control.Lock();
  _Klm_LLBRKEvolve_RKFBase54_ComputeEnergyChange::offset = 0;
  _Klm_LLBRKEvolve_RKFBase54_ComputeEnergyChange::job_control.Unlock();

  OC_INT4m ibs = (vecsize + MaxThreadCount - 1)/MaxThreadCount;
  ibs = (ibs + 7)/8;
  if(ibs%8) { ibs += 8 - (ibs%8); }
  if(ibs> static_cast<OC_INT4m>((vecsize*3+3)/4)) ibs = vecsize;

  OC_INT4m ithread;
  Oxs_ThreadTree threadtree;
  for(ithread=0;ithread<MaxThreadCount;++ithread) {
    thread_data[ithread].mesh = mesh;
    thread_data[ithread].current_energy_ptr = &current_energy;
    thread_data[ithread].candidate_energy_ptr = &candidate_energy;
    thread_data[ithread].vecsize= vecsize;
    thread_data[ithread].block_size= ibs;
    if(ithread!=0) threadtree.Launch(thread_data[ithread],0);
  }
  threadtree.LaunchRoot(thread_data[0],0);

  // Accumulate results
  dE = thread_data[0].dE;
  var_dE = thread_data[0].var_dE;
  total_E = thread_data[0].total_E;
  for(ithread=1;ithread<MaxThreadCount;++ithread) {
    dE += thread_data[ithread].dE;
    var_dE += thread_data[ithread].var_dE;
    total_E += thread_data[ithread].total_E;
  }
}

#endif


OC_BOOL
Klm_LLB_RKEvolve::InitNewStage
(const Klm_TimeDriver* /* driver */,
 Oxs_ConstKey<Oxs_SimState> state,
 Oxs_ConstKey<Oxs_SimState> prevstate)
{
  // Update derived data in state.
  const Oxs_SimState& cstate = state.GetReadReference();
  const Oxs_SimState* pstate_ptr = prevstate.GetPtr();
  UpdateDerivedOutputs(cstate,pstate_ptr);

  // Note 1: state is a copy-by-value import, so its read lock
  //         will be released on exit.
  // Note 2: pstate_ptr will be NULL if prevstate has
  //         "INVALID" status.

  return 1;
}

OC_BOOL
Klm_LLB_RKEvolve::Step(const Klm_TimeDriver* driver,
                      Oxs_ConstKey<Oxs_SimState> current_state_key,
                      const Oxs_DriverStepInfo& step_info,
                      Oxs_Key<Oxs_SimState>& next_state_key)
{
#if REPORT_TIME
  steponlytime.Start();
#endif // REPORT_TIME

  const OC_REAL8m bad_energy_cut_ratio = 0.75;
  const OC_REAL8m bad_energy_step_increase = 1.3;

  const OC_REAL8m previous_next_timestep = next_timestep;

  const Oxs_SimState& cstate = current_state_key.GetReadReference();

  CheckCache(cstate);

  // Note if start_dM or start_dt is being used
  OC_BOOL start_cond_active=0;
  if(next_timestep<=0.0 ||
     (cstate.stage_iteration_count<1
      && step_info.current_attempt_count==0)) {
    if(cstate.stage_number==0
       || stage_init_step_control == SISC_START_COND) {
      start_cond_active = 1;
    } else if(stage_init_step_control == SISC_CONTINUOUS) {
      start_cond_active = 0;  // Safety
    } else if(stage_init_step_control == SISC_AUTO) {
      // Automatic detection based on energy values across
      // stage boundary.
      OC_REAL8m total_E,E_diff;
      if(cstate.GetDerivedData("Total E",total_E) &&
         cstate.GetDerivedData("Delta E",E_diff)  &&
         fabs(E_diff) <= 256*OC_REAL8_EPSILON*fabs(total_E) ) {
        // The factor of 256 in the preceding line is a fudge factor,
        // selected with no particular justification.
        start_cond_active = 0;  // Continuous case
      } else {
        start_cond_active = 1;  // Assume discontinuous
      }
    } else {
      throw Oxs_ExtError(this,
           "Klm_LLB_RKEvolve::Step; Programming error:"
           " unrecognized stage_init_step_control value");
    }
  }

  // Negotiate timestep, and also initialize both next_state and
  // temp_state structures.
  Oxs_SimState* work_state = &(next_state_key.GetWriteReference());
  OC_BOOL force_step=0,driver_set_step=0;
  NegotiateTimeStep(driver,cstate,*work_state,next_timestep,
                    start_cond_active,force_step,driver_set_step);
  OC_REAL8m stepsize = work_state->last_timestep;
  work_state=NULL; // Safety: disable pointer

  // Step
  OC_REAL8m error_estimate,norm_error;
  OC_REAL8m global_error_order;
  OC_BOOL new_energy_and_dmdt_computed;
  OC_BOOL reject_step=0;
  (this->*rkstep_ptr)(stepsize,current_state_key,
                      dM_dt_output.cache.value,
                      next_state_key,
                      error_estimate,global_error_order,norm_error,
                      new_energy_and_dmdt_computed);
  const Oxs_SimState& nstate = next_state_key.GetReadReference();
  driver->FillStateDerivedData(cstate,nstate);

  OC_REAL8m max_dM_dt;
  cstate.GetDerivedData("Max dM/dt",max_dM_dt);
  OC_REAL8m reference_stepsize = stepsize;
  if(driver_set_step) reference_stepsize = previous_next_timestep;
  OC_BOOL good_step = CheckError(global_error_order,error_estimate,
                              stepsize,reference_stepsize,
                              max_dM_dt,next_timestep);
  /// Note: Might want to use average or larger of max_dM_dt
  /// and new_max_dM_dt (computed below.)

  if(!good_step && !force_step) {
    // Bad step; The only justfication to do energy and dM_dt
    // computation would be to get an energy-based stepsize
    // adjustment estimate, which we only need to try if 
    // next_timestep is larger than cut applied by energy
    // rejection code (i.e., bad_energy_cut_ratio).
    if(next_timestep<=stepsize*bad_energy_cut_ratio) {
      AdjustStepHeadroom(1);
#if REPORT_TIME
      steponlytime.Stop();
#endif // REPORT_TIME
      return 0; // Don't bother with energy calculation
    }
    reject_step=1; // Otherwise, mark step rejected and see what energy
    /// info suggests for next stepsize
  }

  if(start_cond_active && !force_step) {
    if(start_dM>=0.0) {
      // Check that no spin has moved by more than start_dM
      OC_REAL8m diff = MaxDiff_mul(cstate.spin,*(cstate.Ms.GetPtr()),
                                nstate.spin,*(nstate.Ms.GetPtr()));
      if(diff>start_dM) {
        next_timestep = step_headroom * stepsize * (start_dM/diff);
        if(next_timestep<=stepsize*bad_energy_cut_ratio) {
          AdjustStepHeadroom(1);
#if REPORT_TIME
          steponlytime.Stop();
#endif // REPORT_TIME
          return 0; // Don't bother with energy calculation
        }
        reject_step=1; // Otherwise, mark step rejected and see what energy
        /// info suggests for next stepsize
      }
    }
  }

#ifdef OLDE_CODE
  if(norm_error>0.0005) {
    fprintf(stderr,
            "Iteration %u passed error check; norm_error=%8.5f\n",
            nstate.iteration_count,norm_error);
  } /**/
#endif // OLDE_CODE

  // Energy timestep control:
  //   The relationship between energy error and stepsize appears to be
  // highly nonlinear, so that estimating appropriate stepsize from energy
  // increase is difficult.  Perhaps it is possible to include energy
  // interpolation into RK step routines, but for the present we just
  // reduce the step by a fixed ratio if we detect energy increase beyond
  // that which can be attributed to numerical errors.  Of course, this
  // doesn't take into account the expected energy decrease (which depends
  // on the damping ratio alpha), which is another reason to try to build
  // it into the high order RK step routines.
  OC_REAL8m pE_pt,new_pE_pt=0.;
  cstate.GetDerivedData("pE/pt",pE_pt);
#if KL_DEBUG
  OC_BOOL kltmpXfilled = FALSE;
#endif // KL_DEBUG
  if(new_energy_and_dmdt_computed) {
    nstate.GetDerivedData("pE/pt",new_pE_pt);
    // KL(m) Here mxH_output is filled, I guess
  } else {
    OC_REAL8m new_total_E;
    GetEnergyDensity(nstate,temp_energy,
                     &mxH_output.cache.value,
                     &kltmpX,  // KL(m)
                     new_pE_pt,new_total_E);
#if KL_DEBUG
    kltmpXfilled = TRUE;
#endif // KL_DEBUG
    mxH_output.cache.state_id=nstate.Id();
    if(!nstate.AddDerivedData("pE/pt",new_pE_pt)) {
      throw Oxs_ExtError(this,
           "Klm_LLB_RKEvolve::Step:"
           " Programming error; data cache (pE/pt) already set.");
    }
    if(!nstate.AddDerivedData("Total E",new_total_E)) {
      throw Oxs_ExtError(this,
           "Klm_LLB_RKEvolve::Step:"
           " Programming error; data cache (Total E) already set.");
    }
  }

#if REPORT_TIME_RKDEVEL
timer[0].Start(); /**/
#endif // REPORT_TIME_RKDEVEL
  OC_REAL8m dE,var_dE,total_E;
  ComputeEnergyChange(nstate.mesh,energy,temp_energy,dE,var_dE,total_E);
  if(!nstate.AddDerivedData("Delta E",dE)) {
    throw Oxs_ExtError(this,
         "Klm_LLB_RKEvolve::Step:"
         " Programming error; data cache (Delta E) already set.");
  }
#if REPORT_TIME_RKDEVEL
timer[0].Stop(); /**/
++timer_counts[0].pass_count;
 timer_counts[0].bytes += (nstate.mesh->Size())*(2*sizeof(OC_REAL8m));
 timer_counts[0].name = "ComputeEnergyChange";
#endif // REPORT_TIME_RKDEVEL

  if(expected_energy_precision>=0.) {
    var_dE *= expected_energy_precision * expected_energy_precision;
    /// Variance, assuming error in each energy[i] term is independent,
    /// uniformly distributed, 0-mean, with range
    ///        +/- expected_energy_precision*energy[i].
    /// It would probably be better to get an error estimate directly
    /// from each energy term.
    OC_REAL8m E_numerror = OC_MAX(fabs(total_E)*expected_energy_precision,
                               2*sqrt(var_dE));
    OC_REAL8m pE_pt_max = OC_MAX(pE_pt,new_pE_pt); // Might want to
    /// change this from a constant to a linear function in timestep.

    OC_REAL8m reject_dE = 2*E_numerror + pE_pt_max * stepsize;
    if(dE>reject_dE) {
      OC_REAL8m teststep = bad_energy_cut_ratio*stepsize;
      if(teststep<next_timestep) {
        next_timestep=teststep;
        max_step_increase = bad_energy_step_increase;
        // Damp the enthusiasm of the RK stepper routine for
        // stepsize growth, for a step or two.
      }
      if(!force_step) reject_step=1;
    }
  }

  // KL(m)
  // Here was |M| length check code.
  // Moved it to AdjustState (max_mn, min_m are now in exchange module).
  // Later again here, but a little bit below
  
  if(!force_step && reject_step) {
    AdjustStepHeadroom(1);
#if REPORT_TIME
    steponlytime.Stop();
#endif // REPORT_TIME
    return 0;
  }

  // Otherwise, we are accepting the new step.

  // KL(m)
  // |M| length check code.  
  const OC_INDEX size = nstate.mesh->Size();
  const Oxs_MeshValue<OC_REAL8m>& Ms = (*(nstate.Ms.GetPtr()));
  const Oxs_MeshValue<OC_REAL8m>& Ms_T0 = (*(nstate.GetPtr_Ms_T0()));
  OC_INDEX i;
  for(i=0;i<size;++i) {
    if(Ms[i]!= 0.0 &&    // avoid "empty" cells
       Ms_T0[i]!= 0.0) { // security: division by zero
      OC_REAL8m m_norm = Ms[i] / Ms_T0[i]; // normalized with Ms_T0
      // max_m handling
      if(expected_m_max_precision>= 0.0) {
        if(m_norm-1.0 > expected_m_max_precision) {
          char buf[1024];
          Oc_Snprintf(buf,sizeof(buf), "|M| exceeds Ms_T0.\n"
            "This is probably due to large parameter Klm_LLB_Term::chi_parallel -"
            " it should small enough.\n"
            "Max(|M|/Ms_T0)-1= %g, while m_max_precision is %g.\n",
            m_norm-1.0, expected_m_max_precision);
          static Oxs_WarningMessage foo(3); // how often to repeat it
          foo.Send(revision_info,OC_STRINGIFY(__LINE__), buf);
          // How react for such situations? Set it to Ms_T0?
        }
      }
      // min_m handling
      if(expected_m_min>= 0.0) {
        if(m_norm < expected_m_min) {
          char buf[1024];
          Oc_Snprintf(buf,sizeof(buf), "|M| is veeeery small.\n"
            "This could be dengerous.\n"
            "Min(|M|/Ms_T0)= %g, while parameter m_min is %g.\n",
            m_norm, expected_m_min);
          static Oxs_WarningMessage foo(3); // how often to repeat it
          foo.Send(revision_info,OC_STRINGIFY(__LINE__), buf);
          // How react for such situations???
        }
      }
    }  
  }
  
  // Calculate dM_dt at new point, and fill in cache.
  if(new_energy_and_dmdt_computed) {
    dM_dt_output.cache.value.Swap(vtmpA);
  } else {
    OC_REAL8m new_max_dM_dt,new_dE_dt,new_timestep_lower_bound;
    if(mxH_output.cache.state_id != nstate.Id()) { // Safety
    throw Oxs_ExtError(this,
         "Klm_LLB_RKEvolve::Step:"
         " Programming error; mxH_output cache improperly filled.");
    }
#if KL_DEBUG
    if(!kltmpXfilled) { // Safety
    throw Oxs_ExtError(this,
         "Klm_LLB_RKEvolve::Step:"
         " Programming error; kltmpX kl-object improperly filled.");
    }
#endif // KL_DEBUG
#if !OOMMF_THREADS
    Calculate_dM_dt(nstate,
                    mxH_output.cache.value,
                    kltmpX, // KL(m)
                    new_pE_pt,
                    dM_dt_output.cache.value,new_max_dM_dt,
                    new_dE_dt,new_timestep_lower_bound);
#else // OOMMF_THREADS
    {
      vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadA> thread_data;
!!! m->M
      Initialize_Threaded_Calculate_dm_dt(nstate,mxH_output.cache.value,
                                          dm_dt_output.cache.value,
                                          thread_data);
      Oxs_ThreadTree threadtree;
      const OC_INT4m thread_count = thread_data.size();
      for(OC_INT4m ithread=1;ithread<thread_count;++ithread) {
        threadtree.Launch(thread_data[ithread],0);
      }
      threadtree.LaunchRoot(thread_data[0],0);
!!! m->M
      Finalize_Threaded_Calculate_dm_dt(thread_data,new_pE_pt,
                                        new_max_dm_dt,new_dE_dt,
                                        new_timestep_lower_bound);
    }
#endif // OOMMF_THREADS
    if(!nstate.AddDerivedData("Timestep lower bound",
                              new_timestep_lower_bound) ||
       !nstate.AddDerivedData("Max dM/dt",new_max_dM_dt) ||
       !nstate.AddDerivedData("dE/dt",new_dE_dt)) {
      throw Oxs_ExtError(this,
                           "Klm_LLB_RKEvolve::Step:"
                           " Programming error; data cache already set.");
    }
  }
  dM_dt_output.cache.state_id = nstate.Id();

  energy.Swap(temp_energy);
  energy_state_id = nstate.Id();

  AdjustStepHeadroom(0);
  if(!force_step && max_step_increase<max_step_increase_limit) {
    max_step_increase *= max_step_increase_adj_ratio;
  }
  if(max_step_increase>max_step_increase_limit) {
    max_step_increase = max_step_increase_limit;
  }

#if REPORT_TIME
  steponlytime.Stop();
#endif // REPORT_TIME
  return 1; // Accept step
}

void Klm_LLB_RKEvolve::UpdateDerivedOutputs(const Oxs_SimState& state,
                                                const Oxs_SimState* prevstate_ptr)
{ // This routine fills all the Klm_LLB_RKEvolve Oxs_ScalarOutput's to
  // the appropriate value based on the import "state", and any of
  // Oxs_VectorOutput's that have CacheRequest enabled are filled.
  // It also makes sure all the expected WOO objects in state are
  // filled.
  max_dM_dt_output.cache.state_id
    = dE_dt_output.cache.state_id
    = delta_E_output.cache.state_id
//    = min_m_output.cache.state_id
//    = max_m_output.cache.state_id
    = 0;  // Mark change in progress

  OC_REAL8m dummy_value;
  if(!state.GetDerivedData("Max dM/dt",max_dM_dt_output.cache.value) ||
     !state.GetDerivedData("dE/dt",dE_dt_output.cache.value) ||
     !state.GetDerivedData("Delta E",delta_E_output.cache.value) ||
//     !state.GetDerivedData("Min m",min_m_output.cache.value) ||
//     !state.GetDerivedData("Max m",max_m_output.cache.value) ||
     !state.GetDerivedData("pE/pt",dummy_value) ||
     !state.GetDerivedData("Total E",dummy_value) ||
     !state.GetDerivedData("Timestep lower bound",dummy_value) ||
     (dM_dt_output.GetCacheRequestCount()>0
      && dM_dt_output.cache.state_id != state.Id()) ||
     (mxH_output.GetCacheRequestCount()>0
      && mxH_output.cache.state_id != state.Id())) {

    // Missing at least some data, so calculate from scratch

    // Calculate H and mxH outputs
    Oxs_MeshValue<ThreeVector>& mxH = mxH_output.cache.value;
    OC_REAL8m pE_pt,total_E;
    GetEnergyDensity(state,energy,&mxH,
                     &kltmpX, // KL(m)
                     pE_pt,total_E);
    energy_state_id=state.Id();
    mxH_output.cache.state_id=state.Id();
    if(!state.GetDerivedData("pE/pt",dummy_value)) {
      state.AddDerivedData("pE/pt",pE_pt);
    }
    if(!state.GetDerivedData("Total E",dummy_value)) {
      state.AddDerivedData("Total E",total_E);
    }

    // Calculate dM/dt, Max dM/dt and dE/dt
    Oxs_MeshValue<ThreeVector>& dM_dt
      = dM_dt_output.cache.value;
    dM_dt_output.cache.state_id=0;
    OC_REAL8m timestep_lower_bound;
    Calculate_dM_dt(state,mxH,kltmpX,pE_pt,dM_dt,
                    max_dM_dt_output.cache.value,
                    dE_dt_output.cache.value,timestep_lower_bound);
    dM_dt_output.cache.state_id=state.Id();

    if(!state.GetDerivedData("Max dM/dt",dummy_value)) {
      state.AddDerivedData("Max dM/dt",max_dM_dt_output.cache.value);
    }

    if(!state.GetDerivedData("dE/dt",dummy_value)) {
      state.AddDerivedData("dE/dt",dE_dt_output.cache.value);
    }

    if(!state.GetDerivedData("Timestep lower bound",dummy_value)) {
      state.AddDerivedData("Timestep lower bound",
                           timestep_lower_bound);
    }

/*  // KL(m)
    const Oxs_MeshValue<OC_REAL8m>& Ms = *(state.Ms);
//    const Oxs_MeshValue<OC_REAL8m>& Ms_T0 = *(state.GetPtr_Ms_T0()); // for tests
    const Oxs_MeshValue<OC_REAL8m>& Me_T = *(state.GetPtr_Me_T()); // for tests
    const OC_INDEX size = state.mesh->Size();
    OC_INDEX i;
    OC_REAL8m min_m_  = FLT_MAX;
    OC_REAL8m max_m_  = FLT_MIN;
    for(i=0;i<size;i++)
      if(Ms[i]!=0.0) {
        if(min_m_>Ms[i]/Me_T[i]) min_m_=Ms[i]/Me_T[i];
        if(max_m_<Ms[i]/Me_T[i]) max_m_=Ms[i]/Me_T[i];
      }
    if(!state.GetDerivedData("Min m",dummy_value)) {
      state.AddDerivedData("Min m",
                           min_m_);
    }
    if(!state.GetDerivedData("Max m",dummy_value)) {
      state.AddDerivedData("Max m",
                           max_m_);
    }
*/
    //

    if(!state.GetDerivedData("Delta E",dummy_value)) {
      if(state.previous_state_id == 0) {
        // No previous state
        dummy_value = 0.0;
      } else if(prevstate_ptr!=NULL
                && state.previous_state_id == prevstate_ptr->Id()) {
        OC_REAL8m old_E;
        if(!prevstate_ptr->GetDerivedData("Total E",old_E)) {
          throw Oxs_ExtError(this,
                             "Klm_LLB_RKEvolve::UpdateDerivedOutputs:"
                             " \"Total E\" not set in previous state.");
        }
        dummy_value = total_E - old_E; // This is less accurate than adding
        /// up the cell-by-cell differences (as is used in the main code),
        /// but w/o the entire energy map for prevstate this is the best
        /// we can do.
      } else {
        throw Oxs_ExtError(this,
           "Klm_LLB_RKEvolve::UpdateDerivedOutputs:"
           " Can't derive Delta E from single state.");
      }
      state.AddDerivedData("Delta E",dummy_value);
    }
    delta_E_output.cache.value=dummy_value;

  }

  // KL(m)
  // Changed units, now: A/s*m
  // max_dM_dt_output.cache.value*=(180e-9/PI);
  /// Convert from radians/second to deg/ns

  max_dM_dt_output.cache.state_id
    = dE_dt_output.cache.state_id
    = delta_E_output.cache.state_id
//    = min_m_output.cache.state_id
//    = max_m_output.cache.state_id
    = state.Id();
}
