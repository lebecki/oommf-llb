/* FILE: kl_llb_util.h            -*-Mode: c++-*-
 *
 * Utility tools for Landau-Lifshitz-Bloch modules
 * 
 * Kristof M. Lebecki,
 * Universität konstanz, 2010
 */
 
#include "director.h"

/* End includes */

// Value used in comparisons, like
// instead:
//   if(small_num==0)
// write:
//   if(fabs(small_num)<SMALL_VAL)
#define SMALL_VAL 1e-10

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

// Returne TRUE if the driver is Klm_TimeDriver
// This troutin has to be rewritten if there are more than one
// driver in OOMMF allowed
OC_BOOL VarMs_driver(const Oxs_Director* director);

// Set by Klm_LLB_Term.
// Checked/cleared by Klm_LLB_RKEvolve.
void Set_Use_LLB_Term();
OC_BOOL Get_Use_LLB_Term();
void Clear_Use_LLB_Term();

