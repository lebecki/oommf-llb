/* FILE: kl_llb_util.cc            -*-Mode: c++-*-
 *
 * Utility tools for Landau-Lifshitz-Bloch modules
 * 
 * Kristof M. Lebecki,
 * Universität konstanz, 2010
 */
 
#include "kl_llb_util.h"

/* End includes */

OC_BOOL VarMs_driver(const Oxs_Director* director) {
  String drv_name;
  director->GetDriverInstanceName(drv_name);
  drv_name = drv_name.substr(0,20);
  if(drv_name.compare("Klm_TimeDriver:")==0)
    return TRUE;
  else
    return FALSE;
}

OC_BOOL LLB_Term_used = FALSE; // Initialization
void Set_Use_LLB_Term() {
  LLB_Term_used = TRUE;
}
OC_BOOL Get_Use_LLB_Term() {
  return LLB_Term_used;
}
void Clear_Use_LLB_Term() {
  LLB_Term_used = FALSE;
}


