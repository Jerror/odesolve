/** @brief Implement rkab_results API (ie., provide methods for deletion)
 * @details I'm providing a C interface; see adaptive_step_rk.h for details.
 * @author Jeremiah O'Neil
 * @copyright GNU Public License. */
#include "rkab.hpp" // templates
#include "adaptive_step_rk.h" // interface

// Instantiate as defined in adaptive_step_rk.h:
MAP_TARGETS_TO(INST_DELETE_RKAB_RESULTS) 
