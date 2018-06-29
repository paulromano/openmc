#ifndef NUCLIDE_H
#define NUCLIDE_H

#include <string>

#include "error.h"

namespace openmc {

inline std::string nuclide_name(int i) {
  // Get nuclide name as a char*
  char* name;
  int err = openmc_nuclide_name(i, &name);
  if (err) {
    openmc::fatal_error(openmc_err_msg);
  }

  // Since char* is not null-terminated, find first blank and return string with
  // that many characters
  size_t j;
  for (j = 0; j < 20; ++j) {
    if (name[j] == ' ') break;
  }
  return {name, j};
}

} // namespace openmc

#endif // NUCLIDE_H
