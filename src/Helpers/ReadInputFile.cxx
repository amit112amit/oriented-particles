#include "HelperFunctions.h"

namespace OPS {
// Borrowed from StackOverflow
InputParameters readKeyValueInput(std::string fileName) {
  InputParameters params;
  std::ifstream input(fileName.c_str());
  assert(input);
  std::string line;
  while (std::getline(input, line)) {
    std::istringstream is_line(line);
    std::string key, value;
    is_line >> key >> value;
    params.insert(std::make_pair(key, value));
  }
  input.close();
  return params;
}
} // namespace OPS
