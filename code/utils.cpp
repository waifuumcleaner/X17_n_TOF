#include "utils.h"

void handle_exception(const std::exception &e) {
    if (dynamic_cast<const std::out_of_range *>(&e)) {
        std::cerr << "Out of range error: " << e.what() << std::endl;
        std::terminate();
    } else if (dynamic_cast<const std::invalid_argument *>(&e)) {
        std::cerr << "Invalid argument error: " << e.what() << std::endl;
        std::terminate();
    } else {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

