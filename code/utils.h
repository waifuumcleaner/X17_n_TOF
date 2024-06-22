#ifndef UTILS_H
#define UTILS_H

#include <exception>
#include <iostream>
#include <stdexcept>

/**
 * @brief
 * Prints an error line, terminates the program for specific exceptions
 *
 * @param e standard C++ exception
 */
void handle_exception(const std::exception &e);

#endif // UTILS_H
