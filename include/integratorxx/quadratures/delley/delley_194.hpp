#pragma once

namespace IntegratorXX {
namespace DelleyGrids {

/**
 *  \brief Delley Quadrature specification for order 23 grid with 194 points
 *
 */
template <typename T>
struct delley_194 {
  static constexpr std::array<cartesian_pt_t<T>, 194> points = {
      0.10000000000000000E+01,  0.00000000000000000E+00,
      0.00000000000000000E+00,  -0.10000000000000000E+01,
      0.00000000000000000E+00,  0.00000000000000000E+00,
      0.00000000000000000E+00,  0.10000000000000000E+01,
      0.00000000000000000E+00,  0.00000000000000000E+00,
      -0.10000000000000000E+01, 0.00000000000000000E+00,
      0.00000000000000000E+00,  0.00000000000000000E+00,
      0.10000000000000000E+01,  0.00000000000000000E+00,
      0.00000000000000000E+00,  -0.10000000000000000E+01,
      0.57735026918962584E+00,  0.57735026918962584E+00,
      0.57735026918962584E+00,  0.57735026918962584E+00,
      0.57735026918962584E+00,  -0.57735026918962584E+00,
      0.57735026918962584E+00,  -0.57735026918962584E+00,
      0.57735026918962584E+00,  0.57735026918962584E+00,
      -0.57735026918962584E+00, -0.57735026918962584E+00,
      -0.57735026918962584E+00, 0.57735026918962584E+00,
      0.57735026918962584E+00,  -0.57735026918962584E+00,
      0.57735026918962584E+00,  -0.57735026918962584E+00,
      -0.57735026918962584E+00, -0.57735026918962584E+00,
      0.57735026918962584E+00,  -0.57735026918962584E+00,
      -0.57735026918962584E+00, -0.57735026918962584E+00,
      0.70710678118654746E+00,  0.70710678118654746E+00,
      0.00000000000000000E+00,  0.00000000000000000E+00,
      0.70710678118654746E+00,  0.70710678118654746E+00,
      0.70710678118654746E+00,  0.00000000000000000E+00,
      0.70710678118654746E+00,  0.70710678118654746E+00,
      -0.70710678118654746E+00, 0.00000000000000000E+00,
      0.00000000000000000E+00,  0.70710678118654746E+00,
      -0.70710678118654746E+00, -0.70710678118654746E+00,
      0.00000000000000000E+00,  0.70710678118654746E+00,
      -0.70710678118654746E+00, 0.70710678118654746E+00,
      0.00000000000000000E+00,  0.00000000000000000E+00,
      -0.70710678118654746E+00, 0.70710678118654746E+00,
      0.70710678118654746E+00,  0.00000000000000000E+00,
      -0.70710678118654746E+00, -0.70710678118654746E+00,
      -0.70710678118654746E+00, 0.00000000000000000E+00,
      0.00000000000000000E+00,  -0.70710678118654746E+00,
      -0.70710678118654746E+00, -0.70710678118654746E+00,
      0.00000000000000000E+00,  -0.70710678118654746E+00,
      0.77749321931476711E+00,  0.44469331787174371E+00,
      0.44469331787174371E+00,  0.44469331787174371E+00,
      0.77749321931476711E+00,  0.44469331787174371E+00,
      0.44469331787174371E+00,  0.44469331787174371E+00,
      0.77749321931476711E+00,  0.77749321931476711E+00,
      0.44469331787174371E+00,  -0.44469331787174371E+00,
      0.44469331787174371E+00,  0.77749321931476711E+00,
      -0.44469331787174371E+00, 0.44469331787174371E+00,
      0.44469331787174371E+00,  -0.77749321931476711E+00,
      0.77749321931476711E+00,  -0.44469331787174371E+00,
      0.44469331787174371E+00,  0.44469331787174371E+00,
      -0.77749321931476711E+00, 0.44469331787174371E+00,
      0.44469331787174371E+00,  -0.44469331787174371E+00,
      0.77749321931476711E+00,  0.77749321931476711E+00,
      -0.44469331787174371E+00, -0.44469331787174371E+00,
      0.44469331787174371E+00,  -0.77749321931476711E+00,
      -0.44469331787174371E+00, 0.44469331787174371E+00,
      -0.44469331787174371E+00, -0.77749321931476711E+00,
      -0.77749321931476711E+00, 0.44469331787174371E+00,
      0.44469331787174371E+00,  -0.44469331787174371E+00,
      0.77749321931476711E+00,  0.44469331787174371E+00,
      -0.44469331787174371E+00, 0.44469331787174371E+00,
      0.77749321931476711E+00,  -0.77749321931476711E+00,
      0.44469331787174371E+00,  -0.44469331787174371E+00,
      -0.44469331787174371E+00, 0.77749321931476711E+00,
      -0.44469331787174371E+00, -0.44469331787174371E+00,
      0.44469331787174371E+00,  -0.77749321931476711E+00,
      -0.77749321931476711E+00, -0.44469331787174371E+00,
      0.44469331787174371E+00,  -0.44469331787174371E+00,
      -0.77749321931476711E+00, 0.44469331787174371E+00,
      -0.44469331787174371E+00, -0.44469331787174371E+00,
      0.77749321931476711E+00,  -0.77749321931476711E+00,
      -0.44469331787174371E+00, -0.44469331787174371E+00,
      -0.44469331787174371E+00, -0.77749321931476711E+00,
      -0.44469331787174371E+00, -0.44469331787174371E+00,
      -0.44469331787174371E+00, -0.77749321931476711E+00,
      0.91250909686747372E+00,  0.28924656275754385E+00,
      0.28924656275754385E+00,  0.28924656275754385E+00,
      0.91250909686747372E+00,  0.28924656275754385E+00,
      0.28924656275754385E+00,  0.28924656275754385E+00,
      0.91250909686747372E+00,  0.91250909686747372E+00,
      0.28924656275754385E+00,  -0.28924656275754385E+00,
      0.28924656275754385E+00,  0.91250909686747372E+00,
      -0.28924656275754385E+00, 0.28924656275754385E+00,
      0.28924656275754385E+00,  -0.91250909686747372E+00,
      0.91250909686747372E+00,  -0.28924656275754385E+00,
      0.28924656275754385E+00,  0.28924656275754385E+00,
      -0.91250909686747372E+00, 0.28924656275754385E+00,
      0.28924656275754385E+00,  -0.28924656275754385E+00,
      0.91250909686747372E+00,  0.91250909686747372E+00,
      -0.28924656275754385E+00, -0.28924656275754385E+00,
      0.28924656275754385E+00,  -0.91250909686747372E+00,
      -0.28924656275754385E+00, 0.28924656275754385E+00,
      -0.28924656275754385E+00, -0.91250909686747372E+00,
      -0.91250909686747372E+00, 0.28924656275754385E+00,
      0.28924656275754385E+00,  -0.28924656275754385E+00,
      0.91250909686747372E+00,  0.28924656275754385E+00,
      -0.28924656275754385E+00, 0.28924656275754385E+00,
      0.91250909686747372E+00,  -0.91250909686747372E+00,
      0.28924656275754385E+00,  -0.28924656275754385E+00,
      -0.28924656275754385E+00, 0.91250909686747372E+00,
      -0.28924656275754385E+00, -0.28924656275754385E+00,
      0.28924656275754385E+00,  -0.91250909686747372E+00,
      -0.91250909686747372E+00, -0.28924656275754385E+00,
      0.28924656275754385E+00,  -0.28924656275754385E+00,
      -0.91250909686747372E+00, 0.28924656275754385E+00,
      -0.28924656275754385E+00, -0.28924656275754385E+00,
      0.91250909686747372E+00,  -0.91250909686747372E+00,
      -0.28924656275754385E+00, -0.28924656275754385E+00,
      -0.28924656275754385E+00, -0.91250909686747372E+00,
      -0.28924656275754385E+00, -0.28924656275754385E+00,
      -0.28924656275754385E+00, -0.91250909686747372E+00,
      0.31419699418258629E+00,  0.67129734426952259E+00,
      0.67129734426952259E+00,  0.67129734426952259E+00,
      0.31419699418258629E+00,  0.67129734426952259E+00,
      0.67129734426952259E+00,  0.67129734426952259E+00,
      0.31419699418258629E+00,  0.31419699418258629E+00,
      0.67129734426952259E+00,  -0.67129734426952259E+00,
      0.67129734426952259E+00,  0.31419699418258629E+00,
      -0.67129734426952259E+00, 0.67129734426952259E+00,
      0.67129734426952259E+00,  -0.31419699418258629E+00,
      0.31419699418258629E+00,  -0.67129734426952259E+00,
      0.67129734426952259E+00,  0.67129734426952259E+00,
      -0.31419699418258629E+00, 0.67129734426952259E+00,
      0.67129734426952259E+00,  -0.67129734426952259E+00,
      0.31419699418258629E+00,  0.31419699418258629E+00,
      -0.67129734426952259E+00, -0.67129734426952259E+00,
      0.67129734426952259E+00,  -0.31419699418258629E+00,
      -0.67129734426952259E+00, 0.67129734426952259E+00,
      -0.67129734426952259E+00, -0.31419699418258629E+00,
      -0.31419699418258629E+00, 0.67129734426952259E+00,
      0.67129734426952259E+00,  -0.67129734426952259E+00,
      0.31419699418258629E+00,  0.67129734426952259E+00,
      -0.67129734426952259E+00, 0.67129734426952259E+00,
      0.31419699418258629E+00,  -0.31419699418258629E+00,
      0.67129734426952259E+00,  -0.67129734426952259E+00,
      -0.67129734426952259E+00, 0.31419699418258629E+00,
      -0.67129734426952259E+00, -0.67129734426952259E+00,
      0.67129734426952259E+00,  -0.31419699418258629E+00,
      -0.31419699418258629E+00, -0.67129734426952259E+00,
      0.67129734426952259E+00,  -0.67129734426952259E+00,
      -0.31419699418258629E+00, 0.67129734426952259E+00,
      -0.67129734426952259E+00, -0.67129734426952259E+00,
      0.31419699418258629E+00,  -0.31419699418258629E+00,
      -0.67129734426952259E+00, -0.67129734426952259E+00,
      -0.67129734426952259E+00, -0.31419699418258629E+00,
      -0.67129734426952259E+00, -0.67129734426952259E+00,
      -0.67129734426952259E+00, -0.31419699418258629E+00,
      0.98297230270725333E+00,  0.12993354476500668E+00,
      0.12993354476500668E+00,  0.12993354476500668E+00,
      0.98297230270725333E+00,  0.12993354476500668E+00,
      0.12993354476500668E+00,  0.12993354476500668E+00,
      0.98297230270725333E+00,  0.98297230270725333E+00,
      0.12993354476500668E+00,  -0.12993354476500668E+00,
      0.12993354476500668E+00,  0.98297230270725333E+00,
      -0.12993354476500668E+00, 0.12993354476500668E+00,
      0.12993354476500668E+00,  -0.98297230270725333E+00,
      0.98297230270725333E+00,  -0.12993354476500668E+00,
      0.12993354476500668E+00,  0.12993354476500668E+00,
      -0.98297230270725333E+00, 0.12993354476500668E+00,
      0.12993354476500668E+00,  -0.12993354476500668E+00,
      0.98297230270725333E+00,  0.98297230270725333E+00,
      -0.12993354476500668E+00, -0.12993354476500668E+00,
      0.12993354476500668E+00,  -0.98297230270725333E+00,
      -0.12993354476500668E+00, 0.12993354476500668E+00,
      -0.12993354476500668E+00, -0.98297230270725333E+00,
      -0.98297230270725333E+00, 0.12993354476500668E+00,
      0.12993354476500668E+00,  -0.12993354476500668E+00,
      0.98297230270725333E+00,  0.12993354476500668E+00,
      -0.12993354476500668E+00, 0.12993354476500668E+00,
      0.98297230270725333E+00,  -0.98297230270725333E+00,
      0.12993354476500668E+00,  -0.12993354476500668E+00,
      -0.12993354476500668E+00, 0.98297230270725333E+00,
      -0.12993354476500668E+00, -0.12993354476500668E+00,
      0.12993354476500668E+00,  -0.98297230270725333E+00,
      -0.98297230270725333E+00, -0.12993354476500668E+00,
      0.12993354476500668E+00,  -0.12993354476500668E+00,
      -0.98297230270725333E+00, 0.12993354476500668E+00,
      -0.12993354476500668E+00, -0.12993354476500668E+00,
      0.98297230270725333E+00,  -0.98297230270725333E+00,
      -0.12993354476500668E+00, -0.12993354476500668E+00,
      -0.12993354476500668E+00, -0.98297230270725333E+00,
      -0.12993354476500668E+00, -0.12993354476500668E+00,
      -0.12993354476500668E+00, -0.98297230270725333E+00,
      0.00000000000000000E+00,  0.34577021976112826E+00,
      0.93831921813759156E+00,  0.93831921813759156E+00,
      0.00000000000000000E+00,  0.34577021976112826E+00,
      0.34577021976112826E+00,  0.93831921813759156E+00,
      0.00000000000000000E+00,  0.00000000000000000E+00,
      0.93831921813759156E+00,  0.34577021976112826E+00,
      0.34577021976112826E+00,  0.00000000000000000E+00,
      0.93831921813759156E+00,  0.93831921813759156E+00,
      0.34577021976112826E+00,  0.00000000000000000E+00,
      0.00000000000000000E+00,  0.34577021976112826E+00,
      -0.93831921813759156E+00, -0.93831921813759156E+00,
      0.00000000000000000E+00,  0.34577021976112826E+00,
      0.34577021976112826E+00,  -0.93831921813759156E+00,
      0.00000000000000000E+00,  0.00000000000000000E+00,
      -0.93831921813759156E+00, 0.34577021976112826E+00,
      0.34577021976112826E+00,  0.00000000000000000E+00,
      -0.93831921813759156E+00, -0.93831921813759156E+00,
      0.34577021976112826E+00,  0.00000000000000000E+00,
      0.00000000000000000E+00,  -0.34577021976112826E+00,
      0.93831921813759156E+00,  0.93831921813759156E+00,
      0.00000000000000000E+00,  -0.34577021976112826E+00,
      -0.34577021976112826E+00, 0.93831921813759156E+00,
      0.00000000000000000E+00,  0.00000000000000000E+00,
      0.93831921813759156E+00,  -0.34577021976112826E+00,
      -0.34577021976112826E+00, 0.00000000000000000E+00,
      0.93831921813759156E+00,  0.93831921813759156E+00,
      -0.34577021976112826E+00, 0.00000000000000000E+00,
      0.00000000000000000E+00,  -0.34577021976112826E+00,
      -0.93831921813759156E+00, -0.93831921813759156E+00,
      0.00000000000000000E+00,  -0.34577021976112826E+00,
      -0.34577021976112826E+00, -0.93831921813759156E+00,
      0.00000000000000000E+00,  0.00000000000000000E+00,
      -0.93831921813759156E+00, -0.34577021976112826E+00,
      -0.34577021976112826E+00, 0.00000000000000000E+00,
      -0.93831921813759156E+00, -0.93831921813759156E+00,
      -0.34577021976112826E+00, 0.00000000000000000E+00,
      0.15904171053835295E+00,  0.52511857244364202E+00,
      0.83603601548245887E+00,  0.15904171053835295E+00,
      0.83603601548245887E+00,  0.52511857244364202E+00,
      0.83603601548245887E+00,  0.15904171053835295E+00,
      0.52511857244364202E+00,  0.52511857244364202E+00,
      0.15904171053835295E+00,  0.83603601548245887E+00,
      0.52511857244364202E+00,  0.83603601548245887E+00,
      0.15904171053835295E+00,  0.83603601548245887E+00,
      0.52511857244364202E+00,  0.15904171053835295E+00,
      0.15904171053835295E+00,  0.52511857244364202E+00,
      -0.83603601548245887E+00, 0.15904171053835295E+00,
      -0.83603601548245887E+00, 0.52511857244364202E+00,
      -0.83603601548245887E+00, 0.15904171053835295E+00,
      0.52511857244364202E+00,  0.52511857244364202E+00,
      0.15904171053835295E+00,  -0.83603601548245887E+00,
      0.52511857244364202E+00,  -0.83603601548245887E+00,
      0.15904171053835295E+00,  -0.83603601548245887E+00,
      0.52511857244364202E+00,  0.15904171053835295E+00,
      0.15904171053835295E+00,  -0.52511857244364202E+00,
      0.83603601548245887E+00,  0.15904171053835295E+00,
      0.83603601548245887E+00,  -0.52511857244364202E+00,
      0.83603601548245887E+00,  0.15904171053835295E+00,
      -0.52511857244364202E+00, -0.52511857244364202E+00,
      0.15904171053835295E+00,  0.83603601548245887E+00,
      -0.52511857244364202E+00, 0.83603601548245887E+00,
      0.15904171053835295E+00,  0.83603601548245887E+00,
      -0.52511857244364202E+00, 0.15904171053835295E+00,
      0.15904171053835295E+00,  -0.52511857244364202E+00,
      -0.83603601548245887E+00, 0.15904171053835295E+00,
      -0.83603601548245887E+00, -0.52511857244364202E+00,
      -0.83603601548245887E+00, 0.15904171053835295E+00,
      -0.52511857244364202E+00, -0.52511857244364202E+00,
      0.15904171053835295E+00,  -0.83603601548245887E+00,
      -0.52511857244364202E+00, -0.83603601548245887E+00,
      0.15904171053835295E+00,  -0.83603601548245887E+00,
      -0.52511857244364202E+00, 0.15904171053835295E+00,
      -0.15904171053835295E+00, 0.52511857244364202E+00,
      0.83603601548245887E+00,  -0.15904171053835295E+00,
      0.83603601548245887E+00,  0.52511857244364202E+00,
      0.83603601548245887E+00,  -0.15904171053835295E+00,
      0.52511857244364202E+00,  0.52511857244364202E+00,
      -0.15904171053835295E+00, 0.83603601548245887E+00,
      0.52511857244364202E+00,  0.83603601548245887E+00,
      -0.15904171053835295E+00, 0.83603601548245887E+00,
      0.52511857244364202E+00,  -0.15904171053835295E+00,
      -0.15904171053835295E+00, 0.52511857244364202E+00,
      -0.83603601548245887E+00, -0.15904171053835295E+00,
      -0.83603601548245887E+00, 0.52511857244364202E+00,
      -0.83603601548245887E+00, -0.15904171053835295E+00,
      0.52511857244364202E+00,  0.52511857244364202E+00,
      -0.15904171053835295E+00, -0.83603601548245887E+00,
      0.52511857244364202E+00,  -0.83603601548245887E+00,
      -0.15904171053835295E+00, -0.83603601548245887E+00,
      0.52511857244364202E+00,  -0.15904171053835295E+00,
      -0.15904171053835295E+00, -0.52511857244364202E+00,
      0.83603601548245887E+00,  -0.15904171053835295E+00,
      0.83603601548245887E+00,  -0.52511857244364202E+00,
      0.83603601548245887E+00,  -0.15904171053835295E+00,
      -0.52511857244364202E+00, -0.52511857244364202E+00,
      -0.15904171053835295E+00, 0.83603601548245887E+00,
      -0.52511857244364202E+00, 0.83603601548245887E+00,
      -0.15904171053835295E+00, 0.83603601548245887E+00,
      -0.52511857244364202E+00, -0.15904171053835295E+00,
      -0.15904171053835295E+00, -0.52511857244364202E+00,
      -0.83603601548245887E+00, -0.15904171053835295E+00,
      -0.83603601548245887E+00, -0.52511857244364202E+00,
      -0.83603601548245887E+00, -0.15904171053835295E+00,
      -0.52511857244364202E+00, -0.52511857244364202E+00,
      -0.15904171053835295E+00, -0.83603601548245887E+00,
      -0.52511857244364202E+00, -0.83603601548245887E+00,
      -0.15904171053835295E+00, -0.83603601548245887E+00,
      -0.52511857244364202E+00, -0.15904171053835295E+00};

  static constexpr std::array<T, 194> weights = {
      0.17823404472446112E-02, 0.17823404472446112E-02, 0.17823404472446112E-02,
      0.17823404472446112E-02, 0.17823404472446112E-02, 0.17823404472446112E-02,
      0.55733831788487382E-02, 0.55733831788487382E-02, 0.55733831788487382E-02,
      0.55733831788487382E-02, 0.55733831788487382E-02, 0.55733831788487382E-02,
      0.55733831788487382E-02, 0.55733831788487382E-02, 0.57169059499771017E-02,
      0.57169059499771017E-02, 0.57169059499771017E-02, 0.57169059499771017E-02,
      0.57169059499771017E-02, 0.57169059499771017E-02, 0.57169059499771017E-02,
      0.57169059499771017E-02, 0.57169059499771017E-02, 0.57169059499771017E-02,
      0.57169059499771017E-02, 0.57169059499771017E-02, 0.55187714672736135E-02,
      0.55187714672736135E-02, 0.55187714672736135E-02, 0.55187714672736135E-02,
      0.55187714672736135E-02, 0.55187714672736135E-02, 0.55187714672736135E-02,
      0.55187714672736135E-02, 0.55187714672736135E-02, 0.55187714672736135E-02,
      0.55187714672736135E-02, 0.55187714672736135E-02, 0.55187714672736135E-02,
      0.55187714672736135E-02, 0.55187714672736135E-02, 0.55187714672736135E-02,
      0.55187714672736135E-02, 0.55187714672736135E-02, 0.55187714672736135E-02,
      0.55187714672736135E-02, 0.55187714672736135E-02, 0.55187714672736135E-02,
      0.55187714672736135E-02, 0.55187714672736135E-02, 0.51582377118053833E-02,
      0.51582377118053833E-02, 0.51582377118053833E-02, 0.51582377118053833E-02,
      0.51582377118053833E-02, 0.51582377118053833E-02, 0.51582377118053833E-02,
      0.51582377118053833E-02, 0.51582377118053833E-02, 0.51582377118053833E-02,
      0.51582377118053833E-02, 0.51582377118053833E-02, 0.51582377118053833E-02,
      0.51582377118053833E-02, 0.51582377118053833E-02, 0.51582377118053833E-02,
      0.51582377118053833E-02, 0.51582377118053833E-02, 0.51582377118053833E-02,
      0.51582377118053833E-02, 0.51582377118053833E-02, 0.51582377118053833E-02,
      0.51582377118053833E-02, 0.51582377118053833E-02, 0.56087040825879972E-02,
      0.56087040825879972E-02, 0.56087040825879972E-02, 0.56087040825879972E-02,
      0.56087040825879972E-02, 0.56087040825879972E-02, 0.56087040825879972E-02,
      0.56087040825879972E-02, 0.56087040825879972E-02, 0.56087040825879972E-02,
      0.56087040825879972E-02, 0.56087040825879972E-02, 0.56087040825879972E-02,
      0.56087040825879972E-02, 0.56087040825879972E-02, 0.56087040825879972E-02,
      0.56087040825879972E-02, 0.56087040825879972E-02, 0.56087040825879972E-02,
      0.56087040825879972E-02, 0.56087040825879972E-02, 0.56087040825879972E-02,
      0.56087040825879972E-02, 0.56087040825879972E-02, 0.41067770281693937E-02,
      0.41067770281693937E-02, 0.41067770281693937E-02, 0.41067770281693937E-02,
      0.41067770281693937E-02, 0.41067770281693937E-02, 0.41067770281693937E-02,
      0.41067770281693937E-02, 0.41067770281693937E-02, 0.41067770281693937E-02,
      0.41067770281693937E-02, 0.41067770281693937E-02, 0.41067770281693937E-02,
      0.41067770281693937E-02, 0.41067770281693937E-02, 0.41067770281693937E-02,
      0.41067770281693937E-02, 0.41067770281693937E-02, 0.41067770281693937E-02,
      0.41067770281693937E-02, 0.41067770281693937E-02, 0.41067770281693937E-02,
      0.41067770281693937E-02, 0.41067770281693937E-02, 0.50518460646148088E-02,
      0.50518460646148088E-02, 0.50518460646148088E-02, 0.50518460646148088E-02,
      0.50518460646148088E-02, 0.50518460646148088E-02, 0.50518460646148088E-02,
      0.50518460646148088E-02, 0.50518460646148088E-02, 0.50518460646148088E-02,
      0.50518460646148088E-02, 0.50518460646148088E-02, 0.50518460646148088E-02,
      0.50518460646148088E-02, 0.50518460646148088E-02, 0.50518460646148088E-02,
      0.50518460646148088E-02, 0.50518460646148088E-02, 0.50518460646148088E-02,
      0.50518460646148088E-02, 0.50518460646148088E-02, 0.50518460646148088E-02,
      0.50518460646148088E-02, 0.50518460646148088E-02, 0.55302489162330935E-02,
      0.55302489162330935E-02, 0.55302489162330935E-02, 0.55302489162330935E-02,
      0.55302489162330935E-02, 0.55302489162330935E-02, 0.55302489162330935E-02,
      0.55302489162330935E-02, 0.55302489162330935E-02, 0.55302489162330935E-02,
      0.55302489162330935E-02, 0.55302489162330935E-02, 0.55302489162330935E-02,
      0.55302489162330935E-02, 0.55302489162330935E-02, 0.55302489162330935E-02,
      0.55302489162330935E-02, 0.55302489162330935E-02, 0.55302489162330935E-02,
      0.55302489162330935E-02, 0.55302489162330935E-02, 0.55302489162330935E-02,
      0.55302489162330935E-02, 0.55302489162330935E-02, 0.55302489162330935E-02,
      0.55302489162330935E-02, 0.55302489162330935E-02, 0.55302489162330935E-02,
      0.55302489162330935E-02, 0.55302489162330935E-02, 0.55302489162330935E-02,
      0.55302489162330935E-02, 0.55302489162330935E-02, 0.55302489162330935E-02,
      0.55302489162330935E-02, 0.55302489162330935E-02, 0.55302489162330935E-02,
      0.55302489162330935E-02, 0.55302489162330935E-02, 0.55302489162330935E-02,
      0.55302489162330935E-02, 0.55302489162330935E-02, 0.55302489162330935E-02,
      0.55302489162330935E-02, 0.55302489162330935E-02, 0.55302489162330935E-02,
      0.55302489162330935E-02, 0.55302489162330935E-02};
};
}  // namespace DelleyGrids
}  // namespace IntegratorXX