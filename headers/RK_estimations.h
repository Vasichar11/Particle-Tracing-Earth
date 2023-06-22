#pragma once

#include "common.h"
#include "constants.h"

//Function to approximate all values of Runge Kutta's block.
void new_values_RK4(real &latitude, real &ppar, real &pper, real &eta, real &alpha, real &aeq, real l1, real l2, real l3, real l4, real m1, real m2, real m3, real m4, real n1, real n2, real n3, real n4, real o1, real o2, real o3, real o4, real p1, real p2, real p3, real p4, real q1, real q2, real q3, real q4);

//Overloaded for noWPI -> without eta,aeq estimation
void new_values_RK4(real &latitude, real &ppar, real &pper, real &alpha, real l1, real l2, real l3, real l4, real m1, real m2, real m3, real m4, real o1, real o2, real o3, real o4, real p1, real p2, real p3, real p4);
