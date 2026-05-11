/* Random number generation for certain distributions.
 * 
 * Created on 15/06/2018
 * Copyright (C) 2018 Attila Szabó <as2372@cam.ac.uk>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the Creative Commons Attribution License (CC-BY),
 * version 4.0 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Please find a copy of the Creative Commons CC-BY License on
 * <https://creativecommons.org/licenses/by/4.0/>.
 */

#ifndef random_hh
#define random_hh

#include <random>
#include <complex>
#include <cmath>



/* Returns a sample from the von Mises distribution
 *     p(x) ~ exp(kappa * cos(x-mu))
 * in the form exp(i*x).
 *
 * Uses the rejection algorithm compared against the wrapped Cauchy
 * distribution suggested by Best and Fisher and documented in
 * Chapter 9 of Luc's Non-Uniform Random Variate Generation.
 *
 * Addition from A.Sz.: very large kappa is handled with a wrapped Gaussian
 * distribution. It is accurate to almost macheps there.
 *
 * Code adapted from the numpy.random library which has the following
 * copyright notice:
 *
 * Copyright 2005 Robert Kern (robert.kern@gmail.com)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
template <class Generator>
double von_mises(Generator& g, double mu, double kappa) {
    using namespace std;

    /* if kappa is very large, the distribution is basically Gaussian */
    if (kappa > 1e14) {
        normal_distribution<double> normal(mu,pow(kappa,-0.5));
        return  normal(g);
    }
    
    double s;
    double U, V, W, Y, Z;
    double result;
    
    uniform_real_distribution<double> uniform(0.0,1.0);

    /* if kappa is very small, the distribution is uniform */
    if (kappa < 1e-8) {
        return 2*M_PI*uniform(g);
    }
    else {
        /* with double precision rho is zero until 1.4e-8 */
        if (kappa < 1e-5) {
            /* second order taylor expansion around kappa = 0
             * precise until relatively large kappas as second order is 0 */
            s = (1./kappa + kappa);
        }
        else {
            double r = 1 + sqrt(1 + 4*kappa*kappa);
            double rho = (r - sqrt(2*r)) / (2*kappa);
            s = (1 + rho*rho)/(2*rho);
        }

        while (1) {
            U = uniform(g);
            Z = cos(M_PI*U);
            W = (1 + s*Z)/(s + Z);
            Y = kappa * (s - W);
            V = uniform(g);
            if ((Y*(2-Y) - V >= 0) || (log(Y/V)+1 - Y >= 0))
                break;
        }

        U = uniform(g);
        result = acos(W);
        if (U < 0.5)
            result = -result;
        return result+mu;
    }
}

#endif
