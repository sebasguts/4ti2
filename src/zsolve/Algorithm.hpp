/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2013 4ti2 team.
Main author(s): Thomas Kahle.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA. 
*/

#ifndef _4ti2_zsolve__Algorithm_
#define _4ti2_zsolve__Algorithm_

#include <map>
#include "zsolve/BitSet.h"
#include "zsolve/LinearSystem.hpp"
#include "zsolve/Controller.hpp"
#include "zsolve/Norms.hpp"
#include "zsolve/Heuristics.hpp"
#include "zsolve/Timer.h"

namespace _4ti2_zsolve_
{

template <typename T> class Controller;


/**
 * \brief Interface to ZSolve algorithms
 *
 * This class contains only pure virtual methods and can not be
 * instantiated for its lack of a constructor.
 */
template <typename T> 
class Algorithm
{
public:
    virtual ~Algorithm () {};

    virtual void init (LinearSystem <T> * system, Controller <T>* controller) = 0;
    virtual void init (Lattice <T> * lattice, Controller <T> * controller) = 0;
    virtual void init (std::ifstream& stream, Controller <T> * controller) = 0;
    
    virtual void compute (int backup_frequency) = 0;
    
    virtual size_t get_result_variables () const = 0;
    virtual void log_maxnorm() = 0;

    virtual void extract_zsolve_results (  VectorArray <T>& inhoms, 
					   VectorArray <T>& homs, 
					   VectorArray <T>& free) = 0;

    virtual void extract_graver_results (  VectorArray <T>& graver) = 0;
    virtual void extract_hilbert_results (  VectorArray <T>& hilbert) = 0;
    virtual T extract_maxnorm_results (VectorArray <T> & maxnorm) = 0;
};


} // namespace _4ti2_zsolve_
  
#endif
