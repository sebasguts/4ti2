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

#ifndef _4ti2_groebner__GraverBasis_
#define _4ti2_groebner__GraverBasis_

#include "groebner/Feasible.h"
#include "groebner/GeneratingSet.h"

namespace _4ti2_
{

class GraverBasis
{
public:
    GraverBasis( Feasible& feasible);
    GraverBasis( GeneratingSet& gs);
    GraverBasis( GraverBasis& gb);
    ~GraverBasis();

    const VectorArray& get_graver_basis();
    Feasible& get_feasible();

protected:
    GraverBasis();

    void compute();
    VectorArray* basis;  ///< Vectors in the Graver basis
    Feasible* feasible;  ///< Problem (to be) solved
};

} // namespace _4ti2_

#endif
