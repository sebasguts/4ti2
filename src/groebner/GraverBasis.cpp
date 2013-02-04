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

#include "GraverBasis.h"
#include "Globals.h"
#include "ParallelGraver.h"

#include <iostream>
#include <iomanip>
#include "VectorStream.h"
#include "VectorArrayStream.h"
#include "BitSetStream.h"

using namespace _4ti2_;

/// Compute Graver basis from a feasible problem
GraverBasis::GraverBasis( Feasible& _feasible)
{
    feasible = &_feasible;
    compute();
}

/// Compute Graver basis from Markov basis et al.
GraverBasis::GraverBasis( GeneratingSet& gs)
{
    feasible = &gs.get_feasible();
    basis = new VectorArray(gs.get_generating_set());
    compute();
}

/// Copy constructor
GraverBasis::GraverBasis( GraverBasis& gb)
{
    feasible = &gb.get_feasible();
    basis = new VectorArray(gb.get_graver_basis());
}

GraverBasis::~GraverBasis()
{
    delete basis;
}

const VectorArray&
GraverBasis::get_graver_basis()
{
    assert(gens != 0);
    return *basis;
}

Feasible&
GraverBasis::get_feasible()
{
    assert(feasible != 0);
    return *feasible;
}


void
GraverBasis::compute()
{
    assert(basis != 0 && feasible != 0);
    ParallelGraver algorithm;
    algorithm.compute(*feasible, *basis);
    basis->sort();
}
