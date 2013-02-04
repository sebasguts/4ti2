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

#include <iostream>
#include <fstream>

#include "groebner/graver_main.h"
#include "groebner/Vector.h"
#include "groebner/VectorStream.h"
#include "groebner/VectorArray.h"
#include "groebner/VectorArrayStream.h"
#include "groebner/BitSet.h"
#include "groebner/BitSetStream.h"
#include "groebner/Feasible.h"
#include "groebner/FeasibleStream.h"
#include "groebner/GeneratingSet.h"
#include "groebner/GraverBasis.h"
#include "groebner/Globals.h"
#include "groebner/Options.h"

#include <string>

using namespace _4ti2_;

int
_4ti2_::graver_main(int argc, char **argv)
{
    Options::instance()->process_options(argc, argv);

    print_banner();

    std::cout << "In Graver main routine !!!" << std::endl;

    // Read in the sets of fibers.
    Feasible* feasible = input_Feasible(Options::instance()->filename.c_str());

    // The problem is correctly initialized now
    std::cout << feasible->get_matrix() << std::endl; 
    std::cout << feasible->get_basis() << std::endl; 

    // Read in the file with the generating set ?
    // std::string gens_filename(Options::instance()->filename + ".mar");
    // VectorArray* gens = input_VectorArray(feasible->get_dimension(), gens_filename.c_str());

    GraverBasis* gb = new GraverBasis ( *feasible );

    // Output the Graver basis.
    std::string graver_filename(Options::instance()->filename + ".gra");
    output(graver_filename.c_str(), gb->get_graver_basis());

    delete gb;
    delete feasible;

    return 0;
}
