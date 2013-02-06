/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2008 4ti2 team.
Main author(s): Peter Malkin.

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

#ifndef _4ti2_zsolve__GraverAPI_
#define _4ti2_zsolve__GraverAPI_

#include "4ti2/4ti2xx.h"
#include "zsolve/ZSolveAPI.hpp"
#include "zsolve/ExtendedPottier.hpp"

namespace _4ti2_zsolve_ {

template <class T>
class GraverAPI : public ZSolveAPI<T> {
public:
    GraverAPI();

    virtual void compute ();

    virtual void write(const char* project);

    virtual _4ti2_matrix* get_matrix(const char* name);

protected:
    virtual void check_consistency();

    // Extract the output after running the algorithm.
    virtual void extract_results(Algorithm <T>* algorithm);
};

template <class T>
GraverAPI<T>::GraverAPI()
{
    ZSolveAPI<T>::free_default = false;
    ZSolveAPI<T>::lower_default = 1;
    ZSolveAPI<T>::upper_default = -1;
}

/** 
 * \brief Run the actual computation
 * 
 * Calling this function kicks off the actual computation by
 * instantiating an Algorithm and calling its compute function.
 */
template <class T>
void
GraverAPI<T>::compute()
{
    check_consistency();
    
    Algorithm <T>* algorithm;
    DefaultController <T> * controller;
    std::ofstream* log_file = 0;
    if (ZSolveAPI<T>::options.loglevel () > 0) {
        std::string log_name = ZSolveAPI<T>::options.project () + ".log";
        log_file = new std::ofstream (log_name.c_str(), ZSolveAPI<T>::options.resume () ? std::ios::out | std::ios::app : std::ios::out);
    }
    controller = new DefaultController <T> (&std::cout, log_file, ZSolveAPI<T>::options);


    if (ZSolveAPI<T>::mat) {
        /// @TODO: transfer rhs, ub, lb, sign and rel.
        T* rhs_vec = create_zero_vector <T> (ZSolveAPI<T>::mat->data.height());
        if (ZSolveAPI<T>::rhs) { 
            for (size_t i = 0; i < ZSolveAPI<T>::rhs->data.width(); ++i) {
                rhs_vec[i] = ZSolveAPI<T>::rhs->data[0][i];
            }
        }
        LinearSystem <T> * system = new LinearSystem <T> (ZSolveAPI<T>::mat->data, rhs_vec,
							  ZSolveAPI<T>::free_default, 
							  ZSolveAPI<T>::lower_default, 
							  ZSolveAPI<T>::upper_default);
        delete_vector(rhs_vec);
        if (ZSolveAPI<T>::sign) {
            for (size_t i = 0; i < ZSolveAPI<T>::sign->data.width(); ++i) {
                switch (ZSolveAPI<T>::sign->data[0][i]) {
		case 0:
		    system->get_variable(i).set(true);
		    break;
		case 1:
		    system->get_variable(i).set(false, 0, -1);
		    break;
		case -1:
		    system->get_variable(i).set(false, 1, 0);
		    break;
		case 2:
		    system->get_variable(i).set(false);
		    break;
		default:
		    /// @TODO: The following error message should be more informative.
		    throw IOException("Unknown sign value.");
                }
            }
        }
        if (ZSolveAPI<T>::rel) {
            for (size_t i = 0; i < ZSolveAPI<T>::rel->data.width(); ++i) {
                switch (ZSolveAPI<T>::rel->data[0][i]) {
		case 0:
		    system->get_relation(i).set(Relation<T> :: Equal);
		    break;
		case 1:
		    system->get_relation(i).set(Relation<T> :: GreaterEqual);
		    break;
		case -1:
		    system->get_relation(i).set(Relation<T> :: LesserEqual);
		    break;
		default:
		    /// @TODO: The following error message should be more informative.
		    throw IOException("Unknown relation value.");
                }
            }
        }
        if (ZSolveAPI<T>::lb) {
            for (size_t i = 0; i < ZSolveAPI<T>::lb->data.width(); ++i) {
                system->get_variable(i).set_bound(true, ZSolveAPI<T>::lb->data[0][i]);
            }
        }
        if (ZSolveAPI<T>::ub) {
            for (size_t i = 0; i < ZSolveAPI<T>::ub->data.width(); ++i) {
                system->get_variable(i).set_bound(false, ZSolveAPI<T>::ub->data[0][i]);
            }
        }

        system->cancel_down();
        algorithm = new ExtendedPottier <T>;
	algorithm->init(system, controller);
        delete system;
    }
    else if (ZSolveAPI<T>::lat) {
        /// @TODO: transfer ub, lb, and sign.
        Lattice <T> * lattice = new Lattice <T> (& ZSolveAPI<T>::lat->data, 
						 ZSolveAPI<T>::free_default, 
						 ZSolveAPI<T>::lower_default, 
						 ZSolveAPI<T>::upper_default);
        if (ZSolveAPI<T>::sign) {
            for (size_t i = 0; i < ZSolveAPI<T>::sign->data.width(); ++i) {
                switch (ZSolveAPI<T>::sign->data[0][i]) {
		case 0:
		    lattice->get_variable(i).set(true);
		    break;
		case 1:
		    lattice->get_variable(i).set(false, 0, -1);
		    break;
		case -1:
		    lattice->get_variable(i).set(false, 1, 0);
		    break;
		case 2:
		    lattice->get_variable(i).set(false);
		    break;
		default:
		    /// @TODO: The following error message should be more informative.
		    throw IOException("Unknown sign value.");
                }
            }
        }
        if (ZSolveAPI<T>::lb) {
            for (size_t i = 0; i < ZSolveAPI<T>::lb->data.width(); ++i) {
                lattice->get_variable(i).set_bound(true, ZSolveAPI<T>::lb->data[0][i]);
            }
        }
        if (ZSolveAPI<T>::ub) {
            for (size_t i = 0; i < ZSolveAPI<T>::ub->data.width(); ++i) {
                lattice->get_variable(i).set_bound(false, ZSolveAPI<T>::ub->data[0][i]);
            }
        }

        lattice->reduce_gaussian();
        algorithm = new ExtendedPottier <T>;
	algorithm->init(lattice, controller);
        delete lattice;
    }
    else {
        throw IOException ("Neither " + ZSolveAPI<T>::options.project () + ".mat, " + 
			   ZSolveAPI<T>::options.project () + ".lat, nor "
			   + ZSolveAPI<T>::options.project () + ".backup found!");
    }

    // Actual computation starts here.
    algorithm->compute (ZSolveAPI<T>::options.backup_frequency ());

    algorithm->log_maxnorm ();

    extract_results(algorithm);

    delete algorithm;
    delete controller;
    if (log_file) { delete log_file; }
}

template <class T>
void
GraverAPI<T>::check_consistency()
{
    ZSolveAPI<T>::check_consistency();

    if (ZSolveAPI<T>::rhs) {
        throw IOException ("No `rhs' allowed for `graver' executable. Use `zsolve' instead!\n");
    }
    if (ZSolveAPI<T>::rel) {
        throw IOException ("No `rel' allowed for `graver' executable. Use `zsolve' instead.");
    }
}

template <class T>
void
GraverAPI<T>::write(const char* project_c_str)
{
    std::string project(project_c_str);   

    if (ZSolveAPI<T>::zhom) { ZSolveAPI<T>::zhom->write((project + ".gra").c_str()); }
    if (ZSolveAPI<T>::zfree && ZSolveAPI<T>::zfree->data.height() > 0) {
        ZSolveAPI<T>::zfree->write((project + ".zfree").c_str());
    }
}

template <class T>
_4ti2_matrix*
GraverAPI<T>::get_matrix(const char* name)
{
    if (!strcmp(name, "gra")) { return ZSolveAPI<T>::zhom; }
    return ZSolveAPI<T>::get_matrix(name);
}

template <class T>
void
GraverAPI<T>::extract_results(Algorithm <T>* algorithm)
{
    delete ZSolveAPI<T>::zhom;
    ZSolveAPI<T>::zhom = new VectorArrayAPI <T> (0, algorithm->get_result_num_variables ());
    algorithm->extract_graver_results (ZSolveAPI<T>::zhom->data);
}

} // namespace _4ti2_zsolve_

#endif
