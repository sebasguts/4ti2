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

#ifndef _4ti2_zsolve__ParallelGraver_
#define _4ti2_zsolve__ParallelGraver_

#include "zsolve/Algorithm.hpp"
#include "zsolve/BitSet.h"
#include "zsolve/LinearSystem.hpp"
#include "zsolve/Controller.hpp"

namespace _4ti2_zsolve_
{

template <typename T> class Controller;

template <typename T> 
class ParallelGraver : public Algorithm <T>
{
private:
    size_t m_num_variables;

    void localInit ();

protected:
    Controller <T> * m_controller; ///< Controller for computation (unused until now)
    Lattice <T> * m_current_gens; ///< Generating set of a lattice.  Holds result after compute()

public:
    ParallelGraver () { };

    void init (LinearSystem <T> * system, Controller <T>* controller)
    {
        m_controller = controller;

        // system
	if (m_controller != NULL)
	    m_controller->log_system (system);

	//std::cout << *system << std::endl;

        // homogenized system
        LinearSystem <T> * homo = homogenize_linear_system (system);
	if (m_controller != NULL)
	    m_controller->log_homogenized_system (homo);

        // lattice
        m_current_gens = generate_lattice (homo);
        delete homo;
	localInit();
    }

    void init (Lattice <T> * lattice, Controller <T> * controller)
    {
        m_controller = controller;
        m_current_gens = new Lattice <T> (*lattice);

	localInit();
    }

    void init (std::ifstream& stream, Controller <T> * controller) { };

    ~ParallelGraver ()
    {
        delete m_current_gens;
    }

    Lattice <T>& lattice () const
    {
        return *m_current_gens;
    }

    void compute (int backup_frequency = 0);

    void extract_zsolve_results (VectorArray <T>& inhoms, VectorArray <T>& homs, VectorArray <T>& free) { };

    void extract_graver_results (VectorArray <T>& graver);
    
    void extract_hilbert_results (VectorArray <T>& hilbert) { };

    T extract_maxnorm_results (VectorArray <T> & maxnorm) { return 0;};

    size_t get_result_num_variables () const { return m_num_variables;};

    void log_maxnorm () { };
};

template <typename T>
void
ParallelGraver<T>::localInit (){
    // Log if desired:
    if (m_controller != NULL)
    {
	//m_controller->save_lattice (m_current_gens);
	m_controller->log_lattice (m_current_gens);
    }
    // Locally store number of variables
    m_num_variables = m_current_gens->num_variables();
}

template <typename T>
void 
ParallelGraver<T>::compute (int backup_frequency)
{
    // The big loop over the variables
    for (size_t varindex = 0; varindex < m_num_variables; varindex++){
	std::cout << "Now doing variable " << varindex << std::endl;
	/// ... 
    }

}

template <typename T>
void
ParallelGraver<T>::extract_graver_results (VectorArray<T> &graver) {
    // This assures that no variable has property -2 ???
    // assert (m_current_gens->get_splitter () == -1);
    // ??
    // assert (m_current_gens->get_result_num_variables () == m_num_variables);

    graver.clear ();

    for (size_t i = 0; i < m_current_gens->num_vectors (); i++)
    {
	T* vector = (*m_current_gens)[i];
	T* result = copy_vector <T> (vector, m_num_variables);

	std::cout << i << std::endl;
	print_vector (std::cout, result, m_num_variables);
	std::cout << std::endl;
	graver.append_vector (result);
    }
    
    if (m_controller != NULL)
	m_controller->log_result (1, m_current_gens->vectors (), 0);
};

} // namespace _4ti2_zsolve_
  
#endif
