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

#include <algorithm>
#include <future>
#include <map>
#include <stdexcept>

#include "zsolve/Algorithm.hpp"
#include "zsolve/BitSet.h"
#include "zsolve/Controller.hpp"
#include "zsolve/GraverJob.hpp"
#include "zsolve/LinearSystem.hpp"
#include "zsolve/Norms.hpp"
#include "zsolve/Vector.hpp"

namespace _4ti2_zsolve_
{

// Binary search tree for current lattice vectors, sorted according to
// (sub)-norms.  This is the way to do it in the future, but only
// supported in gcc-4.7 and up.  
// template <typename T> using NormBST = std::map<T, VectorArray<T>* >;

template <typename T> class Controller;

template <typename T> 
class ParallelGraver : public Algorithm <T>
{
    typedef std::map <T, VectorArray<T>* > NormBST;

private:
    size_t m_num_variables; ///< Number of variables in the Graver problem

    void localInit (); ///< Common initialization of local data
    NormBST* createNormBST (size_t dims);

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

/** 
 * \brief Update the binary search tree
 * 
 * @param dims How many dimensions to consider
 * 
 * @return new binary search tree for norms.
 */
template <typename T>
std::map< T, VectorArray<T>* >* 
ParallelGraver<T>::createNormBST (size_t dims) {
    /// @TODO: Parallelize this computation by splitting the todolist
    typedef std::map <T, VectorArray<T>* > NormBST;
    NormBST *result = new NormBST ();

    for (size_t i = 0; i < m_current_gens->num_vectors(); i++){
	T current_norm = norm_vector <T> ( (*m_current_gens)[i], dims);
	VectorArray<T> *current_vectors = NULL;
	// The map container gives a strong guarentee of not being
	// changed
	try {
	    current_vectors = result->at(current_norm);
	    current_vectors->append_vector( copy_vector<T> ( (*m_current_gens)[i], 
							     m_num_variables) );
	}
	catch (std::out_of_range& oor) {
	    current_vectors = new VectorArray<T> ();
	    current_vectors->append_vector( copy_vector<T> ( (*m_current_gens)[i], 
							     m_num_variables) );
	    result->insert( std::pair <T, VectorArray<T>* > 
			    (current_norm, current_vectors) );
	};
    }
    return result;
}

template <typename T>
void 
ParallelGraver<T>::compute (int backup_frequency)
{
    typedef std::map <T, VectorArray<T>* > NormBST;

    // First step: Add all negatives to basis.  This makes for
    // sign-consistent representations on the first variable:
    m_current_gens->append_negatives();

    // Now comes the big loop over the remaining variables
    for (size_t varindex = 1; varindex < m_num_variables; varindex++){
	std::cout << "Now doing variable " << varindex+1 << std::endl;

	NormBST *s = createNormBST(1);
	std::cout << "Number of different norms: " << s->size() << std::endl;

	T max_norm = s->rbegin()->first;
	std::cout << "Maximum norm: " << max_norm << std::endl;

	std::vector< NormPair<T> > jobs;
	// If there are degree one vectors, then we have our first job: combine them
	if (s->cbegin()->first == 1) 
	    jobs.push_back( NormPair<T>(1,1) );
	T completed_norm = 1;
	while (completed_norm < 2*max_norm) {
	    // Check if there are any not done jobs at the current level
	    std::vector < std::future < VectorArray<T>* > > m_futures;
	    for (auto it = jobs.begin(); it != jobs.end(); it++) {
		if ((*it).sum <= completed_norm+1) {
			// do it!
			std::cout << "DenkDenkDenk... " << (*it).sum << std::endl;
			// A job with the lowest sum is always in the beginning of the
			// jobs vector: (Check?)

			// Check if there are vectors with this norm at all:
			if (s->find(it->first) != s->end() && s->find(it->second) != s->end()) {
			    // Now we call the future.  This is tricky.  Since graverJob is a
			    // function template it can not be passed to async.  Instead we
			    // create an anonymous function that captures all relevant data
			    // and calls graverJob.  The compiler will then instantiate
			    // graverJob.
			    std::future < VectorArray<T>* > fut = std::async(
				std::launch::async, // Compiler can decide launch order?
				[s,m_current_gens, varindex, it] () {
				    std::cout << "Creating job"<< std::endl;
				    return graverJob (*(s->at(it->first)), 
						      *(s->at(it->second)), 
						      *(m_current_gens), 
						      varindex);});
			    // Wait for everybody to do their jobs NOW?
                            // fut.wait();
			    m_futures.push_back (std::move(fut)); // Lose ownership of the future
			}
		    }
	    }
	    // Synchronize
	    std::cout << "Waiting for results ...";
	    for (auto it = m_futures.begin(); it != m_futures.end(); ++it)
		it->wait();
	    std::cout << "Done." << std::endl;
	    // Retrieve results
	    for (auto it = m_futures.begin(); it != m_futures.end(); ++it){
		VectorArray<T> *result = (*it).get();
		if (result->num_vectors() == 0) // nothing new
		    continue;
		T norm = norm_vector<T> ((*result)[0], varindex); // Workaround, norm should be part of result
		std::cout << "Current Norm : " << norm << std::endl;
		// Check if maximum norm has increased
		if (max_norm < norm) {
		    std::cout << "New norm bound : " << norm << std::endl;
		    max_norm = norm;
		}
		// Join the result with what is already saved in s under this norm
		if (s->find(norm) == s->end()){
		    // not found
		    s->insert( std::pair <T, VectorArray<T>* > 
			       (norm, result));
		    // Ownership transferred to s
		} 
		else { 
		    // found this norm, append vectors.
		    /// @TODO Change this copy to a move!
		    for (size_t vc = 0; vc < result->num_vectors(); vc++)
			s->at(norm)->append_vector ( copy_vector<T> ((*result)[vc], m_num_variables));
		    delete result; // Must then not delete the vectors
		}
	    }
	    // Done with futures, clean up:
	    m_futures.clear();

	    // Now all jobs with norm sum <= completed_norm+1 are done.
	    // Schedule new jobs:
	    completed_norm++;
	    // Add new jobs for each norm pair (i, completed_norm), i=
	    // 1..completed_norm such that there are moves in the
	    // respective degrees.
	    if (s->find(completed_norm) != s->end())
		for (T i = 1; i <= completed_norm; i++)
		    if (s->find(i) != s->end())
			jobs.push_back ( NormPair<T> (i, completed_norm));
	    // Keep jobs sorted according to total norm
	    std::sort(jobs.begin(), jobs.end());
	}
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

	// In the output we remove symmetric copies
	bool has_symmetric = true;
	// Check if there is a symmetric vector
	for (size_t j = 0; j < m_num_variables; j++)
	    if (!m_current_gens->get_variable (j).check_bounds (-vector[j]))
		has_symmetric = false;
	// Check if this vector or its negative should be removed:
	int lex_cmp = lex_compare_vector_with_negative (vector, m_num_variables);
	if (!has_symmetric || lex_cmp > 0)
	    graver.append_vector (result);
    }
    
    if (m_controller != NULL)
	m_controller->log_result (1, m_current_gens->vectors (), 0);
};

} // namespace _4ti2_zsolve_
  
#endif
