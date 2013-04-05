/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2013 4ti2 team.
Main author(s): Mathias Walter, Thomas Kahle.

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

#ifndef _4ti2_zsolve__ParallelPottier_
#define _4ti2_zsolve__ParallelPottier_

#define NUMJOBS 20

#include <algorithm>
#include <chrono>
#include <map>
#include <mutex>
#include <future>

#include "zsolve/Algorithm.hpp"
#include "zsolve/BitSet.h"
#include "zsolve/LinearSystem.hpp"
#include "zsolve/Controller.hpp"
#include "zsolve/Norms.hpp"
#include "zsolve/Heuristics.hpp"
#include "zsolve/Timer.h"
#include "zsolve/ValueTree.hpp"
#include "zsolve/UniqueVectorsHash.hpp"

namespace _4ti2_zsolve_
{

template <typename T>
struct threeTempVectors {
    T *first;
    T *second;
    T *sum;
};

template <typename T>
struct enum_second_job {
    threeTempVectors<T> temp_vectors;
    std::future <void> fut;
    std::vector <T*> *result_vector;
};

template <typename T> class Controller;

template <typename T> 
class ParallelPottier : public Algorithm <T>
{
protected:
    typedef std::map <T, ValueTree <T> *> RootMap;
    typedef std::map <NormPair <T>, bool> NormMap;

    Controller <T> * m_controller; ///< Controller for computation
    Lattice <T> * m_lattice; ///< Integer lattice that contains the solutions

    T m_maxnorm; ///< current maximum norm
    size_t m_current_variable; ///< current component
    size_t m_variables; ///< components

    NormMap m_norms;
    RootMap m_roots; ///< roots of valuetrees

    std::mutex m_resultMapMutex;
    std::map <NormPair <T>, std::vector <T*> > m_resultMap; // Locations to store result vectors

    enum_second_job<T> m_jobs[NUMJOBS];
    bool m_jobsFree[NUMJOBS];

    Timer m_backup_timer; ///< timer
    int m_backup_frequency; ///< frequency for backup

protected:

    /**
     * \brief ???
     *
     * 
     */
    void insert_tree (ValueTree <T> *& tree, size_t vid, bool split_recursive)
    {
        if (tree->level < 0)
        {
//            std::cout << "insert_tree::insert " << vid << std::endl;
            tree->vector_indices.push_back (vid);
            if (split_recursive)
                split_tree (tree);
        }
        else
        {
            T value = (*m_lattice)[vid][tree->level];
//            std::cout << "insert_tree::search on level " << tree->level << " with value " << value << std::endl;
            if (value > 0)
            {
		auto iter = tree->pos.begin();
                for (; iter != tree->pos.end (); iter++)
                    if (value <= (*iter)->value)
                        break;
                if (iter == tree->pos.end() || value != (*iter)->value)
                    tree->pos.insert (iter, new ValueTreeNode <T> (value, vid));
                else
                    insert_tree ((*iter)->sub_tree, vid, split_recursive);
            }
            else if (value < 0)
            {
		auto iter = tree->neg.begin();
                for (; iter != tree->neg.end (); iter++)
                    if (value >= (*iter)->value)
                        break;
                if (iter == tree->neg.end() || value != (*iter)->value)
                    tree->neg.insert (iter, new ValueTreeNode <T> (value, vid));
                else
                    insert_tree ((*iter)->sub_tree, vid, split_recursive);
                
            }
            else
            {
                if (tree->zero == NULL)
                    tree->zero = new ValueTree <T> ();
                insert_tree (tree->zero, vid, split_recursive);
            }
        }
    }

    void split_tree (ValueTree <T> * tree, int start = -1)
    {
        int compo = start < 0 ? m_current_variable : start;
        bool has_pos, has_neg, has_zero;

        if (tree->level >= 0)
            return;

        for (; start < (int) m_current_variable; start++)
        {
            compo = start < 0 ? m_current_variable : start;
            has_pos = has_neg = has_zero = false;
            for (size_t i = 0; i < tree->vector_indices.size (); i++)
            {
                T value = (*m_lattice)[tree->vector_indices[i]][compo];
                if (value > 0)
                    has_pos = true;
                else if (value < 0)
                    has_neg = true;
                else
                    has_zero = true;
                if (has_pos && has_neg)
                    break;
            }
            if (has_pos && has_neg)
                break;
        }

        if ((start < (int) m_current_variable) && (tree->vector_indices.size () >= 1))
        {
//            std::cout << "Splitting on " << compo << std::endl;
            tree->level = compo;
            for (size_t i = 0; i < tree->vector_indices.size (); i++)
                insert_tree (tree, tree->vector_indices[i], false);
            start++;
            if (tree->zero != NULL)
                split_tree (tree->zero, start);
            for (size_t i = 0; i < tree->pos.size (); i++)
                split_tree (tree->pos[i]->sub_tree, start);
            for (size_t i = 0; i < tree->neg.size (); i++)
                split_tree (tree->neg[i]->sub_tree, start);
        }
    }

    void insert_trees (T* vector, T norm)
    {
        int vid = m_lattice->append_vector (copy_vector (vector, m_variables));

        if (m_roots.find (norm) == m_roots.end ())
        {
            m_roots[norm] = new ValueTree <T> ();
            
            for (auto i = m_roots.begin (); i != m_roots.end (); i++)
            {
                NormPair <T> pair (i->first, norm);
                m_norms[pair] = true;
            }
        }

        insert_tree (m_roots[norm], vid, true);
    }

    int wait_a_little (std::vector<T*>& results) { // returns job id.
	int empty_job_id = -1;
	while (empty_job_id == -1) {
	    for (uint i = 0; i < NUMJOBS; i++){
		if (m_jobsFree[i]) // good, just use i
		    return i;
		// Contrary to the standard, wait_for returns a bool.  Let's assume it means ready or not ready.
		// std::future_status status = m_jobs[i].fut.wait_for(std::chrono::seconds(1));
		bool status_ready = m_jobs[i].fut.wait_for(std::chrono::seconds(0));
		if (status_ready) {
		    // This job is ready
		    for (auto it = m_jobs[i].result_vector->begin(); it != m_jobs[i].result_vector->end(); it++){
			results.push_back(*it);
		    }
		    m_jobs[i].result_vector->clear();
		    m_jobsFree[i] = true;
		    empty_job_id = i;
		}
	    }
	    if (empty_job_id==-1) 	// No job was ready
		std::this_thread::sleep_for(std::chrono::milliseconds(2));
	}
	return empty_job_id;
    }

    void enum_first (ValueTree <T> * tree, const NormPair<T>& norms, std::vector<T*>& results)
    {
        if (tree->level < 0)
        {
	    // enum_second_job<T> *job;

// 	    std::vector<T*> *current_result_vector;
// 	    std::vector < std::vector<T*> * > resultDumpDump;
// 	    std::vector < T* > vectors_to_clean;
// 	    std::vector < std::future <void> > futures;
            for (size_t i = 0; i < tree->vector_indices.size(); i++)
            {
		T *first = (*m_lattice)[tree->vector_indices[i]];
                //std::cout << "enum_first enumerated [";
                //print_vector (std::cout, m_first_vector, m_variables);
                //std::cout << "]" << std::endl;
                if (first[m_current_variable] > 0) {
		    int job_slot = wait_a_little(results);
		    if (job_slot < 0 || job_slot >= NUMJOBS) {
			std::cout << "an error occurred\n";
			exit(1);
		    }
		    enum_second_job<T> *job = &(m_jobs[job_slot]);
		    job->temp_vectors.first = first;
		    std::future<void> fut;
		    // We do start many new threads in this program
		    // and once in a while this fails.  In this case,
		    // just wait a little and try again.
		    bool start_success = false;
		    while (!start_success) {
			try {
			    fut = std::async(
				std::launch::async,
				&ParallelPottier<T>::enum_second,
				this,
				m_roots[norms.second], 
				job->temp_vectors,
				norms,
				std::ref(*(job->result_vector)));
			    start_success = true;
			}
			catch (std::system_error e){
			    std::cout << "Failed to launch a new thread, will keep trying ..." << std::endl;
			    start_success = false;
			}
		    }
		    job->fut = std::move (fut);
		    m_jobsFree[job_slot]=false;
		}
            }
        }
        else
        {
            if (tree->zero != NULL)
                enum_first (tree->zero, norms, results);
            for (size_t i = 0; i < tree->pos.size (); i++)
                enum_first (tree->pos[i]->sub_tree, norms, results);
            for (size_t i = 0; i < tree->neg.size (); i++)
                enum_first (tree->neg[i]->sub_tree, norms, results);
        }
    }

    void enum_second (ValueTree <T> * tree, threeTempVectors<T>& tmp, const NormPair<T>& norms, std::vector<T*>& resultDump)
    {
        if (tree->level < 0) // arrived at leaf
        {
            for (size_t i = 0; i < tree->vector_indices.size(); i++)
            {
                tmp.second = (*m_lattice)[tree->vector_indices[i]];
                // std::cout << "enum_second enumerated [";
                // print_vector (std::cout, tmp.second, m_variables);
                // std::cout << "]" << std::endl;
                build_sum (tmp, norms, resultDump);
            }
            if (tree->level >= 0)
                enum_second (tree, tmp, norms, resultDump);
        }
        else if (tree->level == (int) m_current_variable)  // Comparing currently being lifted column
        {
            T value = tmp.first[tree->level];
            if (value <= 0)
                for (size_t i = 0; i < tree->pos.size(); i++)
                    enum_second (tree->pos[i]->sub_tree, tmp, norms, resultDump);
            if (value >= 0)
                for (size_t i = 0; i < tree->neg.size(); i++)
                    enum_second (tree->neg[i]->sub_tree, tmp, norms, resultDump);
        }
        else  // Comparing any other column
        {
            T value = tmp.first[tree->level];
            if (tree->zero != NULL)
                enum_second (tree->zero,tmp, norms, resultDump);
            if (value >= 0)
                for (size_t i = 0; i < tree->pos.size(); i++)
                    enum_second (tree->pos[i]->sub_tree,tmp, norms, resultDump);
            if (value <= 0)
                for (size_t i = 0; i < tree->neg.size(); i++)
                    enum_second (tree->neg[i]->sub_tree,tmp, norms, resultDump);
        }
    }

    bool enum_reducer (ValueTree <T> * tree, threeTempVectors<T>& tmp)
    {
        if (tree->level < 0)
        {
            for (int i = tree->vector_indices.size ()-1; i >= 0; i--)
            {
                T* reducer = (*m_lattice)[tree->vector_indices[i]];
                bool flag = true;
                for (size_t j = 0; j <= m_current_variable; j++)
                {
                    if (reducer[j] < 0)
                    {
                        if (tmp.sum[j] >= 0 || abs (reducer[j]) > abs (tmp.sum[j]))
                        {
                            flag = false;
                            break;
                        }
                    }
                    else if (reducer[j] > 0)
                    {
                        if (tmp.sum[j] <= 0 || abs (reducer[j]) > abs (tmp.sum[j]))
                        {
                            flag = false;
                            break;
                        }
                    }
                }
                if (flag)
                {
                    //std::cout << "Reduced (";
                    //print_vector (std::cout, tmp.sum, m_current_variable+1);
                    //std::cout << ") by (";
                    //print_vector (std::cout, reducer, m_current_variable+1);
                    //std::cout << ")" << std::endl;
                    //count_reduces++;
                    return true;
                }
            }
        }
        else
        {
            T value = tmp.sum[tree->level];
            if (value > 0)
            {
                for (auto iter = tree->pos.begin (); iter != tree->pos.end(); iter++)
                {
                    if ((*iter)->value > value)
                        break;
                    if (enum_reducer ((*iter)->sub_tree, tmp))
                        return true;
                }
            }
            else if (value < 0)
            {
                for (auto iter = tree->neg.begin (); iter != tree->neg.end(); iter++)
                {
                    if ((*iter)->value < value)
                        break;
                    if (enum_reducer ((*iter)->sub_tree, tmp))
                        return true;
                }
                
            }
            if (tree->zero != NULL)
                if (enum_reducer (tree->zero, tmp))
                    return true;
        }
        return false;
    }

    void build_sum (threeTempVectors<T>& tmp, const NormPair<T>& norms, std::vector<T*>& resultDump)
    {
        // std::cout << "buildSum (";
        // print_vector (std::cout, tmp.first, m_variables);
        // std::cout << ") + (";
        // print_vector (std::cout, tmp.second, m_variables);
        // std::cout << ")" << std::endl;

        if (tmp.first == tmp.second)
            return;

        //count_builds++;

        if (tmp.first[m_current_variable] <= 0 && tmp.second[m_current_variable] <= 0)
            return;
        if (tmp.first[m_current_variable] >= 0 && tmp.second[m_current_variable] >= 0)
            return;

        for (size_t i = 0; i < m_current_variable; i++)
        {
            if (tmp.first[i] < 0 && tmp.second[i] > 0)
                return;
            if (tmp.first[i] > 0 && tmp.second[i] < 0)
                return;
        }

	// TODO: Quickcheck goes here

        //bool flag = false;
        for (size_t i = 0; i < m_variables; i++)
        {
            tmp.sum[i] = tmp.first[i] + tmp.second[i];
            //if (tmp.second[i] != 0)
            //    flag = true;

//	    int current_precision = calcPrecision (tmp.sum[i]);
//	    int max_precision = maxPrecision (tmp.sum[i]);
//	    if (max_precision > 0 && current_precision + 1 >= max_precision)
//	    {
//		throw PrecisionException (max_precision);
//	    }
        }
        //if (!flag)
        //    std::cerr << "\n\nGENERATED ZERO VECTOR!!!!\n\n" << std::endl;

	// This is equal to norms.sum, no need to recompute
        // T norm = norm_vector (tmp.sum, m_current_variable);

//        if (norm == 0)
//            return;

        // TODO: norm / 2 ??
        for (auto iter = m_roots.begin (); iter != m_roots.end() && iter->first <= norms.sum/2; iter++)
        {
            if (enum_reducer (iter->second, tmp))
                return;
        }

	// Delegate clean-up to receiver of the resultDump
	resultDump.push_back(copy_vector<T> (tmp.sum, m_variables));
    }

/*
    void dump_tree (ValueTree <T> * tree)
    {
        if (tree->level < 0)
        {
            std::cout << "dump_tree::leaf:\n";
            for (size_t i = 0; i < tree->vector_indices.size(); i++)
            {
                std::cout << "  [" << tree->vector_indices[i] << "] = ";
                print_vector (std::cout, (*m_lattice)[tree->vector_indices[i]], m_variables);
                std::cout << std::endl;
            }
        }
        else
        {
            std::cout << "dump_tree::node at level " << tree->level << "\n";
            if (tree->zero != NULL)
            {
                std::cout << "dump_tree::zero\n";
                dump_tree (tree->zero);
            }
            for (size_t i = 0; i < tree->pos.size (); i++)
            {
                std::cout << "dump_tree::pos (" << tree->pos[i]->value << ")" << std::endl;
                dump_tree (tree->pos[i]->sub_tree);
            }
            for (size_t i = 0; i < tree->neg.size (); i++)
            {
                std::cout << "dump_tree::neg (" << tree->neg[i]->value << ")" << std::endl;
                dump_tree (tree->neg[i]->sub_tree);
            }
        }
    }

    void dump_trees ()
    {
        std::cout << "============================ DUMP ==========================" << std::endl;
        for (auto iter = m_roots.begin (); iter != m_roots.end (); iter++)
        {
            std::cout << "dump_trees::root[" << iter->first << "]" << std::endl;
            dump_tree (iter->second);
            std::cout << std::endl;
        }
        std::cout << "============================ /DUMP ==========================" << std::endl;
    }
*/
    
    
    /** 
     * \brief Initialize reduction trees
     * 
     */
    void create_trees ()
    {
        m_maxnorm = -1;
        for (size_t i = 0; i < m_lattice->vectors (); i++)
        {
            T* vector = (*m_lattice)[i];
            T norm = norm_vector (vector, m_current_variable);
            if (norm == 0 && vector[m_current_variable] == 0)
                continue;
            m_maxnorm = max (m_maxnorm, norm);
            if (m_roots.find (norm) == m_roots.end ())
            {
                m_roots[norm] = new ValueTree <T> ();
            }
            //std::cout << "Inserting vector " << i << " into roots[" << norm << "]" << std::endl;
            insert_tree (m_roots[norm], i, false);
        }

//        std::cout << "Before SPLIT:\n\n";
//        dump_trees ();

        for (auto iter = m_roots.begin (); iter != m_roots.end (); iter++)
        {
            split_tree (iter->second);
        }

//        std::cout << "After SPLIT:\n\n";
//        dump_trees ();
    }

    void delete_trees ()
    {
        for (auto iter = m_roots.begin (); iter != m_roots.end (); iter++)
        {
            delete iter->second;
        }
        m_roots.clear ();
        m_maxnorm = -1;
    }

    int chooseNextVariable ()
    {
        //std::cout << "chooseNextVariable " << std::endl;
        BitSet allowed (m_variables, true);
        BitSet tempset (m_variables, false);

        // handled, free, range
        T best_range = 0;
        int best_infinity_count = 3;
        for (size_t i = 0; i < m_variables; i++)
        {
            VariableProperty <T> & var = m_lattice->get_variable (i);
            if ( i < m_current_variable || var.free ())
            {
                allowed.unset (i);
                continue;
            }

            int infinity_count = var.count_infinity ();
            T range = var.get_range ();

            if (infinity_count < best_infinity_count || (infinity_count == best_infinity_count && range < best_range))
            {
                best_range = range;
                best_infinity_count = infinity_count;
                tempset.zero ();
                tempset.set (i);
            }
            else if (infinity_count == best_infinity_count && range == best_range)
                tempset.set (i);
        }
        allowed.set_intersection (tempset);

        // best gcd
        T best_gcd = -1;
        tempset.zero ();


        for (size_t i = 0; i < m_variables; i++)
        {
            if (!allowed[i])
                continue;
            T gcd = m_lattice->gcd_column (i, 0, m_lattice->vectors ());
            if (best_gcd < 0 || gcd < best_gcd)
            {
                best_gcd = gcd;
                tempset.zero ();
                tempset.set (i);
            }
            else if (gcd == best_gcd)
                tempset.set (i);
        }
        allowed.set_intersection (tempset);

        return Heuristics <T> :: chooseNextVariable (*m_lattice, allowed);
    }
    
    /**
     * \brief Do some mysterious preprocessing
     *
     * This seems to be a row reduction.
     */
    void preprocess ()
    {
        T* last_reducer = NULL;
        bool repeat;
        do
        {
            repeat = false;
            for (size_t reducer_row = 0; reducer_row < m_lattice->vectors (); reducer_row++)
            {
                T* reducer = (*m_lattice)[reducer_row];
                T norm = norm_vector (reducer, m_current_variable);
                if (norm != 0 || reducer[m_current_variable] == 0)
                    continue;

                last_reducer = reducer;

                for (size_t row = 0; row < m_lattice->vectors (); row++)
                {
                    if (row == reducer_row)
                        continue;
                    T* vec = (*m_lattice)[row];
                    if (abs(vec[m_current_variable]) >= abs(reducer[m_current_variable]))
                    {
                        T factor = abs(vec[m_current_variable]) / abs(reducer[m_current_variable]);
                        if (factor != 0)
                        {
                            repeat = true;
                            if (vec[m_current_variable] * reducer[m_current_variable] > 0)
                                factor = -factor;
                            m_lattice->combine_rows (row, factor, reducer_row);
                        }
                    }
                }
            }
        }
        while (repeat);
        if (last_reducer != NULL)
        {
            last_reducer = copy_vector (last_reducer, m_lattice->variables ());
            negate_vector (last_reducer, m_lattice->variables ());
            m_lattice->append_vector (last_reducer);
        }
    }

    void complete (const NormPair<T>& norms)
    {
	if (m_controller != NULL)
	    m_controller->log_status (m_current_variable+1, norms.sum, m_maxnorm, norms.first, m_lattice->vectors (), m_backup_frequency, m_backup_timer);
        //std::cout << "Variable: " << m_current_variable+1 << ", Sum: " << norms.first << " + " << norms.second << " = " << norms.sum << ", Solutions = " << m_lattice->vectors () << std::endl;

        //dump_trees ();

        if ((m_roots.find (norms.first) != m_roots.end ()) && (m_roots.find (norms.second) != m_roots.end ()))
        {
	    for (int i = 0; i < NUMJOBS; i++){
		m_jobsFree[i] = true;
		m_jobs[i].temp_vectors.sum = create_vector <T> (m_variables);
		m_jobs[i].result_vector = new std::vector<T*>;
	    }
            //std::cout << "enum_first (roots[" << norms.first << "])" << std::endl;
	    std::vector<T*> result;
	    // The following call starts a lot of jobs:
            enum_first (m_roots[norms.first], norms, result);
            //std::cout << "enum_first finished." << std::endl;

	    // All jobs are scheduled now.  Wait for all of them to finish.
	    for (int i = 0; i < NUMJOBS; i++){
		if (m_jobsFree[i])
		    continue;
		m_jobs[i].fut.wait();
		for (auto it = m_jobs[i].result_vector->begin(); it != m_jobs[i].result_vector->end(); it++){
		    result.push_back (*it);
		}
	    }

	    m_resultMapMutex.lock();
	    m_resultMap[norms] = std::move(result);
	    m_resultMapMutex.unlock();
	    for (int i = 0; i < NUMJOBS; i++){
		delete_vector<T> (m_jobs[i].temp_vectors.sum);
		delete m_jobs[i].result_vector;
	    }
        }
    }

public:
    ParallelPottier () { };

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
        m_lattice = generate_lattice (homo);
        delete homo;
	if (m_controller != NULL)
	{
	    //m_controller->save_lattice (m_lattice);
	    m_controller->log_lattice (m_lattice);
	}

        m_maxnorm = -1;
        m_current_variable = 0;
        m_variables = m_lattice->variables ();
    }

    void init (Lattice <T> * lattice, Controller <T> * controller)
    {
        m_controller = controller;
        m_lattice = new Lattice <T> (*lattice);

        m_controller->log_lattice (m_lattice);

        m_maxnorm = -1;
        m_current_variable = 0;
        m_variables = m_lattice->variables ();
    }

    void init (std::ifstream& stream, Controller <T> * controller)
    {
        m_controller = controller;
        controller->read_backup (stream);
	// Backup functionality is broken in parallel mode!
        stream >> m_current_variable;
        int vectors;
        stream >> vectors >> m_variables;

        m_maxnorm = -1;

        VariableProperties <T> * properties = new VariableProperties <T> (m_variables, false, 0, 0);
        for (size_t i = 0; i < m_variables; i++)
        {
            int column;
            bool free;
            T lower, upper;
            stream >>  column >> free >> lower >> upper;
            VariableProperty <T> & var = properties->get_variable (i);
            var.set (column, free, lower, upper);
        }

        m_lattice = new Lattice <T> (properties);

        delete properties;

        for (int i = 0; i < vectors; i++)
        {
            T* vector = read_vector <T> (stream, m_variables);
            m_lattice->append_vector (vector);
        }

	// The two zero arguments represent the current sum and current first norm that 
	// are not applicable for backups of the parallel algorithm.
        m_controller->log_resume (m_variables, m_current_variable+1, 0,0, vectors);
    }

    ~ParallelPottier ()
    {
        delete m_lattice;
    }


    /**
     * What does compute do?
     *
     * At time of start, some things have been initialized:
     *    m_maxnorm = -1;
     *    m_current_variable = 0;
     *    m_variables = m_lattice->variables ();
     *     
     * The linear system and a controller have been set up.
     * 
     */
    void compute (int backup_frequency = 0)
    {
        m_norms.clear ();

        m_backup_frequency = backup_frequency;
        m_backup_timer.reset ();

        // m_sum_vector = create_vector <T> (m_variables);

	//std::cout << "\n\nSTARTED ON " << m_current_variable << ", " << m_sum_norm << ", " << m_first_norm << std::endl;

	// This for loop iterates of the variables of the system.  It
	// assumes that the system is already solved on the variables
	// that are done and solves it for the next variable.
        for (; m_current_variable < m_variables; m_current_variable++)
        {
	    //std::cout << "\n\n\nLattice before choose:\n" << (*m_lattice) << "\nvectors = " << m_lattice->vectors () << std::endl;
	    int next = chooseNextVariable ();
	    if (next < 0)
		break;
	    
	    //std::cout << "\n-------------------------------------------" << std::endl;
	    //std::cout << "\nChoice for next component is: " << next << std::endl;
	    
	    if (m_controller != NULL)
		m_controller->log_variable_start (m_current_variable+1, m_lattice->vectors ());
	    
	    m_lattice->swap_columns (m_current_variable, next);
                
	    //std::cout << "Lattice after choosing new component:\n" << (*m_lattice) << "\nMaxnorm = " << m_maxnorm << std::endl;
	    preprocess ();
	    
	    create_trees ();

            // create norms
            m_norms.clear ();
            for (auto i = m_roots.begin (); i != m_roots.end (); i++)
            {
                for (auto j = m_roots.begin (); j != m_roots.end (); j++)
                {
                    NormPair <T> pair (i->first, j->first);
                    m_norms[pair] = true;
                }
            }

	    std::cout << "Norm pairs to start with: ";
	    for (auto it = m_norms.begin(); it != m_norms.end(); it++) {
		std::cout << "(" << it->first.first << "," << it->first.second << "), ";
	    }
	    std::cout << std::endl;

	    m_maxnorm = m_norms.rbegin()->first.sum;
 	    std::cout << "Maximum norm on first components:" << m_maxnorm << "\n";

            // norm pairs
            // T old_sum = -1;

 	    T current_norm = 0;
 	    // The big job-loop
 	    while (current_norm < m_maxnorm) {
 		current_norm++;
 		std::cout << "Now doing norm "<< current_norm << " of " << m_maxnorm << " for variable " << m_current_variable+1 << " of " << m_variables << std::endl;
 		std::vector < std::future <void> > m_futures;

		for (auto it = m_norms.cbegin(); it != m_norms.cend(); it++){
		    if (it->first.sum == current_norm) {
			std::cout << "Starting job :" << it->first.first << "," << it->first.second << "\n";
			complete (it->first);
//  			std::future<void> fut = std::async(
//  			    std::launch::async, // This directive makes it launch a new thread for each job (not good if there are many!)
//  			    &ParallelPottier<T>::complete,
//			    this,
//  			    it->first);
//			// This parallelization must be disabled until the enum_second_job pointer array is norm_job specific.
//			fut.wait();
//			m_futures.push_back (std::move(fut));
		    }
		}

		// Synchronize:
		std::cout << "Waiting for sync ... ";
		std::cout.flush();
		for (auto it = m_futures.begin(); it != m_futures.end(); ++it)
		    it->wait();
		std::cout << "Done.\n";

		if (m_controller != NULL)
		    m_controller->log_status (m_current_variable+1, current_norm, m_maxnorm, 0, m_lattice->vectors (), m_backup_frequency, m_backup_timer);

		// collecting results
		UniqueVectorsHash<T> *unique_res = new UniqueVectorsHash<T> (m_variables);
		T *neg = create_vector<T> (m_variables);
		for (auto it = m_norms.cbegin(); it != m_norms.cend(); it++) {
		    if (it->first.sum == current_norm) {
			if (m_resultMap[it->first].size() > 0){
			    std::cout << "Inserting: " << m_resultMap[it->first].size() << " results.\n";
			    m_maxnorm = it->first.sum * 2;
			}
			for (auto jt = m_resultMap[it->first].begin(); jt != m_resultMap[it->first].end(); jt++ ){
			    if (!unique_res->is_present(*jt)) {
				// jt need not be deleted in this branch, can be reused
				for (size_t i = 0; i < m_variables; i++){
				    neg[i] = -(*jt)[i];
				}
				// This copies *jt
				insert_trees(*jt, it->first.sum);
				// This pilfers *jt.
				unique_res->insert(*jt);
				insert_trees(neg, it->first.sum);
				unique_res->insert_copy(neg);
			    }
			    else {
				delete_vector (*jt);
			    }
			}
		    }
		}
		delete_vector (neg);
		delete unique_res;
	    } // looping the current norm
		
            delete_trees ();
	    if (m_controller != NULL)
		m_controller->log_variable_end (m_current_variable+1, m_lattice->vectors ());
        }

        //std::cout << "Finished with " << m_lattice->vectors () << " solutions." << std::endl;
        m_lattice->sort_columns ();
    }

    Lattice <T>& lattice () const
    {
        return *m_lattice;
    }

    void extract_zsolve_results (VectorArray <T>& inhoms, VectorArray <T>& homs, VectorArray <T>& free)
    {
	std::cout << "Error.  I am not a zsolve implementation!\n";
	exit (1);
    }

    void extract_graver_results (VectorArray <T>& graver)
    {
        assert (m_lattice->get_splitter () == -1);
        assert (m_lattice->get_result_num_variables () == m_variables);

        graver.clear ();

        for (size_t i = 0; i < m_lattice->vectors (); i++)
        {
            T* vector = (*m_lattice)[i];
            T* result = copy_vector <T> (vector, m_variables);

            bool has_symmetric = true;
            for (size_t j = 0; j < m_variables; j++)
                if (!m_lattice->get_variable (j).check_bounds (-vector[j]))
                    has_symmetric = false;
            
            int lex_cmp = lex_compare_vector_with_negative (vector, m_variables);

            if (!has_symmetric || lex_cmp > 0)
                graver.append_vector (result);
        }

	    if (m_controller != NULL)
	        m_controller->log_result (1, m_lattice->vectors (), 0);
    }
    
    void extract_hilbert_results (VectorArray <T>& hilbert)
    {
	std::cout << "Error, I am not a hilbert implementation!\n";
	exit(1);
    }

    T extract_maxnorm_results (VectorArray <T> & maxnorm)
    {
        int result_variables = m_lattice->get_result_num_variables ();

        maxnorm.clear ();

        m_maxnorm = -1;
        for (size_t i = 0; i < m_lattice->vectors (); i++)
        {
            T* vector = (*m_lattice)[i];
            T norm = norm_vector (vector, result_variables);
            if (norm > m_maxnorm)
            {
                m_maxnorm = norm;
                maxnorm.clear ();
            }
            if (norm == m_maxnorm)
                maxnorm.append_vector (copy_vector <T> (vector, result_variables));
        }
        return m_maxnorm;
    }

    size_t get_result_num_variables () const
    {
        return m_lattice->get_result_num_variables ();
    }

    void log_maxnorm ()
    {
	if (m_controller != NULL)
	    m_controller->log_maxnorm (this, true);
    }

    //int count_builds;
    //int count_reduces;
    //int count_insertions;
    //int complete_calls;
};

} // namespace _4ti2_zsolve_
  
#endif
