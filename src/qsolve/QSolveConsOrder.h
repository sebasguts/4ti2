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

#ifndef _4ti2_qsolve__QSolveConsOrder_
#define _4ti2_qsolve__QSolveConsOrder_

#include "qsolve/Cone.h"
#include "qsolve/VectorArray.h"

namespace _4ti2_
{

enum QSolveConsOrder { MAXINTER, MININDEX, MAXCUTOFF, MINCUTOFF};

class ConsOrder {
public:
    //enum Ordering { MAXINTER, MININDEX, MAXCUTOFF, MINCUTOFF };
    ConsOrder();
    ConsOrder(QSolveConsOrder o);

    void set_constraint_order(QSolveConsOrder o);
    QSolveConsOrder get_constraint_order() const;

    bool (*compare)(Index next_pos_count, Index next_neg_count, Index next_zero_count,
                    Index pos_count, Index neg_count, Index zero_count);
    bool (*circuit_compare)(Index next_pos_count, Index next_neg_count, Index next_zero_count,
                    Index pos_count, Index neg_count, Index zero_count);

private:
    QSolveConsOrder order;

    static bool maxinter(
                    Index next_pos_count, Index next_neg_count, Index next_zero_count,
                    Index pos_count, Index neg_count, Index zero_count);
    static bool maxcutoff(
                    Index next_pos_count, Index next_neg_count, Index next_zero_count,
                    Index pos_count, Index neg_count, Index zero_count);
    static bool mincutoff(
                    Index next_pos_count, Index next_neg_count, Index next_zero_count,
                    Index pos_count, Index neg_count, Index zero_count);
    static bool minindex(
                    Index next_pos_count, Index next_neg_count, Index next_zero_count,
                    Index pos_count, Index neg_count, Index zero_count);
};

} // namespace 4ti2

#endif
