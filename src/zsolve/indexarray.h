/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.

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

#ifndef _INDEXARRAY_H
#define _INDEXARRAY_H

#include <stdio.h>

typedef struct indexarray_t
{
	int Size;
	int *Data;
} indexarray_t;
typedef indexarray_t *IndexArray;

IndexArray createIndexArray();
/* create an indexarray */

void deleteIndexArray(IndexArray);
/* delete an indexarray (array) */

void appendToIndexArray(IndexArray, int);
/* append an index at the end of the array (array, index) */

void fprintIndexArray(FILE *, IndexArray);
/* print indexarray to a stream (stream, array) */

#endif
