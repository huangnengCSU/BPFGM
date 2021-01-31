/*
    Kalign - a multiple sequence alignment program

    Copyright 2006, 2019 Timo Lassmann

    This file is part of kalign.

    Kalign is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef WEAVE_ALIGNMENT_H
#define WEAVE_ALIGNMENT_H

#include "global.h"
#include "aln_task.h"
#include "msa.h"

//extern int weave(struct msa* msa, int** map, int* tree);

extern int weave(struct msa* msa,struct aln_tasks*t);

extern int clean_aln(struct msa* msa);

#endif
