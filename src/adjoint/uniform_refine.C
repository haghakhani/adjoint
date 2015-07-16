/*******************************************************************
 * Copyright (C) 2003 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 *******************************************************************
 */
//Jan 30, 2015
//haghakha
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "../header/hpfem.h"

#define KEY0   3777862041
#define KEY1   2576980374
//#define DEBUG

void uinform_refine(MeshCTX* meshctx, PropCTX* propctx, int numprocs, int myid) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;

	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	Element* Curr_El = NULL;
	int rescomp = 1;

	//for debugging perpose
	unsigned key[2] = { KEY0, KEY1 };
	double max=0;

#ifdef DEBUG
	double dummyv_star = 0.0;
	int adjflag = 1;
	tecplotter(El_Table, NodeTable, matprops_ptr, timeprops_ptr, mapname_ptr, 0.,adjflag);

	int nonz1 = num_nonzero_elem(El_Table);

	cout << "number of elements before refinement  " << nonz1 << endl;

	int *dbgvec = new int[nonz1];
	int *pass = new int[nonz1];
	for (int i = 0; i < nonz1; i++) {
		dbgvec[i] = 0;
		pass[i] = 0;
	}

#endif

	htflush(El_Table, NodeTable, 1);
	move_data(numprocs, myid, El_Table, NodeTable, timeprops_ptr);


	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() >= NOTRECADAPTED) {
					Curr_El->put_adapted_flag(NOTRECADAPTED);
				}
				currentPtr = currentPtr->next;
			}
		}
	}

	ElemPtrList RefinedList(num_nonzero_elem(El_Table));

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				currentPtr = currentPtr->next;
				if (Curr_El->get_adapted_flag() == NOTRECADAPTED) {

					refinewrapper(El_Table, NodeTable, matprops_ptr, &RefinedList,
							Curr_El, rescomp);
				}

			}
		}
	}

//	cout << "here is the problem   "<<checkElement(El_Table, &max, key) << endl;

#ifdef DEBUG

	cout << "number of elements -7   " << num_nonzero_elem(El_Table, -7) << endl
	<< "number of elements -6   " << num_nonzero_elem(El_Table, -6) << endl
	<< "number of elements  0   " << num_nonzero_elem(El_Table, 0) << endl
	<< "number of elements  1   " << num_nonzero_elem(El_Table, 1) << endl
	<< "number of elements  2   " << num_nonzero_elem(El_Table, 2) << endl
	<< "number of elements  3   " << num_nonzero_elem(El_Table, 3) << endl
	<< "number of elements  4   " << num_nonzero_elem(El_Table, 4) << endl
	<< "number of elements  5   " << num_nonzero_elem(El_Table, 5) << endl;
#endif

	bilinear_interp(El_Table);	//this function reconstruct linear interpolation

	refine_neigh_update(El_Table, NodeTable, numprocs, myid, (void*) &RefinedList,
			timeprops_ptr);	//this function delete old father elements

//	cout << "here is the problem   "<<checkElement(El_Table, &max, key) << endl;

	RefinedList.trashlist();

	move_data(numprocs, myid, El_Table, NodeTable, timeprops_ptr);

	return;
}
