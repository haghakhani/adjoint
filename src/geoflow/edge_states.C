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
 * Author: 
 * Description: 
 *
 *******************************************************************
 * $Id: edge_states.C 128 2007-06-07 19:51:52Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"

#define KEY0   310842289
#define KEY1   991146299
#define ITER   10

/*! calc_edge_states() cycles through the element Hashtable (listing of all 
 *  elements) and for each element (that has not been refined this iteration 
 *  and is not a ghost_element) calls Element member function 
 *  Element::calc_edge_states() (which calculates the Riemann fluxes across 
 *  the element's boundaries), and adds local boundary-element outflow to 
 *  GIS map's cummulative outflow (defined as the mass flow off of the
 *  GIS map).  Also, the elements are checked for multiple pile-height values 
 */
void calc_edge_states(HashTable* El_Table, HashTable* NodeTable, MatProps* matprops_ptr,
    TimeProps* timeprops_ptr, int myid, int* order_flag, double *outflow, ResFlag resflag) {
	int i, j, k, counter, iter, keys1, keys2;
	double tiny = GEOFLOW_TINY;
	int el_counter = 0;
	double evalue = 1;

	MapNames mapnames;
	char *b, *c, *d;
	char a[5] = "abs"; // ,b[5],c[5],d[5];
	b = c = d = a;
	int ce = 0;
	mapnames.assign(a, b, c, d, ce);

	//-------------------go through all the elements of the subdomain and
	//-------------------find the edge states

	HashEntryPtr* buck = El_Table->getbucketptr();
	double localoutflow;
	*outflow = 0.0;
	/* mdj 2007-04 */
	double localoutflow_sum = 0.0;
	HashEntryPtr currentPtr;
	Element* Curr_El;
#pragma omp parallel for private(currentPtr,Curr_El) reduction(+:localoutflow_sum)
	for (i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {
					//if this element doesn't belong on this processor don't involve

					if (*(Curr_El->pass_key()) == KEY0 && *(Curr_El->pass_key() + 1) == KEY1
					    && timeprops_ptr->iter == ITER) {
						int ddd, aa = 0;
						int gg = ddd;
					}

					Curr_El->calc_edge_states(El_Table, NodeTable, matprops_ptr, myid, timeprops_ptr->dtime,
					    order_flag, &localoutflow, resflag, resflag);
					localoutflow_sum += localoutflow;
				}
				currentPtr = currentPtr->next;
			}
		}
	*outflow = localoutflow_sum;
	return;
}

void calc_flux(MeshCTX* meshctx, PropCTX* propctx, int myid, ResFlag resflag) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;

	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	Element* Curr_El;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {
					//if this element doesn't belong on this processor don't involve

					if (*(Curr_El->pass_key()) == KEY0 && *(Curr_El->pass_key() + 1) == KEY1
					    && timeprops_ptr->iter == ITER) {
						int ddd, aa = 0;
						int gg = ddd;
					}

					Curr_El->calc_fluxes(El_Table, NodeTable, myid, resflag, resflag);

				}
				currentPtr = currentPtr->next;
			}
		}

	return;
}
