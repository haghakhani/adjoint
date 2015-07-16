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

#define KEY0 3777862041
#define KEY1 2576980374
#define DEBUG

void bilinear_interp(HashTable* El_Table) {

	HashEntryPtr currentPtr;
	Element *Curr_El, *father, *elem11, *elem12, *elem21, *elem22;
	HashEntryPtr *buck = El_Table->getbucketptr();
	int which_son;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() == NEWSON) {
					father = (Element*) El_Table->lookup(Curr_El->getfather());
					assert(father->get_adapted_flag()==OLDFATHER);

					which_son = Curr_El->get_which_son();
#ifdef DEBUG
					double aa, bb = .1;
					if (*(Curr_El->pass_key()) == KEY0 && *(Curr_El->pass_key() + 1) == KEY1)
						aa = bb;
#endif
					switch (which_son) {

						case 0: {
							elem22 = father;
							elem12 = father->get_side_neighbor(El_Table, 2);
							elem21 = father->get_side_neighbor(El_Table, 3);
							elem11 = father->get_side_neighbor(El_Table, 6);
							bilinear_interp_elem(elem11, elem21, elem12, elem22, Curr_El);
						}
							break;
						case 1: {
							elem22 = father->get_side_neighbor(El_Table, 0);
							elem12 = father;
							elem21 = father->get_side_neighbor(El_Table, 7);
							elem11 = father->get_side_neighbor(El_Table, 3);
							bilinear_interp_elem(elem11, elem21, elem12, elem22, Curr_El);

						}
							break;
						case 2: {
							elem22 = father->get_side_neighbor(El_Table, 4);
							elem12 = father->get_side_neighbor(El_Table, 1);
							elem21 = father->get_side_neighbor(El_Table, 0);
							elem11 = father;
							bilinear_interp_elem(elem11, elem21, elem12, elem22, Curr_El);

						}
							break;
						case 3: {
							elem22 = father->get_side_neighbor(El_Table, 1);
							elem12 = father->get_side_neighbor(El_Table, 5);
							elem21 = father;
							elem11 = father->get_side_neighbor(El_Table, 2);
							bilinear_interp_elem(elem11, elem21, elem12, elem22, Curr_El);

						}
							break;
						default:
							cout << "incorrect son please check me" << endl;

					}

				}
				currentPtr = currentPtr->next;
			}
		}
	return;
}

double bilinear_interp_value(double x1, double x2, double y1, double y2, double f11, double f21,
    double f12, double f22, double xinterp, double yinterp, int type) {

//	xmin = x1
//	xmax = x2
//	ymin = y1
//	ymax = y2
//
//	f12--------f22
//	|           |
//	|           |
//	|  interp   |
//	|           |
//	|           |
//	f11---------f21

	double interp = 0.0;
	switch (type) {
		case 0:			//for bilinear interpolation
			interp = (f11 * (x2 - xinterp) * (y2 - yinterp) + f21 * (xinterp - x1) * (y2 - yinterp)
			    + f12 * (x2 - xinterp) * (yinterp - y1) + f22 * (xinterp - x1) * (yinterp - y1))
			    / ((x2 - x1) * (y2 - y1));
			break;
		case 1:
			//interpolation in y
			//no other modification on bilinear interpolation is required,
			//since we call the function such that put zero for those elements that do not exist
			//we just removed (x2 - x1) from denominator
			interp = (f11 * (x2 - xinterp) * (y2 - yinterp) + f21 * (xinterp - x1) * (y2 - yinterp)
			    + f12 * (x2 - xinterp) * (yinterp - y1) + f22 * (xinterp - x1) * (yinterp - y1))
			    / (y2 - y1);
			break;
		case 2:
			//interpolation in x
			//no other modification on bilinear interpolation is required,
			//since we call the function such that put zero for those elements that do not exist
			//we just removed (y2 - y1) from denominator
			interp = (f11 * (x2 - xinterp) * (y2 - yinterp) + f21 * (xinterp - x1) * (y2 - yinterp)
			    + f12 * (x2 - xinterp) * (yinterp - y1) + f22 * (xinterp - x1) * (yinterp - y1))
			    / (x2 - x1);
			break;

		default:
			cout << "not a valid type in interp_value function " << endl;
	}
	if (isnan(interp) || isinf(interp))
		cout << "it is so sad that I found you" << endl;
	return (interp);
}

void bilinear_interp_elem(Element *elem11, Element *elem21, Element *elem12, Element *elem22,
    Element *Curr_El) {

	double *state_vars, *elem11_state, *elem12_state, *elem21_state, *elem22_state;
	double *prev_state_vars, *elem11_prev_state, *elem12_prev_state, *elem21_prev_state,
	    *elem22_prev_state;
	double *adjoint, *elem11_adjoint, *elem12_adjoint, *elem21_adjoint, *elem22_adjoint;
	double *coord, *elem11_coord, *elem12_coord, *elem21_coord, *elem22_coord;

	int type = 0;			//this is just a flag that indicates the type of element

	state_vars = Curr_El->get_state_vars();
	prev_state_vars = Curr_El->get_prev_state_vars();
	coord = Curr_El->get_coord();
	adjoint = Curr_El->get_adjoint();

	elem11_state = elem11->get_state_vars();
	elem12_state = elem12->get_state_vars();
	elem21_state = elem21->get_state_vars();
	elem22_state = elem22->get_state_vars();

	elem11_prev_state = elem11->get_prev_state_vars();
	elem12_prev_state = elem12->get_prev_state_vars();
	elem21_prev_state = elem21->get_prev_state_vars();
	elem22_prev_state = elem22->get_prev_state_vars();

	elem11_coord = elem11->get_coord();
	elem12_coord = elem12->get_coord();
	elem21_coord = elem21->get_coord();
	elem22_coord = elem22->get_coord();

	elem11_adjoint = elem11->get_adjoint();
	elem12_adjoint = elem12->get_adjoint();
	elem21_adjoint = elem21->get_adjoint();
	elem22_adjoint = elem22->get_adjoint();

	if (elem11 && elem12 && elem21 && elem22) {
		//this is an ordinary case for an element inside the domain
		//type = 0; we initialized type=0
		double aaa, bbb = .1;
		for (int j = 0; j < NUM_STATE_VARS; j++)
			if ( isnan(
			    elem11_state[j]) || isnan(elem12_state[j]) || isnan(elem21_state[j]) || isnan(elem22_state[j]) ||
			    isinf(elem11_state[j]) || isinf(elem12_state[j]) || isinf(elem21_state[j]) || isinf(elem22_state[j]))
				aaa = bbb;

		for (int j = 0; j < NUM_STATE_VARS; j++) {
			state_vars[j] = bilinear_interp_value(elem11_coord[0], elem21_coord[0], elem21_coord[1],
			    elem22_coord[1], elem11_state[j], elem21_state[j], elem12_state[j], elem22_state[j],
			    coord[0], coord[1], type);

			prev_state_vars[j] = bilinear_interp_value(elem11_coord[0], elem21_coord[0], elem21_coord[1],
			    elem22_coord[1], elem11_prev_state[j], elem21_prev_state[j], elem12_prev_state[j],
			    elem22_prev_state[j], coord[0], coord[1], type);

			adjoint[j] = bilinear_interp_value(elem11_coord[0], elem21_coord[0], elem21_coord[1],
			    elem22_coord[1], elem11_adjoint[j], elem21_adjoint[j], elem12_adjoint[j],
			    elem22_adjoint[j], coord[0], coord[1], type);
		}
	} else if ((!elem11 && !elem12 && elem21 && elem22)	//interpolation only in y
	|| (elem11 && elem12 && !elem21 && !elem22)	//left or right side of father is boundary
	    ) {
		type = 1;
		if (elem11) {	// in this case elem21 & elem22 do not exist, so we replace their value with zero
			for (int j = 0; j < NUM_STATE_VARS; j++) {
				state_vars[j] = bilinear_interp_value(0,
				    0,	//interpolation is in y, so x position is not important
				    elem11_coord[1], elem12_coord[1], elem11_state[j], 0, elem12_state[j], 0, coord[0],
				    coord[1], type);

				prev_state_vars[j] = bilinear_interp_value(0,
				    0,	//interpolation is in y, so x position is not important
				    elem11_coord[1], elem12_coord[1], elem11_prev_state[j], 0, elem12_prev_state[j], 0,
				    coord[0], coord[1], type);

				adjoint[j] = bilinear_interp_value(0,
				    0,	//interpolation is in y, so x position is not important
				    elem11_coord[1], elem12_coord[1], elem11_adjoint[j], 0, elem12_adjoint[j], 0, coord[0],
				    coord[1], type);
			}

		} else {	// in this case elem11 & elem12 do not exist, so we replace their value with zero

			for (int j = 0; j < NUM_STATE_VARS; j++) {
				state_vars[j] = bilinear_interp_value(0,
				    0,	//interpolation is in y, so x position is not important
				    elem21_coord[1], elem22_coord[1], 0, elem21_state[j], 0, elem22_state[j], coord[0],
				    coord[1], type);

				prev_state_vars[j] = bilinear_interp_value(0,
				    0,	//interpolation is in y, so x position is not important
				    elem21_coord[1], elem22_coord[1], 0, elem21_prev_state[j], 0, elem22_prev_state[j],
				    coord[0], coord[1], type);

				adjoint[j] = bilinear_interp_value(0,
				    0,	//interpolation is in y, so x position is not important
				    elem21_coord[1], elem22_coord[1], 0, elem21_adjoint[j], 0, elem22_adjoint[j], coord[0],
				    coord[1], type);
			}
		}
	} else if ((!elem11 && elem12 && !elem21 && elem22)	//interpolation only in x
	|| (elem11 && !elem12 && elem21 && !elem22)	//top or bottom side of father is boundary
	    ) {
		type = 2;

		if (elem11) {	// in this case elem12 & elem22 do not exist, so we replace their value with zero
			for (int j = 0; j < NUM_STATE_VARS; j++) {
				state_vars[j] = bilinear_interp_value(elem11_coord[0], elem21_coord[0], 0, 0,	//interpolation is in x, so y position is not important
				    elem11_state[j], elem21_state[j], 0, 0, coord[0], coord[1], type);

				prev_state_vars[j] = bilinear_interp_value(elem11_coord[0], elem21_coord[0], 0, 0,//interpolation is in x, so y position is not important
				    elem11_prev_state[j], elem21_prev_state[j], 0, 0, coord[0], coord[1], type);

				adjoint[j] = bilinear_interp_value(elem11_coord[0], elem21_coord[0], 0, 0,//interpolation is in x, so y position is not important
				    elem11_adjoint[j], elem21_adjoint[j], 0, 0, coord[0], coord[1], type);
			}

		} else {	// in this case elem11 & elem21 do not exist, so we replace their value with zero

			for (int j = 0; j < NUM_STATE_VARS; j++) {
				state_vars[j] = bilinear_interp_value(elem12_coord[0], elem22_coord[0], 0, 0,	//interpolation is in x, so y position is not important
				    0, 0, elem12_state[j], elem22_state[j], coord[0], coord[1], type);

				prev_state_vars[j] = bilinear_interp_value(elem12_coord[0], elem22_coord[0], 0, 0,//interpolation is in x, so y position is not important
				    0, 0, elem12_prev_state[j], elem22_prev_state[j], coord[0], coord[1], type);

				adjoint[j] = bilinear_interp_value(elem12_coord[0], elem22_coord[0], 0, 0,//interpolation is in x, so y position is not important
				    0, 0, elem12_adjoint[j], elem22_adjoint[j], coord[0], coord[1], type);
			}
		}

	} else if ((elem11 && !elem12 && !elem21 && !elem22)//father is in corner so there is no element for interp
	|| (!elem11 && elem12 && !elem21 && !elem22)	//we do not do any extrapolation and leave as it is,
	    || (!elem11 && !elem12 && elem21 && !elem22)//which in refinement constructor should be the value of father element
	    || (!elem11 && !elem12 && !elem21 && elem22)) {
		//do not do anything

	} else {
		cout << "something is wrong in this configuration" << endl;
		*(Curr_El->get_state_vars()) = 50;
		return;
	}

	return;

}
