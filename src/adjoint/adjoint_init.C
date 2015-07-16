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
 * $Id: calc_jacobian.C 164 2014-05-22 15:27:22Z haghakhani $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"

void adjoint_init(HashTable *El_Table, HashTable *Node_Table) {

	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	Element* Curr_El = NULL;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0 && *(Curr_El->get_state_vars()) > 0) {
					*(Curr_El->get_state_vars() + 6) = *(Curr_El->get_dx()) * *(Curr_El->get_dx() + 1);
					if (isinf(*(Curr_El->get_state_vars() + 6)))
						cout << "there is something wrong in initialization" << endl;

					*(Curr_El->get_state_vars() + 7) = 0.;
					*(Curr_El->get_state_vars() + 8) = 0.;
				}
				currentPtr = currentPtr->next;
			}
		}
	}

//  for(int i=0; i<El_Table->get_no_of_buckets(); i++){
//    if(*(buck+i)) {
//      currentPtr = *(buck+i);
//      while(currentPtr){
//	Curr_El=(Element*)(currentPtr->value);
////	  if (isinf(*(Curr_El->get_state_vars()+6))) 
//	    cout<<"there is something wrong in initialization   "<<   *(Curr_El->get_state_vars()+6)<<endl;
//
//	currentPtr=currentPtr->next;
//      }
//    }
//  }
//

	return;
}

void compute_funcsens(Element* element, TimeProps* timep, double* sensitivity_curr,
    double* sensitivity_prev) {

	double *pos = element->get_dx();
	double dt = timep->dt.back();

	double height_curr = *(element->get_state_vars());
	double height_prev = *(element->get_prev_state_vars());

	// in trapezoidal time integrator that we use for functional,
	// we have the effect of current and previous time step

	*(sensitivity_curr) = 0.5 * dt * pos[0] * pos[1] * height_curr;
	*(sensitivity_curr + 1) = 0.;
	*(sensitivity_curr + 2) = 0.;

	*(sensitivity_prev) = 0.5 * dt * pos[0] * pos[1] * height_prev;
	*(sensitivity_prev + 1) = 0.;
	*(sensitivity_prev + 2) = 0.;

	return;
}

void compute_funcsens(Element* element, double dt, double* func_sens) {

	double *dx = element->get_dx();

	double height = *(element->get_state_vars());

	// in trapezoidal time integrator that we use for functional,
	// we have the effect of current and previous time step

	func_sens[0] = 0.5 * dt * dx[0] * dx[1] * height;
	func_sens[1] = 0.;
	func_sens[2] = 0.;

//	func_sens->sensitivity = sensitivity;

	return;
}
