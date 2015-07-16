/*
 * adjoint_test.C
 *
 *  Created on: May 22, 2015
 *      Author: haghakha
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "../header/hpfem.h"

void perturbU(HashTable* El_Table, PertElemInfo* pelinf, int iter) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	const double epsilon = INCREMENT;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {
					if (dabs(*(Curr_El->get_coord()) - pelinf->elempos[0]) < epsilon
					    && dabs(*(Curr_El->get_coord() + 1) - pelinf->elempos[1]) < epsilon) {

						// perturbing U
						*(Curr_El->get_state_vars()) += epsilon;

						// to break the for loop
						i = El_Table->get_no_of_buckets();

						// break the while loop
						break;
					}

				}
				currentPtr = currentPtr->next;
			}
		}

	return;
}

void fill_pertelem_info(HashTable* El_Table, PertElemInfo* eleminfo) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();
	Jacobian *jacobian, *neighjac;
	Element *neigh_elem;
	double*** jacobianmat;
	double *curr_adj_ptr, *prev_adj_ptr;

	const double epsilon = INCREMENT;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {
					if (dabs(*(Curr_El->get_coord()) - eleminfo->elempos[0]) < epsilon
					    && dabs(*(Curr_El->get_coord() + 1) - eleminfo->elempos[1]) < epsilon) {

						eleminfo->h = *(Curr_El->get_state_vars());

						eleminfo->elem_size[0] = *(Curr_El->get_dx());

						eleminfo->elem_size[1] = *(Curr_El->get_dx() + 1);

//						eleminfo->

						jacobianmat = Curr_El->get_jacobian();
						eleminfo->func_sens = *(Curr_El->get_func_sens());

						for (int effelement = 0; effelement < 5; effelement++) { //0 for the element itself, and the rest id for neighbour elements

							if (effelement == 0) {		    //this part of code

								curr_adj_ptr = (Curr_El->get_adjoint());
								prev_adj_ptr = (Curr_El->get_prev_adjoint());

								for (int j = 0; j < 3; ++j) {
									eleminfo->neigh_jac[effelement].curr_adj[j] = curr_adj_ptr[j];
									eleminfo->neigh_jac[effelement].prev_adj[j] = prev_adj_ptr[j];
								}

								jacobianmat = Curr_El->get_jacobian();

								for (int k = 0; k < 3; ++k)
									for (int l = 0; l < 3; ++l)
										eleminfo->neigh_jac[effelement].jacobianMat[k][l] =
										    jacobianmat[effelement][k][l];

							} else {

								neigh_elem = Curr_El->get_side_neighbor(El_Table, effelement - 1);//basically we are checking all neighbor elements, and start from xp neighbor
								if (neigh_elem) {

									curr_adj_ptr = (neigh_elem->get_state_vars() + 6);
									prev_adj_ptr = (neigh_elem->get_prev_state_vars() + 6);

									for (int j = 0; j < 3; ++j) {
										eleminfo->neigh_jac[effelement].curr_adj[j] = curr_adj_ptr[j];
										eleminfo->neigh_jac[effelement].prev_adj[j] = prev_adj_ptr[j];
									}

									jacobianmat = neigh_elem->get_jacobian();

									int jacind;

									switch (effelement) {
										case 1:	//in xp neighbor I have to read jacobian of xm, because position of curr_el is in xm side of that neighbor
											jacind = 3;
											break;
										case 2:		    //for yp return ym
											jacind = 4;
											break;
										case 3:		    //for xm return xp
											jacind = 1;
											break;
										case 4:		    //for ym return yp
											jacind = 2;
											break;
										default:
											cout << "invalid neighbor position" << endl;
									}

									for (int k = 0; k < 3; ++k)
										for (int l = 0; l < 3; ++l)
											eleminfo->neigh_jac[effelement].jacobianMat[k][l] = jacobianmat[jacind][k][l];

								}
							}
						}

						// to break the for loop
						i = El_Table->get_no_of_buckets();

						// break the while loop
						break;
					}
				}
				currentPtr = currentPtr->next;
			}
		}

	return;
}
void find_test_elem(HashTable* El_Table, PertElemInfo** pelinf, int iter) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	double maxh = 0.;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {
					if (*(Curr_El->get_state_vars()) > maxh)
						maxh = *(Curr_El->get_state_vars());

				}
				currentPtr = currentPtr->next;
			}
		}

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {
					if (*(Curr_El->get_state_vars()) == maxh) {

						double* pos = Curr_El->get_coord();

						(*pelinf) = new PertElemInfo(pos, maxh, iter);

						// to break the for loop
						i = El_Table->get_no_of_buckets();

						// break the while loop
						break;
					}

				}
				currentPtr = currentPtr->next;
			}
		}

	return;
}

int check_elem_exist(HashTable *El_Table, unsigned *key) {
	int num = 0;			//myid
	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0)
					if (*(Curr_El->pass_key()) == key[0] && *(Curr_El->pass_key() + 1) == key[1])
						return (1);

				currentPtr = currentPtr->next;
			}
		}

	return (0);
}

int checkElement(HashTable *El_Table, HashTable *NodeTable, double *max, unsigned *key) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	int gg = 0;

	double fluxold[4][NUM_STATE_VARS];

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {
					if (*(Curr_El->pass_key()) == key[0] && *(Curr_El->pass_key() + 1) == key[1]) {
						int xp = Curr_El->get_positive_x_side(); //finding the direction of element
						int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;

						Node* nxp = (Node*) NodeTable->lookup(Curr_El->getNode() + (xp + 4) * 2);

						Node* nyp = (Node*) NodeTable->lookup(Curr_El->getNode() + (yp + 4) * 2);

						Node* nxm = (Node*) NodeTable->lookup(Curr_El->getNode() + (xm + 4) * 2);

						Node* nym = (Node*) NodeTable->lookup(Curr_El->getNode() + (ym + 4) * 2);

						for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
							fluxold[0][ivar] = nxp->flux[ivar];
							fluxold[1][ivar] = nyp->flux[ivar];
							fluxold[2][ivar] = nxm->flux[ivar];
							fluxold[3][ivar] = nym->flux[ivar];
						}
						gg = 1;
					}
				}
				currentPtr = currentPtr->next;
			}
		}

	cout << " \n flux xp: \n";
	for (int i = 0; i < NUM_STATE_VARS; ++i)
		cout << " " << fluxold[0][i];

	cout << " \n flux yp: \n";
	for (int i = 0; i < NUM_STATE_VARS; ++i)
		cout << " " << fluxold[1][i];

	cout << " \n flux xm: \n";
	for (int i = 0; i < NUM_STATE_VARS; ++i)
		cout << " " << fluxold[2][i];

	cout << " \n flux ym: \n";
	for (int i = 0; i < NUM_STATE_VARS; ++i)
		cout << " " << fluxold[3][i];

	return (gg);
}

