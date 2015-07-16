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
 * Author: Hossein Aghakhani
 * Date: Jan 9 2015
 *******************************************************************
 */
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "../header/hpfem.h"

#define	DEBUG

void residual(double* residual, double *state_vars, double *prev_state_vars, //3
    double *fluxxp, double *fluxyp, double *fluxxm, double *fluxym, double dtdx, //5
    double dtdy, double dt, double *d_state_vars_x, double *d_state_vars_y, //4
    double *curvature, double intfrictang, double bedfrict, double *gravity, //4
    double *dgdx, double kactxyelem, double fric_tiny, double* orgSrcSgn, //4
    double increment, double epsilon, int* check_stop_crit, int srcflag, int org_res_flag) {

	double velocity[DIMENSION];
	double kactxy[DIMENSION];
	double tmp[NUM_STATE_VARS];
	//double bedfrict;

	if (prev_state_vars[0] > GEOFLOW_TINY) {
		for (int k = 0; k < DIMENSION; k++)
			kactxy[k] = kactxyelem;

		if ((prev_state_vars[1] == 0. && prev_state_vars[2] == increment)
		    || (prev_state_vars[2] == 0. && prev_state_vars[1] == increment)) {
			velocity[0] = 0.;
			velocity[1] = 0.;
		} else {

			// fluid velocities
			velocity[0] = prev_state_vars[1] / prev_state_vars[0];
			velocity[1] = prev_state_vars[2] / prev_state_vars[0];
		}

	} else {
		for (int k = 0; k < DIMENSION; k++) {
			kactxy[k] = epsilon;
			velocity[k] = 0.;
		}
		//bedfrict = bedfrictin;
	}

	for (int i = 0; i < NUM_STATE_VARS; i++) {
		residual[i] = 0.0;
		tmp[i] = 0.;
	}

	for (int i = 0; i < NUM_STATE_VARS; i++)
		tmp[i] = prev_state_vars[i] - dtdx * (fluxxp[i] - fluxxm[i]) - dtdy * (fluxyp[i] - fluxym[i]);

	if (prev_state_vars[0] > GEOFLOW_TINY && srcflag) {

		double unitvx = 0., unitvy = 0., h_inv = 0., speed = 0.;

		speed = sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1]);

		if (speed > 0.) {
			unitvx = velocity[0] / speed;
			unitvy = velocity[1] / speed;
		}

		//x dir
		double s1 = gravity[0] * prev_state_vars[0];

		double sin_int_fric = sin(intfrictang);
		double s2 = orgSrcSgn[0] * prev_state_vars[0] * kactxy[0]
		    * (gravity[2] * d_state_vars_y[0] + dgdx[1] * prev_state_vars[0]) * sin_int_fric;

		double tan_bed_fric = tan(bedfrict);
		double s3 = unitvx
		    * max(gravity[2] * prev_state_vars[0] + velocity[0] * prev_state_vars[1] * curvature[0],
		        0.0) * tan_bed_fric;
		if (prev_state_vars[1]==increment)
			s3=0;

			if (dabs(tmp[1] + dt * s1) > dabs(dt * (s2 + s3)) && !check_stop_crit[0])
				tmp[1] += dt * (s1 - s2 - s3);
			else {
				tmp[1] = 0.;
				if (org_res_flag)
					check_stop_crit[0] = 1;
			}

		//y dir

		s1 = gravity[1] * prev_state_vars[0];

		s2 = orgSrcSgn[1] * prev_state_vars[0] * kactxy[0]
		    * (gravity[2] * d_state_vars_x[0] + dgdx[0] * prev_state_vars[0]) * sin_int_fric;

		s3 = unitvy
		    * max(gravity[2] * prev_state_vars[0] + velocity[1] * prev_state_vars[2] * curvature[1],
		        0.0) * tan_bed_fric;

		if (prev_state_vars[2]==increment)
			s3=0;

		if (dabs(tmp[2] + dt * s1) > dabs(dt * (s2 + s3)) && !check_stop_crit[1])
			tmp[2] += dt * (s1 - s2 - s3);
		else {
			tmp[2] = 0.;
			if (org_res_flag)
				check_stop_crit[1] = 1;
		}

	}

	if (org_res_flag)
		for (int i = 0; i < NUM_STATE_VARS; i++)
			state_vars[i] = tmp[i];


	for (int i = 0; i < NUM_STATE_VARS; i++)
		residual[i] = state_vars[i] - tmp[i];

#ifdef DEBUG

//	for (int k = 0; k < 3; k++)
//		if (residual[k] > 1e-5) {
//			cout << "something that has to be checked" << endl << flush;
//			exit(-2);
//		}

	for (int k = 0; k < 3; k++)
		if (isnan(residual[k])) {
			cout << "exit for NAN in residual" << endl << flush;
			exit(-1);
		}

	for (int k = 0; k < 3; k++)
		if (isinf(residual[k])) {
			cout << "exit for Inf in residual" << endl << flush;
			exit(-2);
		}
#endif

	return;
}

