#pragma once

#include "Header.h"

/* Find search direction */
std::vector<double> findSearchDirection(std::vector<double> states, int method)
{
	std::vector<double> dir(2);
	std::vector<double> gradient(2);
	std::vector<std::vector<double>> hessianInverse;
	switch (method)
	{
	case STEEPEST_DESCENT:
		/* normalized */
		//dir[0] = df_dx1_at_x(states[0], states[1]) / abs(df_dx1_at_x(states[0], states[1]));
		//dir[1] = df_dx2_at_x(states[0], states[1]) / abs(df_dx2_at_x(states[0], states[1]));
		/* non normalized */
		dir[0] = df_dx1_at_x(states[0], states[1]);
		dir[1] = df_dx2_at_x(states[0], states[1]);
		break;
	case NEWTON:
		// search direction is given by (Hessian_f_at_x)^-1 * grad_f_at_x
		gradient[0] = df_dx1_at_x(states[0], states[1]);
		gradient[1] = df_dx2_at_x(states[0], states[1]);
		hessianInverse = invert_2by2(hessian_f_at_x(states));
		// Compute (Hessian_f_at_x)^-1 * grad_f_at_x
		dir[0] = hessianInverse[0][0] * gradient[0] + hessianInverse[0][1] * gradient[1];
		dir[1] = hessianInverse[1][0] * gradient[0] + hessianInverse[1][1] * gradient[1];
		/* normalize */
		//dir[0] = dir[0] / abs(dir[0]);
		//dir[1] = dir[1] / abs(dir[1]);
		break;
	default:
		std::cout << "Unknown method\n";
		break;
	}
	//multiply by -1 to convert gradient from ascent to descent
	dir[0] = -dir[0];
	dir[1] = -dir[1];

#if DEBUG_DIRECTION_STEP_BY_STEP
	std::cout << dir[0] << "\t" << dir[1] << std::endl;
#endif
	return dir;
}
