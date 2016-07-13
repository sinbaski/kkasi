#include <iostream>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <random>
#include <algorithm>
#include <armadillo>
#include <string.h>
#include <array>

using namespace std;
using namespace arma;
#define SAMPLE_SIZE 400000

random_device gen;
// mat<double> pool;

struct tail_fun_param
{
    double *coef;
    double *pool;
};

struct prod_tail_fun_param
{
    const double *A;
    const double *B;
    const double (*pool)[SAMPLE_SIZE];
};

double fun(double xi, void *param)
{
    struct tail_fun_param *pack = (struct tail_fun_param *)param;
    double *coef = pack->coef;
    double *pool = pack->pool;
    double r;
    r = accumulate(pool, pool + SAMPLE_SIZE, (double)0,
		   [xi, coef](double s, double x) {
		       return s + pow(coef[0] * x + coef[1], xi)/SAMPLE_SIZE;
		   }
	);
    return r - 1;
}

double prod_tail_fun(double xi, void *param)
{
    struct prod_tail_fun_param *pack = (struct prod_tail_fun_param *)param;
    const double *A = pack->A;
    const double *B = pack->B;
    const double (*pool)[SAMPLE_SIZE] = pack->pool;
    double r = 0;
#pragma omp parallel for    
    for (int i = 0; i < SAMPLE_SIZE; i++) {
	#pragma omp atomic
	r += pow(
	    (A[0] * pow(pool[0][i], 2) + A[1])*
	    (B[0] * pow(pool[1][i], 2) + B[1]),
	    xi) / SAMPLE_SIZE;
    }
    return r - 1;
    
}

double tail_index(double coef[2])
{
    gsl_root_fsolver *solver;
    gsl_function F;
    int iter = 0;
    int status = 0;
    int max_iter = 100;
    double lb, ub = 30, xi = -1;

    chi_squared_distribution<double> dist(1.0);
    static double pool[SAMPLE_SIZE];
    static bool flag = true;

    if (flag) {
#pragma omp parallel for
	for (int i = 0; i < SAMPLE_SIZE; i++) {
	    pool[i] = dist(gen);
	}
	flag = false;
    }

    double a;
    double bounds[2];
    struct tail_fun_param param = {coef, pool};
    for (a = 2; a > 0 && fun(a, &param) > 0; a -= 1);
    if (a > 0) {
    	bounds[0] = a;
    } else {
    	printf("%s %.2f.\n", "lower bound less than ", a);
    	return -1;
    }
    for (a = 2; a < ub && fun(a, &param) < 0; a += 1);
    if (a < ub) {
    	bounds[1] = a;
    } else {
    	printf("%s %.2f.\n", "Upper bound larger than ", (double)ub);
    	return -1;
    }

    F.function = fun;
    F.params = &param;

    solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(solver, &F, bounds[0], bounds[1]);

    do {
	iter++;
	status = gsl_root_fsolver_iterate (solver);
	xi = gsl_root_fsolver_root (solver);
	lb = gsl_root_fsolver_x_lower(solver);
	ub = gsl_root_fsolver_x_upper(solver);
	status = gsl_root_test_interval(lb, ub, 1.0e-6, 0);
	// if (status == GSL_SUCCESS)
	//     cout << "Tail index found: xi = " << xi << endl;
    } while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free(solver);

    if (status != GSL_SUCCESS) {
	// cout << "The Brent algorithm did not converge after " << max_iter
	//      << " iterations." << endl;
	xi = -1;
    }
    return xi;
}

double prod_tail_index(double coef1[2], double coef2[2], double cov)
{
    normal_distribution<double> dist;
    double pool[2][SAMPLE_SIZE];
#pragma omp parallel for
    for (int i = 0; i < SAMPLE_SIZE; i++) {
	pool[0][i] = dist(gen);
	pool[1][i] = pool[0][i] * sqrt(cov) + sqrt(1 - cov) * dist(gen);
    }

    gsl_root_fsolver *solver;
    gsl_function F;
    int iter = 0;
    int status = 0;
    int max_iter = 100;
    double lb, ub = 30, xi = -1;

    double a = 2;
    double bounds[2];
    struct prod_tail_fun_param param = {coef1, coef2, pool};
    if (prod_tail_fun(a, &param) > 0) {
	while (prod_tail_fun(a, &param) > 0) a -= 1;
	bounds[0] = a;
	bounds[1] = a + 1;
    } else {
	while (prod_tail_fun(a, &param) < 0) a += 1;
	bounds[0] = a - 1;
	bounds[1] = a;
    }

    F.function = prod_tail_fun;
    F.params = &param;

    solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(solver, &F, bounds[0], bounds[1]);

    do {
	iter++;
	status = gsl_root_fsolver_iterate (solver);
	xi = gsl_root_fsolver_root (solver);
	lb = gsl_root_fsolver_x_lower(solver);
	ub = gsl_root_fsolver_x_upper(solver);
	status = gsl_root_test_interval(lb, ub, 1.0e-6, 0);
	// if (status == GSL_SUCCESS)
	//     cout << "Tail index found: xi = " << xi << endl;
    } while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free(solver);

    if (status != GSL_SUCCESS) {
	// cout << "The Brent algorithm did not converge after " << max_iter
	//      << " iterations." << endl;
	xi = -1;
    }
    return xi;
}

int main(int argc, char *argv[])
{
    // double coef[2];
    // coef[0] = 1.0e-7;
    // coef[1] = stod(argv[1]);
    // coef[2] = stod(argv[2]);

    double coef[][2] = {
	{0.0333595465136058,  0.935462645685956},
	{0.0129814902539519,  0.982156632926710},
	{0.0355747855877685,  0.959357802762939},
	{0.0160964734982745,  0.978232215797576},
	{0.0691075334956582,  0.881676778393785},
	{0.0685060663688339,  0.883690620885534},
	{0.0697364292173322,  0.892126609515694},
	{0.0393300789355985,  0.953107657508420},
	{0.0398855744760083,  0.948950588639140},
	{0.0592371608362950,  0.932987967998614},
	{0.0608751942941075,  0.900850053580394},
	{0.0397208752411525,  0.951096271431723},
	{0.0276614424135275,  0.968677095189546},
	{0.0310147328057462,  0.961939617257637},
	{0.0198512714753665,  0.973700216994867},
	{0.0537300110625453,  0.924474900212370},
	{0.0392929313584584,  0.953560592370058}
    };
    const int d = sizeof(coef)/sizeof(coef[0]);
    double varcov [][17] = {
	{1, 0.669373380195363, 0.464786070633786, 0.232241044065783, 0.301931053094307, 0.300819121771269, 0.471859859133159, 0.461783943939229, 0.241426921444126, 0.364352744551913, 0.632216615688306, 0.355932142424239, 0.582497296688298, 0.315136952413253, 0.761249405223686, 0.648734819278494, 0.455067503601205},
	{0.669373380195363, 1, 0.682216456857965, 0.28564933040021, 0.403200607937382, 0.399262088083759, 0.614613578696547, 0.696521094381235, 0.167223689064189, 0.510527261775464, 0.664674951486154, 0.484201903632382, 0.615133642738711, 0.358333442794494, 0.598801542195879, 0.755826752595378, 0.695850629017868},
	{0.464786070633786, 0.682216456857965, 1, 0.295472569461776, 0.541543103975116, 0.537845544947713, 0.697881335786559, 0.952053242275239, 0.0377028078812318, 0.737624153269522, 0.658324580069529, 0.664439796807126, 0.507937558031846, 0.244289590175636, 0.422932455154974, 0.858972696964364, 0.95460279437536},
	{0.232241044065783, 0.28564933040021, 0.295472569461776, 1, 0.758327675455036, 0.758273490533149, 0.372337785020873, 0.299518256024636, 0.519623144014893, 0.287890848079932, 0.273157434632663, 0.594729950127318, 0.233585460525581, 0.248962749035133, 0.244416584538928, 0.388891006500512, 0.299217604005703},
	{0.301931053094307, 0.403200607937382, 0.541543103975116, 0.758327675455036, 1, 0.998939766976158, 0.571038246597074, 0.549541706370983, 0.448215524917665, 0.520996546755114, 0.405134498547661, 0.834780773836922, 0.291946382845743, 0.314505759477557, 0.321298839914128, 0.58724331538277, 0.552375431152225},
	{0.300819121771269, 0.399262088083759, 0.537845544947713, 0.758273490533149, 0.998939766976158, 1, 0.568271666805865, 0.546154865513154, 0.449564394216971, 0.519007852852421, 0.402333886291538, 0.834336832181381, 0.292050889632662, 0.314662158606586, 0.320709306169859, 0.584895400268014, 0.549222113475339},
	{0.471859859133159, 0.614613578696547, 0.697881335786559, 0.372337785020873, 0.571038246597074, 0.568271666805865, 1, 0.723292920441249, 0.184474639012553, 0.588148471959225, 0.559564101013708, 0.605447900238591, 0.446770912825022, 0.298373972055242, 0.457433612232791, 0.727955872155998, 0.725044220915399},
	{0.461783943939229, 0.696521094381235, 0.952053242275239, 0.299518256024636, 0.549541706370983, 0.546154865513154, 0.723292920441249, 1, 0.0344817307458453, 0.764809766063511, 0.659197075807358, 0.675696589595487, 0.515948265791794, 0.245194647881277, 0.426121169381712, 0.877170785742397, 0.995980032523156},
	{0.241426921444126, 0.167223689064189, 0.0377028078812318, 0.519623144014893, 0.448215524917665, 0.449564394216971, 0.184474639012553, 0.0344817307458453, 1, 0.00253043824878582, 0.202339858787924, 0.316790873190302, 0.259167677878366, 0.206521853362169, 0.239865593910444, 0.191375297169156, 0.0309863915237204},
	{0.364352744551913, 0.510527261775464, 0.737624153269522, 0.287890848079932, 0.520996546755114, 0.519007852852421, 0.588148471959225, 0.764809766063511, 0.00253043824878582, 1, 0.499626263230065, 0.589350485565211, 0.314546071670067, 0.191504689630016, 0.352370235584457, 0.732387061978369, 0.766888496266103},
	{0.632216615688306, 0.664674951486154, 0.658324580069529, 0.273157434632663, 0.405134498547661, 0.402333886291538, 0.559564101013708, 0.659197075807358, 0.202339858787924, 0.499626263230065, 1, 0.492875034021309, 0.55487358505922, 0.30473078275775, 0.552201229768552, 0.790467752060754, 0.653727301011015},
	{0.355932142424239, 0.484201903632382, 0.664439796807126, 0.594729950127318, 0.834780773836922, 0.834336832181381, 0.605447900238591, 0.675696589595487, 0.316790873190302, 0.589350485565211, 0.492875034021309, 1, 0.351095316357393, 0.281533156379539, 0.338765043170621, 0.670760988236642, 0.678074331216711},
	{0.582497296688298, 0.615133642738711, 0.507937558031846, 0.233585460525581, 0.291946382845743, 0.292050889632662, 0.446770912825022, 0.515948265791794, 0.259167677878366, 0.314546071670067, 0.55487358505922, 0.351095316357393, 1, 0.275527806143073, 0.503940837541148, 0.612521333396119, 0.511868452195206},
	{0.315136952413253, 0.358333442794494, 0.244289590175636, 0.248962749035133, 0.314505759477557, 0.314662158606586, 0.298373972055242, 0.245194647881277, 0.206521853362169, 0.191504689630016, 0.30473078275775, 0.281533156379539, 0.275527806143073, 1, 0.289331685108988, 0.326692632695095, 0.242891575572914},
	{0.761249405223686, 0.598801542195879, 0.422932455154974, 0.244416584538928, 0.321298839914128, 0.320709306169859, 0.457433612232791, 0.426121169381712, 0.239865593910444, 0.352370235584457, 0.552201229768552, 0.338765043170621, 0.503940837541148, 0.289331685108988, 1, 0.588431344934906, 0.420920439011483},
	{0.648734819278494, 0.755826752595378, 0.858972696964364, 0.388891006500512, 0.58724331538277, 0.584895400268014, 0.727955872155998, 0.877170785742397, 0.191375297169156, 0.732387061978369, 0.790467752060754, 0.670760988236642, 0.612521333396119, 0.326692632695095, 0.588431344934906, 1, 0.87776751952863},
	{0.455067503601205, 0.695850629017868, 0.95460279437536, 0.299217604005703, 0.552375431152225, 0.549222113475339, 0.725044220915399, 0.995980032523156, 0.0309863915237204, 0.766888496266103, 0.653727301011015, 0.678074331216711, 0.511868452195206, 0.242891575572914, 0.420920439011483, 0.87776751952863, 1}
    };
    array<array<double, d>, d> indices;
    for (int i = 0; i < d; i++) {
    	for (int j = 0; j <= i; j++) {
	    if (j == i) {
		indices[i][i] = tail_index(coef[i]);
		printf("%-8.3f", indices[i][i]);
	    } else {
		indices[i][j] = prod_tail_index(coef[i], coef[j], varcov[i][j]);
		printf("%-8.3f", indices[i][j]);
	    }
    	}
	printf("\n");
    }
    // double I[d];
    // for (int i = 0; i < d; i++) {
    // 	if (i == 7)
    // 	    I[i] = tail_index(coef[7]);
    // 	else
    // 	    I[i] = prod_tail_index(coef[i], coef[7], varcov[i][7]);
    // 	printf("% 5.3f\n", I[i]);
    // }

    // for (a = 2; a > 0 && fun(a, coef) > 0; a -= 1);
    // if (a > 0) {
    // 	bounds[0] = a;
    // } else {
    // 	printf("%s %.2f.\n", "lower bound less than ", a);
    // 	return 0;
    // }
    // for (a = 2; a < ub && fun(a, coef) < 0; a += 1);
    // if (a < ub) {
    // 	bounds[1] = a;
    // } else {
    // 	printf("%s %.2f.\n", "Upper bound larger than ", (double)ub);
    // 	return 0;
    // }
    // double xi = tail_index(coef);
    // cout << xi << endl;
}
