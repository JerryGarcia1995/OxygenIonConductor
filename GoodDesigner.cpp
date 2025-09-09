#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cctype>
#include <sstream>
#include <vector>
#include <cmath>
#include <complex>
#include <set>
#include <chrono>

#define PI 3.14159265359

typedef struct {

	std::string name;
	long double min, ave, max, stddv;

}DESCRIPTOR;

std::vector<DESCRIPTOR> descriptor;

typedef struct {

	std::string atomname;
	int iatom;

}NEIGHBOR;

typedef struct {

	std::string atomname;
	std::vector<long double> v;
	std::vector<long double> zv;
	std::vector<NEIGHBOR> neighbor;

}ATOMDAT;

std::vector<ATOMDAT> atomdat;

typedef struct {

	std::string name;
	std::vector<std::string> sitename;
	std::vector<long double> v;
	std::vector<std::string> strv;

}SDESCRIPTOR;

typedef struct {

	std::string strname;
	std::vector<std::string> sitename;
	std::vector<SDESCRIPTOR> sdescriptor;

}STRDAT;

std::vector<STRDAT> strdat;

typedef struct {

	std::string name, ASB, postfix;
	int ides;
	int** idesB;
	int UType_numerator_S0A1E2N3;
	int UType_denominator_S0A1E2N3;
	std::vector<std::string> UType_numerator_AtomKinds;
	std::vector<std::string> UType_denominator_AtomKinds;
	long double UType_numerator_value;
	long double UType_denominator_value;
	int UType_m0d1;

}ANADESCRIPTOR;

std::vector<ANADESCRIPTOR> anadescriptor;

typedef struct {

	std::string AtomName;
	bool face;
	std::vector<long double> vFeature;

}FACE;

std::vector<FACE> face;

typedef struct {

	std::string name;
	std::vector<long double> coeff_0;
	std::vector<long double> power_0;
	std::vector<long double> coeff_1;
	std::vector<long double> power_1;
	std::string nonlinear_funcname;
	bool sum_only;
	std::vector<long double> coeff_2;
	std::vector<long double> power_2;
	std::vector<int> ianades_2;
	std::vector<int> imamssk_2;
	std::vector<int> ilogln_2;
	std::vector<int> ipow10exp_2;
	std::vector<int> ierf_2;
	bool added;
	std::vector<long double> coeff_3;
	std::vector<long double> power_3;
	std::vector<int> ianades_3;
	std::vector<int> imamssk_3;
	std::vector<int> ilogln_3;
	std::vector<int> ipow10exp_3;
	std::vector<int> ierf_3;

	std::vector<long double> v_list;
	long double v_ave, v_stddv;

}FUNCTYPE;

std::vector<FUNCTYPE> functype;

typedef struct {

	std::vector<std::string> userdefined_name;
	std::string FName;
	std::vector<long double> coeff;

}ORTHFUNCTYPE;

typedef struct {

	bool identified;
	std::string name;
	int type;
	std::string CationName;
	int i0anades;
	int i0mamssk;
	int i1anades[2];
	int i1mamssk[2];
	int i1m0d1;
	int i2whichlastfunc;
	int i3whichlastfunc[2];
	int i3m0d1;

}FEATURE;

typedef struct {

	bool identified, feature_identified;
	std::string name;
	long double coeff;
	int whichorthfunctype;
	int whichfeature;

}FUNC;

typedef struct {

	std::string filename;
	std::vector<FEATURE> feature;
	std::vector<FUNC> func;
	long double intercept;
	long double score_inverse_slope, score_inverse_intercept;
	long double revise_slope, revise_intercept;
	bool allow_beta;

}FUNCCHAIN;

typedef struct {

	std::string name;
	std::string beta_target;
	int ibeta_target;
	std::vector<ORTHFUNCTYPE> orthfunctype;
	std::vector<FUNCCHAIN> funcchain;

	std::vector<long double> v_list;
	long double v_ave, v_stddv;

}REGTYPE;

std::vector<REGTYPE> regtype;

typedef struct {
	
	std::string name, composition;
	std::string strname;
	int istr;
	std::vector<std::string> atomname;
	std::vector<int> iatom;
	std::vector<std::string> sitename;
	std::vector<int> isite;
	std::vector<long double> content;
	long double anion_content;
	long double OM;
	std::vector<long double> CO;

}MAT;

typedef struct {

	std::vector<std::string> atomname;
	std::vector<long double> content;

}SMAT;

typedef struct {

	std::vector<SMAT> smat_min; 
	std::vector<SMAT> smat_max;

}SMATLIST;

std::vector<SMATLIST> smatlist;

std::string AtomDatListFilename;
std::string FaceFilename;
std::string StructureListFilename;
int UseStrName;
std::string MaterialsFilename;
int HeadExists;
std::vector<std::string> PostTagNames;
std::vector<std::string> RegModelFilename;
std::string MapOutSpaceFilename;
std::vector<long double> TargetVReg;
std::vector<long double> ZDistWeight;
std::vector<long double> ModulateContentFactor;
int UseOffset;
int MaxIteration;
std::vector<std::string> NeighborSelectSpace;
long double NeighborZCrit;
std::vector<std::string> ExcludeNeighbors;
std::vector<std::string> DoNotFindNeighbors;
std::string AnionName;
std::vector<std::string> CationNames;
long double AnionValence;
long double VInfinite;
std::string OutputFilename;
std::string AnalysisOutputFilename;
std::vector<std::string> DoNotAnalyze;

bool mapout = false;
int mapout_preNRandom, mapout_NRandom, mapout_Verbose;
std::vector<std::string> sdescriptorname;
std::vector<std::string> sitename;
std::vector<std::string> FaceFeatureName;
int Nmat;
long double** xbox;
long double xbox_OM;
long double* xbox_CO;
long double* xbox_user;
std::vector<long double> zTargetVReg;
int CValence = -1;

std::string sbuf;
long double ldbuf;
std::ifstream readpara, readatomdat, readstrdat, readmatdat, readreg, readregchain, readbpregressor, readfuncchain, readmapout, readface;
std::ofstream writeoutput, writeanalysis;

long double folder(long double v) {

	long double vv;
	if (v < 0.0) {
		while (true) {
			v += 1.0;
			if (v >= 0.0 && v < 1.0) {
				break;
			}
		}
		vv = v;
	}
	else if (v >= 0.0 && v < 1.0) {
		vv = v;
	}
	else {
		while (true) {
			v -= 1.0;
			if (v >= 0.0 && v < 1.0) {
				break;
			}
		}
		vv = v;
	}

	return vv;

}

long double erfinv(long double* x) {

	long double w, p;
	w = -logl((1.0 - *x) * (1.0 + *x));
	if (w < 5.000000) {
		w = w - 2.500000;
		p = 2.81022636e-08;
		p = 3.43273939e-07 + p * w;
		p = -3.5233877e-06 + p * w;
		p = -4.39150654e-06 + p * w;
		p = 0.00021858087 + p * w;
		p = -0.00125372503 + p * w;
		p = -0.00417768164 + p * w;
		p = 0.246640727 + p * w;
		p = 1.50140941 + p * w;
	}
	else {
		w = sqrtl(w) - 3.000000;
		p = -0.000200214257;
		p = 0.000100950558 + p * w;
		p = 0.00134934322 + p * w;
		p = -0.00367342844 + p * w;
		p = 0.00573950773 + p * w;
		p = -0.0076224613 + p * w;
		p = 0.00943887047 + p * w;
		p = 1.00167406 + p * w;
		p = 2.83297682 + p * w;
	}

	return p * *x;

}

long double get_v_special(std::string* ftype, long double* x) {

	if (*ftype == "logit" || *ftype == "isigmoid") {
		return(logl(*x / (1.0 - *x)));
	}
	else if (*ftype == "ilogit" || *ftype == "sigmoid") {
		return(1.0 / (1.0 + exp(-*x)));
	}
	else if (*ftype == "probit" || *ftype == "icdf") {
		long double newx = 2.0 * *x - 1.0;
		return(sqrtl(2.0) * erfinv(&newx));
	}
	else if (*ftype == "iprobit" || *ftype == "cdf") {
		return(0.5 * (1.0 + erfl(*x / sqrtl(2.0))));
	}
	else if (*ftype == "cloglog" || *ftype == "icexpexp") {
		return(logl(-logl(1.0 - *x)));
	}
	else if (*ftype == "icloglog" || *ftype == "cexpexp") {
		return(1.0 - expl(-expl(*x)));
	}
	else if (*ftype == "cauchit" || *ftype == "itangenttype") {
		return(tanl(PI * (*x - 0.5)));
	}
	else if (*ftype == "icauchit" || *ftype == "tangenttype") {
		return(atanl(*x) / PI + 0.5);
	}
	else if (*ftype == "nloglog" || *ftype == "loglog" || *ftype == "iexpexp") {
		return(-logl(-logl(*x)));
	}
	else if (*ftype == "inloglog" || *ftype == "iloglog" || *ftype == "expexp") {
		return(expl(-expl(-*x)));
	}
	else {
		return 0.0;
	}

}

long double get_v(ORTHFUNCTYPE* functype, long double* x) {

	if (functype->FName == "pow") {
		return powl(*x + functype->coeff[0], functype->coeff[1]);
	}
	else if (functype->FName == "log") {
		return powl(logl(*x + functype->coeff[0]) / logl(functype->coeff[1]), functype->coeff[2]);
	}
	else if (functype->FName == "sin") {
		return powl(sinl(powl(functype->coeff[1] * (*x + functype->coeff[0]), functype->coeff[2])), functype->coeff[3]);
	}
	else if (functype->FName == "cos") {
		return powl(cosl(powl(functype->coeff[1] * (*x + functype->coeff[0]), functype->coeff[2])), functype->coeff[3]);
	}
	else if (functype->FName == "tan") {
		return powl(tanl(powl(functype->coeff[1] * (*x + functype->coeff[0]), functype->coeff[2])), functype->coeff[3]);
	}
	else if (functype->FName == "sinh") {
		return powl(sinhl(powl(functype->coeff[1] * (*x + functype->coeff[0]), functype->coeff[2])), functype->coeff[3]);
	}
	else if (functype->FName == "cosh") {
		return powl(coshl(powl(functype->coeff[1] * (*x + functype->coeff[0]), functype->coeff[2])), functype->coeff[3]);
	}
	else if (functype->FName == "tanh") {
		return powl(tanhl(powl(functype->coeff[1] * (*x + functype->coeff[0]), functype->coeff[2])), functype->coeff[3]);
	}
	else if (functype->FName == "exp") {
		return powl(powl(functype->coeff[3], functype->coeff[2] * powl(*x + functype->coeff[0], functype->coeff[1])), functype->coeff[4]);
	}
	else if (functype->FName == "erf") {
		return powl(erfl(functype->coeff[1] * (*x + functype->coeff[0])), functype->coeff[2]);
	}
	else if (functype->FName == "erfc") {
		return powl(erfcl(functype->coeff[1] * (*x + functype->coeff[0])), functype->coeff[2]);
	}
	else if (functype->FName == "abs") {
		return fabsl(powl(*x + functype->coeff[0], functype->coeff[1]));
	}
	else {
		return 0;
	}

}

bool isNumber(const std::string& s) {

	if (s.empty()) { 
		return false; 
	}

	std::istringstream iss(s);
	double d;
	char c;

	if (!(iss >> d)) {
		return false;
	}

	return !(iss >> c);

}

void get_dv(std::vector<long double>* dv, MAT* mat, std::vector<ATOMDAT>* atomdat, std::vector<STRDAT>* strdat, std::vector<ANADESCRIPTOR>* anadescriptor, int* WhichAnaDescriptor) {

	if ((*anadescriptor)[*WhichAnaDescriptor].ASB == "A") {
		for (int k = 0; k < (signed)mat->iatom.size(); k++) {
			dv->push_back((*atomdat)[mat->iatom[k]].v[(*anadescriptor)[*WhichAnaDescriptor].ides]);
		}
	}
	else if ((*anadescriptor)[*WhichAnaDescriptor].ASB == "S") {
		for (int k = 0; k < (signed)mat->isite.size(); k++) {
			dv->push_back((*strdat)[mat->istr].sdescriptor[(*anadescriptor)[*WhichAnaDescriptor].ides].v[mat->isite[k]]);
		}
	}
	else if ((*anadescriptor)[*WhichAnaDescriptor].ASB == "B") {
		for (int k = 0; k < (signed)mat->iatom.size(); k++) {
			dv->push_back((*atomdat)[mat->iatom[k]].v[(*anadescriptor)[*WhichAnaDescriptor].idesB[mat->istr][mat->isite[k]]]);
		}
	}

}

void get_mamssk(std::vector<long double>* mamssk, std::vector<long double>* dv, long double* VInfinite, MAT* mat, std::vector<std::string>* DoNotAnalyze) {

	std::vector<long double>().swap(*mamssk);

	std::vector<long double> finite_dv;
	std::vector<long double> finite_content;
	for (int i = 0; i < (signed)dv->size(); i++) {
		bool analyze = true;
		for (int j = 0; j < (signed)DoNotAnalyze->size(); j++) {
			if (mat->atomname[i] == (*DoNotAnalyze)[j]) {
				analyze = false;
				break;
			}
		}
		if (analyze && (*dv)[i] != 0.0 && mat->content[i] != 0.0) {
			finite_dv.push_back((*dv)[i]);
			finite_content.push_back(mat->content[i]);
		}
	}

	long double min = *VInfinite;
	long double ave = 0.0;
	long double max = -*VInfinite;

	long double content_sum = 0.0;
	for (int i = 0; i < (signed)finite_content.size(); i++) {
		content_sum += finite_content[i];
	}

	for (int i = 0; i < (signed)finite_dv.size(); i++) {
		if (finite_dv[i] < min) {
			min = finite_dv[i];
		}
		ave += finite_content[i] * finite_dv[i];
		if (finite_dv[i] > max) {
			max = finite_dv[i];
		}
	}
	ave /= content_sum;

	long double stddv = 0.0;
	for (int i = 0; i < (signed)finite_dv.size(); i++) {
		stddv += finite_content[i] * powl(finite_dv[i] - ave, 2.0);
	}
	stddv = sqrtl(stddv / content_sum);

	long double skew = 0.0;
	long double kurto = 0.0;
	if (fabsl(stddv) > 1 / *VInfinite) {
		for (int i = 0; i < (signed)finite_dv.size(); i++) {
			skew += finite_content[i] * powl((finite_dv[i] - ave) / stddv, 3.0);
			kurto += finite_content[i] * powl((finite_dv[i] - ave) / stddv, 4.0);
		}
		skew /= content_sum;
		kurto /= content_sum;
	}
	else {
		skew = 0.0;
		kurto = *VInfinite;
	}

	mamssk->push_back(min);
	mamssk->push_back(ave);
	mamssk->push_back(max);
	mamssk->push_back(stddv);
	mamssk->push_back(skew);
	mamssk->push_back(kurto);

	std::vector<long double>().swap(finite_dv);
	std::vector<long double>().swap(finite_content);

}

long double calc_reg(FUNCTYPE* functype) {

	/*xbox_CO!!!*/

	long double value = 0.0;
	for (int i = 0; i < (signed)functype->coeff_2.size(); i++) {
		long double sub_value = functype->coeff_2[i];
		if (functype->power_2[i] != 0.0) {
			if (functype->ianades_2[i] == -1) {
				if (functype->ilogln_2[i] == 1) {
					sub_value = sub_value * powl(log10l(xbox_OM), functype->power_2[i]);
				}
				else if (functype->ilogln_2[i] == 2) {
					sub_value = sub_value * powl(logl(xbox_OM), functype->power_2[i]);
				}
				else if (functype->ipow10exp_2[i] == 1) {
					sub_value = sub_value * powl(powl(10.0, xbox_OM), functype->power_2[i]);
				}
				else if (functype->ipow10exp_2[i] == 2) {
					sub_value = sub_value * powl(powl(2.71828, xbox_OM), functype->power_2[i]);
				}
				else if (functype->ierf_2[i] != -1000) {
					sub_value = sub_value * powl(erfl(powl(10.0, (long double)functype->ierf_2[i]) * xbox_OM), functype->power_2[i]);
				}
				else {
					sub_value = sub_value * powl(xbox_OM, functype->power_2[i]);
				}
			}
			else { 
				if (functype->ilogln_2[i] == 1) {
					sub_value = sub_value * powl(log10l(xbox[functype->ianades_2[i]][functype->imamssk_2[i]]), functype->power_2[i]);
				}
				else if (functype->ilogln_2[i] == 2) {
					sub_value = sub_value * powl(logl(xbox[functype->ianades_2[i]][functype->imamssk_2[i]]), functype->power_2[i]);
				}
				else if (functype->ipow10exp_2[i] == 1) {
					sub_value = sub_value * powl(powl(10.0, xbox[functype->ianades_2[i]][functype->imamssk_2[i]]), functype->power_2[i]);
				}
				else if (functype->ipow10exp_2[i] == 2) {
					sub_value = sub_value * powl(powl(2.71828, xbox[functype->ianades_2[i]][functype->imamssk_2[i]]), functype->power_2[i]);
				}
				else if (functype->ierf_2[i] != -1000) {
					sub_value = sub_value * powl(erfl(powl(10.0, (long double)functype->ierf_2[i]) * xbox[functype->ianades_2[i]][functype->imamssk_2[i]]), functype->power_2[i]);
				}
				else {
					sub_value = sub_value * powl(xbox[functype->ianades_2[i]][functype->imamssk_2[i]], functype->power_2[i]);
				}
			}
		}
		value += sub_value;
	}

	long double value_sofar = value;
	if (functype->nonlinear_funcname != "") {
		long double super_value = 0.0;
		for (int i = 0; i < (signed)functype->coeff_1.size(); i++) {
			long double sub_value = functype->coeff_1[i];
			if (functype->power_1[i] != 0.0) {
				sub_value = sub_value * powl(get_v_special(&functype->nonlinear_funcname, &value), functype->power_1[i]);
			}
			super_value += sub_value;
		}
		long double ssuper_value = 0.0;
		for (int i = 0; i < (signed)functype->coeff_0.size(); i++) {
			long double sub_value = functype->coeff_0[i];
			if (functype->power_0[i] != 0.0) {
				sub_value = sub_value * powl(super_value, functype->power_0[i]);
			}
			ssuper_value += sub_value;
		}
		value_sofar = ssuper_value;
	}

	if (functype->added) {
		long double addedv = 0.0;
		for (int i = 0; i < (signed)functype->coeff_3.size(); i++) {
			long double sub_value = functype->coeff_3[i];
			if (functype->power_3[i] != 0.0) {
				if (functype->ianades_3[i] == -1) {
					if (functype->ilogln_3[i] == 1) {
						sub_value = sub_value * powl(log10l(xbox_OM), functype->power_3[i]);
					}
					else if (functype->ilogln_3[i] == 2) {
						sub_value = sub_value * powl(logl(xbox_OM), functype->power_3[i]);
					}
					else if (functype->ipow10exp_3[i] == 1) {
						sub_value = sub_value * powl(powl(10.0, xbox_OM), functype->power_3[i]);
					}
					else if (functype->ipow10exp_3[i] == 2) {
						sub_value = sub_value * powl(powl(2.71828, xbox_OM), functype->power_3[i]);
					}
					else if (functype->ierf_3[i] != -1000) {
						sub_value = sub_value * powl(erfl(powl(10.0, (long double)functype->ierf_3[i]) * xbox_OM), functype->power_3[i]);
					}
					else {
						sub_value = sub_value * powl(xbox_OM, functype->power_3[i]);
					}
				}
				else {
					if (functype->ilogln_3[i] == 1) {
						sub_value = sub_value * powl(log10l(xbox[functype->ianades_3[i]][functype->imamssk_3[i]]), functype->power_3[i]);
					}
					else if (functype->ilogln_3[i] == 2) {
						sub_value = sub_value * powl(logl(xbox[functype->ianades_3[i]][functype->imamssk_3[i]]), functype->power_3[i]);
					}
					else if (functype->ipow10exp_3[i] == 1) {
						sub_value = sub_value * powl(powl(10.0, xbox[functype->ianades_3[i]][functype->imamssk_3[i]]), functype->power_3[i]);
					}
					else if (functype->ipow10exp_3[i] == 2) {
						sub_value = sub_value * powl(powl(2.71828, xbox[functype->ianades_3[i]][functype->imamssk_3[i]]), functype->power_3[i]);
					}
					else if (functype->ierf_3[i] != -1000) {
						sub_value = sub_value * powl(erfl(powl(10.0, (long double)functype->ierf_3[i]) * xbox[functype->ianades_3[i]][functype->imamssk_3[i]]), functype->power_3[i]);
					}
					else {
						sub_value = sub_value * powl(xbox[functype->ianades_3[i]][functype->imamssk_3[i]], functype->power_3[i]);
					}
				}
			}
			addedv += sub_value;
		}
		value_sofar += addedv;
	}

	return value_sofar;

}

long double calc_reg_chain(REGTYPE* regtype, MAT* objMAT, std::vector<ATOMDAT>* atomdat, std::vector<STRDAT>* strdat, std::vector<ANADESCRIPTOR>* anadescriptor, long double* VInfinite, std::vector<std::string>* DoNotAnalyze, std::vector<std::string>* CationNames) {

	std::vector<long double> funcele, funcele_last;
	for (int i = 0; i < (signed)regtype->funcchain.size(); i++) {
		funcele.clear();
		for (int j = 0; j < regtype->funcchain[i].func.size(); j++) {
			long double v = 0.0;
			if (regtype->funcchain[i].func[j].identified && regtype->funcchain[i].func[j].feature_identified) {
				int ift = regtype->funcchain[i].func[j].whichfeature;
				if (regtype->funcchain[i].feature[ift].type == 0) {
					std::vector<long double> dv, mamssk;
					get_dv(&dv, objMAT, atomdat, strdat, anadescriptor, &regtype->funcchain[i].feature[ift].i0anades);
					get_mamssk(&mamssk, &dv, VInfinite, objMAT, DoNotAnalyze);
					v = mamssk[regtype->funcchain[i].feature[ift].i0mamssk];
					std::vector<long double>().swap(dv);
					std::vector<long double>().swap(mamssk);
				}
				else if (regtype->funcchain[i].feature[ift].type == 1) {
					long double v0 = 0.0;
					long double v1 = 0.0;
					if (regtype->funcchain[i].feature[ift].i1mamssk[0] >= 0) {
						std::vector<long double> dv, mamssk;
						get_dv(&dv, objMAT, atomdat, strdat, anadescriptor, &regtype->funcchain[i].feature[ift].i1anades[0]);
						get_mamssk(&mamssk, &dv, VInfinite, objMAT, DoNotAnalyze);
						v0 = mamssk[regtype->funcchain[i].feature[ift].i1mamssk[0]];
						std::vector<long double>().swap(dv);
						std::vector<long double>().swap(mamssk);
					}
					else if (regtype->funcchain[i].feature[ift].i1mamssk[0] == -1) {
						long double cation_sum = 0.0;
						for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
							cation_sum += objMAT->content[ia];
						}
						v0 = objMAT->anion_content / cation_sum;
					}
					else if (regtype->funcchain[i].feature[ift].i1mamssk[0] == -2) {
						long double cation_sum = 0.0;
						for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
							if (objMAT->atomname[ia] == (*CationNames)[regtype->funcchain[i].feature[ift].i1anades[0]]) {
								cation_sum += objMAT->content[ia];
							}
						}
						v0 = cation_sum / objMAT->anion_content;
					}
					else if (regtype->funcchain[i].feature[ift].i1mamssk[0] == -3) {
						int ianades = regtype->funcchain[i].feature[ift].i1anades[0];
						long double vv0 = 0.0;
						if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 0) {
							for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
								for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_numerator_AtomKinds.size(); ja++) {
									if (objMAT->atomname[ia] == (*anadescriptor)[ianades].UType_numerator_AtomKinds[ja]) {
										vv0 += objMAT->content[ia];
										break;
									}
								}
							}
						}
						else if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 1) {
							for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
								vv0 += objMAT->content[ia];
							}
						}
						else if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 2) {
							for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
								for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_numerator_AtomKinds.size(); ja++) {
									if (objMAT->atomname[ia] != (*anadescriptor)[ianades].UType_numerator_AtomKinds[ja]) {
										vv0 += objMAT->content[ia];
										break;
									}
								}
							}
						}
						else if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 3) {
							vv0 = (*anadescriptor)[ianades].UType_numerator_value;
						}
						long double vv1 = 0.0;
						if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 0) {
							for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
								for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_denominator_AtomKinds.size(); ja++) {
									if (objMAT->atomname[ia] == (*anadescriptor)[ianades].UType_denominator_AtomKinds[ja]) {
										vv1 += objMAT->content[ia];
										break;
									}
								}
							}
						}
						else if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 1) {
							for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
								vv1 += objMAT->content[ia];
							}
						}
						else if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 2) {
							for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
								for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_denominator_AtomKinds.size(); ja++) {
									if (objMAT->atomname[ia] != (*anadescriptor)[ianades].UType_denominator_AtomKinds[ja]) {
										vv1 += objMAT->content[ia];
										break;
									}
								}
							}
						}
						else if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 3) {
							vv1 = (*anadescriptor)[ianades].UType_denominator_value;
						}
						if ((*anadescriptor)[ianades].UType_m0d1 == 0) {
							v0 = vv0 * vv1;
						}
						else {
							v0 = vv0 / vv1;
						}
					}
					if (regtype->funcchain[i].feature[ift].i1mamssk[1] >= 0) {
						std::vector<long double> dv, mamssk;
						get_dv(&dv, objMAT, atomdat, strdat, anadescriptor, &regtype->funcchain[i].feature[ift].i1anades[1]);
						get_mamssk(&mamssk, &dv, VInfinite, objMAT, DoNotAnalyze);
						v1 = mamssk[regtype->funcchain[i].feature[ift].i1mamssk[1]];
						std::vector<long double>().swap(dv);
						std::vector<long double>().swap(mamssk);
					}
					else if (regtype->funcchain[i].feature[ift].i1mamssk[1] == -1) {
						long double cation_sum = 0.0;
						for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
							cation_sum += objMAT->content[ia];
						}
						v1 = objMAT->anion_content / cation_sum;
					}
					else if (regtype->funcchain[i].feature[ift].i1mamssk[1] == -2) {
						long double cation_sum = 0.0;
						for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
							if (objMAT->atomname[ia] == (*CationNames)[regtype->funcchain[i].feature[ift].i1anades[1]]) {
								cation_sum += objMAT->content[ia];
							}
						}
						v1 = cation_sum / objMAT->anion_content;
					}
					else if (regtype->funcchain[i].feature[ift].i1mamssk[1] == -3) {
						int ianades = regtype->funcchain[i].feature[ift].i1anades[1];
						long double vv0 = 0.0;
						if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 0) {
							for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
								for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_numerator_AtomKinds.size(); ja++) {
									if (objMAT->atomname[ia] == (*anadescriptor)[ianades].UType_numerator_AtomKinds[ja]) {
										vv0 += objMAT->content[ia];
										break;
									}
								}
							}
						}
						else if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 1) {
							for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
								vv0 += objMAT->content[ia];
							}
						}
						else if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 2) {
							for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
								for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_numerator_AtomKinds.size(); ja++) {
									if (objMAT->atomname[ia] != (*anadescriptor)[ianades].UType_numerator_AtomKinds[ja]) {
										vv0 += objMAT->content[ia];
										break;
									}
								}
							}
						}
						else if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 3) {
							vv0 = (*anadescriptor)[ianades].UType_numerator_value;
						}
						long double vv1 = 0.0;
						if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 0) {
							for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
								for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_denominator_AtomKinds.size(); ja++) {
									if (objMAT->atomname[ia] == (*anadescriptor)[ianades].UType_denominator_AtomKinds[ja]) {
										vv1 += objMAT->content[ia];
										break;
									}
								}
							}
						}
						else if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 1) {
							for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
								vv1 += objMAT->content[ia];
							}
						}
						else if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 2) {
							for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
								for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_denominator_AtomKinds.size(); ja++) {
									if (objMAT->atomname[ia] != (*anadescriptor)[ianades].UType_denominator_AtomKinds[ja]) {
										vv1 += objMAT->content[ia];
										break;
									}
								}
							}
						}
						else if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 3) {
							vv1 = (*anadescriptor)[ianades].UType_denominator_value;
						}
						if ((*anadescriptor)[ianades].UType_m0d1 == 0) {
							v1 = vv0 * vv1;
						}
						else {
							v1 = vv0 / vv1;
						}
					}
					if (regtype->funcchain[i].feature[ift].i1m0d1 == 0) {
						v = v0 * v1;
					}
					else {
						v = v0 / v1;
					}
				}
				else if (regtype->funcchain[i].feature[ift].type == 2) {
					v = funcele_last[regtype->funcchain[i].feature[ift].i2whichlastfunc];
				}
				else if (regtype->funcchain[i].feature[ift].type == 3) {
					long double v0 = funcele_last[regtype->funcchain[i].feature[ift].i3whichlastfunc[0]];
					long double v1 = funcele_last[regtype->funcchain[i].feature[ift].i3whichlastfunc[1]];
					if (regtype->funcchain[i].feature[ift].i3m0d1 == 0) {
						v = v0 * v1;
					}
					else {
						v = v0 / v1;
					}
				}
				else if (regtype->funcchain[i].feature[ift].type == -1) {
					long double cation_sum = 0.0;
					for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
						cation_sum += objMAT->content[ia];
					}
					v = objMAT->anion_content / cation_sum;
				}
				else if (regtype->funcchain[i].feature[ift].type == -2) {
					long double cation_sum = 0.0;
					for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
						if (objMAT->atomname[ia] == regtype->funcchain[i].feature[ift].CationName) {
							cation_sum += objMAT->content[ia];
						}
					}
					v = cation_sum / objMAT->anion_content;
				}
				else if (regtype->funcchain[i].feature[ift].type == -3) {
					int ianades = regtype->funcchain[i].feature[ift].i0anades;
					long double v0 = 0.0;
					if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 0) {
						for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
							for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_numerator_AtomKinds.size(); ja++) {
								if (objMAT->atomname[ia] == (*anadescriptor)[ianades].UType_numerator_AtomKinds[ja]) {
									v0 += objMAT->content[ia];
									break;
								}
							}
						}
					}
					else if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 1) {
						for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
							v0 += objMAT->content[ia];
						}
					}
					else if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 2) {
						for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
							for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_numerator_AtomKinds.size(); ja++) {
								if (objMAT->atomname[ia] != (*anadescriptor)[ianades].UType_numerator_AtomKinds[ja]) {
									v0 += objMAT->content[ia];
									break;
								}
							}
						}
					}
					else if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 3) {
						v0 = (*anadescriptor)[ianades].UType_numerator_value;
					}
					long double v1 = 0.0;
					if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 0) {
						for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
							for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_denominator_AtomKinds.size(); ja++) {
								if (objMAT->atomname[ia] == (*anadescriptor)[ianades].UType_denominator_AtomKinds[ja]) {
									v1 += objMAT->content[ia];
									break;
								}
							}
						}
					}
					else if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 1) {
						for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
							v1 += objMAT->content[ia];
						}
					}
					else if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 2) {
						for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
							for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_denominator_AtomKinds.size(); ja++) {
								if (objMAT->atomname[ia] != (*anadescriptor)[ianades].UType_denominator_AtomKinds[ja]) {
									v1 += objMAT->content[ia];
									break;
								}
							}
						}
					}
					else if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 3) {
						v1 = (*anadescriptor)[ianades].UType_denominator_value;
					}
					if ((*anadescriptor)[ianades].UType_m0d1 == 0) {
						v = v0 * v1;
					}
					else {
						v = v0 / v1;
					}
				}
				v = get_v(&(regtype->orthfunctype)[regtype->funcchain[i].func[j].whichorthfunctype], &v);
			}
			funcele.push_back(v);
		}
		funcele_last.clear();
		funcele_last = funcele;
	}

	long double value = 0.0;
	for (int j = 0; j < regtype->funcchain[regtype->funcchain.size() - 1].func.size(); j++) {
		value += regtype->funcchain[regtype->funcchain.size() - 1].func[j].coeff * funcele[j];
	}
	value += regtype->funcchain[regtype->funcchain.size() - 1].intercept;
	
	std::vector<long double>().swap(funcele);
	std::vector<long double>().swap(funcele_last);

	if (regtype->funcchain[regtype->funcchain.size() - 1].allow_beta) {
		std::string ftype = "i" + regtype->beta_target;
		value = regtype->funcchain[regtype->funcchain.size() - 1].score_inverse_slope * get_v_special(&ftype, &value) + regtype->funcchain[regtype->funcchain.size() - 1].score_inverse_intercept;
		value = regtype->funcchain[regtype->funcchain.size() - 1].revise_slope * value + regtype->funcchain[regtype->funcchain.size() - 1].revise_intercept;
	}

	return value;

}

void read_para() {

	readpara.open("basic_properties_GoodDesigner.txt");

	readpara >> sbuf >> AtomDatListFilename;
	readpara >> sbuf >> FaceFilename;
	readpara >> sbuf >> StructureListFilename >> UseStrName;
	readpara >> sbuf >> MaterialsFilename >> HeadExists;
	std::getline(readpara, sbuf);
	std::stringstream lcstr;
	lcstr << sbuf;
	while (true) {
		lcstr >> sbuf;
		if (sbuf != "") {
			PostTagNames.push_back(sbuf);
		}
		if (lcstr.eof()) {
			break;
		}
	}

	readpara >> sbuf;
	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	while (true) {
		lcstr >> sbuf;
		if (sbuf == "-") {
			break;
		}
		RegModelFilename.push_back(sbuf);
		if (lcstr.eof()) {
			break;
		}
	}

	readpara >> sbuf >> MapOutSpaceFilename;
	if (MapOutSpaceFilename != "-") {
		mapout = true;
	}

	readpara >> sbuf;
	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	for (int i = 0; i < (signed)RegModelFilename.size(); i++) {
		lcstr >> ldbuf;
		TargetVReg.push_back(ldbuf);
		if (lcstr.eof()) {
			break;
		}
	}

	readpara >> sbuf;
	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	for (int i = 0; i < (signed)RegModelFilename.size(); i++) {
		lcstr >> ldbuf;
		ZDistWeight.push_back(ldbuf);
		if (lcstr.eof()) {
			break;
		}
	}

	readpara >> sbuf;
	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	while (true) {
		lcstr >> ldbuf;
		ModulateContentFactor.push_back(ldbuf);
		if (lcstr.eof()) {
			break;
		}
	}
	
	readpara >> sbuf >> UseOffset;
	readpara >> sbuf >> MaxIteration;

	readpara >> sbuf;
	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	while (true) {
		lcstr >> sbuf;
		NeighborSelectSpace.push_back(sbuf);
		if (lcstr.eof()) {
			break;
		}
	}

	readpara >> sbuf >> NeighborZCrit;
	readpara >> sbuf;
	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	if (sbuf.find("-") == std::string::npos) {
		while (true) {
			lcstr >> sbuf;
			ExcludeNeighbors.push_back(sbuf);
			if (lcstr.eof()) {
				break;
			}
		}
	}

	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	lcstr >> sbuf;
	while (true) {
		lcstr >> sbuf;
		if (sbuf != "-") {
			DoNotFindNeighbors.push_back(sbuf);
		}
		if (lcstr.eof()) {
			break;
		}
	}
	
	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	lcstr >> sbuf >> AnionName;
	if (!lcstr.eof()) {
		while (true) {
			lcstr >> sbuf;
			CationNames.push_back(sbuf);
			if (lcstr.eof()) {
				break;
			}
		}
	}
	
	readpara >> sbuf >> AnionValence;
	readpara >> sbuf >> VInfinite;
	readpara >> sbuf >> OutputFilename;
	readpara >> sbuf >> AnalysisOutputFilename;
	std::getline(readpara, sbuf);

	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	lcstr >> sbuf;
	if (!lcstr.eof()) {
		while (true) {
			lcstr >> sbuf;
			DoNotAnalyze.push_back(sbuf);
			if (lcstr.eof()) {
				break;
			}
		}
	}

	readpara >> sbuf;
	int Ucount = 0;
	while (true) {

		ANADESCRIPTOR objANADESCRIPTOR;
		readpara >> objANADESCRIPTOR.name;
		if (objANADESCRIPTOR.name == "STOP") {
			break;
		}

		readpara >> objANADESCRIPTOR.ASB;
		if (objANADESCRIPTOR.ASB == "B") {
			readpara >> objANADESCRIPTOR.postfix;
		}
		else if (objANADESCRIPTOR.ASB == "U") {
			Ucount++;
			objANADESCRIPTOR.UType_numerator_S0A1E2N3 = 0;
			objANADESCRIPTOR.UType_denominator_S0A1E2N3 = 0;
			while (true) {
				readpara >> sbuf;
				if (sbuf == "*" || sbuf == "/") {
					if (sbuf == "*") {
						objANADESCRIPTOR.UType_m0d1 = 0;
					}
					else {
						objANADESCRIPTOR.UType_m0d1 = 1;
					}
					break;
				}
				else {
					if (sbuf == "All") {
						objANADESCRIPTOR.UType_numerator_S0A1E2N3 = 1;
					}
					else if (sbuf == "!") {
						objANADESCRIPTOR.UType_numerator_S0A1E2N3 = 2;
					}
					else {
						if (isNumber(sbuf)) {
							objANADESCRIPTOR.UType_numerator_S0A1E2N3 = 3;
							std::stringstream ldbuf;
							ldbuf.str("");
							ldbuf.clear();
							ldbuf << sbuf;
							ldbuf >> objANADESCRIPTOR.UType_numerator_value;
						}
						else {
							objANADESCRIPTOR.UType_numerator_AtomKinds.push_back(sbuf);
						}
					}
				}
			}
			std::getline(readpara, sbuf);
			std::stringstream decompbuf;
			decompbuf.str("");
			decompbuf.clear();
			decompbuf << sbuf;
			while (true) {
				decompbuf >> sbuf;
				if (sbuf == "All") {
					objANADESCRIPTOR.UType_denominator_S0A1E2N3 = 1;
				}
				else if (sbuf == "!") {
					objANADESCRIPTOR.UType_denominator_S0A1E2N3 = 2;
				}
				else {
					if (isNumber(sbuf)) {
						objANADESCRIPTOR.UType_denominator_S0A1E2N3 = 3;
						std::stringstream ldbuf;
						ldbuf.str("");
						ldbuf.clear();
						ldbuf << sbuf;
						ldbuf >> objANADESCRIPTOR.UType_denominator_value;
					}
					else {
						objANADESCRIPTOR.UType_denominator_AtomKinds.push_back(sbuf);
					}
				}
				if (decompbuf.eof()) {
					break;
				}
			}
		}

		if (!mapout) {
			anadescriptor.push_back(objANADESCRIPTOR);
		}
		else if (objANADESCRIPTOR.ASB == "A" || objANADESCRIPTOR.ASB == "U") {
			anadescriptor.push_back(objANADESCRIPTOR);
		}

		if (readpara.eof()) {
			break;
		}

	}

	xbox = new long double *[anadescriptor.size()];
	for (int i = 0; i < (signed)anadescriptor.size(); i++) {
		xbox[i] = new long double[6];
		for (int j = 0; j < 6; j++) {
			xbox[i][j] = 0.0;
		}
	}
	xbox_OM = 0.0;
	xbox_CO = new long double[CationNames.size()];
	for (int i = 0; i < (signed)CationNames.size(); i++) {
		xbox_CO[i] = 0.0;
	}
	xbox_user = new long double[Ucount];
	for (int i = 0; i < Ucount; i++) {
		xbox_user[i] = 0.0;
	}

	readpara.close();

}

void read_atomdat(std::vector<ATOMDAT>* atomdat, std::vector<DESCRIPTOR> * descriptor, std::vector<ANADESCRIPTOR> * anadescriptor, std::string * AtomDatListFilename) {

	readatomdat.open(*AtomDatListFilename);

	std::getline(readatomdat, sbuf);
	std::stringstream lcstr;
	lcstr << sbuf;
	lcstr >> sbuf;
	while (true) {
		DESCRIPTOR objDESCRIPTOR;
		lcstr >> objDESCRIPTOR.name;
		if (objDESCRIPTOR.name == "Valence") {
			CValence = (signed)descriptor->size();
		}
		descriptor->push_back(objDESCRIPTOR);
		if (lcstr.eof()) {
			break;
		}
	}

	for (int ii = 0; ii < (signed)anadescriptor->size(); ii++) {
		if ((*anadescriptor)[ii].ASB == "A") {
			for (int i = 0; i < (signed)descriptor->size(); i++) {
				if ((*descriptor)[i].name == (*anadescriptor)[ii].name) {
					(*anadescriptor)[ii].ides = i;
					break;
				}
			}
		}
	}
	
	while (true) {
		ATOMDAT objATOMDAT;
		readatomdat >> objATOMDAT.atomname;
		if (objATOMDAT.atomname == "STOP") {
			break;
		}
		for(int i = 0; i < (signed)descriptor->size(); i++){
			readatomdat >> ldbuf;
			objATOMDAT.v.push_back(ldbuf);
		}
		atomdat->push_back(objATOMDAT);
		if (readatomdat.eof()) {
			break;
		}
	}

	readatomdat.close();

	for (int j = 0; j < (signed)descriptor->size(); j++) {
		(*descriptor)[j].min = 10000000000000000000000.0;
		(*descriptor)[j].ave = 0.0;
		(*descriptor)[j].max = -10000000000000000000000.0;
		for (int i = 0; i < (signed)atomdat->size(); i++) {
			if ((*atomdat)[i].v[j] < (*descriptor)[j].min) {
				(*descriptor)[j].min = (*atomdat)[i].v[j];
			}
			(*descriptor)[j].ave += (*atomdat)[i].v[j];
			if ((*atomdat)[i].v[j] > (*descriptor)[j].max) {
				(*descriptor)[j].max = (*atomdat)[i].v[j];
			}
		}
		(*descriptor)[j].ave /= (long double)(signed)atomdat->size();
		(*descriptor)[j].stddv = 0.0;
		for (int i = 0; i < (signed)atomdat->size(); i++) {
			(*descriptor)[j].stddv += powl((*atomdat)[i].v[j] - (*descriptor)[j].ave, 2.0);
		}
		(*descriptor)[j].stddv /= (long double)(signed)atomdat->size();
		(*descriptor)[j].stddv = sqrtl((*descriptor)[j].stddv);
		for (int i = 0; i < (signed)atomdat->size(); i++) {
			(*atomdat)[i].zv.push_back(((*atomdat)[i].v[j] - (*descriptor)[j].ave) / (*descriptor)[j].stddv);
		}
	}

	writeoutput << "Atom\t";
	writeoutput.flush();
	for (int j = 0; j < (signed)descriptor->size(); j++) {
		writeoutput << "z" << (*descriptor)[j].name;
		writeoutput.flush();
		if (j != (signed)descriptor->size() - 1) {
			writeoutput << "\t";
			writeoutput.flush();
		}
		else {
			writeoutput << "\n";
			writeoutput.flush();
		}
	}
	for (int i = 0; i < (signed)atomdat->size(); i++) {
		writeoutput << (*atomdat)[i].atomname << "\t";
		writeoutput.flush();
		for (int j = 0; j < (signed)descriptor->size(); j++) {
			writeoutput << (*atomdat)[i].zv[j];
			writeoutput.flush();
			if (j != (signed)descriptor->size() - 1) {
				writeoutput << "\t";
				writeoutput.flush();
			}
			else {
				writeoutput << "\n";
				writeoutput.flush();
			}
		}
	}

	writeoutput << "=================================================================================\n";
	writeoutput.flush();
	
}

void read_face(std::vector<FACE>* face, std::vector<std::string>* FaceFeatureName, std::string* FaceFilename) {

	readface.open(*FaceFilename);

	std::getline(readface, sbuf);
	std::stringstream lcstream;
	lcstream << sbuf;
	lcstream >> sbuf >> sbuf;
	while (true) {
		lcstream >> sbuf;
		FaceFeatureName->push_back(sbuf);
		if (lcstream.eof()) {
			break;
		}
	}

	while (true) {
		std::getline(readface, sbuf);
		if (readface.eof() || sbuf.find("STOP") != std::string::npos) {
			break;
		}
		else {
			lcstream.str("");
			lcstream.clear();
			lcstream << sbuf;
			FACE objFACE;
			lcstream >> objFACE.AtomName;
			lcstream >> sbuf;
			if (sbuf == "-") {
				objFACE.face = false;
			}
			else {
				objFACE.face = true;
			}
			long double ldbuf;
			for (int i = 0; i < (signed)FaceFeatureName->size(); i++) {
				lcstream >> ldbuf;
				objFACE.vFeature.push_back(ldbuf);
			}
			face->push_back(objFACE);
		}
	}

	readface.close();

}

void pick_neighbor(std::vector<ATOMDAT>* atomdat, std::vector<std::string>* NeighborSelectSpace, long double * NeighborZCrit, std::vector<std::string>* ExcludeNeighbors, std::vector<std::string>* DoNotFindNeighbors, std::vector<DESCRIPTOR> * descriptor) {

	std::vector<int> ides;
	for (int i = 0; i < (signed)NeighborSelectSpace->size(); i++) {
		for (int j = 0; j < (signed)descriptor->size(); j++) {
			if ((*NeighborSelectSpace)[i] == (*descriptor)[j].name) {
				ides.push_back(j);
				break;
			}
		}
	}

	for (int ii = 0; ii < (signed)atomdat->size(); ii++) {
		bool needtofind = true;
		for (int jj = 0; jj < (signed)DoNotFindNeighbors->size(); jj++) {
			if ((*atomdat)[ii].atomname == (*DoNotFindNeighbors)[jj]) {
				needtofind = false;
				break;
			}
		}
		if (needtofind) {
			for (int jj = 0; jj < (signed)atomdat->size(); jj++) {
				bool exclude = false;
				for (int kk = 0; kk < (signed)ExcludeNeighbors->size(); kk++) {
					if ((*atomdat)[jj].atomname == (*ExcludeNeighbors)[kk]) {
						exclude = true;
						break;
					}
				}
				if (ii != jj && !exclude) {
					long double zdist = 0.0;
					for (int i = 0; i < (signed)ides.size(); i++) {
						zdist += powl((*atomdat)[jj].zv[ides[i]] - (*atomdat)[ii].zv[ides[i]], 2.0);
					}
					zdist = sqrtl(zdist);
					if (zdist <= *NeighborZCrit) {
						NEIGHBOR objNEIGHBOR;
						objNEIGHBOR.atomname = (*atomdat)[jj].atomname;
						objNEIGHBOR.iatom = jj;
						(*atomdat)[ii].neighbor.push_back(objNEIGHBOR);
					}
				}
			}
		}
	}
	
	for (int ii = 0; ii < (signed)atomdat->size(); ii++) {
		writeoutput << (*atomdat)[ii].atomname << "\t...NEIGHBOR...";
		writeoutput.flush();
		if ((*atomdat)[ii].neighbor.size() > 0) {
			writeoutput << "\t";
			writeoutput.flush();
			for (int jj = 0; jj < (signed)(*atomdat)[ii].neighbor.size(); jj++) {
				writeoutput << (*atomdat)[ii].neighbor[jj].atomname;
				writeoutput.flush();
				if (jj != (signed)(*atomdat)[ii].neighbor.size() - 1) {
					writeoutput << "\t";
					writeoutput.flush();
				}
				else {
					writeoutput << "\n";
					writeoutput.flush();
				}
			}
		}
		else {
			writeoutput << "\n";
			writeoutput.flush();
		}
	}

	writeoutput << "=================================================================================\n";
	writeoutput.flush();
	
}

void read_strdat(std::vector<STRDAT> * strdat, std::vector<ANADESCRIPTOR>* anadescriptor, std::vector<std::string> * sdescriptorname, std::vector<std::string>* sitename, std::string * StructureListFilename, std::vector<DESCRIPTOR>* descriptor, std::vector<ATOMDAT>* atomdat, bool *usestrname) {

	readstrdat.open(*StructureListFilename);

	readstrdat >> sbuf;
	std::getline(readstrdat, sbuf);
	std::stringstream lcstr;
	lcstr << sbuf;
	while (true) {
		lcstr >> sbuf;
		sdescriptorname->push_back(sbuf);
		if (lcstr.eof()) {
			break;
		}
	}

	readstrdat >> sbuf;
	std::getline(readstrdat, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	while (true) {
		lcstr >> sbuf;
		sitename->push_back(sbuf);
		if (lcstr.eof()) {
			break;
		}
	}

	std::getline(readstrdat, sbuf);

	while (true) {
		readstrdat >> sbuf;
		if (sbuf == "STOP" || readstrdat.eof()) {
			break;
		}
		else {
			STRDAT objSTRDAT;
			if (*usestrname) {
				objSTRDAT.strname = sbuf;
			}
			else {
				objSTRDAT.strname = "None";
			}
			for (int i = 0; i < (signed)sdescriptorname->size(); i++) {
				SDESCRIPTOR objSDESCRIPTOR;
				objSDESCRIPTOR.name = (*sdescriptorname)[i];
				for (int j = 0; j < (signed)sitename->size(); j++) {
					readstrdat >> sbuf;
					if (sbuf != "X") {
						objSDESCRIPTOR.sitename.push_back((*sitename)[j]);
						objSDESCRIPTOR.strv.push_back(sbuf);
						std::stringstream lcstr;
						lcstr << sbuf;
						lcstr >> ldbuf;
						objSDESCRIPTOR.v.push_back(ldbuf);
					}
				}
				objSTRDAT.sdescriptor.push_back(objSDESCRIPTOR);
			}
			strdat->push_back(objSTRDAT);
		}
		if (readstrdat.eof()) {
			break;
		}
	}

	for (int i = 0; i < (signed)strdat->size(); i++) {
		for (int j = 0; j < (signed)(*strdat)[i].sdescriptor[0].sitename.size(); j++) {
			(*strdat)[i].sitename.push_back((*strdat)[i].sdescriptor[0].sitename[j]);
		}
	}

	readstrdat.close();

	for (int ii = 0; ii < (signed)anadescriptor->size(); ii++) {
		if ((*anadescriptor)[ii].ASB == "S") {
			for (int i = 0; i < (signed)sdescriptorname->size(); i++) {
				if ((*sdescriptorname)[i] == (*anadescriptor)[ii].name) {
					(*anadescriptor)[ii].ides = i;
					break;
				}
			}
		}
		else if ((*anadescriptor)[ii].ASB == "B") {
			(*anadescriptor)[ii].idesB = new int* [strdat->size()];
			for (int i = 0; i < (signed)strdat->size(); i++) {
				(*anadescriptor)[ii].idesB[i] = new int[sitename->size()];
				for (int j = 0; j < (signed)sitename->size(); j++) {
					(*anadescriptor)[ii].idesB[i][j] = -1;
					std::string target_name = (*anadescriptor)[ii].name;
					int target_ipostfix = -1;
					for (int k = 0; k < (signed)sdescriptorname->size(); k++) {
						if ((*sdescriptorname)[k] == (*anadescriptor)[ii].postfix) {
							target_ipostfix = k;
							break;
						}
					}
					if (target_ipostfix != -1 && j < (*strdat)[i].sdescriptor[target_ipostfix].strv.size()) {
						target_name += (*strdat)[i].sdescriptor[target_ipostfix].strv[j];
					}
					for (int k = 0; k < (signed)descriptor->size(); k++) {
						if ((*descriptor)[k].name == target_name) {
							(*anadescriptor)[ii].idesB[i][j] = k;
							break;
						}
					}
				}
			}
		}
	}

}

void find_nonlinear(std::string * functype, std::string * termline) {

	*functype = "";

	if ((termline->find("logit") != std::string::npos && termline->find("inverse") == std::string::npos) || (termline->find("sigmoid") != std::string::npos && termline->find("inverse") != std::string::npos)) {
		*functype = "logit";
	}
	else if ((termline->find("sigmoid") != std::string::npos && termline->find("inverse") == std::string::npos) || (termline->find("logit") != std::string::npos && termline->find("inverse") != std::string::npos)) {
		*functype = "sigmoid";
	}
	else if ((termline->find("probit") != std::string::npos && termline->find("inverse") == std::string::npos) || (termline->find("cdf") != std::string::npos && termline->find("inverse") != std::string::npos)) {
		*functype = "probit";
	}
	else if ((termline->find("cdf") != std::string::npos && termline->find("inverse") == std::string::npos) || (termline->find("probit") != std::string::npos && termline->find("inverse") != std::string::npos)) {
		*functype = "cdf";
	}
	else if ((termline->find("cloglog") != std::string::npos && termline->find("inverse") == std::string::npos) || (termline->find("cexpexp") != std::string::npos && termline->find("inverse") != std::string::npos)) {
		*functype = "cloglog";
	}
	else if ((termline->find("cexpexp") != std::string::npos && termline->find("inverse") == std::string::npos) || (termline->find("cloglog") != std::string::npos && termline->find("inverse") != std::string::npos)) {
		*functype = "cexpexp";
	}
	else if ((termline->find("cauchit") != std::string::npos && termline->find("inverse") == std::string::npos) || (termline->find("tangenttype") != std::string::npos && termline->find("inverse") != std::string::npos)) {
		*functype = "cauchit";
	}
	else if ((termline->find("tangenttype") != std::string::npos && termline->find("inverse") == std::string::npos) || (termline->find("cauchit") != std::string::npos && termline->find("inverse") != std::string::npos)) {
		*functype = "tangenttype";
	}
	else if ((termline->find("nloglog") != std::string::npos && termline->find("inverse") == std::string::npos) || (termline->find("expexp") != std::string::npos && termline->find("inverse") != std::string::npos)) {
		*functype = "nloglog";
	}
	else if ((termline->find("expexp") != std::string::npos && termline->find("inverse") == std::string::npos) || (termline->find("nloglog") != std::string::npos && termline->find("inverse") != std::string::npos)) {
		*functype = "expexp";
	}

}

void analyze_term(long double* coeff, long double* power, std::vector<std::string>* term) {

	if (term->size() == 0) {
		*coeff = 0.0;
		*power = 0.0;
	}
	else if (term->size() == 1) {
		std::stringstream lcstr;
		lcstr << (*term)[0];
		lcstr >> ldbuf;
		*coeff = ldbuf;
		*power = 0.0;
	}
	else {
		if ((*term)[0] == (*term)[1]) {
			(*term)[0] = "1.0";
		}
		std::stringstream lcstr;
		lcstr << (*term)[0];
		lcstr >> ldbuf;
		*coeff = ldbuf;
		if ((*term)[term->size() - 1].find("^") == std::string::npos) {
			*power = 1.0;
		}
		else {
			sbuf = (*term)[term->size() - 1].substr((*term)[term->size() - 1].find("^") + 1);
			lcstr.str("");
			lcstr.clear();
			lcstr << sbuf;
			lcstr >> ldbuf;
			*power = ldbuf;
		}
	}

}

void identify_anades(int* ianades, int* imamssk, int* ilogln, int* ipow10exp, int* ierf, std::string * term, std::vector<ANADESCRIPTOR> * anadescriptor) {

	*ianades = -1; // OM
	size_t max_length = 0;
	for (int i = 0; i < (signed)anadescriptor->size(); i++) {
		if((term->find("^") != std::string::npos && term->find((*anadescriptor)[i].name + "^") != std::string::npos) || (term->find("^") == std::string::npos && term->find((*anadescriptor)[i].name) != std::string::npos)) {
			if ((*anadescriptor)[i].name.length() > max_length) {
				*ianades = i;
				max_length = (*anadescriptor)[i].name.length();
			}
		}
	}

	*ipow10exp = -1;
	if (term->find("Min") != std::string::npos) {
		*imamssk = 0;
		if (term->find("10Min") != std::string::npos && term->find("log10Min") == std::string::npos && term->find("Log10Min") == std::string::npos) {
			*ipow10exp = 1;
		}
		else if (term->find("expMin") != std::string::npos || term->find("ExpMin") != std::string::npos) {
			*ipow10exp = 2;
		}
	}
	else if (term->find("Ave") != std::string::npos) {
		*imamssk = 1;
		if (term->find("10Ave") != std::string::npos && term->find("log10Ave") == std::string::npos && term->find("Log10Ave") == std::string::npos) {
			*ipow10exp = 1;
		}
		else if (term->find("expAve") != std::string::npos || term->find("ExpAve") != std::string::npos) {
			*ipow10exp = 2;
		}
	}
	else if (term->find("Max") != std::string::npos) {
		*imamssk = 2;
		if (term->find("10Max") != std::string::npos && term->find("log10Max") == std::string::npos && term->find("Log10Max") == std::string::npos) {
			*ipow10exp = 1;
		}
		else if (term->find("expMax") != std::string::npos || term->find("ExpMax") != std::string::npos) {
			*ipow10exp = 2;
		}
	}
	else if (term->find("StdDv") != std::string::npos) {
		*imamssk = 3;
		if (term->find("10StdDv") != std::string::npos && term->find("log10StdDv") == std::string::npos && term->find("Log10StdDv") == std::string::npos) {
			*ipow10exp = 1;
		}
		else if (term->find("expStdDv") != std::string::npos || term->find("ExpStdDv") != std::string::npos) {
			*ipow10exp = 2;
		}
	}
	else if (term->find("Skew") != std::string::npos) {
		*imamssk = 4;
		if (term->find("10Skew") != std::string::npos && term->find("log10Skew") == std::string::npos && term->find("Log10Skew") == std::string::npos) {
			*ipow10exp = 1;
		}
		else if (term->find("expSkew") != std::string::npos || term->find("ExpSkew") != std::string::npos) {
			*ipow10exp = 2;
		}
	}
	else if (term->find("Kurto") != std::string::npos) {
		*imamssk = 5;
		if (term->find("10Kurto") != std::string::npos && term->find("log10Kurto") == std::string::npos && term->find("Log10Kurto") == std::string::npos) {
			*ipow10exp = 1;
		}
		else if (term->find("expKurto") != std::string::npos || term->find("ExpKurto") != std::string::npos) {
			*ipow10exp = 2;
		}
	}
	else {
		*imamssk = -1; //intercept
	}

	if (term->find("Log10") != std::string::npos || term->find("log10") != std::string::npos) {
		*ilogln = 1;
	}
	else if (term->find("Log") != std::string::npos || term->find("log") != std::string::npos || term->find("Ln") != std::string::npos || term->find("ln") != std::string::npos) {
		*ilogln = 2;
	}
	else {
		*ilogln = -1;
	}

	*ierf = -1000;
	if (term->find("erf") != std::string::npos || term->find("Erf") != std::string::npos) {
		std::string terminology = "erf";
		if (term->find("Erf") != std::string::npos) {
			terminology = "Erf";
		}
		*ierf = 0;
		for (int ierf_find = -9; ierf_find <= 9; ierf_find++) {
			if (term->find(terminology + std::to_string(ierf_find)) != std::string::npos) {
				*ierf = ierf_find;
				break;
			}
		}
	}

}

bool identify_anades_chain(int* itype, int* ianades, int* imamssk, std::string* CationName, std::string* term, std::vector<ANADESCRIPTOR>* anadescriptor, std::string* AnionName, std::vector<std::string>* CationNames) {

	bool found = false;

	if (*term == *AnionName + "M") {
		*imamssk = -1;
		*ianades = -1;
		*itype = -1;
		found = true;
	}
	else {
		for (int i = 0; i < (signed)anadescriptor->size(); i++) {
			if (*term == "Min" + (*anadescriptor)[i].name) {
				*imamssk = 0;
				*ianades = i;
				*itype = 0;
				found = true;
			}
			else if (*term == "Ave" + (*anadescriptor)[i].name) {
				*imamssk = 1;
				*ianades = i;
				*itype = 0;
				found = true;
			}
			else if (*term == "Max" + (*anadescriptor)[i].name) {
				*imamssk = 2;
				*ianades = i;
				*itype = 0;
				found = true;
			}
			else if (*term == "StdDv" + (*anadescriptor)[i].name) {
				*imamssk = 3;
				*ianades = i;
				*itype = 0;
				found = true;
			}
			else if (*term == "Skew" + (*anadescriptor)[i].name) {
				*imamssk = 4;
				*ianades = i;
				*itype = 0;
				found = true;
			}
			else if (*term == "Kurto" + (*anadescriptor)[i].name) {
				*imamssk = 5;
				*ianades = i;
				*itype = 0;
				found = true;
			}
			else if (*term == (*anadescriptor)[i].name && (*anadescriptor)[i].ASB == "U") {
				*imamssk = -3;
				*ianades = i;
				*itype = -3;
				found = true;
			}
			if (found) {
				break;
			}
		}
	}

	if (!found) {
		for (int i = 0; i < CationNames->size(); i++) {
			if (*term == (*CationNames)[i] + *AnionName) {
				*imamssk = -2;
				*ianades = i;
				*itype = -2;
				*CationName = (*CationNames)[i];
				found = true;
				break;
			}
		}
	}

	return found;

}


void read_reg(std::vector<FUNCTYPE> * functype, std::string * RegFilename, std::vector<ANADESCRIPTOR>* anadescriptor) {

	readreg.open(*RegFilename);

	if (readreg.is_open()) {
		FUNCTYPE objFUNCTYPE;
		objFUNCTYPE.added = false;
		readreg >> sbuf >> objFUNCTYPE.name;
		std::getline(readreg, sbuf);
		std::getline(readreg, sbuf);
		std::vector<std::string> slinebuf;
		while (true) {
			std::getline(readreg, sbuf);
			if (sbuf.find("=====") == std::string::npos) {
				slinebuf.push_back(sbuf);
			}
			else {
				break;
			}
		}
		if (slinebuf.size() >= 2) {
			objFUNCTYPE.sum_only = false;
			find_nonlinear(&objFUNCTYPE.nonlinear_funcname, &slinebuf[1]);
			for (int iline = 0; iline < 2; iline++) {
				std::stringstream lcstr;
				lcstr << slinebuf[iline];
				while (true) {
					lcstr >> sbuf;
					if (sbuf == "=") {
						break;
					}
				}
				while (true) {
					std::vector<std::string> term;
					while (true) {
						lcstr >> sbuf;
						if (sbuf == "+") {
							break;
						}
						else {
							if (sbuf != "*") {
								term.push_back(sbuf);
							}
						}
						if (lcstr.eof()) {
							break;
						}
					}
					long double coeff, power;
					analyze_term(&coeff, &power, &term);
					if (iline == 0) {
						objFUNCTYPE.coeff_0.push_back(coeff);
						objFUNCTYPE.power_0.push_back(power);
					}
					else {
						objFUNCTYPE.coeff_1.push_back(coeff);
						objFUNCTYPE.power_1.push_back(power);
					}
					std::vector<std::string>().swap(term);
					if (lcstr.eof()) {
						break;
					}
				}
			}
		}
		else {
			objFUNCTYPE.sum_only = true;
		}
		std::vector<std::string>().swap(slinebuf);
		std::getline(readreg, sbuf);
		bool addfeatures_exist = false;
		while (true) {
			std::vector<std::string> term;
			std::getline(readreg, sbuf);
			if (sbuf.find("STOP") != std::string::npos || sbuf.find("*****") != std::string::npos || readreg.eof()) {
				if (sbuf.find("*****") != std::string::npos) {
					addfeatures_exist = true;
				}
				break;
			}
			std::stringstream lcstr;
			lcstr << sbuf;
			lcstr >> sbuf;
			if (sbuf != "intercept") {
				term.push_back(sbuf);
				term.push_back(sbuf);
			}
			lcstr >> sbuf;
			if (term.size() != 0) {
				term[0] = sbuf;
			}
			else {
				term.push_back(sbuf);
			}
			long double coeff, power;
			analyze_term(&coeff, &power, &term);
			objFUNCTYPE.coeff_2.push_back(coeff);
			objFUNCTYPE.power_2.push_back(power);
			if (term.size() == 1) {
				objFUNCTYPE.ianades_2.push_back(-1);
				objFUNCTYPE.imamssk_2.push_back(-1);
				objFUNCTYPE.ilogln_2.push_back(-1);
				objFUNCTYPE.ipow10exp_2.push_back(-1);
				objFUNCTYPE.ierf_2.push_back(-1000);
			}
			else {
				int ianades, imamssk, ilogln, ipow10exp, ierf;
				identify_anades(&ianades, &imamssk, &ilogln, &ipow10exp, &ierf, &term[1], anadescriptor);
				objFUNCTYPE.ianades_2.push_back(ianades);
				objFUNCTYPE.imamssk_2.push_back(imamssk);
				objFUNCTYPE.ilogln_2.push_back(ilogln);
				objFUNCTYPE.ipow10exp_2.push_back(ipow10exp);
				objFUNCTYPE.ierf_2.push_back(ierf);
			}
			std::vector<std::string>().swap(term);
		}
		if (addfeatures_exist) {
			objFUNCTYPE.added = true;
			std::getline(readreg, sbuf);
			while (true) {
				std::vector<std::string> term;
				std::getline(readreg, sbuf);
				if (sbuf.find("STOP") != std::string::npos || readreg.eof()) {
					break;
				}
				std::stringstream lcstr;
				lcstr << sbuf;
				lcstr >> sbuf;
				if (sbuf != "intercept") {
					term.push_back(sbuf);
					term.push_back(sbuf);
				}
				lcstr >> sbuf;
				if (term.size() != 0) {
					term[0] = sbuf;
				}
				else {
					term.push_back(sbuf);
				}
				long double coeff, power;
				analyze_term(&coeff, &power, &term);
				objFUNCTYPE.coeff_3.push_back(coeff);
				objFUNCTYPE.power_3.push_back(power);
				if (term.size() == 1) {
					objFUNCTYPE.ianades_3.push_back(-1);
					objFUNCTYPE.imamssk_3.push_back(-1);
					objFUNCTYPE.ilogln_3.push_back(-1);
					objFUNCTYPE.ipow10exp_3.push_back(-1);
					objFUNCTYPE.ierf_3.push_back(-1000);
				}
				else {
					int ianades, imamssk, ilogln, ipow10exp, ierf;
					identify_anades(&ianades, &imamssk, &ilogln, &ipow10exp, &ierf, &term[1], anadescriptor);
					objFUNCTYPE.ianades_3.push_back(ianades);
					objFUNCTYPE.imamssk_3.push_back(imamssk);
					objFUNCTYPE.ilogln_3.push_back(ilogln);
					objFUNCTYPE.ipow10exp_3.push_back(ipow10exp);
					objFUNCTYPE.ierf_3.push_back(ierf);
				}
				std::vector<std::string>().swap(term);
			}
		}
		functype->push_back(objFUNCTYPE);
	}

	readreg.close();

}

REGTYPE read_reg_chain(std::string* RegFilename, std::vector<ANADESCRIPTOR>* anadescriptor, std::string* AnionName, std::vector<std::string>* CationNames) {

	readregchain.open(*RegFilename);
	REGTYPE objREGTYPE;
	std::stringstream localstr;

	if (readregchain.is_open()) {

		readregchain >> sbuf >> objREGTYPE.name;
		readregchain >> sbuf >> objREGTYPE.beta_target;
		
		std::string bpregressor_filename;
		readregchain >> sbuf >> bpregressor_filename;
		readbpregressor.open(bpregressor_filename);
		if (readbpregressor.is_open()) {
			while (true) {
				readbpregressor >> sbuf;
				if (sbuf == "Functypes") {
					break;
				}
			}
			while (true) {
				std::getline(readbpregressor, sbuf);
				if (sbuf.find("BetaRegression") == std::string::npos) {
					ORTHFUNCTYPE objORTHFUNCTYPE;
					localstr.str("");
					localstr.clear();
					localstr << sbuf;
					localstr >> sbuf >> sbuf;
					std::stringstream getname(sbuf);
					std::string seg;
					while (std::getline(getname, seg, '^')) {
						objORTHFUNCTYPE.userdefined_name.push_back(seg);
					}
					localstr >> objORTHFUNCTYPE.FName;
					int NCoeff = 0;
					if (objORTHFUNCTYPE.FName == "pow" || objORTHFUNCTYPE.FName == "abs") {
						NCoeff = 2;
					}
					else if (objORTHFUNCTYPE.FName == "log" || objORTHFUNCTYPE.FName == "erf" || objORTHFUNCTYPE.FName == "erfc") {
						NCoeff = 3;
					}
					else if (objORTHFUNCTYPE.FName == "sin" || objORTHFUNCTYPE.FName == "cos" || objORTHFUNCTYPE.FName == "tan" || objORTHFUNCTYPE.FName == "sinh" || objORTHFUNCTYPE.FName == "cosh" || objORTHFUNCTYPE.FName == "tanh") {
						NCoeff = 4;
					}
					else if (objORTHFUNCTYPE.FName == "exp") {
						NCoeff = 5;
					}
					for (int j = 0; j < NCoeff; j++) {
						long double ldbuf;
						localstr >> sbuf >> ldbuf;
						objORTHFUNCTYPE.coeff.push_back(ldbuf);
					}
					objREGTYPE.orthfunctype.push_back(objORTHFUNCTYPE);
				}
				else {
					break;
				}
			}
			localstr.str("");
			localstr.clear();
			localstr << sbuf;
			localstr >> sbuf;
			objREGTYPE.ibeta_target = 0;
			if (objREGTYPE.beta_target != "none") {
				objREGTYPE.ibeta_target = 1;
				while (true) {
					localstr >> sbuf;
					if (sbuf == objREGTYPE.beta_target || localstr.eof()) {
						break;
					}
					else {
						objREGTYPE.ibeta_target++;
					}
				}
			}
		}
		readbpregressor.close();

		readregchain >> sbuf;
		std::getline(readregchain, sbuf);
		while (true) {
			std::getline(readregchain, sbuf);
			if (sbuf == "STOP" || readregchain.eof()) {
				break;
			}
			else {
				FUNCCHAIN objFUNCCHAIN;
				localstr.str("");
				localstr.clear();
				localstr << sbuf;
				localstr >> objFUNCCHAIN.filename;
				objREGTYPE.funcchain.push_back(objFUNCCHAIN);
			}
		}

	}

	for (int ichain = 0; ichain < (signed)objREGTYPE.funcchain.size(); ichain++) {

		std::vector<std::string> dataset_firstcolumn;
		std::vector<std::string> beta_score_inverse;
		std::vector<std::string> beta_revised;
		readfuncchain.open(objREGTYPE.funcchain[ichain].filename);
		beta_revised.push_back("");
		if (readfuncchain.is_open()) {
			while (true) {
				readfuncchain >> sbuf;
				if (readfuncchain.eof()) {
					break;
				}
				else if (sbuf == "******************************DATASET******************************") {
					std::getline(readfuncchain, sbuf);
					readfuncchain >> sbuf;
					if (sbuf.find("Name") != std::string::npos) {
						readfuncchain >> sbuf;
					}
					std::getline(readfuncchain, sbuf);
					localstr.str("");
					localstr.clear();
					localstr << sbuf;
					std::string dataset_firstcolumn_element = "";
					while (true) {
						localstr >> sbuf;
						dataset_firstcolumn_element += sbuf;
						if (localstr.eof()) {
							break;
						}
						else {
							dataset_firstcolumn_element += " ";
						}
					}
					dataset_firstcolumn.push_back(dataset_firstcolumn_element);
				}
				else if (sbuf == "ModelDF") {
					readfuncchain >> sbuf;
					std::getline(readfuncchain, sbuf);
					beta_score_inverse.push_back(sbuf);
				}
				else if (sbuf == "det_r2_revised") {
					readfuncchain >> sbuf >> sbuf;
					if (sbuf == ":") {
						std::getline(readfuncchain, sbuf);
						beta_revised.push_back(sbuf);
					}
				}
			}
			
		}
		readfuncchain.close();
		
		std::string target_dataset_firstcolumn = "";
		std::string target_beta_score_inverse = "";
		std::string target_beta_revised = "";
		int ibeta_t = objREGTYPE.ibeta_target;
		if (objREGTYPE.ibeta_target >= (signed)dataset_firstcolumn.size()) {
			ibeta_t = 0;
		}
		target_dataset_firstcolumn = dataset_firstcolumn[ibeta_t];
		ibeta_t = objREGTYPE.ibeta_target;
		if (objREGTYPE.ibeta_target >= (signed)beta_score_inverse.size()) {
			ibeta_t = 0;
		}
		target_beta_score_inverse = beta_score_inverse[ibeta_t];
		ibeta_t = objREGTYPE.ibeta_target;
		if (objREGTYPE.ibeta_target >= (signed)beta_revised.size()) {
			ibeta_t = 0;
		}
		target_beta_revised = beta_revised[ibeta_t];

		objREGTYPE.funcchain[ichain].allow_beta = false;
		if (target_beta_score_inverse != "" && target_beta_revised != "") {
			objREGTYPE.funcchain[ichain].allow_beta = true;
		}

		if (target_beta_score_inverse != "") {
			localstr.str("");
			localstr.clear();
			localstr << target_beta_score_inverse;
			localstr >> sbuf >> sbuf >> objREGTYPE.funcchain[ichain].score_inverse_slope >> sbuf >> sbuf >> sbuf >> objREGTYPE.funcchain[ichain].score_inverse_intercept;
		}
		if (target_beta_revised != "") {
			localstr.str("");
			localstr.clear();
			localstr << target_beta_revised;
			localstr >> sbuf >> sbuf >> objREGTYPE.funcchain[ichain].revise_intercept >> sbuf >> objREGTYPE.funcchain[ichain].revise_slope >> sbuf;
		}

		std::vector<std::string>().swap(dataset_firstcolumn);
		std::vector<std::string>().swap(beta_score_inverse);
		std::vector<std::string>().swap(beta_revised);
		localstr.str("");
		localstr.clear();
		localstr << target_dataset_firstcolumn;
		while (true) {
			FEATURE objFEATURE;
			localstr >> objFEATURE.name;
			objFEATURE.identified = identify_anades_chain(&objFEATURE.type, &objFEATURE.i0anades, &objFEATURE.i0mamssk, &objFEATURE.CationName, &objFEATURE.name, anadescriptor, AnionName, CationNames);
			objREGTYPE.funcchain[ichain].feature.push_back(objFEATURE);
			if (localstr.eof()) {
				break;
			}
		}

		for (int ift = 0; ift < (signed)objREGTYPE.funcchain[ichain].feature.size(); ift++) {
			if (!objREGTYPE.funcchain[ichain].feature[ift].identified) {
				std::cout << objREGTYPE.name << " -- identifying funcchain " << ichain << " : " << ift << "\n";
				std::cout.flush();
				for (int jft = 0; jft < (signed)objREGTYPE.funcchain[ichain].feature.size(); jft++) {
					if (objREGTYPE.funcchain[ichain].feature[jft].identified) {
						for (int kft = 0; kft < (signed)objREGTYPE.funcchain[ichain].feature.size(); kft++) {
							if (objREGTYPE.funcchain[ichain].feature[kft].identified) {
								if (objREGTYPE.funcchain[ichain].feature[ift].name == "(" + objREGTYPE.funcchain[ichain].feature[jft].name + "*" + objREGTYPE.funcchain[ichain].feature[kft].name + ")") {
									objREGTYPE.funcchain[ichain].feature[ift].identified = true;
									objREGTYPE.funcchain[ichain].feature[ift].type = 1;
									objREGTYPE.funcchain[ichain].feature[ift].i1mamssk[0] = objREGTYPE.funcchain[ichain].feature[jft].i0mamssk;
									objREGTYPE.funcchain[ichain].feature[ift].i1anades[0] = objREGTYPE.funcchain[ichain].feature[jft].i0anades;
									objREGTYPE.funcchain[ichain].feature[ift].i1mamssk[1] = objREGTYPE.funcchain[ichain].feature[kft].i0mamssk;
									objREGTYPE.funcchain[ichain].feature[ift].i1anades[1] = objREGTYPE.funcchain[ichain].feature[kft].i0anades;
									objREGTYPE.funcchain[ichain].feature[ift].i1m0d1 = 0;
								}
								else if (objREGTYPE.funcchain[ichain].feature[ift].name == "(" + objREGTYPE.funcchain[ichain].feature[jft].name + "/" + objREGTYPE.funcchain[ichain].feature[kft].name + ")") {
									objREGTYPE.funcchain[ichain].feature[ift].identified = true;
									objREGTYPE.funcchain[ichain].feature[ift].type = 1;
									objREGTYPE.funcchain[ichain].feature[ift].i1mamssk[0] = objREGTYPE.funcchain[ichain].feature[jft].i0mamssk;
									objREGTYPE.funcchain[ichain].feature[ift].i1anades[0] = objREGTYPE.funcchain[ichain].feature[jft].i0anades;
									objREGTYPE.funcchain[ichain].feature[ift].i1mamssk[1] = objREGTYPE.funcchain[ichain].feature[kft].i0mamssk;
									objREGTYPE.funcchain[ichain].feature[ift].i1anades[1] = objREGTYPE.funcchain[ichain].feature[kft].i0anades;
									objREGTYPE.funcchain[ichain].feature[ift].i1m0d1 = 1;
								}
							}
							if (objREGTYPE.funcchain[ichain].feature[ift].identified) {
								break;
							}
						}
					}
					if (objREGTYPE.funcchain[ichain].feature[ift].identified) {
						break;
					}
				}
			}
		}

		if (ichain > 0) {

			for (int ift = 0; ift < (signed)objREGTYPE.funcchain[ichain].feature.size(); ift++) {
				if (!objREGTYPE.funcchain[ichain].feature[ift].identified) {
					std::cout << objREGTYPE.name << " -- identifying funcchain " << ichain << " : " << ift << "\n";
					std::cout.flush();
					for (int ifc = 0; ifc < (signed)objREGTYPE.funcchain[ichain - 1].func.size(); ifc++) {
						if (objREGTYPE.funcchain[ichain - 1].func[ifc].identified && objREGTYPE.funcchain[ichain - 1].func[ifc].feature_identified) {
							if (objREGTYPE.funcchain[ichain].feature[ift].name == "(" + objREGTYPE.funcchain[ichain - 1].func[ifc].name + ")") {
								objREGTYPE.funcchain[ichain].feature[ift].identified = true;
								objREGTYPE.funcchain[ichain].feature[ift].type = 2;
								objREGTYPE.funcchain[ichain].feature[ift].i2whichlastfunc = ifc;
								break;
							}
						}
					}
				}
			}

			for (int ift = 0; ift < (signed)objREGTYPE.funcchain[ichain].feature.size(); ift++) {
				if (!objREGTYPE.funcchain[ichain].feature[ift].identified) {
					std::cout << objREGTYPE.name << " -- identifying funcchain " << ichain << " : " << ift << "\n";
					std::cout.flush();
					for (int ifc = 0; ifc < (signed)objREGTYPE.funcchain[ichain - 1].func.size(); ifc++) {
						if (objREGTYPE.funcchain[ichain - 1].func[ifc].identified && objREGTYPE.funcchain[ichain - 1].func[ifc].feature_identified) {
							for (int jfc = 0; jfc < (signed)objREGTYPE.funcchain[ichain - 1].func.size(); jfc++) {
								if (objREGTYPE.funcchain[ichain - 1].func[jfc].identified && objREGTYPE.funcchain[ichain - 1].func[jfc].feature_identified) {
									if (objREGTYPE.funcchain[ichain].feature[ift].name == "(" + objREGTYPE.funcchain[ichain - 1].func[ifc].name + "*" + objREGTYPE.funcchain[ichain - 1].func[jfc].name + ")") {
										objREGTYPE.funcchain[ichain].feature[ift].identified = true;
										objREGTYPE.funcchain[ichain].feature[ift].type = 3;
										objREGTYPE.funcchain[ichain].feature[ift].i3whichlastfunc[0] = ifc;
										objREGTYPE.funcchain[ichain].feature[ift].i3whichlastfunc[1] = jfc;
										objREGTYPE.funcchain[ichain].feature[ift].i3m0d1 = 0;
										break;
									}
									else if (objREGTYPE.funcchain[ichain].feature[ift].name == "(" + objREGTYPE.funcchain[ichain - 1].func[ifc].name + "/" + objREGTYPE.funcchain[ichain - 1].func[jfc].name + ")") {
										objREGTYPE.funcchain[ichain].feature[ift].identified = true;
										objREGTYPE.funcchain[ichain].feature[ift].type = 3;
										objREGTYPE.funcchain[ichain].feature[ift].i3whichlastfunc[0] = ifc;
										objREGTYPE.funcchain[ichain].feature[ift].i3whichlastfunc[1] = jfc;
										objREGTYPE.funcchain[ichain].feature[ift].i3m0d1 = 1;
										break;
									}
								}
							}
						}
						if (objREGTYPE.funcchain[ichain].feature[ift].identified) {
							break;
						}
					}
				}
			}

		}

		for (int ift = 0; ift < (signed)objREGTYPE.funcchain[ichain].feature.size(); ift++) {
			if (!objREGTYPE.funcchain[ichain].feature[ift].identified) {
				std::cout << "[ERROR] " << objREGTYPE.funcchain[ichain].filename << "\n";
				std::cout.flush();
				std::cout << "[ERROR] " << objREGTYPE.funcchain[ichain].feature[ift].name << " cannot be identified!!\n";
				std::cout.flush();
			}
		}

		readfuncchain.open(objREGTYPE.funcchain[ichain].filename);
		if (readfuncchain.is_open()) {
			int beta_count = -1;
			while (true) {
				readfuncchain >> sbuf;
				if (readfuncchain.eof()) {
					break;
				}
				else if (sbuf == "fit_variable") {
					beta_count++;
					if (beta_count == objREGTYPE.ibeta_target) {
						break;
					}
				}
			}
			std::getline(readfuncchain, sbuf);
			while (true) {
				std::getline(readfuncchain, sbuf);
				if (sbuf.find("dataset.size()") != std::string::npos) {
					break;
				}
				else if (sbuf.find("intercept") == std::string::npos) {
					std::stringstream lcstr;
					lcstr.str("");
					lcstr.clear();
					lcstr << sbuf;
					FUNC objFUNC;
					lcstr >> objFUNC.name >> objFUNC.coeff;
					std::cout << objREGTYPE.name << " -- identifying funcchain " << ichain << " : " << objFUNC.name << "\n";
					std::cout.flush();
					objFUNC.identified = false;
					for (int ift = 0; ift < (signed)objREGTYPE.funcchain[ichain].feature.size(); ift++) {
						for (int jfc = 0; jfc < (signed)objREGTYPE.orthfunctype.size(); jfc++) {
							if (objREGTYPE.orthfunctype[jfc].userdefined_name.size() > 0) {
								std::string referred_name = objREGTYPE.orthfunctype[jfc].userdefined_name[0] + objREGTYPE.funcchain[ichain].feature[ift].name;
								if (objREGTYPE.orthfunctype[jfc].userdefined_name.size() > 1) {
									for (int iname = 1; iname < (signed)objREGTYPE.orthfunctype[jfc].userdefined_name.size(); iname++) {
										referred_name += "^" + objREGTYPE.orthfunctype[jfc].userdefined_name[iname];
									}
								}
								if (objFUNC.name == referred_name) {
									objFUNC.identified = true;
									objFUNC.whichfeature = ift;
									objFUNC.whichorthfunctype = jfc;
									objFUNC.feature_identified = objREGTYPE.funcchain[ichain].feature[ift].identified;
									break;
								}
							}
						}
						if (objFUNC.identified) {
							break;
						}
					}
					objREGTYPE.funcchain[ichain].func.push_back(objFUNC);
				}
				else {
					std::stringstream lcstr;
					lcstr.str("");
					lcstr.clear();
					lcstr << sbuf;
					lcstr >> sbuf >> objREGTYPE.funcchain[ichain].intercept;
				}
			}
		}
		readfuncchain.close();

	}

	readregchain.close();
	return objREGTYPE;

}

void set_par(std::string* MaterialsFilename, bool *headexists) {

	readmatdat.open(*MaterialsFilename);

	if (readmatdat.is_open()) {

		if (*headexists) {
			std::getline(readmatdat, sbuf);
		}

		Nmat = 0;
		while (true) {
			std::getline(readmatdat, sbuf);
			if (sbuf.find("STOP") != std::string::npos) {
				break;
			}
			else {
				Nmat++;
			}
			if (readmatdat.eof()) {
				break;
			}
		}

		readmatdat.close();

		if (sbuf == "") {
			Nmat--;
		}

		if (Nmat <= 1) {
			writeoutput << Nmat << " material is found...\n";
			writeoutput.flush();
		}
		else {
			writeoutput << Nmat << " materials are found...\n";
			writeoutput.flush();
		}
		writeoutput << "=================================================================================\n";
		writeoutput.flush();

	}

}

void designmate_iatom(MAT * mat, std::vector<ATOMDAT> * atomdat) {

	std::vector<int>().swap(mat->iatom);

	for (int i = 0; i < mat->atomname.size(); i++) {
		int iatom = -1;
		for (int l = 0; l < atomdat->size(); l++) {
			if ((*atomdat)[l].atomname == mat->atomname[i]) {
				iatom = l;
				break;
			}
		}
		mat->iatom.push_back(iatom);
	}
	
}

void designmate_isite(MAT * mat, std::vector<std::string> * sitename) {

	for (int i = 0; i < mat->sitename.size(); i++) {
		int isite = -1;
		for (int l = 0; l < sitename->size(); l++) {
			if ((*sitename)[l] == mat->sitename[i]) {
				isite = l;
				break;
			}
		}
		mat->isite.push_back(isite);
	}

}

void write_zdist(MAT * mat, long double * zdist, std::string* AnionName) {

	writeoutput << mat->composition << "\t";
	writeoutput.flush();
	for (int i = 0; i < (signed)mat->atomname.size(); i++) {
		writeoutput << mat->atomname[i] << "\t" << mat->content[i] << "\t";
		writeoutput.flush();
	}
	if (*AnionName != "-") {
		writeoutput << *AnionName << "\t" << mat->anion_content << "\tzDist=\t" << *zdist;
		writeoutput.flush();
	}

}

void process_matdat(std::vector<long double>* zTargetVReg, std::vector<long double>* ZDistWeight, std::vector<long double>* TargetVReg, std::vector<long double>* ModulateContentFactor, std::string* MaterialsFilename, std::vector<REGTYPE>* regtype, std::vector<ATOMDAT>* atomdat, std::string* FaceFilename, std::vector<std::string>* FaceFeatureName, std::vector<FACE>* face, std::vector<STRDAT>* strdat, std::vector<std::string>* sitename, std::vector<ANADESCRIPTOR>* anadescriptor, std::string* AnionName, std::vector<std::string>* CationNames, long double* AnionValence, long double* VInfinite, bool* collect_only, bool* UseOffset, int* MaxIteration, bool* headexists, bool* usestrname, std::vector<std::string>* PostTagNames, std::vector<std::string>* DoNotAnalyze) {

	readmatdat.open(*MaterialsFilename);
	
	if (readmatdat.is_open()) {

		if (*headexists) {
			std::getline(readmatdat, sbuf);
		}

		MAT gMAT;
		std::vector<long double> gmonitor_value;
		if (!*collect_only) {
			for (int k = 0; k < (signed)regtype->size(); k++) {
				gmonitor_value.push_back(0.0);
			}
		}
		long double global_min_zdist = *VInfinite;
		std::string gstr, pgstr;

		for (int i = 0; i < Nmat; i++) {

			std::getline(readmatdat, sbuf);
			std::string linebuf = sbuf;

			if (i >= 0 && i <= Nmat - 1) {

				std::vector<std::string> PostTags;
				if(PostTagNames->size() > 0){ 
					for (int ipost = 0; ipost < (signed)PostTagNames->size(); ipost++) {
						PostTags.push_back("");
					}
					std::stringstream linestream;
					linestream << linebuf;
					std::string subbuf;
					while (true) {
						linestream >> subbuf;
						if (PostTags.size() > 1) {
							for (unsigned int ii = 1; ii < PostTags.size(); ii++) {
								PostTags[ii - 1] = PostTags[ii];
							}
						}
						PostTags[PostTags.size() - 1] = subbuf;
						if (linestream.eof()) {
							break;
						}
					}
				}

				MAT objMAT;
				std::stringstream lcstr;
				lcstr << sbuf;
				lcstr >> objMAT.name;
				lcstr >> sbuf;
				lcstr >> objMAT.strname;
				if (!*usestrname) {
					objMAT.strname = "None";
				}
				objMAT.istr = -1;
				for (int j = 0; j < (signed)strdat->size(); j++) {
					if (objMAT.strname == (*strdat)[j].strname) {
						objMAT.istr = j;
						break;
					}
				}
				while (true) {
					lcstr >> sbuf;
					if (sbuf == "||") {
						break;
					}
				}
				objMAT.composition = "";
				long double cation_sum = 0.0;
				long double* scation_sum;
				scation_sum = new long double[CationNames->size()];
				for (int j = 0; j < (signed)CationNames->size(); j++) {
					scation_sum[j] = 0.0;
				}
				for (int j = 0; j < (signed)sitename->size(); j++) {
					while (true) {
						lcstr >> sbuf;
						if (sbuf != "||") {
							objMAT.atomname.push_back(sbuf);
							objMAT.sitename.push_back((*sitename)[j]);
							lcstr >> ldbuf;
							objMAT.content.push_back(ldbuf);
							cation_sum += ldbuf;
							for (int k = 0; k < (signed)CationNames->size(); k++) {
								if (objMAT.atomname[objMAT.atomname.size() - 1] == (*CationNames)[k]) {
									scation_sum[k] += ldbuf;
									break;
								}
							}
							objMAT.composition += sbuf;
							if (ldbuf != 1.0) {
								objMAT.composition += std::to_string(ldbuf).substr(0, 5);
							}
						}
						else {
							break;
						}
					}
				}
				designmate_iatom(&objMAT, atomdat);
				designmate_isite(&objMAT, sitename);
				if (*AnionName != "-") {
					lcstr >> sbuf >> objMAT.anion_content;
					objMAT.composition += sbuf;
					if (objMAT.anion_content != 1.0) {
						objMAT.composition += std::to_string(objMAT.anion_content).substr(0, 5);
					}
					objMAT.OM = objMAT.anion_content / cation_sum;
					for (int k = 0; k < (signed)CationNames->size(); k++) {
						objMAT.CO.push_back(scation_sum[k] / objMAT.anion_content);
					}
				}
				else {
					objMAT.OM = 0.0;
				}
				delete[] scation_sum;

				writeoutput << i + 1 << "\t...\t" << objMAT.name << "\t...\t" << objMAT.composition << "\t...\t" << objMAT.strname << "\n";
				writeoutput.flush();
				if (*collect_only) {
					writeanalysis << i + 1 << "\t" << objMAT.name << "\t" << objMAT.composition << "\t" << objMAT.strname << "\t";
					writeanalysis.flush();
					if (*AnionName != "-") {
						writeanalysis << objMAT.OM << "\t";
						writeanalysis.flush();
						for (int k = 0; k < (signed)CationNames->size(); k++) {
							writeanalysis << objMAT.CO[k] << "\t";
							writeanalysis.flush();
						}
					}
				}

				if (*collect_only) {	

					writeoutput << *AnionName << "M\t" << objMAT.OM << "\n";
					writeoutput.flush();
					for (int k = 0; k < (signed)CationNames->size(); k++) {
						writeoutput << (*CationNames)[k] << *AnionName << "\t" << objMAT.CO[k];
						writeoutput.flush();
						if (k != (signed)CationNames->size() - 1) {
							writeoutput << "\t";
							writeoutput.flush();
						}
						else {
							writeoutput << "\n";
							writeoutput.flush();
						}
					}
					writeoutput << "Atomname\t";
					writeoutput.flush();
					for (int j = 0; j < (signed)objMAT.atomname.size(); j++) {
						writeoutput << objMAT.atomname[j];
						writeoutput.flush();
						if (j != (signed)objMAT.atomname.size() - 1) {
							writeoutput << "\t";
							writeoutput.flush();
						}
						else {
							writeoutput << "\n";
							writeoutput.flush();
						}
					}
					writeoutput << "Content\t";
					writeoutput.flush();
					for (int j = 0; j < (signed)objMAT.content.size(); j++) {
						writeoutput << objMAT.content[j];
						writeoutput.flush();
						if (j != (signed)objMAT.content.size() - 1) {
							writeoutput << "\t";
							writeoutput.flush();
						}
						else {
							writeoutput << "\n";
							writeoutput.flush();
						}
					}
				}

				if (*FaceFilename != "") {
					
					std::string found_face = "";
					std::vector<long double> vface;
					for (int j = 0; j < (signed)FaceFeatureName->size(); j++) {
						vface.push_back(0.0);
					}

					long double find_max_content = -1000000000000000.0;
					for (int jj = 0; jj < (signed)objMAT.atomname.size(); jj++) {
						for (int j = 0; j < (signed)face->size(); j++) {
							if (objMAT.atomname[jj] == (*face)[j].AtomName && objMAT.content[jj] > find_max_content && (*face)[j].face) {
								found_face = (*face)[j].AtomName;
								find_max_content = objMAT.content[jj];
								for (int iv = 0; iv < (signed)FaceFeatureName->size(); iv++) {
									vface[iv] = (*face)[j].vFeature[iv];
								}
								break;
							}
						}
					}
					
					if (found_face != "") {

						writeoutput << "FACE_at\t" << found_face << "\n";
						writeoutput.flush();
						for (int iv = 0; iv < (signed)FaceFeatureName->size(); iv++) {
							writeoutput << (*FaceFeatureName)[iv];
							writeoutput.flush();
							if (iv != (signed)FaceFeatureName->size() - 1) {
								writeoutput << "\t";
								writeoutput.flush();
							}
							else {
								writeoutput << "\n";
								writeoutput.flush();
							}
						}
						for (int iv = 0; iv < (signed)vface.size(); iv++) {
							writeoutput << vface[iv] << "\t";
							writeoutput.flush();
						}

					}

					for (int iv = 0; iv < (signed)vface.size(); iv++) {
						writeanalysis << vface[iv] << "\t";
						writeanalysis.flush();
					}

					std::vector<long double>().swap(vface);

				}

				if (objMAT.istr != -1) {
					int notUlastcolumn = (signed)anadescriptor->size() - 1;
					for (int j = 0; j < (signed)anadescriptor->size(); j++) {
						if ((*anadescriptor)[j].ASB != "U") {
							notUlastcolumn = j;
						}
					}
					for (int j = 0; j < (signed)anadescriptor->size(); j++) {
						if ((*anadescriptor)[j].ASB != "U") {
							std::vector<long double> dv, mamssk;
							get_dv(&dv, &objMAT, atomdat, strdat, anadescriptor, &j);
							get_mamssk(&mamssk, &dv, VInfinite, &objMAT, DoNotAnalyze);
							if (*collect_only) {
								writeoutput << (*anadescriptor)[j].name << "\t";
								writeoutput.flush();
								for (int k = 0; k < (signed)dv.size(); k++) {
									writeoutput << dv[k];
									writeoutput.flush();
									if (k != (signed)dv.size() - 1) {
										writeoutput << "\t";
										writeoutput.flush();
									}
									else {
										writeoutput << "\t***\tMin=\t" << mamssk[0] << "\tAve=\t" << mamssk[1] << "\tMax=\t" << mamssk[2] << "\tStdDv=\t" << mamssk[3] << "\tSkew=\t" << mamssk[4] << "\tKurto=\t" << mamssk[5];
										writeoutput.flush();
										writeanalysis << mamssk[0] << "\t" << mamssk[1] << "\t" << mamssk[2] << "\t" << mamssk[3] << "\t" << mamssk[4] << "\t" << mamssk[5];
										writeanalysis.flush();
										if (j != notUlastcolumn) {
											writeanalysis << "\t";
											writeanalysis.flush();
										}
									}
								}
							}
							for (int k = 0; k < 6; k++) {
								xbox[j][k] = mamssk[k];
							}
							std::vector<long double>().swap(dv);
							std::vector<long double>().swap(mamssk);
						}
						else {
							long double v = 0.0;
							long double v0 = 0.0;
							if ((*anadescriptor)[j].UType_numerator_S0A1E2N3 == 0) {
								for (int ia = 0; ia < (signed)objMAT.content.size(); ia++) {
									for (int ja = 0; ja < (signed)(*anadescriptor)[j].UType_numerator_AtomKinds.size(); ja++) {
										if (objMAT.atomname[ia] == (*anadescriptor)[j].UType_numerator_AtomKinds[ja]) {
											v0 += objMAT.content[ia];
											break;
										}
									}
								}
							}
							else if ((*anadescriptor)[j].UType_numerator_S0A1E2N3 == 1) {
								for (int ia = 0; ia < (signed)objMAT.content.size(); ia++) {
									v0 += objMAT.content[ia];
								}
							}
							else if ((*anadescriptor)[j].UType_numerator_S0A1E2N3 == 2) {
								for (int ia = 0; ia < (signed)objMAT.content.size(); ia++) {
									for (int ja = 0; ja < (signed)(*anadescriptor)[j].UType_numerator_AtomKinds.size(); ja++) {
										if (objMAT.atomname[ia] != (*anadescriptor)[j].UType_numerator_AtomKinds[ja]) {
											v0 += objMAT.content[ia];
											break;
										}
									}
								}
							}
							else if ((*anadescriptor)[j].UType_numerator_S0A1E2N3 == 3) {
								v0 = (*anadescriptor)[j].UType_numerator_value;
							}
							long double v1 = 0.0;
							if ((*anadescriptor)[j].UType_denominator_S0A1E2N3 == 0) {
								for (int ia = 0; ia < (signed)objMAT.content.size(); ia++) {
									for (int ja = 0; ja < (signed)(*anadescriptor)[j].UType_denominator_AtomKinds.size(); ja++) {
										if (objMAT.atomname[ia] == (*anadescriptor)[j].UType_denominator_AtomKinds[ja]) {
											v1 += objMAT.content[ia];
											break;
										}
									}
								}
							}
							else if ((*anadescriptor)[j].UType_denominator_S0A1E2N3 == 1) {
								for (int ia = 0; ia < (signed)objMAT.content.size(); ia++) {
									v1 += objMAT.content[ia];
								}
							}
							else if ((*anadescriptor)[j].UType_denominator_S0A1E2N3 == 2) {
								for (int ia = 0; ia < (signed)objMAT.content.size(); ia++) {
									for (int ja = 0; ja < (signed)(*anadescriptor)[j].UType_denominator_AtomKinds.size(); ja++) {
										if (objMAT.atomname[ia] != (*anadescriptor)[j].UType_denominator_AtomKinds[ja]) {
											v1 += objMAT.content[ia];
											break;
										}
									}
								}
							}
							else if ((*anadescriptor)[j].UType_denominator_S0A1E2N3 == 3) {
								v1 = (*anadescriptor)[j].UType_denominator_value;
							}
							if ((*anadescriptor)[j].UType_m0d1 == 0) {
								v = v0 * v1;
							}
							else {
								v = v0 / v1;
							}
							if (*collect_only) {
								writeoutput << (*anadescriptor)[j].name << "\t" << v;
								writeoutput.flush();
								writeanalysis << "\t" << v;
								writeanalysis.flush();
							}
						}
						if (*collect_only) {
							writeoutput << "\n";
							writeoutput.flush();
						}
					}

					if (*collect_only) {
						if (PostTags.size() > 0) {
							for (int ipost = 0; ipost < (signed)PostTags.size(); ipost++) {
								writeanalysis << "\t" << PostTags[ipost];
								writeanalysis.flush();
							}
						}
						writeanalysis << "\n";
						writeanalysis.flush();
					}

					xbox_OM = 0.0;
					if (*AnionName != "-") {
						xbox_OM = objMAT.OM;
						for (int k = 0; k < (signed)CationNames->size(); k++) {
							xbox_CO[k] = objMAT.CO[k];
						}
					}
					if ((signed)regtype->size() > 0) {
						if (*collect_only) {

							writeoutput << "!!!\t";
							writeoutput.flush();
							if (*UseOffset) {
								lcstr >> sbuf;
							}
							for (int k = 0; k < (signed)regtype->size(); k++) {
								long double value = 0.0; 
								if (!*UseOffset) {
									value = calc_reg_chain(&(*regtype)[k], &objMAT, atomdat, strdat, anadescriptor, VInfinite, DoNotAnalyze, CationNames);
								}
								else {
									lcstr >> value;
									std::string target_ref_v_str = "";
									bool target_ref_v_found = false;
									for (int ipost = 0; ipost < (signed)PostTagNames->size(); ipost++) {
										if ((*PostTagNames)[ipost] == (*regtype)[k].name) {
											target_ref_v_str = PostTags[ipost];
											target_ref_v_found = true;
											break;
										}
									}
									if (target_ref_v_found) {
										std::stringstream glcstr;
										glcstr.str("");
										glcstr.clear();
										glcstr << target_ref_v_str;
										glcstr >> value;
									}
								}
								(*regtype)[k].v_list.push_back(value);
								writeoutput << (*regtype)[k].name << "=\t" << value;
								writeoutput.flush();
								if (k != (signed)regtype->size() - 1) {
									writeoutput << "\t";
									writeoutput.flush();
								}
								else {
									writeoutput << "\n";
									writeoutput.flush();
								}
							}

						}
						else {

							std::vector<long double> init_zvalue;
							std::vector<long double> monitor_value;
							std::vector<long double> init_monitor_value;
							std::vector<long double> offset;
							long double init_zdist = 0.0;
							writeoutput << "Init:\t";
							writeoutput.flush();
							pgstr = "Init:\t";
							if (*UseOffset) {
								lcstr >> sbuf;
							}
							for (int k = 0; k < (signed)regtype->size(); k++) {
								long double init_value = calc_reg_chain(&(*regtype)[k], &objMAT, atomdat, strdat, anadescriptor, VInfinite, DoNotAnalyze, CationNames);
								long double offset_value = 0.0;
								if (*UseOffset) {
									lcstr >> ldbuf;
									std::string target_ref_v_str = "";
									bool target_ref_v_found = false;
									for (int ipost = 0; ipost < (signed)PostTagNames->size(); ipost++) {
										if ((*PostTagNames)[ipost] == (*regtype)[k].name) {
											target_ref_v_str = PostTags[ipost];
											target_ref_v_found = true;
											break;
										}
									}
									if (target_ref_v_found) {
										std::stringstream glcstr;
										glcstr.str("");
										glcstr.clear();
										glcstr << target_ref_v_str;
										glcstr >> ldbuf;
									}
									offset_value = ldbuf - init_value;
								}
								long double inz = (init_value + offset_value - (*regtype)[k].v_ave) / (*regtype)[k].v_stddv;
								init_zvalue.push_back(inz);
								monitor_value.push_back(init_value + offset_value);
								init_monitor_value.push_back(init_value + offset_value);
								offset.push_back(offset_value);
								init_zdist += (*ZDistWeight)[k] * powl(inz - (*zTargetVReg)[k], 2.0);
							}
							init_zdist = sqrtl(init_zdist);
							writeoutput << "zDist=\t" << init_zdist << "\t";
							writeoutput.flush();
							pgstr += "zDist=\t" + std::to_string(init_zdist) + "\t";
							for (int k = 0; k < (signed)init_zvalue.size(); k++) {
								writeoutput << (*regtype)[k].name << "=\t" << monitor_value[k];
								writeoutput.flush();
								pgstr += (*regtype)[k].name + "=\t" + std::to_string(monitor_value[k]);
								if (k != (signed)init_zvalue.size() - 1) {
									writeoutput << "\t";
									writeoutput.flush();
									pgstr += "\t";
								}
								else {
									writeoutput << "\n";
									writeoutput.flush();
								}
							}

							int super_monitor = 0;
							MAT updMAT, startMAT;
							bool update = false;
							startMAT = objMAT;
							
							while (true) {

								for (int itest = 0; itest < 2; itest++) {

									int monitor = 0;

									while (true) {

										int target_t = -1;
										int target_tt = -1;
										long double min_upd_zdist = *VInfinite;
										bool success = false;
										MAT targetMAT = objMAT;

										for (int t = 0; t < (signed)objMAT.iatom.size(); t++) {

											int tt_size = 0;
											if (itest == 0) {
												tt_size = (signed)(*atomdat)[startMAT.iatom[t]].neighbor.size();
											}
											else {
												tt_size = (signed)ModulateContentFactor->size();
											}

											for (int tt = 0; tt < tt_size; tt++) {

												updMAT = objMAT;
												if (itest == 0) {
													updMAT.atomname[t] = (*atomdat)[startMAT.iatom[t]].neighbor[tt].atomname;
													updMAT.iatom[t] = (*atomdat)[startMAT.iatom[t]].neighbor[tt].iatom;
												}
												else {
													updMAT.content[t] = (*ModulateContentFactor)[tt] * startMAT.content[t];
												}
												designmate_iatom(&objMAT, atomdat);
												long double upd_CationSum = 0.0;
												long double* upd_scation_sum;
												upd_scation_sum = new long double[CationNames->size()];
												for (int upd_j = 0; upd_j < (signed)CationNames->size(); upd_j++) {
													upd_scation_sum[upd_j] = 0.0;
												}
												if (*AnionName != "-") {
													long double CationValenceTotal = 0.0;
													for (int ttt = 0; ttt < (signed)updMAT.iatom.size(); ttt++) {
														CationValenceTotal += updMAT.content[ttt] * (*atomdat)[updMAT.iatom[ttt]].v[CValence];
														upd_CationSum += updMAT.content[ttt];
													}
													updMAT.anion_content = -CationValenceTotal / *AnionValence;
													for (int ttt = 0; ttt < (signed)updMAT.atomname.size(); ttt++) {
														for (int upd_k = 0; upd_k < (signed)CationNames->size(); upd_k++) {
															if (updMAT.atomname[ttt] == (*CationNames)[upd_k]) {
																upd_scation_sum[upd_k] += updMAT.content[ttt];
																break;
															}
														}
													}
												}
												else {
													updMAT.OM = 0.0;
												}

												updMAT.composition = "";
												for (int ttt = 0; ttt < (signed)updMAT.atomname.size(); ttt++) {
													updMAT.composition += updMAT.atomname[ttt];
													if (updMAT.content[ttt] != 1.0) {
														updMAT.composition += std::to_string(updMAT.content[ttt]).substr(0, 5);
													}
												}
												if (*AnionName != "-") {
													updMAT.composition += *AnionName;
													if (updMAT.anion_content != 1.0) {
														updMAT.composition += std::to_string(updMAT.anion_content).substr(0, 5);
													}
													updMAT.OM = updMAT.anion_content / upd_CationSum;
													updMAT.CO.clear();
													for (int upd_k = 0; upd_k < (signed)CationNames->size(); upd_k++) {
														updMAT.CO.push_back(upd_scation_sum[upd_k] / updMAT.anion_content);
													}
												}
												else {
													updMAT.OM = 0.0;
												}

												std::vector<long double> dv, mamssk;
												for (int j = 0; j < (signed)anadescriptor->size(); j++) {
													get_dv(&dv, &updMAT, atomdat, strdat, anadescriptor, &j);
													get_mamssk(&mamssk, &dv, VInfinite, &updMAT, DoNotAnalyze);
													for (int k = 0; k < 6; k++) {
														xbox[j][k] = mamssk[k];
													}
													std::vector<long double>().swap(dv);
													std::vector<long double>().swap(mamssk);
												}
												xbox_OM = 0.0;
												if (*AnionName != "-") {
													xbox_OM = updMAT.OM;
													for (int k = 0; k < (signed)CationNames->size(); k++) {
														xbox_CO[k] = updMAT.CO[k];
													}
												}
												std::vector<long double> upd_zvalue;
												long double upd_zdist = 0.0;
												for (int k = 0; k < (signed)regtype->size(); k++) {
													long double upd_value = calc_reg_chain(&(*regtype)[k], &updMAT, atomdat, strdat, anadescriptor, VInfinite, DoNotAnalyze, CationNames);
													monitor_value[k] = upd_value + offset[k];
													long double upz = (upd_value + offset[k] - (*regtype)[k].v_ave) / (*regtype)[k].v_stddv;
													upd_zvalue.push_back(upz);
													upd_zdist += (*ZDistWeight)[k] * powl(upz - (*zTargetVReg)[k], 2.0);
												}
												std::vector<long double>().swap(upd_zvalue);
												upd_zdist = sqrtl(upd_zdist);
												if (upd_zdist < min_upd_zdist && upd_zdist < init_zdist) {
													min_upd_zdist = upd_zdist;
													target_t = t;
													target_tt = tt;
													targetMAT = updMAT;
													success = true;
													update = true;
												}
											}
										}

										if (success) {
											objMAT = targetMAT;
											init_zdist = min_upd_zdist;
											if (min_upd_zdist < global_min_zdist) {
												global_min_zdist = min_upd_zdist;
												gMAT = targetMAT;
												for (int k = 0; k < (signed)regtype->size(); k++) {
													gmonitor_value[k] = monitor_value[k];
												}
												gstr = std::to_string(i + 1) + "\t...\t" + startMAT.name + "\t...\t" + startMAT.composition + "\t...\t" + startMAT.strname;
												gstr += "\n" + pgstr;
											}
											write_zdist(&objMAT, &min_upd_zdist, AnionName);
											writeoutput << "\t";
											writeoutput.flush();
											for (int k = 0; k < (signed)regtype->size(); k++) {
												writeoutput << (*regtype)[k].name << "=\t" << monitor_value[k] << "\t";
												writeoutput.flush();
												if (k != (signed)regtype->size() - 1) {
													writeoutput << "\t";
													writeoutput.flush();
												}
												else {
													writeoutput << "\n";
													writeoutput.flush();
												}
											}
										}

										monitor++;
										if (monitor == *MaxIteration || !success) {
											break;
										}

									}

								}

								super_monitor++;
								if (super_monitor == *MaxIteration || !update) {
									break;
								}

							}
							
							if (update) {
								writeoutput << "@@@\t";
								writeoutput.flush();
								write_zdist(&objMAT, &init_zdist, AnionName);
								writeoutput << "\t";
								writeoutput.flush();
								for (int k = 0; k < (signed)regtype->size(); k++) {
									writeoutput << (*regtype)[k].name << "=\t" << monitor_value[k] << "\t";
									writeoutput.flush();
									if (k != (signed)regtype->size() - 1) {
										writeoutput << "\t";
										writeoutput.flush();
									}
									else {
										writeoutput << "\n";
										writeoutput.flush();
									}
								}
							}
							else {
								if (init_zdist < global_min_zdist) {
									global_min_zdist = init_zdist;
									gMAT = startMAT;
									for (int k = 0; k < (signed)regtype->size(); k++) {
										gmonitor_value[k] = init_monitor_value[k];
									}
									gstr = std::to_string(i + 1) + "\t...\t" + startMAT.name + "\t...\t" + startMAT.composition + "\t...\t" + startMAT.strname;
									gstr += "\n" + pgstr;
								}
							}

							std::vector<long double>().swap(init_zvalue);
							std::vector<long double>().swap(monitor_value);
							std::vector<long double>().swap(init_monitor_value);
							std::vector<long double>().swap(offset);

						}
					}
					writeoutput << "=================================================================================\n";
					writeoutput.flush();
				}
				else {
					std::cout << "UseStrName was not identified ... Check StructureListFilename or UseStrName...\n";
					std::cout.flush();
				}

				std::vector<std::string>().swap(PostTags);

			}

			if (i == Nmat - 1) {
				break;
			}

		}

		if (!*collect_only) {
			writeoutput << gstr << "\n@@@\t";
			writeoutput.flush();
			write_zdist(&gMAT, &global_min_zdist, AnionName);
			writeoutput << "\t";
			writeoutput.flush();
			for (int k = 0; k < (signed)regtype->size(); k++) {
				writeoutput << (*regtype)[k].name << "=\t" << gmonitor_value[k] << "\t";
				writeoutput.flush();
				if (k != (signed)regtype->size() - 1) {
					writeoutput << "\t";
					writeoutput.flush();
				}
				else {
					writeoutput << "\n";
					writeoutput.flush();
				}
			}
		}

		readmatdat.close();
	
	}

	if (*collect_only && (signed)regtype->size() > 0) {
		for (int k = 0; k < (signed)regtype->size(); k++) {
			(*regtype)[k].v_ave = 0.0;
			(*regtype)[k].v_stddv = 0.0;
			int nancount = 0;
			for (int kk = 0; kk < (signed)(*regtype)[k].v_list.size(); kk++) {
				if (!std::isnan((*regtype)[k].v_list[kk])) {
					(*regtype)[k].v_ave += (*regtype)[k].v_list[kk];
				}
				else {
					nancount++;
				}
			}
			int Nsample = (signed)(*regtype)[k].v_list.size() - nancount;
			(*regtype)[k].v_ave /= (long double)Nsample;
			for (int kk = 0; kk < (signed)(*regtype)[k].v_list.size(); kk++) {
				if (!std::isnan((*regtype)[k].v_list[kk])) {
					(*regtype)[k].v_stddv += powl((*regtype)[k].v_list[kk] - (*regtype)[k].v_ave, 2.0);
				}
			}
			(*regtype)[k].v_stddv = sqrtl((*regtype)[k].v_stddv / (long double)Nsample);
			zTargetVReg->push_back(((*TargetVReg)[k] - (*regtype)[k].v_ave) / (*regtype)[k].v_stddv);
			writeoutput << "___\t" << (*regtype)[k].name << "\t:\tAve=\t" << (*regtype)[k].v_ave << "\tStdDv=\t" << (*regtype)[k].v_stddv << "\tzTarget=\t" << (*zTargetVReg)[k] << "\n";
			writeoutput.flush();
		}
		writeoutput << "=================================================================================\n";
		writeoutput.flush();
	}

}

void process_mapout(std::vector<long double>* zTargetVReg, std::vector<long double>* ZDistWeight, std::vector<long double>* TargetVReg, std::string* MapoutFilename, std::vector<REGTYPE>* regtype, std::vector<ATOMDAT>* atomdat, std::vector<STRDAT>* strdat, std::vector<ANADESCRIPTOR>* anadescriptor, long double* VInfinite, std::vector<std::string>* DoNotAnalyze, std::vector<std::string>* CationNames) {

	readmapout.open(*MapoutFilename);

	if (readmapout.is_open()) {

		readmapout >> sbuf >> mapout_preNRandom >> mapout_NRandom;
		readmapout >> sbuf >> mapout_Verbose;
		std::getline(readmapout, sbuf);
		std::getline(readmapout, sbuf);
		int system_count = 0;
		long double global_min;
		global_min = (*VInfinite);
		std::string global_min_str;
		SMATLIST objSMATLIST;
		std::vector<long double> ave, stdev;
		for (int ifunc = 0; ifunc < (signed)regtype->size(); ifunc++) {
			ave.push_back(0.0);
			stdev.push_back(0.0);
		}
		int max_comp = -1;

		int sample_count = 0;
		writeoutput << "Gathering ave vector...\n";
		writeoutput.flush();
		while (true) {

			std::getline(readmapout, sbuf);
			if (sbuf == "STOP" || readmapout.eof()) {
				break;
			}

			MAT objMAT;
			std::stringstream lcstr;
			lcstr << sbuf;
			if (sbuf.find("*") == std::string::npos) {
				while (true) {
					lcstr >> sbuf;
					objMAT.atomname.push_back(sbuf);
					objMAT.content.push_back(0.0);
					if (lcstr.eof()) {
						break;
					}
				}
			}
			else {
				lcstr >> sbuf;
				while (true) {
					long double ldbuf;
					lcstr >> sbuf >> ldbuf;
					objMAT.atomname.push_back(sbuf);
					objMAT.content.push_back(ldbuf);
					if (lcstr.eof()) {
						break;
					}
				}
			}
			
			system_count++;
			designmate_iatom(&objMAT, atomdat);
			if ((signed)objMAT.content.size() > max_comp) {
				max_comp = (signed)objMAT.content.size();
			}

			for (int i = 0; i < (signed)objMAT.content.size() + (signed)objMAT.content.size() * mapout_preNRandom + mapout_preNRandom; i++) {

				if (i < objMAT.content.size()) {
					for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
						if (ii == i) {
							objMAT.content[ii] = 1.0;
						}
						else {
							objMAT.content[ii] = 0.0;
						}
					}
				}
				else if (objMAT.content.size() > 2 && ((signed)objMAT.content.size() + (signed)objMAT.content.size() * mapout_preNRandom)) {
					int isetzero = (i - (signed)objMAT.content.size()) / mapout_preNRandom;
					if (isetzero < 0) {
						isetzero = 0;
					}
					else if (isetzero >= (signed)objMAT.content.size()) {
						isetzero = (signed)objMAT.content.size() - 1;
					}
					long double norm = 0.0;
					for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
						if (ii != isetzero) {
							objMAT.content[ii] = folder((long double)rand() / RAND_MAX);
						}
						else {
							objMAT.content[ii] = 0.0;
						}
						norm += objMAT.content[ii];
					}
					for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
						objMAT.content[ii] /= norm;
					}
				}
				else {
					long double norm = 0.0;
					for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
						objMAT.content[ii] = folder((long double)rand() / RAND_MAX);
						norm += objMAT.content[ii];
					}
					for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
						objMAT.content[ii] /= norm;
					}
				}

				for (int j = 0; j < (signed)anadescriptor->size(); j++) {
					std::vector<long double> dv, mamssk;
					get_dv(&dv, &objMAT, atomdat, strdat, anadescriptor, &j);
					get_mamssk(&mamssk, &dv, VInfinite, &objMAT, DoNotAnalyze);
					for (int k = 0; k < 6; k++) {
						xbox[j][k] = mamssk[k];
					}
					xbox_OM = 0.0;
					std::vector<long double>().swap(dv);
					std::vector<long double>().swap(mamssk);
				}

				sample_count++;
				for (int k = 0; k < (signed)regtype->size(); k++) {
					ave[k] += calc_reg_chain(&(*regtype)[k], &objMAT, atomdat, strdat, anadescriptor, VInfinite, DoNotAnalyze, CationNames);
				}

			}

		}

		for (int ifunc = 0; ifunc < (signed)regtype->size(); ifunc++) {
			ave[ifunc] /= (long double)sample_count;
			writeoutput << (*regtype)[ifunc].name << "\t[AVE]\t" << ave[ifunc] << "\n";
			writeoutput.flush();
		}

		readmapout.close();
		readmapout.open(*MapoutFilename);
		for (int iline = 0; iline < 3; iline++) {
			std::getline(readmapout, sbuf);
		}

		writeoutput << "Gathering stdev vector...\n";
		writeoutput.flush();
		for (int isystem = 0; isystem < system_count; isystem++) {

			std::getline(readmapout, sbuf);
			MAT objMAT;
			std::stringstream lcstr;
			lcstr << sbuf;
			if (sbuf.find("*") == std::string::npos) {
				while (true) {
					lcstr >> sbuf;
					objMAT.atomname.push_back(sbuf);
					objMAT.content.push_back(0.0);
					if (lcstr.eof()) {
						break;
					}
				}
			}
			else {
				lcstr >> sbuf;
				while (true) {
					long double ldbuf;
					lcstr >> sbuf >> ldbuf;
					objMAT.atomname.push_back(sbuf);
					objMAT.content.push_back(ldbuf);
					if (lcstr.eof()) {
						break;
					}
				}
			}

			designmate_iatom(&objMAT, atomdat);

			for (int i = 0; i < (signed)objMAT.content.size() + (signed)objMAT.content.size() * mapout_preNRandom + mapout_preNRandom; i++) {

				if (i < objMAT.content.size()) {
					for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
						if (ii == i) {
							objMAT.content[ii] = 1.0;
						}
						else {
							objMAT.content[ii] = 0.0;
						}
					}
				}
				else if (objMAT.content.size() > 2 && ((signed)objMAT.content.size() + (signed)objMAT.content.size() * mapout_preNRandom)) {
					int isetzero = (i - (signed)objMAT.content.size()) / mapout_preNRandom;
					if (isetzero < 0) {
						isetzero = 0;
					}
					else if (isetzero >= (signed)objMAT.content.size()) {
						isetzero = (signed)objMAT.content.size() - 1;
					}
					long double norm = 0.0;
					for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
						if (ii != isetzero) {
							objMAT.content[ii] = folder((long double)rand() / RAND_MAX);
						}
						else {
							objMAT.content[ii] = 0.0;
						}
						norm += objMAT.content[ii];
					}
					for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
						objMAT.content[ii] /= norm;
					}
				}
				else {
					long double norm = 0.0;
					for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
						objMAT.content[ii] = folder((long double)rand() / RAND_MAX);
						norm += objMAT.content[ii];
					}
					for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
						objMAT.content[ii] /= norm;
					}
				}

				for (int j = 0; j < (signed)anadescriptor->size(); j++) {
					std::vector<long double> dv, mamssk;
					get_dv(&dv, &objMAT, atomdat, strdat, anadescriptor, &j);
					get_mamssk(&mamssk, &dv, VInfinite, &objMAT, DoNotAnalyze);
					for (int k = 0; k < 6; k++) {
						xbox[j][k] = mamssk[k];
					}
					xbox_OM = 0.0;
					std::vector<long double>().swap(dv);
					std::vector<long double>().swap(mamssk);
				}

				for (int k = 0; k < (signed)regtype->size(); k++) {
					stdev[k] += powl(calc_reg_chain(&(*regtype)[k], &objMAT, atomdat, strdat, anadescriptor, VInfinite, DoNotAnalyze, CationNames) - ave[k], 2.0);
				}

			}

		}

		for (int ifunc = 0; ifunc < (signed)regtype->size(); ifunc++) {
			stdev[ifunc] /= (long double)sample_count;
			stdev[ifunc] = sqrtl(stdev[ifunc]);
			zTargetVReg->push_back(((*TargetVReg)[ifunc] - ave[ifunc])/stdev[ifunc]);
			writeoutput << (*regtype)[ifunc].name << "\t[STDEV]\t" << stdev[ifunc] << "\t[TargetVReg]\t" << (*TargetVReg)[ifunc] << "\t[zTargetVReg]\t" << (*zTargetVReg)[ifunc] << "\n";
			writeoutput.flush();
		}

		readmapout.close();
		readmapout.open(*MapoutFilename);
		for (int iline = 0; iline < 3; iline++) {
			std::getline(readmapout, sbuf);
		}

		for (int isystem = 0; isystem < system_count; isystem++) {

			bool specific_content = false;
			std::getline(readmapout, sbuf);
			MAT objMAT;
			std::stringstream lcstr;
			lcstr << sbuf;
			if (sbuf.find("*") == std::string::npos) {
				while (true) {
					lcstr >> sbuf;
					objMAT.atomname.push_back(sbuf);
					objMAT.content.push_back(0.0);
					if (lcstr.eof()) {
						break;
					}
				}
			}
			else {
				specific_content = true;
				lcstr >> sbuf;
				while (true) {
					long double ldbuf;
					lcstr >> sbuf >> ldbuf;
					objMAT.atomname.push_back(sbuf);
					objMAT.content.push_back(ldbuf);
					if (lcstr.eof()) {
						break;
					}
				}
			}

			designmate_iatom(&objMAT, atomdat);

			long double system_min = (*VInfinite);
			long double* v_at_system_min;
			long double* z_at_system_min;
			v_at_system_min = new long double[regtype->size()];
			z_at_system_min = new long double[regtype->size()];
			for (int j = 0; j < (signed)regtype->size(); j++) {
				v_at_system_min[j] = 0.0;
				z_at_system_min[j] = 0.0;
			}
			long double* min_content;
			min_content = new long double [(signed)objMAT.content.size()];
			for (int j = 0; j < (signed)objMAT.content.size(); j++) {
				min_content[j] = 0.0;
			}

			int upto = 1;
			if (!specific_content) {
				upto = (signed)objMAT.content.size() + (signed)objMAT.content.size() * mapout_NRandom + mapout_NRandom;
			}

			for (int i = 0; i < upto; i++) {

				if (!specific_content) {
					if (i < objMAT.content.size()) {
						for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
							if (ii == i) {
								objMAT.content[ii] = 1.0;
							}
							else {
								objMAT.content[ii] = 0.0;
							}
						}
					}
					else if (objMAT.content.size() > 2 && (i < (signed)objMAT.content.size() + (signed)objMAT.content.size() * mapout_NRandom)) {
						int isetzero = (i - (signed)objMAT.content.size()) / mapout_NRandom;
						if (isetzero < 0) {
							isetzero = 0;
						}
						else if (isetzero >= (signed)objMAT.content.size()) {
							isetzero = (signed)objMAT.content.size() - 1;
						}
						long double norm = 0.0;
						for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
							if (ii != isetzero) {
								objMAT.content[ii] = folder((long double)rand() / RAND_MAX);
							}
							else {
								objMAT.content[ii] = 0.0;
							}
							norm += objMAT.content[ii];
						}
						for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
							objMAT.content[ii] /= norm;
						}
					}
					else {
						long double norm = 0.0;
						for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
							objMAT.content[ii] = folder((long double)rand() / RAND_MAX);
							norm += objMAT.content[ii];
						}
						for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
							objMAT.content[ii] /= norm;
						}
					}
				}
				
				for (int j = 0; j < (signed)anadescriptor->size(); j++) {
					std::vector<long double> dv, mamssk;
					get_dv(&dv, &objMAT, atomdat, strdat, anadescriptor, &j);
					get_mamssk(&mamssk, &dv, VInfinite, &objMAT, DoNotAnalyze);
					for (int k = 0; k < 6; k++) {
						xbox[j][k] = mamssk[k];
					}
					xbox_OM = 0.0;
					std::vector<long double>().swap(dv);
					std::vector<long double>().swap(mamssk);
				}

				std::vector<long double> v, z;
				long double z_dist = 0.0;
				for (int k = 0; k < (signed)regtype->size(); k++) {
					long double v_loc = calc_reg_chain(&(*regtype)[k], &objMAT, atomdat, strdat, anadescriptor, VInfinite, DoNotAnalyze, CationNames);
					long double z_loc = (v_loc - ave[k]) / stdev[k];
					v.push_back(v_loc);
					z.push_back(z_loc);
					z_dist += (*ZDistWeight)[k] * powl(z_loc - (*zTargetVReg)[k], 2.0);
				}

				if (mapout_Verbose == 1) {
					writeoutput << isystem + 1 << "\t" << z_dist;
					writeoutput.flush();
					int finite_comp = 0;
					for (int iatom = 0; iatom < (signed)objMAT.content.size(); iatom++) {
						if (objMAT.content[iatom] != 0.0) {
							finite_comp++;
							writeoutput << "\t" << objMAT.atomname[iatom] << "\t" << objMAT.content[iatom];
							writeoutput.flush();
						}
					}
					for (int vatom = 0; vatom < max_comp - finite_comp; vatom++) {
						writeoutput << "\t\t";
						writeoutput.flush();
					}
					writeoutput << "\t...\t";
					writeoutput.flush();
					for (int k = 0; k < (signed)regtype->size(); k++) {
						writeoutput << (*regtype)[k].name << "\t[v]\t" << v[k] << "\t[z]\t" << z[k];
						writeoutput.flush();
						if (k != (signed)regtype->size() - 1) {
							writeoutput << "\t";
							writeoutput.flush();
						}
						else {
							writeoutput << "\n";
							writeoutput.flush();
						}
					}
				}

				if (z_dist < system_min) {
					system_min = z_dist;
					for (int iatom = 0; iatom < (signed)objMAT.content.size(); iatom++) {
						min_content[iatom] = objMAT.content[iatom];
					}
					for (int k = 0; k < (signed)regtype->size(); k++) {
						v_at_system_min[k] = v[k];
						z_at_system_min[k] = z[k];
					}
				}

				std::vector<long double>().swap(v);
				std::vector<long double>().swap(z);

			}

			
			SMAT objSMATmin, objSMATmax;
			for (int iatom = 0; iatom < (signed)objMAT.content.size(); iatom++) {
				if (min_content[iatom] != 0.0) {
					objSMATmin.atomname.push_back(objMAT.atomname[iatom]);
					objSMATmin.content.push_back(min_content[iatom]);
				}
			}
			bool overlapmin = false;
			for (int ilist = (signed)objSMATLIST.smat_min.size() - 1; ilist >= 0; ilist--) {
				bool is_same = true;
				if (objSMATLIST.smat_min[ilist].atomname.size() == objSMATmin.atomname.size()) {
					for (int iele = 0; iele < (signed)objSMATmin.atomname.size(); iele++) {
						if (objSMATLIST.smat_min[ilist].atomname[iele] != objSMATmin.atomname[iele]) {
							is_same = false;
							break;
						}
					}
				}
				else {
					is_same = false;
				}
				if (is_same) {
					overlapmin = true;
					break;
				}
			}
			if (!overlapmin) {
				objSMATLIST.smat_min.push_back(objSMATmin);
				writeoutput << isystem + 1 << "\t" << system_min;
				writeoutput.flush();
				int finite_comp = 0;
				for (int iatom = 0; iatom < (signed)objMAT.content.size(); iatom++) {
					if (min_content[iatom] != 0.0) {
						finite_comp++;
						writeoutput << "\t" << objMAT.atomname[iatom] << "\t" << min_content[iatom];
						writeoutput.flush();
					}
				}
				for (int vatom = 0; vatom < max_comp - finite_comp; vatom++) {
					writeoutput << "\t\t";
					writeoutput.flush();
				}
				writeoutput << "\t...\t";
				writeoutput.flush();
				for (int k = 0; k < (signed)regtype->size(); k++) {
					writeoutput << (*regtype)[k].name << "\t[v]\t" << v_at_system_min[k] << "\t[z]\t" << z_at_system_min[k];
					writeoutput.flush();
					if (k != (signed)regtype->size() - 1) {
						writeoutput << "\t";
						writeoutput.flush();
					}
					else {
						writeoutput << "\n";
						writeoutput.flush();
					}
				}
			}
			
			if (system_min < global_min) {
				global_min = system_min;
				std::stringstream globalstr;
				globalstr << "GlobalMIN:\t" << system_min << "\n(GlobalMIN_at):\t";
				for (int iatom = 0; iatom < (signed)objMAT.content.size(); iatom++) {
					if (min_content[iatom] != 0.0) {
						globalstr << "\t" << objMAT.atomname[iatom] << "\t" << min_content[iatom] << "\t";
					}
				}
				globalstr << "...\t";
				for (int k = 0; k < (signed)regtype->size(); k++) {
					globalstr << (*regtype)[k].name << "\t[v]\t" << v_at_system_min[k] << "\t[z]\t" << z_at_system_min[k];
					if (k != (signed)regtype->size() - 1) {
						globalstr << "\t";
					}
				}
				global_min_str = globalstr.str();

			}
			
			std::vector<std::string>().swap(objMAT.atomname);
			std::vector<long double>().swap(objMAT.content);

			delete[] v_at_system_min;
			delete[] z_at_system_min;
			delete[] min_content;
			v_at_system_min = nullptr;
			z_at_system_min = nullptr;
			min_content = nullptr;

		}

		writeoutput << "****************************************************************\n";
		writeoutput.flush();
		writeoutput << global_min_str << "\n";
		writeoutput.flush();
		writeoutput << "****************************************************************\n";
		writeoutput.flush();

		std::vector<long double>().swap(ave);
		std::vector<long double>().swap(stdev);

	}

	readmapout.close();

}

int main() {

	read_para();
	if (FaceFilename != "-") {
		read_face(&face, &FaceFeatureName, &FaceFilename);
	}

	writeoutput.open(OutputFilename, std::ios::binary);
	writeoutput << std::fixed << std::setprecision(20);
	writeanalysis.open(AnalysisOutputFilename, std::ios::binary);
	writeanalysis << std::fixed << std::setprecision(20);
	writeanalysis << "No.\tName\tComposition\tStrcture\t";
	writeanalysis.flush();
	if (AnionName != "-") {
		writeanalysis << AnionName << "M\t";
		writeanalysis.flush();
		for (int k = 0; k < (signed)CationNames.size(); k++) {
			writeanalysis << CationNames[k] << AnionName << "\t";
			writeanalysis.flush();
		}
	}
	if (FaceFilename != "-") {
		for (int j = 0; j < (signed)FaceFeatureName.size(); j++) {
			writeanalysis << FaceFeatureName[j] << "\t";
			writeanalysis.flush();
		}
	}
	for (int j = 0; j < (signed)anadescriptor.size(); j++) {
		if (anadescriptor[j].ASB != "U") {
			writeanalysis << "Min" << anadescriptor[j].name << "\t";
			writeanalysis.flush();
			writeanalysis << "Ave" << anadescriptor[j].name << "\t";
			writeanalysis.flush();
			writeanalysis << "Max" << anadescriptor[j].name << "\t";
			writeanalysis.flush();
			writeanalysis << "StdDv" << anadescriptor[j].name << "\t";
			writeanalysis.flush();
			writeanalysis << "Skew" << anadescriptor[j].name << "\t";
			writeanalysis.flush();
			writeanalysis << "Kurto" << anadescriptor[j].name;
			writeanalysis.flush();
		}
		else {
			writeanalysis << anadescriptor[j].name;
			writeanalysis.flush();
		}
		if(j != (signed)anadescriptor.size() - 1) {
			writeanalysis << "\t";
			writeanalysis.flush();
		}
	}
	if (PostTagNames.size() > 0) {
		for (int ipost = 0; ipost < (signed)PostTagNames.size(); ipost++) {
			writeanalysis << "\t" << PostTagNames[ipost];
			writeanalysis.flush();
		}
	}
	writeanalysis << "\n";
	writeanalysis.flush();

	read_atomdat(&atomdat, &descriptor, &anadescriptor, &AtomDatListFilename);
	pick_neighbor(&atomdat, &NeighborSelectSpace, &NeighborZCrit, &ExcludeNeighbors, &DoNotFindNeighbors, &descriptor);
	bool usestrname = true;
	if (UseStrName == 0) {
		usestrname = false;
	}
	read_strdat(&strdat, &anadescriptor, &sdescriptorname, &sitename, &StructureListFilename, &descriptor, &atomdat, &usestrname);
	for (int i = 0; i < (signed)RegModelFilename.size(); i++) {
		//read_reg(&functype, &RegModelFilename[i], &anadescriptor);
		regtype.push_back(read_reg_chain(&RegModelFilename[i], &anadescriptor, &AnionName, &CationNames));
	}

	bool bUseOffset = false;
	if (UseOffset == 1) {
		bUseOffset = true;
	}

	if (!mapout) {
		bool headexists = true;
		if (HeadExists == 0) {
			headexists = false;
		}
		set_par(&MaterialsFilename, &headexists);
		bool collect_only = true;
		process_matdat(&zTargetVReg, &ZDistWeight, &TargetVReg, &ModulateContentFactor, &MaterialsFilename, &regtype, &atomdat, &FaceFilename, &FaceFeatureName, &face, &strdat, &sitename, &anadescriptor, &AnionName, &CationNames, &AnionValence, &VInfinite, &collect_only, &bUseOffset, &MaxIteration, &headexists, &usestrname, &PostTagNames, &DoNotAnalyze);
		if ((signed)regtype.size() > 0) {
			collect_only = false;
			process_matdat(&zTargetVReg, &ZDistWeight, &TargetVReg, &ModulateContentFactor, &MaterialsFilename, &regtype, &atomdat, &FaceFilename, &FaceFeatureName, &face, &strdat, &sitename, &anadescriptor, &AnionName, &CationNames, &AnionValence, &VInfinite, &collect_only, &bUseOffset, &MaxIteration, &headexists, &usestrname, &PostTagNames, &DoNotAnalyze);
		}
	}
	else {
		process_mapout(&zTargetVReg, &ZDistWeight, &TargetVReg, &MapOutSpaceFilename, &regtype, &atomdat, &strdat, &anadescriptor, &VInfinite, &DoNotAnalyze, &CationNames);
	}
	
	writeoutput.close();
	writeanalysis.close();

	return 0;

}