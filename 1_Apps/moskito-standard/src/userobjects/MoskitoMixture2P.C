/**************************************************************************/
/*  MOSKITO - Multiphysics cOupled Simulator toolKIT for wellbOres        */
/*                                                                        */
/*  Copyright (C) 2017 by Maziar Gholami Korzani                          */
/*  Karlsruhe Institute of Technology, Institute of Applied Geosciences   */
/*  Division of Geothermal Research                                       */
/*                                                                        */
/*  This file is part of MOSKITO App                                      */
/*                                                                        */
/*  This program is free software: you can redistribute it and/or modify  */
/*  it under the terms of the GNU General Public License as published by  */
/*  the Free Software Foundation, either version 3 of the License, or     */
/*  (at your option) any later version.                                   */
/*                                                                        */
/*  This program is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          */
/*  GNU General Public License for more details.                          */
/*                                                                        */
/*  You should have received a copy of the GNU General Public License     */
/*  along with this program.  If not, see <http://www.gnu.org/licenses/>  */
/**************************************************************************/

#include "MoskitoMixture2P.h"

registerMooseObject("MoskitoApp", MoskitoMixture2P);

// template <>
// InputParameters
// validParams<MoskitoMixture2P>()
// {
//   InputParameters params = validParams<MoskitoEOS2P>();

InputParameters
MoskitoMixture2P::validParams()
{
  InputParameters params = MoskitoEOS2P::validParams();

  return params;
}

MoskitoMixture2P::MoskitoMixture2P(const InputParameters & parameters)
  : MoskitoEOS2P(parameters)
{
}

void
MoskitoMixture2P::VMFrac_T_from_p_h(
  const Real & press, const Real & enthalpy, Real & vmfrac, Real & temp, Real & phase) const
{
  vmfrac = 0.2;
  temp = 333.15;
  phase = 1.0;

}

void
MoskitoMixture2P::rho_g_from_p_T(
  const Real &  press, const Real &  temp, Real & rho, Real & drho_dp, Real & drho_dT, const Real & phase) const
{
  rho = 1.0;
  drho_dp = 0.000001;
  drho_dT = 0.000001;
}

void
MoskitoMixture2P::rho_l_from_p_T(
  const Real &  press, const Real &  temp, Real & rho, Real & drho_dp, Real & drho_dT, const Real & phase) const
{
  rho = 1000.0;
  drho_dp = 0.000001;
  drho_dT = 0.000001;

}

Real
MoskitoMixture2P::rho_g_from_p_T(const Real & press, const Real & temp, const Real & phase) const
{
  Real rho = 1.0;

  return rho;
}

Real
MoskitoMixture2P::rho_l_from_p_T(const Real &  press, const Real &  temp, const Real & phase) const
{
  Real rho = 1000.0;

  return rho;
}

Real
MoskitoMixture2P::cp_m_from_p_T(
      const Real & press, const Real & temp, const Real & vmfrac, const Real & phase) const
{
  Real cp = 4000.0;

  return cp;
}

Real
MoskitoMixture2P::rho_m_from_p_h(const Real & press, const Real & enthalpy) const
{
  Real rho,temp,phase,vmfrac;

  return 800.0;
}

void
MoskitoMixture2P::h_lat(const Real & press, Real & hlat, Real & hsatl, Real & hsatg) const
{

}
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////

//Real  t;// omf[10], m_mas[10], omf_sum, tmf[6], sm[5], k[6] ;

//t = 333.15;

/*
omf[1] = 0.7;
omf[2] = 0.0;
omf[3] = 0.0;
omf[4] = 0.0;
omf[5] = 0.3;
omf[6] = 0.0;
omf[7] = 0.0;
omf[8] = 0.0;
omf[9] = 0.0;

omf_sum = 0.0;
for(kk=1;kk<=5;kk++)
{
  omf_sum += omf[kk];
}


for(kk=1;kk<=5;kk++)
{
  tmf[kk] = omf[kk]/omf_sum;
}

m_mas[1] = 18.01528e-3;
m_mas[2] = 16.04e-3;
m_mas[3] = 44.01e-3;
m_mas[4] = 28.0134e-3;
m_mas[5] = 34.1e-3;
m_mas[6] = 58.44e-3;
m_mas[7] = 74.5513e-3;
m_mas[8] = 110.98e-3;
m_mas[9] = 95.211e-3;

sm[1] = 0.0;
sm[2] = 0.0;
sm[3] = 0.0;
sm[4] = 0.1;

k[1] = 0.0015;
k[2] = 650.0;
k[3] = 40.0;
k[4] = 1250.0;
k[5] = 10.0;
*/

void
MoskitoMixture2P::PENG_PURE(const Real & press, const Real & temp, Real tc[], Real pc[], Real m[], Real a[], Real b[]) const
{
  int i;
	double w[6], tr[6], squred_part[6], R;

  R = 8.31446;

	w[1] = 0.345;
	w[2] = 0.011;
	w[3] = 0.228;
	w[4] = 0.0403;
	w[5] = 0.1081;

	pc[1] = 22055000;
	pc[2] = 4604000;
	pc[3] = 7382000;
	pc[4] = 3400000;
	pc[5] = 8940000;

	tc[1] = 647.13;
	tc[2] = 190.58;
	tc[3] = 304.19;
	tc[4] = 126.1;
	tc[5] = 373.2;


for(i=1 ; i<=5 ; i++)
{
  m[i] = 0.37464 + 1.54226*w[i] - 0.26992*pow(w[i],2.0);
  b[i] = 0.07780*R*tc[i]/pc[i];
  tr[i] = temp/tc[i];
  squred_part[i] = 1.0+m[i]*(1.0-sqrt(tr[i]));
  a[i] = 0.45724*pow(R,2.0)*pow(tc[i],2.0)/pc[i]*pow(squred_part[i],2.0);
}

}

void
MoskitoMixture2P::PENG_MIX(const Real & press, const Real & temp, Real y[], const Real a[], const Real b[],
 Real & capital_a, Real & capital_b, Real & sum_b, Real & sum_a, Real binary_a[][6], Real kappa[][6]) const
{
  double yy[6], R;
	int i,j;

  R = 8.31446;

	sum_a = 0.0;
	sum_b = 0.0;

	kappa[1][1] = 0.0;
	kappa[1][2] = 0.47893;
	kappa[1][3] = 0.19014;
	kappa[1][4] = 0.32547;
	kappa[1][5] = 0.105;

	kappa[2][1] = 0.47893;
	kappa[2][2] = 0.0;
	kappa[2][3] = 0.100;
	kappa[2][4] = 0.0;
	kappa[2][5] = 0.084;

	kappa[3][1] = 0.19014;
	kappa[3][2] = 0.100;
	kappa[3][3] = 0.0;
	kappa[3][4] = -0.007;
	kappa[3][5] = 0.099;

	kappa[4][1] = 0.32547;
	kappa[4][2] = 0.0;
	kappa[4][3] = -0.007;
	kappa[4][4] = 0.0;
	kappa[4][5] = 0.0;

	kappa[5][1] = 0.105;
	kappa[5][2] = 0.084;
	kappa[5][3] = 0.099;
	kappa[5][4] = 0.0;
	kappa[5][5] = 0.0;


	yy[1] = 0.0;
	for(i=2;i<=5;i++)
	{
		yy[i] = y[i];
	}


    for(i=1;i<=5;i++)
	{
		sum_b += yy[i]*b[i];
    	for(j=1;j<=5;j++)
		{
    		binary_a[i][j]=sqrt(a[i]*a[j])*(1.0-kappa[i][j]);
    		sum_a += yy[i]*yy[j]*binary_a[i][j];
		}
	}


	capital_a = sum_a*press/(pow(R*temp,2));
	capital_b = sum_b*press/R/temp;

}

void
MoskitoMixture2P::ROOT_FINDER(const Real & capital_a, const Real & capital_b, Real & z) const
{
  double ra, rb, rc, rq, rr, sign_root, teta, x_one, x_two, x_three, s_prime, t_prime, PI;
  	int i,j;

    PI = 3.14159;

  	ra = -(1.0 - capital_b);
  	rb = capital_a - 2.0*capital_b - 3.0*pow(capital_b, 2.0);
  	rc = -(capital_a*capital_b) + pow(capital_b, 2.0) + pow(capital_b, 3.0);

  	rq = (ra*ra - 3.0*rb)/9.0;
  	rr = (2.0*pow(ra,3.0) - 9.0*ra*rb + 27.0*rc)/54.0;
  	sign_root = pow(rr,2) - pow(rq,3);


      if(sign_root<=0.0){
      	teta = acos(rr/pow(rq,3.0/2.0)) ;
          x_one = -2.0*sqrt(rq)*cos(teta/3.0) - ra/3.0;
          x_two = -2.0*sqrt(rq)*cos((teta+2*PI)/3.0) - ra/3.0;
          x_three = -2.0*sqrt(rq)*cos((teta-2*PI)/3.0) - ra/3.0;
          z = (x_one>x_two) ? ((x_one>x_three) ? x_one : x_three) : ((x_two>x_three) ? x_two : x_three);
  	} else {
  		s_prime = -rr/fabs(rr)*pow((fabs(rr)+sqrt(sign_root)),1.0/3.0);
  	    t_prime = rq/s_prime;
  	    z = s_prime + t_prime - ra/3.0;
  	}

}

void
MoskitoMixture2P::CONCENTRATION_CALCULATOR(const Real k[], const Real tmf[], Real & non_aqueous_mf, Real x[], Real y[]) const
{
  double u[30][6], u_prime1[30][6], u_prime2[30][6], sum_u[30], sum_u_prime[30], nna[30];
	double zk, zok;
	int i,j;

	nna[1] = 0.1;
	zk = 0.0;
	zok = 0.0;


	for(i=1; i<=5; i++)
	{
		zk += tmf[i]*(k[i]-1.0);
		zok += tmf[i]*(k[i]-1.0)/k[i];
	}


	for(j=1; j<=29; j++)
	{
		sum_u[j] = 0.0;
	    sum_u_prime[j] = 0.0;
		for(i=1; i<=5; i++)
		{
			u[j][i] = tmf[i]*(k[i]-1.0)/(1.0+(k[i]-1.0)*nna[j]);
    	    sum_u[j] += u[j][i];
    	    u_prime1[j][i] = tmf[i]*pow((k[i]-1.0), 2.0)/pow(nna[j]*(k[i]-1.0)+1.0, 2.0);
    	    sum_u_prime[j] += u_prime1[j][i];
		}
		nna[j+1] = nna[j] + sum_u[j]/sum_u_prime[j];
		if(sum_u[j]>=-0.001 && sum_u[j]<=0.001)
		{
			non_aqueous_mf = nna[j+1];
		    break;
		}

	}


	for(i=1;i<=5;i++)
	{
		x[i]=tmf[i]/(1+(k[i]-1)*non_aqueous_mf);
		y[i]=tmf[i]*k[i]/(1+(k[i]-1)*non_aqueous_mf);
	}

}

void
MoskitoMixture2P::FRAC_CONVERTER(const Real omf[], const Real tmf[], const Real & non_aqueous_mf, const Real x[],
   const Real y[], const Real m_mas[], Real s_molal[], Real s_molef[], Real s_massf[], Real molef[], Real massf[],
    Real & non_aqueous_massf, Real s_molality[], Real s_molefrac[], Real & mf_h2o) const

{
  double m[10], zigx, zigxm, ntaux, somf, somfm, ym;
	int i;

	m[1] = 18.01528e-3;
	m[2] = 16.04e-3;
	m[3] = 44.01e-3;
	m[4] = 28.0134e-3;
	m[5] = 34.1e-3;
	m[6] = 58.44e-3;
	m[7] = 74.5513e-3;
	m[8] = 110.98e-3;
	m[9] = 95.211e-3;

	mf_h2o = y[1]*m[1]/(y[1]*m[1]+y[2]*m[2]+y[3]*m[3]+y[4]*m[4]+y[5]*m[5]);


	zigx = 0.0;
	zigxm = 0.0;
	ntaux = 0.0;

	for(i=1;i<=5;i++)
	{
		zigx += x[i];
		zigxm += x[i]*m[i];
		ntaux += omf[i];
	}

	somf = 0.0;
	somfm = 0.0;

	for(i=6;i<=9;i++)
	{
		somf += omf[i];
		somfm += omf[i]*m[i];
	}

	for(i=1;i<=4;i++)
	{
		s_molality[i] = omf[i+5]/(1.0 - non_aqueous_mf)/x[1]/m[1]/ntaux;
		s_molal[i] = omf[i+5]/(1.0 - non_aqueous_mf)/zigxm/ntaux;
		s_molef[i] = omf[i+5]/(somf + ntaux*(1.0 - non_aqueous_mf)*zigx);
		s_molefrac[i] = omf[i+5]/(somf + ntaux*(1.0 - non_aqueous_mf)*x[1]);
		s_massf[i] = omf[i+5]*m[i+5]/(somfm + ntaux*(1.0 - non_aqueous_mf)*zigxm);
	}

	ym = 0.0;

	for(i=1;i<=5;i++)
	{
		molef[i] = ntaux*(1.0 - non_aqueous_mf)*x[i]/(somf + ntaux*(1.0 - non_aqueous_mf)*zigx);
		massf[i] = ntaux*(1.0 - non_aqueous_mf)*x[i]*m[i]/(somfm + ntaux*(1.0 - non_aqueous_mf)*zigxm);
		ym += y[i]*m[i];
	}

	//cout<<"zigx = "<<zigx<<" zigxm = "<<zigxm<<" ntaux = "<<ntaux<<" somf = "<<somf<<" somfm = "<<somfm<<" ym = "<<ym<<endl;

	non_aqueous_massf = non_aqueous_mf*ntaux*ym/(somfm + ntaux*(1.0 - non_aqueous_mf)*zigxm + non_aqueous_mf*ntaux*ym);

}

void
MoskitoMixture2P::ACTIVITY(const Real & press, const Real & t, const Real x[],
   const Real tmf[], const Real & non_aqueous_mf, Real gama[], const Real s_molality[]) const
{
  double lambda[5], exi[5];
	double pp, mna, mk, mca, mmg, mcl;
	int i;

	pp = press/100000.0 ;
	//s_mf[1] = 6.0;

	lambda[1] = -5.7066455e-1 + 7.2997588e-4*t + 1.5176903e2/t + 3.1927112e-5*pp - 1.642651e-5*pp/t ;
	lambda[2] = -0.0652869 + 1.6790636e-4*t + 40.838951/t - 3.9266518e-2*pp/t + 2.1157167e-2*pp/(630.0 - t)
	 + 6.5486487e-6*t*log(pp);
	lambda[3] = -2.0939363 + 3.1445269E-03*t + 3.91E+02/t - 2.9973977E-07*pp - 1.5918098E-05*pp/t;
	lambda[4] = 1.03658689 - 1.1784797E-03*t - 1.7754826E+02/t - 4.5313285E-04*pp + 0.47751650E+02*pp/t/t;

	exi[1] = -2.999084e-3;
	exi[2] = -1.14462e-2 + 2.8274958e-5*t + 1.3980876e-2*pp/t - 1.4349005e-2*pp/(630.0 - t);
	exi[3] = -6.3981858E-03;
	exi[4] = -0.010274152;
	/*
    mna = s_mf[1];
	mca = s_mf[3];
	mk = s_mf[2];
	mmg = s_mf[4];
	mcl = s_mf[1];
	*/
	mna = s_molality[1];
	mca = s_molality[3];
	mk = s_molality[2];
	mmg = s_molality[4];
	mcl = s_molality[1];

	//cout<<"mna = "<<mna<<"  mcl = "<<mcl<<"  mca = "<<mca<<"  mk = "<<mk<<"  mmg = "<<mmg<<endl;

	for(i=1; i<=4; i++)
	{
		gama[i] = exp(2.0*lambda[i]*(mna + mk + 2.0*mca + 2.0*mmg) + mcl*(mna + mk + mca + mmg)*exi[i]);
    //cout<<gama[i]<<endl;
	}

}

void
MoskitoMixture2P::HENRY(const Real & press, const Real & temp, Real hen[]) const
{
  double ktav[5], kbet[5], etah[5], dbeta[5];
	double tet, tet2, tet3, tet4, tet5, ah1, ah2, bh, surv, makhv, v0, v, p, t;
	double ai1, ai2, ai3, ai4, ai5, ai6, tav, nam, ps, fu, tcw, pcw, patm;

  p = press ;
  t = temp;

	  tet = temp - 273.15;
    tet2 = tet*tet;
    tet3 = tet2*tet;
    tet4 = tet3*tet;
    tet5 = tet4*tet;

    ah1 = 3.2891 - 2.3910e-3*tet + 2.8446e-4*tet2 - 2.8200e-6*tet3 + 8.477e-9*tet4;
    ah2 = 6.245e-5 - 3.913e-6*tet - 3.499e-8*tet2 + 7.942e-10*tet3 - 3.299e-12*tet4;
    bh = 19654.320 + 147.037*tet - 2.21554*tet2 + 1.0478e-2*tet3 - 2.2789e-5*tet4;

    surv = 1.0 + 18.159725e-3*tet;
    makhv = 0.9998396 + 18.224944e-3*tet - 7.922210e-6*tet2 - 55.44846e-9*tet3 + 149.7562e-12*tet4 - 393.2952e-15*tet5;

    patm = p/100000.0;

    v0 =  (1.0 + 18.159725e-3*tet)/(0.9998396 + 18.224944e-3*tet - 7.922210e-6*tet2 - 55.44846e-9*tet3 +
	 149.7562e-12*tet4 - 393.2952e-15*tet5);
    v = v0 - v0*patm/(bh + ah1*patm + ah2*patm*patm);

    ai1 = -7.85951783;
    ai2 = 1.84408259;
    ai3 = -11.7866497;
    ai4 = 22.6807411;
    ai5 = -15.9618719;
    ai6 = 1.80122502;

    tcw = 647.096;
    pcw = 220.64;

    tav = 1.0 - t/tcw;
    nam = tcw/t*(ai1*tav + ai2*pow(tav,1.5) + ai3*pow(tav,3.0) + ai4*pow(tav,3.5) + ai5*pow(tav,4.0) + ai6*pow(tav,7.5));
    ps = pcw*exp(nam);
    fu = ps*exp((patm-ps)*v*18.0152/83.1447/t);

    ktav[1] = -5.779280;
    ktav[2] = -5.279063;
    ktav[3] = -5.175337;
    ktav[4] = 0.27049433;

    kbet[1] = 7.262730;
    kbet[2] = 6.187967;
    kbet[3] = 6.906469;
    kbet[4] = 0.27543436;

    etah[1] = -0.092248;
    etah[2] = -0.114535;
    etah[3] = -0.008194;
    etah[4] = 0.77357854;

	for(int i=1; i<=4; i++)
	{
		dbeta[i] = ktav[i] + kbet[i]*sqrt(1000./t);
		hen[i] = exp((1.0-etah[i])*log(fu) + etah[i]*log(83.1447*t/v/18.0152) + 2./v*dbeta[i])*100000.0;
    //hen[i] = tet;
	}

}

void
MoskitoMixture2P::FUG(const Real & capital_a, const Real & capital_b, const Real & z, const Real y[], const Real b[],
   const Real & sum_a, const Real & sum_b, const Real binary_a[][6], Real phi[]) const
{
  double c1, c2, c3[6], c4[6], c5[6], yy[6];
	int i,j;

	c1 = log(z-capital_b);
	c2 = log((z + 2.414*capital_b)/(z - 0.414*capital_b));

	yy[1] = 0.0;
	for(i=2;i<=5;i++)
	{
		yy[i] = y[i];
	}


	for(i=1;i<=5;i++)
	{
		c3[i] = 0.0;
		for(j=1;j<=5;j++)
		{
			c3[i] += yy[j]*binary_a[j][i];
		}
		c4[i] = b[i]/sum_b*(z - 1.0);
		c5[i] = capital_a/2.828/capital_b*(b[i]/sum_b - 2.0*c3[i]/sum_a)*c2;
		phi[i] = exp(c4[i] - c1 + c5[i]);
	}

}

void
MoskitoMixture2P::K_CONVERGENCE(const Real & press, const Real & temp, const Real hen[], const Real gama[],
   const Real phi[], Real k[], Real k_old[]) const
{
  int i,j;
	double tkel, logk0, k0, y2, pbar;

	pbar = press/100000.0;

	k_old[1] = k[1];

	tkel = temp - 273.15;
	logk0 = -2.209 + 3.0972e-2*tkel - 1.098e-4*pow(tkel, 2.0) + 2.048e-7*pow(tkel, 3.0);
	k0 = pow(10.0, logk0);

    k[1] = k0*exp((press/100000.0 - 1.0)*18.18/83.14472/temp)/phi[1]/pbar;


	for(i=2; i<=5; i++)
	{
		k_old[i] = k[i];
		k[i] = hen[i-1]*gama[i-1]/press/phi[i];
	}

}

void
MoskitoMixture2P::TH_DRIESNER(const Real & press, const Real & temp, const Real s_molef[], const Real molef[], Real & th) const
{
  double q1, q2, q10, q11, q12, q20, q21, q22, q23, q1x1, q2x1, pbar, pbar2, xn;

	xn =s_molef[1]/(s_molef[1] + molef[1]);
	//cout<<"xn = "<<xn<<" s_molef[1] = "<<s_molef[1]<<" molef[1] = "<<molef[1]<<endl;
	pbar = press*1.0e-5;
    pbar2 = pbar * pbar;
	q11 = - 32.1724 + 0.0621255 * pbar;
    q21 = - 1.69513 - 4.52781e-4 * pbar - 6.04279e-8 * pbar2;
    q22 = 0.0612567 + 1.88082e-5 * pbar;

    q1x1 = 47.9048 - 9.36994e-3 * pbar + 6.51059e-6 * pbar2;
    q2x1 = 0.241022 + 3.45087e-5 * pbar - 4.28356e-9 * pbar2;

    q12 = - q11 - q1x1;
    q10 = q1x1;

    q20 = 1.0 - q21 * sqrt(q22);
    q23 = q2x1 - q20 - q21 * sqrt(1.0 + q22);

    q1 = q10 + q11 * (1.0 - xn) + q12 * (1.0 - xn) * (1.0 - xn);
    q2 = q20 + q21 * sqrt(xn + q22) + q23 * xn;

    th = q1 + q2 * (temp - 273.15) + 273.15;

}

Real
MoskitoMixture2P::SATURATED_H_WATER(const Real & th) const
{
  double d1, d2, d3, d4, d5, d6, shw, tc;

  tc = th - 273.15;

  d1 = -2.844699e-2;
  d2 = 4.211925;
  d3 = -1.017034e-3;
  d4 = 1.311054e-5;
  d5 = -6.756469e-8;
  d6 = 1.724481e-10;

  shw = d1 + d2*tc + d3*pow(tc,2.0) + d4*pow(tc,3.0) + d5*pow(tc,4.0) + d6*pow(tc,5.0) ;

return shw;

}

Real
MoskitoMixture2P::SPECIFIC_V_WATER(const Real & press, const Real & th) const
{
  double tet, tet2, tet3, tet4, tet5, a1, a2, b, surv, makhv, v0, v, pbar;

    	pbar = press*1.0e-5;
    	tet = th - 273.15;
      tet2 = tet*tet;
      tet3 = tet2*tet;
      tet4 = tet3*tet;
      tet5 = tet4*tet;

      a1 = 3.2891 - 2.3910e-3*tet + 2.8446e-4*tet2 - 2.8200e-6*tet3 + 8.477e-9*tet4;
      a2 = 6.245e-5 - 3.913e-6*tet - 3.499e-8*tet2 + 7.942e-10*tet3 - 3.299e-12*tet4;
      b = 19654.320 + 147.037*tet - 2.21554*tet2 + 1.0478e-2*tet3 - 2.2789e-5*tet4;

      surv = 1.0 + 18.159725e-3*tet;
      makhv = 0.9998396 + 18.224944e-3*tet - 7.922210e-6*tet2 - 55.44846e-9*tet3 + 149.7562e-12*tet4 - 393.2952e-15*tet5;

      v0 = (1.0 + 18.159725e-3*tet)/(0.9998396 + 18.224944e-3*tet - 7.922210e-6*tet2 - 55.44846e-9*tet3 + 149.7562e-12*tet4 - 393.2952e-15*tet5);

      v = (v0 - v0*pbar/(b + a1*pbar + a2*pbar*pbar))/1000.0;

  	return v;

}

Real
MoskitoMixture2P::T_EXPANSION_WATER(const Real & th) const
{
  double at, bt, ct, dt, beta, tc;

  	tc = th - 273.15;

	  at = -6.8785895e-5;
    bt = 2.1687942e-5;
    ct = -2.1236686e-6;
    dt = 7.7200882e-8;

    beta = at + bt*tc + ct*pow(tc, 1.5) + dt*pow(tc, 2.0);

	  return beta;

}

Real
MoskitoMixture2P::P_SAT_WATER(const Real & th) const
{
  double a1, a2, a3, a4, a5, a6, tc, pc, tav, nam, ps;

	  a1 = -7.85951783;
    a2 = 1.84408259;
    a3 = -11.7866497;
    a4 = 22.6807411;
    a5 = -15.9618719;
    a6 = 1.80122502;

    tc = 647.096;
    pc = 220.64;

    tav = 1.0 - th/tc;

    nam = tc/th*(a1*tav + a2*pow(tav,1.5) + a3*pow(tav,3.0) + a4*pow(tav,3.5) + a5*pow(tav,4.0) + a6*pow(tav,7.5));

    ps = pc*exp(nam);

	 return ps*100000.0;

}

void
MoskitoMixture2P::SALT_ENTHALPY(const Real & temp, const Real s_molality[], Real sh[]) const
{
  int i,j;
	double a1[4], a2[4], a3[4], a4[4], b1[4], b2[4], b3[4], b4[4], b5[4], mm[4], t;

  t = temp;

	//KCl
    a1[1] = 0 ;
    a2[1] = 3.9580e3 ;
    a3[1] = -1.5990e2 ;
    a4[1] = 2.1150e-1 ;
    b1[1] = 4.7190e-1 ;
    b2[1] = -2.7380 ;
    b3[1] = -5.2720e1 ;
    b4[1] = 2.1460e4 ;
    b5[1] = 6.2710e2 ;

    mm[1] = 74.5513e-3;

    //CaCl2
    a1[2] = 9.0590E+4 ;
    a2[2] =  0;
    a3[2] =  0;
    a4[2] =  0;
    b1[2] =  5.1170E-1;
    b2[2] =  -2.2890E+0;
    b3[2] =  1.3030E+0;
    b4[2] =  1.9550E+4;
    b5[2] =  6.2030E+2;

    mm[2] = 110.98e-3;

    //MgCl2
    a1[3] = 0 ;
    a2[3] = -4.4950E+4 ;
    a3[3] = 5.9170E+2 ;
    a4[3] = -7.9050E-1 ;
    b1[3] = 5.0200E-1 ;
    b2[3] = -1.5620E+0 ;
    b3[3] = 9.5460E+1 ;
    b4[3] = 9.1770E+3 ;
    b5[3] = 6.0230E+2 ;

    mm[3] = 95.211e-3;

    for(i=1;i<=3;i++)
    {
    	sh[i] = a1[i] + a2[i]*s_molality[i+1] + (a3[i] + a4[i]*t)*t + (pow(s_molality[i+1], b1[i]) + b2[i])*(b3[i]*t - b4[i]*log(1.0 - t/b5[i]));
    	sh[i] /= mm[i];
	}

}

void
MoskitoMixture2P::LIQUID_ENTHALPY(const Real s_massf[], const Real massf[], const Real sh[],
   Real & liq_h, const Real & water_h, const Real gas_ha[], const Real molarv_deriv[]) const
{
  int i;
	double mmg[6];

	mmg[1] = 18.01528e-3;
	mmg[2] = 16.04e-3;
	mmg[3] = 44.01e-3;
	mmg[4] = 28.0134e-3;
	mmg[5] = 34.1e-3;

	liq_h = 0.0;

	for(i=1;i<=3;i++)
	{
		liq_h += sh[i]*s_massf[i+1];
	}

	for(i=2;i<=5;i++)
	{
		liq_h += (gas_ha[i] + molarv_deriv[i-1]/mmg[i])*massf[i];
		//cout<<i<<" = "<<molarv_deriv[i-1]<<"    "<<molarv_deriv[i-1]/mmg[i]<<endl;
	}

	liq_h += (massf[1]+s_massf[1])*water_h;

}

void
MoskitoMixture2P::GAS_H(const Real & t, const Real y[], const Real & z, const Real & sum_a, const Real & sum_b,
   const Real & capital_b, const Real a[], const Real kappa[][6], Real & gas_h, const Real m[], const Real tc[],
   const Real pc[], const Real & mf_h2o, const Real & hh2o) const
{
  int i,j;
	double ah0[6], ah1[6], ah2[6], ah3[6], ah4[6], h1[6], dadt[6], mmg[6], reft[6], refh[6];
	double t1[6], t2[6], t3[6], t4[6], t5[6], r, h1g, h22g, dadt_mix, h21g, idh_mix, h_depart, mmix;

	r = 8.31446;

	mmg[1] = 18.01528e-3;
	mmg[2] = 16.04e-3;
	mmg[3] = 44.01e-3;
	mmg[4] = 28.0134e-3;
	mmg[5] = 34.1e-3;

	mmix = 0.0;

	for(i=2;i<=5;i++)
	{
		mmix += y[i]*mmg[i];
	}

	reft[1] = 273.15;
	reft[2] = 111.667;
	reft[3] = 273.15;
	reft[4] = 273.15;
	reft[5] = 212.85;

	refh[1] = 0.0;
	refh[2] = 8440;
	refh[3] = 21340;
	refh[4] = 8080;
	refh[5] = 18700.0;

	ah0[1] = 4.395;
    ah1[1] = -4.186e-3;
    ah2[1] = 1.405e-5;
    ah3[1] = -1.564e-8;
    ah4[1] = 0.632e-11;

    ah0[2] = 4.568;
    ah1[2] = -8.975e-3;
    ah2[2] = 3.631e-5;
    ah3[2] = -3.407e-8;
    ah4[2] = 1.091e-11;

    ah0[3] = 3.259;
    ah1[3] = 1.356e-3;
    ah2[3] = 1.502e-5;
    ah3[3] = -2.374e-8;
    ah4[3] = 1.056e-11;

/////////////////////////////////////////

    ah0[4] = 3.539;
    ah1[4] = -0.261e-3;
    ah2[4] = 0.007e-5;
    ah3[4] = 0.157e-8;
    ah4[4] = -0.099e-11;

    ah0[5] = 4.266;
    ah1[5] = -3.438e-3;
    ah2[5] = 1.319e-5;
    ah3[5] = -1.331e-8;
    ah4[5] = 0.488e-11;


    for(i=2;i<=5;i++)
    {
    	t1[i] = t - reft[i];
        t2[i] = pow(t,2.0) - pow(reft[i],2.0);
        t3[i] = pow(t,3.0) - pow(reft[i],3.0);
        t4[i] = pow(t,4.0) - pow(reft[i],4.0);
        t5[i] = pow(t,5.0) - pow(reft[i],5.0);
	}

    for(i=2;i<=5;i++)
    {
    	h1[i] = (ah0[i]*t1[i] + ah1[i]/2.0*t2[i] + ah2[i]/3.0*t3[i] + ah3[i]/4.0*t4[i] + ah4[i]/5.0*t5[i])*r + refh[i] ;
	}

	idh_mix = 0.0;

	for(i=2;i<=5;i++)
	{
		idh_mix += y[i]*h1[i];
	}

    //cout<<h1[1]<<' '<<h1[2]<<' '<<h1[3]<<endl;


    h1g = r*t*(z - 1.0);

    h22g = log((z + 2.414*capital_b)/(z - 0.414*capital_b));

    //cout<<h22g<<endl;


    for(i=2;i<=5;i++)
    {
    	dadt[i] = -0.45724*pow(r*tc[i], 2.0)/pc[i]*(1.0 + m[i]*(1.0 - sqrt(t/tc[i])))*(m[i]/sqrt(t*tc[i]));
	}

	//cout<<dadt[1]<<' '<<dadt[2]<<' '<<dadt[3]<<endl;


	dadt_mix = 0.0;


	for(i=2;i<=5;i++)
	{
		for(j=2;j<=5;j++)
		{
			dadt_mix += 0.5*y[i]*y[j]*(1.0 - kappa[i][j])*sqrt(a[i]*a[j])*(1.0/a[i]*dadt[i] + 1.0/a[j]*dadt[j]);
		}
	}

	//cout<<dadt_mix<<endl;

	h21g = (t*dadt_mix - sum_a)/(2.0*sqrt(2.0)*sum_b);

	//cout<<h1g<<' '<<h21g<<' '<<h22g<<endl;

	h_depart = h1g + h21g*h22g;

	//h_depart = h21g*h22g ;

	//cout<<h_depart<<endl;

	//gas_h = (idh_mix + h_depart)/mmix;
	//cout<<"  z1 = "<<gas_h<<endl;
	gas_h = mf_h2o*hh2o + (1.0 - mf_h2o)*((idh_mix + h_depart)/mmix);
  //cout<<gas_h<<endl;
	//cout<<"  z2 = "<<gas_h<<endl;


	//cout<<h1[3]<<endl;


}

void
MoskitoMixture2P::GAS_RO(const Real & p, const Real & t, const Real & z, const Real y[],
   const Real tc[], const Real pc[], Real & g_ro, const Real & roh2o, const Real & mf_h2o) const
{

  int i;
	double sum_mmg, zg, vg, tr, pr, r, c_mix;
	double mmg[6], c[6];

	r = 8.31446e6;

	tr = t/tc[3];
	pr = p/pc[3];

	mmg[1] = 18.01528;
	mmg[2] = 16.04;
	mmg[3] = 44.01;
	mmg[4] = 28.0134;
	mmg[5] = 34.1;

	c[1] = 0.0;
	c[2] = 3.2;
	c[3] = pow(tr, 6.27844)*(-8.18682*pow(pr, -0.94163) + 0.94612) + 2.49642;
	c[4] = 3.9;
	c[5] = 2.17;

	//cout<<c[3]<<endl;

	sum_mmg = 0.0;
	for(i=2;i<=5;i++)
	{
		sum_mmg +=  y[i]*mmg[i];
	}

	c_mix = 0.0;
	for(i=2;i<=5;i++)
	{
		c_mix +=  y[i]*c[i];
	}

	vg = z*r*t/p;

//	cout<<"  vg = "<<vg<<endl;

	//g_ro = sum_mmg/(vg + c_mix)*1000.0;
	//cout<<"1 = "<<g_ro<<"1 = "<<roh2o<<endl;
	g_ro = mf_h2o*roh2o + (1.0 - mf_h2o)*(sum_mmg/(vg + c_mix)*1000.0);
	//cout<<"2 = "<<g_ro<<"1 = "<<mf_h2o<<endl;
//    g_ro = sum_mmg/(vg);

	//cout<<vg + c_mix<<"vg"<<endl;
	//cout<<g_ro<<'r'<<endl;

}

void
MoskitoMixture2P::TV_DRIESNER(const Real & press, const Real & temp, const Real tmf[], Real & tv) const
{
  double n1, n2, n11, n12, n20, n21, n22, n23, n1x1, n2x1, pbar, pbar2, pbar3, xn, m_nacl, m_water, m_mix;

	//xn = tmf[4];
	xn = 0.1;
	m_nacl = 58.44;
	m_water = 18.015;
	m_mix = xn*m_nacl + (1.0 - xn)*m_water;

    pbar = press*1.0e-5;
    pbar2 = pbar * pbar;
    pbar3 = pbar2 * pbar;

    n11 = -54.2958 - 45.7623*exp(-9.44785e-4*pbar);
    n21 = -2.6142 - 2.39092e-4*pbar;
    n22 = 0.0356828 + 4.37235e-6*pbar + 2.0566e-9*pbar2;
    n1x1 = 330.47 + 0.942876*sqrt(pbar) + 0.0817193*pbar - 2.47556e-8*pbar2 + 3.45052e-10*pbar3;
    n2x1 = -0.0370751 + 0.00237723*sqrt(pbar) + 5.42049e-5*pbar + 5.84709e-9*pbar2 - 5.99373e-13*pbar3;
    n12 = -n1x1 - n11;
    n20 = 1.0 - n21*sqrt(n22);
    n23 = n2x1 - n20 - n21*sqrt(1.0 + n22);

    n1 = n1x1 + n11*(1.0 - xn) + n12*(1.0 - xn)*(1.0 - xn);
    n2 = n20 + n21*sqrt(xn + n22) + n23*xn;
    tv = n1 + n2 * (temp - 273.15) + 273.15;

}

void
MoskitoMixture2P::LIQUID_RHO(const Real & p, const Real & t, Real & l_rho, const Real s_molality[],
   const Real s_massf[], const Real massf[], const Real molarv_vol[]) const
{
  int i,j ;
	double rhos, tc, bi, ps;
	double cb[4], b_cap[4], mi[4], sum1[4], sum2[4], rhoref[4], xs[4], m_m[4], m_g[6];


	double a3d[4][3][5] = {
		                     { { 2863.158, -46844.356, 120760.118, -116867.722, 40285.426 }, {-2000.028, 34013.518, -88557.123, 86351.784, -29910.216 }, { 413.046, -7125.857, 18640.780, -18244.074, 6335.275 } },
		                     { { 2332.802, -39637.418, 104801.288, -104266.828, 37030.556 }, {-1287.572, 23543.994, -63846.097, 65023.561, -23586.370 }, { 206.032, -4003.757, 11128.162, -11595.475, 4295.498 } },
		                     { { 2546.760, -39884.946, 102056.957, -98403.334 , 33976.048 }, {-1362.157, 22785.572, -59216.108, 57894.824, -20222.898 }, { 217.778, -3770.645, 9908.135 , -9793.484 , 3455.587 } },
		                     { { 2385.823, -38428.112, 99526.269 , -97041.399 , 33841.139 }, {-1254.938, 21606.295, -56988.274, 56465.943, -19934.064 }, { 192.534, -3480.374, 9345.908 , -9408.904 , 3364.018 } }
	                      };


	double bet3d[4][3][5] = {
		                     {{-1622.40, 9383.80, -14893.80, 7309.10},{241.57, -980.97, 1482.31, -750.98}},
		                     {{-1622.40, 9383.80, -14893.80, 7309.10},{211.49, -888.16, 1400.09, -732.79}},
		                     {{-1622.40, 9383.80, -14893.80, 7309.10},{307.24,-1259.10, 2034.03,-1084.94}},
		                     {{-1622.40, 9383.80, -14893.80, 7309.10},{358.00,-1597.10, 2609.47,-1383.91}}
	                      };


	tc = t-273.15;
	bi = s_molality[1] + s_molality[2] + s_molality[3] + s_molality[4];

    for(i=0; i<=3; i++)
    {
    	if(bi == 0.0)
    	{
    		xs[i] = 0.0;
    		xs[0] = 1.0;
		}
		else
		{
			xs[i] = s_molality[i+1]/bi;
		}
	}

	//cout<<"bi = "<<bi<<","<<xs[0]<<","<<xs[1]<<","<<xs[2]<<","<<xs[3]<<endl;

	//xs[0] = 1.0;
	//xs[1] = 0.0;
	//xs[2] = 0.0;
	//xs[3] = 0.0;


	m_m[0] = 58.44-3;
	m_m[1] = 74.5513-3;
	m_m[2] = 110.98e-3;
	m_m[3] = 95.211e-3;

	m_g[1] = 18.01528e-3;
	m_g[2] = 16.04e-3;
	m_g[3] = 44.01e-3;
	m_g[4] = 28.0134e-3;
	m_g[5] = 34.1e-3;

	rhos = 999.79684 + 0.068317355*tc - 0.010740248*pow(tc, 2.0) + 0.00082140905*pow(tc, 2.5) - 2.3030988e-5*pow(tc, 3.0);

	int snum;

    for(snum=0;snum<=3;snum++)
    {
    	sum1[snum] = 0.0;
	    sum2[snum] = 0.0;
    	for(i=0;i<=2;i++)
	    {
		    sum1[snum] += a3d[snum][i][0]*pow(bi, (i*1.0 + 2.0)/2.0) ;
		     //cout<<sum1<<endl;

		    for(j=1;j<=4;j++)
		    {
			sum2[snum] += a3d[snum][i][j]*pow(bi, (i*1.0 + 2.0)/2.0)*pow(t/647.10, (j*1.0 + 1.0)/2.0);
			//cout<<i<<'i'<<j<<'j'<<sum2<<endl;
	  	    }
	    }
	}

	for(snum=0;snum<=3;snum++)
    {
    	rhoref[snum] = rhos + sum1[snum] + sum2[snum];
	}


	//cout<<rhos<<endl;
	//cout<<sum1<<' '<<sum2<<endl;



	//cout<<rhoref<<endl;

	cb[0] = 0.11725 - 0.00134*bi + 0.00056*pow(bi, 1.5);
	cb[1] = 0.11725 - 0.00170*bi + 0.00083*pow(bi, 1.5);
	cb[2] = 0.11725 - 0.00493*bi + 0.00231*pow(bi, 1.5);
	cb[3] = 0.11725 - 0.00789*bi + 0.00142*pow(bi, 1.5);

	//cout<<cb<<endl;

    for(snum=0;snum<=3;snum++)
    {
    	b_cap[snum] = 0.0;
	    mi[snum] = 0.0;
	    for(i=0;i<=1;i++)
	    {
		    for(j=0;j<=3;j++)
		    {
			    b_cap[snum] += bet3d[snum][i][j]*pow(bi, i*1.0)*pow(t/647.10, j*1.0);
			    mi[snum] += bet3d[snum][i][j];
		   }
	    }
	}



	//cout<<mi<<endl;

	ps = P_SAT_WATER(t);

	//cout<<ps<<' '<<p<<' '<<b_cap<<endl;

	//l_rho = rhoref/ (1.0 - cb*log((b_cap+p/1000000.0)/(b_cap+ps/1000000.0)));

	double auxden[4], ft1, ft2;

	ft1 = 0.0;
	ft2 = 0.0;

	for(snum=0;snum<=3;snum++)
    {
    	auxden[snum] = rhoref[snum]/ (1.0 - cb[snum]*log((b_cap[snum]+p/1000000.0)/(b_cap[snum]+ps/1000000.0)));
    	ft1 += xs[snum]*(1.0 + bi*m_m[snum]);
    	ft2 += xs[snum]*(1.0 + bi*m_m[snum])/auxden[snum];
	}


	//l_rho = rhoref/ (1.0 - cb*log((b_cap+p/1000000.0)/(b_cap+ps/1000000.0)));
	double l_rhoo, w1;

	l_rhoo = ft1/ft2;

	w1 = massf[1] + s_massf[1] + s_massf[2] + s_massf[3] + s_massf[4];

	l_rho = w1/l_rhoo;

	for(i=2; i<=5; i++)
  {
    if (molarv_vol[i+1]!= 0.0)
    {
      l_rho +=  massf[i]/m_g[i]*molarv_vol[i-1];
    }
	}

	l_rho = pow(l_rho, -1);
	//l_rho = l_rhoo;


	//cout<<l_rho<<'r'<<endl;

}
/*
void
MoskitoMixture2P::GAS_VIS(const Real & press, const Real & temp, const Real tc[], const Real pc[],
   const Real y[], const Real & g_ro, Real & g_vis) const
{
  int i;
	double a[5], c[10] ;
	double ts, sumg, gs, eta0, deta, d11, d21, d64, d81, d82, gcc, ts2, gamas, omega, eta_di, be, roz, ux, g_vis2;

	a[0] = 0.235156;
	a[1] = -0.491266;
	a[2] = 5.211155e-2;
	a[3] = 5.347906e-2;
	a[4] = -1.537102e-2;

	d11 = 0.4071119e-2;
	d21 = 0.7198037e-4;
	d64 = 0.2411697e-16;
	d81 = 0.2971072e-22;
	d82 = -0.1627888e-22;

	ts = temp/252.196;

	sumg = 0.0;

	for(i=0;i<=4;i++)
	{
		sumg += a[i]*pow(log(ts), i);
	}

	gs = exp(sumg);

	eta0 = 1.0069*sqrt(temp)/gs;

	deta = d11*g_ro + d21*pow(g_ro, 2.0) + d64*pow(g_ro, 6.0)/pow(ts, 3.0) + d81*pow(g_ro, 8.0) + d82*pow(g_ro, 8.0)/ts ;

	g_vis = eta0 + deta;

	//cout<<g_vis<<'v'<<endl;

//////////////////////////////////////////

	c[1] = -3.0328138281;
	c[2] = 16.918880086;
	c[3] = -37.189364917;
	c[4] = 41.288861858;
	c[5] = -24.615921140;
	c[6] = 8.9488430959;
	c[7] = -1.8739245042;
	c[8] = 0.20966101390;
	c[9] = -9.6570437074e-3;

	gcc = g_ro/1000.0 ;

	ts2 = temp/174.0;

	gamas = 0.0;

	for(i=1;i<=9;i++)
	{
		gamas += c[i]*pow(ts2, ((i*1.0 - 1.0)/3.0 - 1.0));
	}

	omega = 1.0/gamas;

	eta_di = 10.5*sqrt(ts2)/omega;

	be = 6.71 + 1.969e-2*temp ;

	roz = 0.799 + 4.055e-3*temp ;

	ux = log(eta_di*roz/((be - 1.0)*gcc + roz)) ;

	g_vis2 = exp(ux + be*gcc/(roz - gcc)) ;

	//cout<<g_vis2*106.0/100.0<<'2'<<endl;

}
*/
void
MoskitoMixture2P::LIQUID_VIS(const Real & t, const Real s_massf[], const Real massf[], Real & l_vis) const
{
  int i;
	double a[5], c[10], gas_massf ;
	double a0, a1, a2, a3, tc, w_vis, v1[5], v2[5], v3[5], v4[5], v5[5], v6[5], muc1[5], muc2[5], muc[5], pre_vis, w[5];

	gas_massf = massf[2] + massf[3] + massf[4] + massf[5];

	w[0] = massf[1]/(1.0 - gas_massf);
	w[1] = s_massf[1]/(1.0 - gas_massf);
	w[2] = s_massf[2]/(1.0 - gas_massf);
	w[3] = s_massf[3]/(1.0 - gas_massf);
	w[4] = s_massf[4]/(1.0 - gas_massf);

	//cout<<t-273<<"  "<<w[3]<<endl;

	a0 = 137.37;
	a1 = 5.2842;
	a2 = 0.05594;
	a3 = 246.0;

	v1[1] = 16.222;
	v2[1] = 1.3229;
	v3[1] = 1.4849;
	v4[1] = 0.0074691;
	v5[1] = 30.78;
	v6[1] = 2.0583;

	v1[2] = 6.4883;
	v2[2] = 1.3175;
	v3[2] = -0.77785;
	v4[2] = 0.092722;
	v5[2] = -1.3;
	v6[2] = 2.0811;

	v1[3] = 32.028;
	v2[3] = 0.78792;
	v3[3] = -1.1495;
	v4[3] = 0.0026995;
	v5[3] = 780860.0;
	v6[3] = 5.8442;

	v1[4] = 24.032;
	v2[4] = 2.2694;
	v3[4] = 3.7108;
	v4[4] = 0.021853;
	v5[4] = -1.1236;
	v6[4] = 0.14474;

	tc = t - 273.15;

	w_vis =  (tc + a3)/(a2*tc*tc + a1*tc + a0);

	for(i=1;i<=4;i++)
	{
		muc1[i] = (v1[i]*pow((1.0 - w[0]), v2[i]) + v3[i])/(v4[i]*tc + 1.0);
	    muc2[i] = v5[i]*pow((1.0 - w[0]), v6[i]) + 1.0 ;
	    muc[i] = exp(muc1[i])/muc2[i];
	}

	l_vis = pow(w_vis, w[0]);

	for(i=1;i<=4;i++)
	{
		l_vis *=  pow(muc[i],w[i]);
	}

	l_vis *= 0.001;

}

void
MoskitoMixture2P::VAR_SEQUENCE(const Real & var, vector<Real> & var_seq) const
{

  var_seq.push_back(var);

}

void
MoskitoMixture2P::VAR_GRADIENT(const Real & press, const Real & temp, const Real tmf[],
   const vector<Real> & var_seq, Real & var_p, Real & var_t, Real & var_z) const
{

  double delta_p, delta_t, delta_z;

	delta_p = 0.001 * press;
	delta_t = 0.001 * temp;
	delta_z = 0.001 * 1.001*tmf[3];

	var_p = (var_seq[1] - var_seq[2])/delta_p/2.0 ;
	var_t = (var_seq[3] - var_seq[4])/delta_t/2.0 ;
	var_z = (var_seq[5] - var_seq[6])/delta_z/2.0 ;

}

void
MoskitoMixture2P::Transport_PROP(const Real & g_ro, const Real & l_ro, const Real & non_aqueous_mf,
   const Real k[], Real & rog_k_n, Real & rol_n) const
{

  rog_k_n = g_ro*non_aqueous_mf*k[3];
	rol_n = l_ro*non_aqueous_mf;

}

//void
//MoskitoMixture2P::INITIALIZER(const Real & pressure, const Real & temperature, Real & press, Real & temp, const Real omf[], Real tmf[], Real m_mass[], Real k[]) const
//{


//}

void
MoskitoMixture2P::mix_PROP(const Real x[], const Real y[], const Real & non_aqueous_mf,
   const Real & liq_prop, const Real & gas_prop, const Real m_mas[], Real & mix_prop, Real & mass_frac) const
{
  double sum1, sum2;
  int i;

  sum1 = 0.0;
  sum2 = 0.0;

	for(i=1;i<=5;i++)
	{
		sum1 += non_aqueous_mf*y[i]*m_mas[i];
    sum2 += non_aqueous_mf*y[i]*m_mas[i] + (1.0 - non_aqueous_mf)*x[i]*m_mas[i];
	}

  mass_frac = sum1/sum2;

  mix_prop = mass_frac*gas_prop + (1.0 - mass_frac)*liq_prop;

  //mix_prop = gas_prop + liq_prop;

}

void
MoskitoMixture2P::Itera(const Real & h, const Real & p, const Real & m,const Real c_vect[], Real & g_ro, Real & l_rho,
   Real & lhx, Real & gas_h, Real & mass_frac, Real & non_aqueous_mf, Real & rho_mix, Real & temper, Real x[],
    Real aux1_arr[], Real aux2_arr[], Real & viscosity, Real & conductivity) const
{

  int kk;
  Real t, capital_a, capital_b, z, sum_a, sum_b, th, water_h, h_mix;
  Real tv, g_vis, l_vis, omf_sum, liq_h, cp_depart;
  Real k[6], tc[6], pc[6], mm[6], a[6], b[6], omf[10], tmf[6], y[6], m_mas[10], gas_ha[6], roi[6], molarv_vol[5];
  Real s_molal[5], s_mf[5], gama[5], hen[5], phi[6], k_old[6], sm[5], sh[4], gih[4], molarv_deriv[5];
  Real binary_a[6][6], kappa[6][6], vis[6];
  Real s_molef[5], s_massf[5], molef[10], massf[10], non_aqueous_massf, s_molality[5], s_molefrac[5], mf_h2o;
  Real roh2o, hh2o, rov, rol, usv, usl, xv, enthalpy, con[5], g_con, l_con;

  k[1] = 0.0015;
  k[2] = 650.0;
  k[3] = 40.0;
  k[4] = 1250.0;
  k[5] = 10.0;


  omf[1] = 1.0 - c_vect[0] - c_vect[1] - c_vect[2] - c_vect[3] - c_vect[4] - c_vect[5] - c_vect[6] - c_vect[7];
	omf[2] = c_vect[1];
	omf[3] = c_vect[0];
	omf[4] = c_vect[2];
	omf[5] = c_vect[3];
	omf[6] = c_vect[4];
	omf[7] = c_vect[5];
	omf[8] = c_vect[6];
	omf[9] = c_vect[7];

  //cout<<"mixture2p c_vect[0] = "<<c_vect[0]<<" c_vect[1] = "<<c_vect[1]<<" c_vect[2] = "<<c_vect[2]<<endl;

  //cout<<"c_vect[2] = "<<c_vect[2]<<endl;

	omf_sum = 0.0;
	for(kk=1;kk<=5;kk++)
	{
		omf_sum += omf[kk];
	}


	for(kk=1;kk<=5;kk++)
	{
		tmf[kk] = omf[kk]/omf_sum;
	}

  m_mas[1] = 18.01528e-3;
	m_mas[2] = 16.04e-3;
	m_mas[3] = 44.01e-3;
	m_mas[4] = 28.0134e-3;
	m_mas[5] = 34.1e-3;
	m_mas[6] = 58.44e-3;
	m_mas[7] = 74.5513e-3;
	m_mas[8] = 110.98e-3;
	m_mas[9] = 95.211e-3;

  sm[1] = 0.0;
	sm[2] = 0.0;
	sm[3] = 0.0;
	sm[4] = 0.1;

////////////////////////////////////////////////////////////////////////////////

int it = 1;
int ts = 1;
Real t_test[19], h_test[19];

t_test[1] = 283.15;
t_test[2] = 363.15;

for(it=1; it<=17; it++)
{
  t = t_test[it];
  for(ts=1; ts<=3; ts++)
  {
    CONCENTRATION_CALCULATOR(k, tmf, non_aqueous_mf, x, y);
    FRAC_CONVERTER(omf, tmf, non_aqueous_mf, x, y, m_mas, s_molal, s_molef, s_massf, molef, massf, non_aqueous_massf, s_molality, s_molefrac, mf_h2o);
    PENG_PURE(p, t, tc, pc, mm, a, b);
    PENG_MIX(p, t, y, a, b, capital_a, capital_b, sum_b, sum_a, binary_a, kappa);
    ROOT_FINDER(capital_a, capital_b, z);
    ACTIVITY(p, t, x, tmf, non_aqueous_mf, gama, s_molality);
    HENRY(p, t, hen);
    FUG(capital_a, capital_b, z, y, b, sum_a, sum_b, binary_a, phi);
    K_CONVERGENCE(p, t, hen, gama, phi, k, k_old);
  }

  GAS_MV_DERIV( t, p, molarv_deriv);
  GAS_IN_LIQ( t, p, gas_ha, roi);
  TH_DRIESNER(p, t, s_molef, molef, th);
  water_h = SATURATED_H_WATER(th) + SPECIFIC_V_WATER(p,th)*
  (1.0 - th*T_EXPANSION_WATER(th))*(p - P_SAT_WATER(th))/1000.0;
  water_h *= 1000.0;
  SALT_ENTHALPY(t, s_molality, sh);
  LIQUID_ENTHALPY(s_massf, massf, sh, liq_h, water_h, gas_ha, molarv_deriv);
  water_v( y, p, t, roh2o, hh2o, rov, rol, usv, usl, xv);
  GAS_RO(p, t, z, y, tc, pc, g_ro, roh2o, mf_h2o);
  GAS_H(t, y, z, sum_a, sum_b, capital_b, a, kappa, gas_h, mm,  tc, pc, mf_h2o, hh2o);
  enti( non_aqueous_massf,  gas_h, liq_h, enthalpy);

  h_test[it] = h - enthalpy;

  if(it >= 2)
  {
    t_test[it+1] = t_test[it] - h_test[it]*(t_test[it] - t_test[it-1])/(h_test[it] - h_test[it-1]);
  }

  if(h_test[it]/h>=-0.001 && h_test[it]/h<=0.001)
  {
    break;
  }

}
lhx = liq_h;
mass_frac = non_aqueous_massf;

TV_DRIESNER(p, t, tmf, tv);
GAS_LD( t, molarv_vol);
LIQUID_RHO(p, t, l_rho, s_molality, s_massf, massf, molarv_vol);
GAS_VIS(p, t, tc, pc, y, roi, g_vis, vis, xv);
LIQUID_VIS(t, s_massf, massf, l_vis);
CONDUCTIVITY_GAS( t,  p, g_con, g_vis, con, vis, y);
LIQ_CONDUCTIVITY( t, l_con, s_massf, massf);


  enti( non_aqueous_massf,  g_ro, l_rho, rho_mix);
  enti( non_aqueous_massf,  g_vis, l_vis, viscosity);
  enti( non_aqueous_massf,  g_con, l_con, conductivity);


for(it=1; it<=5; it++)
{
  aux1_arr[it] = g_ro*non_aqueous_mf*k[it]/(1.0 + (k[it] - 1.0)*non_aqueous_mf);
  aux2_arr[it] = l_rho*non_aqueous_mf/(1.0 + (k[it] - 1.0)*non_aqueous_mf);
}

for(it=6; it<=9; it++)
{
  aux1_arr[it] = 0.0;
  aux2_arr[it] = l_rho*non_aqueous_mf/(1.0 + (0.0 - 1.0)*non_aqueous_mf);
}

temper = t;

if((p == 1e7) && (m == 8.0) && (c_vect[0] == 0.3) && (c_vect[1] == 0.0))
{
  //std::cout<<"k[1] =" <<y[1]<<" k[2] ="<<y[2]<<" k[3] ="<<y[3]<<" k[4] ="<<y[4]<<" k[5] ="<<y[5]<<" non_aqueous_mf ="<<non_aqueous_mf<<endl;
  //std::cout<<"k[1] =" <<roi[2]<<" k[2] ="<<roi[3]<<" k[3] ="<<roi[4]<<" k[4] ="<<roi[5]<<" non_aqueous_mf ="<<non_aqueous_mf<<endl;
  //std::cout<<"k[1] =" <<con[1]<<" k[2] ="<<con[2]<<" k[3] ="<<con[3]<<" k[4] ="<<con[4]<< " non_aqueous_mf ="<<non_aqueous_mf<<endl;
  //std::cout<<" g_con ="<<g_con<<endl;
  //std::cout<<" l_con ="<<temper<<endl;
  //cout<<g_ro<<" "<<l_rho<<" "<<lhx<<" "<<gas_h<<" "<<mass_frac<<" "<<non_aqueous_mf<<" "<<rho_mix<<" "<<temper<<" "<<aux1_arr[4]<<" "<<aux2_arr[4]<<" "<<endl;
  //std::cout<<" non_aqueous_massf ="<<non_aqueous_massf<<" g_ro ="<<g_ro<<" l_rho ="<<l_rho<<" rho_mix ="<<rho_mix<<endl;
  //std::cout<<" non_aqueous_massf ="<<non_aqueous_massf<<" g_vis ="<<g_vis<<" l_vis ="<<l_vis<<" viscosity ="<<viscosity<<endl;
  //std::cout<<" non_aqueous_massf ="<<non_aqueous_massf<<" g_con ="<<g_con<<" l_con ="<<l_con<<" conductivity ="<<conductivity<<endl;
  //std::cout<<" non_aqueous_massf ="<<non_aqueous_massf<<" gas_h ="<<gas_h<<" liq_h ="<<liq_h<<" enthalpy ="<<enthalpy<<endl;
}
//g_ro = 100.0;
//l_rho = 500.0;
//lhx = 300000;
//gas_h = 500000;
//mass_frac = 0.5;
//non_aqueous_mf = 0.5;
//rho_mix = 300.0;
//temper = 320.0;
//cout<<g_ro<<" "<<l_rho<<" "<<lhx<<" "<<gas_h<<" "<<mass_frac<<" "<<non_aqueous_mf<<" "<<rho_mix<<" "<<temper<<" "<<aux1_arr[4]<<" "<<aux2_arr[4]<<" "<<endl;
//cout<<aux2_arr[1]<<" "<<aux2_arr[2]<<" "<<aux2_arr[3]<<" "<<aux2_arr[4]<<" "<<aux2_arr[5]<<" "<<endl;
//cout<<l_rho<<" "<<endl;
}
///////^^^^^^^^^^^^^^
void
MoskitoMixture2P::water_v(const Real y[], const Real & p, const Real & t, Real & roh2o,
   Real & hh2o, Real & rov, Real & rol, Real & usv, Real & usl, Real & xv) const
{

  double p_others, psat, tc, rhos, rhosv, tr, aux1, aux2;

	double p_h2o, bi, b_cap, tav, tcc;
	double bet[2][4] = {{-1622.40, 9383.80, -14893.80, 7309.10},{241.57, -980.97, 1482.31, -750.98}};
	int i,j;

	p_h2o = y[1]*p;
	p_others = p - p_h2o;

	psat = P_SAT_WATER(t);

	xv = psat/p_h2o;

	tc = t-273.15;
	bi = 0.0;

	rhos = 999.79684 + 0.068317355*tc - 0.010740248*pow(tc, 2.0) + 0.00082140905*pow(tc, 2.5) - 2.3030988e-5*pow(tc, 3.0);

	tcc = 647.096;

    tav = 1.0 - t/tcc;

	rhosv = (1.0/322.0)*exp(2.03150240*pow(tav, 1.0/3.0) + 2.68302940*pow(tav, 2.0/3.0) + 5.38626492*pow(tav, 4.0/3.0) + 17.2991605*pow(tav, 3.0) + 44.7586581*pow(tav, 37.0/6.0) + 63.9201063*pow(tav, 71.0/6.0));
	rhosv = 1.0/rhosv;

	rov = p/p_h2o*rhosv;

	b_cap = 0.0;

	for(j=0;j<=3;j++)
	{
		b_cap += bet[0][j]*pow(t/647.10, j*1.0);
	}

	rol = rhos/ (1.0 - 0.11725*log((b_cap+p/1000000.0)/(b_cap+psat/1000000.0)));


	roh2o = rov + pow(1.0 - xv, 1.8)*rol;

    SATURATED_H_WATER( t);
    tr = t/tcc;
    aux1 = 64.87678 + 11.76476*pow(log(1.0/tr), 0.35) - 11.94431*pow(tr, -2.0) + 6.29015*pow(tr, -3.0) - 0.99893*pow(tr,-4.0);
    aux2 = sqrt(aux1);
    usv = exp(aux2)*1000.0;

    usl = (SATURATED_H_WATER(t) + SPECIFIC_V_WATER(p,t)*(1.0 - t*T_EXPANSION_WATER(t))*(p - P_SAT_WATER(t))/1000.0)*1000.0;

    hh2o = xv*usv + (1.0 - xv)*usl + p_h2o/roh2o;

    //(const Real y[], const Real & p, const Real & t, Real & roh2o,
      // Real & hh2o, Real & rov, Real & rol, Real & usv, Real & usl, Real & xv) const

    //cout<<roh2o<<" "<<hh2o<<" "<<rov<<" "<<rol<<" "<<usv<<" "<<usl<<" "<<xv<<" "<<endl;

}

///////^^^^^^^^^^^^^^
void
MoskitoMixture2P::GAS_VIS(const Real & p, const Real & t, const Real tc[], const Real pc[],
   const Real y[], const Real ro[], Real & g_vis, Real vis[], const Real & xv) const
{

  int i,j;
	double a[5], c[10] ;
	double ts, sumg, gs, eta0, deta, d11, d21, d64, d81, d82, gcc, ts2, gamas, omega, eta_di, be, roz, ux, g_vis2;

	a[0] = 0.235156;
	a[1] = -0.491266;
	a[2] = 5.211155e-2;
	a[3] = 5.347906e-2;
	a[4] = -1.537102e-2;

	d11 = 0.4071119e-2;
	d21 = 0.7198037e-4;
	d64 = 0.2411697e-16;
	d81 = 0.2971072e-22;
	d82 = -0.1627888e-22;

	ts = t/252.196;

	sumg = 0.0;

	for(i=0;i<=4;i++)
	{
		sumg += a[i]*pow(log(ts), i);
	}

	gs = exp(sumg);

	eta0 = 1.0069*sqrt(t)/gs;

	deta = d11*ro[3] + d21*pow(ro[3], 2.0) + d64*pow(ro[3], 6.0)/pow(ts, 3.0) + d81*pow(ro[3], 8.0) + d82*pow(ro[3], 8.0)/ts ;

	vis[3] = 1.03*(eta0 + deta);
	vis[3] *= 0.000001;

//////////////////////////////////////////

double bcap[5][5], sumb;

   bcap[0][0] = 5.00444e2;
   bcap[0][1] = -5.18209;
   bcap[0][2] = 2.30117e-2;
   bcap[0][3] = -4.19261e-5;
   bcap[0][4] = 2.8e-8;

   bcap[1][0] = -9.03675e-1;
   bcap[1][1] = 4.94166e-3;
   bcap[1][2] = -6.07085e-6;
   bcap[1][3] = 0.0;
   bcap[1][4] = 0.0;

   bcap[2][0] = 5.39163e-2;
   bcap[2][1] = -3.33832e-4;
   bcap[2][2] = 6.91861e-7;
   bcap[2][3] = -4.75340e-10;
   bcap[2][4] = 0.0;

   bcap[3][0] = -1.29414e-4 ;
   bcap[3][1] = 7.61842e-7;
   bcap[3][2] = -1.45795e-9;
   bcap[3][3] = 8.93402e-13;
   bcap[3][4] = 0.0;

   bcap[4][0] = 7.06309e-8 ;
   bcap[4][1] = -3.26629e-10;
   bcap[4][2] = 3.75837e-13;
   bcap[4][3] = 0.0;
   bcap[4][4] = 0.0;

   sumb = 0.0;

   for(i=0 ; i<=4 ; i++)
	{
		for(j=0 ; j<=4 ; j++)
		{
			sumb += bcap[j][i]*pow(t,i)*pow(p/100000.0, j);
		}
	}

	vis[2] = sumb*0.0000001;


//////////////////////////////////////////

    vis[5] = 5.448325 + 0.022148*t*exp((0.002784 + 0.007225/t + 72.73416/t/t)*ro[5]);
    vis[5] *= 0.000001;

//////////////////////////////////////////
    double ai[5], c1, c2, boltz, eps, sig, omegi, dv023, ci[5], ronorm, tpnorm, rho, sum, dv23;

    ai[0] = 0.46649;
    ai[1] = -0.57015;
    ai[2] = 0.19164;
    ai[3] = -0.03708;
    ai[4] = 0.00241;

    c1 = 0.3125e6;
    c2 = 2.0442e-49;
    boltz = 1.38062e-23;
    eps = 138.08483e-23;
    sig = 0.36502496e-9;

    omegi = 0.0;

	for(i=0 ; i<=4 ; i++)
	{
		omegi += ai[i]*pow(log(t*boltz/eps), i);
	}

	dv023 = c1*sqrt(c2*t)/(sig*sig*exp(omegi));

	ci[0] = -20.099970;
    ci[1] = 3.4376416;
    ci[2] = -1.4470051;
    ci[3] = -0.27766561e-01;
    ci[4] = -0.21662362;

    ronorm = 314.0;
    tpnorm = 14.0;

    rho = ro[4]/ronorm;
    sum = ci[0]/(rho-ci[1]) + ci[0]/ci[1];

    for(i=2 ; i<=4 ; i++)
	{
		sum += ci[i]*pow(rho, i-1);
	}

	dv23 = sum*tpnorm;

    vis[4] = dv023 + dv23;
    vis[4] *= 0.000001;

////////////////////////////////////////////////////////////
double m_g[6], phi1[6][6], phi2[6][6], phi[6][6], musv;

	double tt, tcel, miu, kmiu, a0, a1, a2, a3;

//water_water_water_water_water_water_water_water_water_water_water_water

	a0 = 137.37;
	a1 = 5.2842;
	a2 = 0.05594;
	a3 = 246.0;

	tcel = t - 273.15;

	miu =  (tcel + a3)/(a2*tcel*tcel + a1*tcel + a0);
	miu *= 0.001;

	musv = (3.6944e-05*tcel*tcel + 0.0294*tcel + 8.9450)*0.000001;

	double mh2o, zvh2o, zlh2o, fh2o, miuh2o;

	mh2o = 18.015;
	zvh2o = xv*sqrt(mh2o)/(xv*sqrt(mh2o) + (1.0 - xv)*sqrt(mh2o));
	zlh2o = 1.0 - zvh2o;
	fh2o = zvh2o*zvh2o/musv + zlh2o*zlh2o/miu + 2.0*zvh2o*zlh2o/sqrt(musv*miu);
	vis[1] = 1.0/fh2o;

//water_water_water_water_water_water_water_water_water_water_water_water

    m_g[1] = 18.015e-3;
	m_g[2] = 16.04e-3;
	m_g[3] = 44.01e-3;
	m_g[4] = 28.0134e-3;
	m_g[5] = 34.1e-3;

	for(i=1 ; i<=5 ; i++)
	{
		for(j=1 ; j<=5 ; j++)
		{
			phi1[i][j] = 1.0 + sqrt(vis[i]/vis[j])*pow(m_g[j]/m_g[i], 0.25);
			phi2[i][j] = pow(phi1[i][j] , 2.0);
			phi[i][j] = phi2[i][j]/sqrt(8.0*(1.0 + m_g[i]/m_g[j]));
			//cout<<phi[i][j]<<endl;
		}
	}

	double makh[6], s, anti[6];

	for(i=1 ; i<=5 ; i++)
	{
		g_vis = 0.0;
		makh[i] = 0.0 ;
		anti[i] = 0.0 ;
		if (y[i] != 0.0)
		{
			for(j=1 ; j<=5 ; j++)
			{
				if (y[j] != 0.0)
				{
					if (j != i)
					{
						makh[i] += y[j]*phi[i][j];
					}
				}
			}
		makh[i] /= y[i];
		anti[i] = vis[i]/(1.0 + makh[i]);
		}
		//cout<<	anti[i] << endl;
	}

	g_vis = 0.0;

	for(i=1 ; i<=5 ; i++)
	{
		g_vis += anti[i];
	}

	//g_vis = musv;

}

///////^^^^^^^^^^^^^^
void
MoskitoMixture2P::GAS_IN_LIQ(const Real & t, const Real & p, Real gas_ha[], Real ro[]) const
{

  int i,j;
	double ah0[6], ah1[6], ah2[6], ah3[6], ah4[6], h1[6], dadt[6], mmg[6], reft[6], refh[6];
	double t1[6], t2[6], t3[6], t4[6], t5[6], r, h1g[6], h22g[6], h21g[6], h_depart[6];

	r = 8.31446;

	mmg[1] = 18.01528e-3;
	mmg[2] = 16.04e-3;
	mmg[3] = 44.01e-3;
	mmg[4] = 28.0134e-3;
	mmg[5] = 34.1e-3;

	reft[1] = 273.15;
	reft[2] = 111.667;
	reft[3] = 273.15;
	reft[4] = 273.15;
	reft[5] = 212.85;

	refh[1] = 0.0;
	refh[2] = 8440;
	refh[3] = 21340;
	refh[4] = 8080;
	refh[5] = 18700.0;

	ah0[1] = 4.395;
    ah1[1] = -4.186e-3;
    ah2[1] = 1.405e-5;
    ah3[1] = -1.564e-8;
    ah4[1] = 0.632e-11;

    ah0[2] = 4.568;
    ah1[2] = -8.975e-3;
    ah2[2] = 3.631e-5;
    ah3[2] = -3.407e-8;
    ah4[2] = 1.091e-11;

    ah0[3] = 3.259;
    ah1[3] = 1.356e-3;
    ah2[3] = 1.502e-5;
    ah3[3] = -2.374e-8;
    ah4[3] = 1.056e-11;

/////////////////////////////////////////

    ah0[4] = 3.539;
    ah1[4] = -0.261e-3;
    ah2[4] = 0.007e-5;
    ah3[4] = 0.157e-8;
    ah4[4] = -0.099e-11;

    ah0[5] = 4.266;
    ah1[5] = -3.438e-3;
    ah2[5] = 1.319e-5;
    ah3[5] = -1.331e-8;
    ah4[5] = 0.488e-11;


    for(i=2;i<=5;i++)
    {
    	t1[i] = t - reft[i];
        t2[i] = pow(t,2.0) - pow(reft[i],2.0);
        t3[i] = pow(t,3.0) - pow(reft[i],3.0);
        t4[i] = pow(t,4.0) - pow(reft[i],4.0);
        t5[i] = pow(t,5.0) - pow(reft[i],5.0);
	}

    for(i=2;i<=5;i++)
    {
    	h1[i] = (ah0[i]*t1[i] + ah1[i]/2.0*t2[i] + ah2[i]/3.0*t3[i] + ah3[i]/4.0*t4[i] + ah4[i]/5.0*t5[i])*r + refh[i] ;
	}

//############################################################################################################
    double wa[6], pca[6], tca[6], ma[6], ba[6], tra[6], squred_parta[6], aa[6], capital_aa[6], capital_ba[6], za[6];
	double  R = 8.31446, tr, pr;
	double ca[6], mmga[6], vg[6];

	wa[1] = 0.345;
	wa[2] = 0.011;
	wa[3] = 0.228;
	wa[4] = 0.0403;
	wa[5] = 0.1081;

	pca[1] = 22055000;
	pca[2] = 4604000;
	pca[3] = 7382000;
	pca[4] = 3400000;
	pca[5] = 8940000;

	tca[1] = 647.13;
	tca[2] = 190.58;
	tca[3] = 304.19;
	tca[4] = 126.1;
	tca[5] = 373.2;


	for(i=2 ; i<=5 ; i++)
	{
		ma[i] = 0.37464 + 1.54226*wa[i] - 0.26992*pow(wa[i],2.0);
		ba[i] = 0.07780*R*tca[i]/pca[i];
		tra[i] = t/tca[i];
		squred_parta[i] = 1.0+ma[i]*(1.0-sqrt(tra[i]));
		aa[i] = 0.45724*pow(R,2.0)*pow(tca[i],2.0)/pca[i]*pow(squred_parta[i],2.0);
		capital_aa[i] = aa[i]*p/(pow(R*t,2));
	    capital_ba[i] = ba[i]*p/R/t;
	}

	ROOT_FINDER(capital_aa[2], capital_ba[2], za[2]);
	ROOT_FINDER(capital_aa[3], capital_ba[3], za[3]);
	ROOT_FINDER(capital_aa[4], capital_ba[4], za[4]);
	ROOT_FINDER(capital_aa[5], capital_ba[5], za[5]);

	tr = t/tca[1];
	pr = p/pca[1];

	mmga[1] = 18.01528;
	mmga[2] = 16.04;
	mmga[3] = 44.01;
	mmga[4] = 28.0134;
	mmga[5] = 34.1;

	ca[1] = 0.0;
	ca[2] = 3.2;
	ca[3] = pow(tr, 6.27844)*(-8.18682*pow(pr, -0.94163) + 0.94612) + 2.49642;
	ca[4] = 3.9;
	ca[5] = 2.17;

	for(i=2 ; i<=5 ; i++)
	{
		vg[i] = za[i]*8.31446e6*t/p;
	    ro[i] = mmga[i]/(vg[i] + ca[i])*1000.0;
	}
//############################################################################################################

    for(i=2;i<=5;i++)
    {
    	h1g[i] = r*t*(za[i] - 1.0);
        h22g[i] = log((za[i] + 2.414*capital_ba[i])/(za[i] - 0.414*capital_ba[i]));
    	dadt[i] = -0.45724*pow(r*tca[i], 2.0)/pca[i]*(1.0 + ma[i]*(1.0 - sqrt(t/tca[i])))*(ma[i]/sqrt(t*tca[i]));
    	h21g[i] = (t*dadt[i] - aa[i])/(2.0*sqrt(2.0)*ba[i]);
	    h_depart[i] = h1g[i] + h21g[i]*h22g[i];
	    gas_ha[i] = (h1[i] + h_depart[i])/mmg[i];
	}

}

///////^^^^^^^^^^^^^^
void
MoskitoMixture2P::CONDUCTIVITY_GAS(const Real & t, const Real & p, Real & g_con, const Real & g_vis,
   Real con[], const Real vis[], const Real y[]) const
{

  double wa[4], pca[4], tca[4], ma[4], ba[4], tra[4], squred_parta[4], aa[4], capital_aa[4], capital_ba[4], za[4];
	double  R = 8.31446, tr, pr;
	double ca[4], mmg[4], vg[4], ro[4];
	int i;

	wa[1] = 0.228;
	wa[2] = 0.011;
	wa[3] = 0.0403;

	pca[1] = 7382000;
	pca[2] = 4604000;
	pca[3] = 3400000;

	tca[1] = 304.19;
	tca[2] = 190.58;
	tca[3] = 126.1;

	for(i=1 ; i<=3 ; i++)
	{
		ma[i] = 0.37464 + 1.54226*wa[i] - 0.26992*pow(wa[i],2.0);
		ba[i] = 0.07780*R*tca[i]/pca[i];
		tra[i] = t/tca[i];
		squred_parta[i] = 1.0+ma[i]*(1.0-sqrt(tra[i]));
		aa[i] = 0.45724*pow(R,2.0)*pow(tca[i],2.0)/pca[i]*pow(squred_parta[i],2.0);
		capital_aa[i] = aa[i]*p/(pow(R*t,2));
	    capital_ba[i] = ba[i]*p/R/t;
	}

	ROOT_FINDER(capital_aa[1], capital_ba[1], za[1]);
	ROOT_FINDER(capital_aa[2], capital_ba[2], za[2]);
	ROOT_FINDER(capital_aa[3], capital_ba[3], za[3]);

	tr = t/tca[1];
	pr = p/pca[1];

	mmg[1] = 44.01;
	mmg[2] = 16.04;
	mmg[3] = 28.0134;

	ca[1] = pow(tr, 6.27844)*(-8.18682*pow(pr, -0.94163) + 0.94612) + 2.49642;
	ca[2] = 3.2;
	ca[3] = 3.9;

	for(i=1 ; i<=3 ; i++)
	{
		vg[i] = za[i]*8.31446e6*t/p;
	    ro[i] = mmg[i]/(vg[i] + ca[i])*1000.0;
	}

	con[2] = -105.161 + 0.9007*ro[1] + 0.0007*pow(ro[1], 2) + 3.50e-15*pow(ro[1]*t, 3) + 3.76e-10*pow(ro[1], 4) + 0.7500*t + 0.0017*pow(t, 2);
	con[2] = con[2]/sqrt(t)*0.001;

	//cout<<"  ro[1] = "<<ro[1]<<endl;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	double  lambda1, trc, prc, ak[4], at, ror, bro, f[10], sami, u, eu1, sami2, cvideal;

	trc = t/tca[2];
	prc = p/pca[2];

	ak[0] = -0.47313668e-2;
	ak[1] = 9.35882659e-2;
	ak[2] = -5.37104001e-2;
	ak[3] = 1.32157337e-2;

	at = 0.0;

	for(i=0 ; i<=3 ; i++)
	{
		at += ak[i]/pow(trc, i);
	}

	lambda1 = sqrt(trc)/at;

	ror =  ro[2]/162.66;

	bro = 15.94004983*ror + 10.20614596*pow(ror, 2) + 0.47065428*pow(ror, 5);

	con[1] = lambda1 + bro ;
	con[1] *= 0.001;
	/////////////////////////////////////////////////////////////////////////////////
	double delta, tau, lambda0, lambdar, rho_critical, t_critical, logtstar, logomega, mu0, mn2;

	double nk[6]={8.862, 31.11, -73.13, 20.03, -0.7096, 0.2672};
	double tk[6]={0.0, 0.03, 0.2, 0.8, 0.6, 1.9};
	double dk[6]={1.0, 2.0, 3.0, 4.0, 8.0, 10.0};
	double lk[6]={0.0, 0.0, 1.0, 2.0, 2.0, 2.0};
	double gammak[6]={0.0, 0.0, 1.0, 1.0, 1.0};
	double bmu[5]={0.431, -0.4623, 0.08406, 0.005341, -0.00331};

	rho_critical = 313.3;
	t_critical = 126.192;

	delta = ro[3] / rho_critical;
	tau = t_critical / t ;
	logtstar = log(t / 98.94);

	logomega = 0.0;

	for(i=0 ; i<=4 ; i++)
	{
		logomega += bmu[i] * pow(logtstar, i);
	}

	mn2 = 28.01348e-3;

	mu0 = 0.0266958 * sqrt(1000.0 * mn2 * t) / (0.3656 * 0.3656 * exp(logomega));
	mu0 *= 0.000001;

	lambda0 = 1.511 * mu0 * 1.0e6 + 2.117 / tau - 3.332 * pow(tau, -0.7);

	lambdar = 0.0;

	for(i=0 ; i<=5 ; i++)
	{
		lambdar += nk[i] * pow(tau, tk[i]) * pow(delta, dk[i]) * exp(-gammak[i] * pow(delta, lk[i]));
	}

	con[3] = (lambda0 + lambdar) * 1.0e-3;

	con[4] =vis[5]*1200;

/////////////////////////////////$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	double ac[5][5], mc[5], yc[5], ya[5], whole[5];
	int j;

	mc[1] = 16.04;
	mc[2] = 44.01;
	mc[3] = 28.0134;
	mc[4] = 34.1;

	for(i=1 ; i<=4 ; i++)
	{
		yc[i] = y[i+1];
	}

	for(i=1 ; i<=4 ; i++)
	{
		for(j=1 ; j<=4 ; j++)
		{
			ac[i][j] = sqrt(mc[j]/mc[i]);
		}
	}

	for(i=1 ; i<=4 ; i++)
	{
		ya[i] = 0.0 ;
		whole[i] = 0.0 ;
		if (yc[i] != 0.0)
		{
			for(j=1 ; j<=4 ; j++)
			{
				if (yc[j] != 0.0)
				{
					ya[i] += yc[j]*ac[i][j];
				}
			}
		whole[i] = yc[i]*con[i]/ya[i];
		}
		//cout<<	anti[i] << endl;
	}

	g_con = 0.0;

	for(i=1 ; i<=4 ; i++)
	{
		g_con += whole[i];
	}

}

void
MoskitoMixture2P::GAS_MV(const Real & t, const Real & p, Real molarv[]) const
{

  double c[5][12], pbar;

	c[1][1] = 0.83143711e1;
	c[1][2] = -0.72772168e-3;
	c[1][3] = 0.21489858e4;
	c[1][4] = -0.14019672e-4;
	c[1][5] = -0.66743449e6;
	c[1][6] = 0.76985890e-2;
	c[1][7] = -0.50253331e-5;
	c[1][8] = -0.30092013e1;
	c[1][9] = 0.48468502e3;

	c[2][1] = 28.9447706;
	c[2][2] = -0.0354581768;
	c[2][3] = -4770.67077;
	c[2][4] = 1.02782768e-5;
	c[2][5] = 33.8126098;
	c[2][6] = 9.04037140e-3;
	c[2][7] = -1.14934031e-3;
	c[2][8] = -0.307405726;
	c[2][9] = -0.0907301486;
	c[2][10] = 9.32713393e-4;

	c[3][1] = -0.23093813E+02;
	c[3][2] = 0.56048525E-01;
	c[3][3] = 0.98808898E+04;
	c[3][4] = -0.51091621E-04;
	c[3][5] = -0.13220298E+07;
	c[3][6] = -0.49542866E-03;
	c[3][7] = 0.12698747E-05;
	c[3][8] = 0.51411144;
	c[3][9] = -0.64733978E-04;

	c[4][1] = 42.564957;
	c[4][2] = -8.6260377E-2;
	c[4][3] = -6084.3775;
	c[4][4] = 6.8714437E-5;
	c[4][5] = -102.76849;
	c[4][6] = 8.4482895E-4;
	c[4][7] = -1.0590768;
	c[4][8] = 3.5665902E-3;

	pbar = p/100000.0 ;

	molarv[1] = c[1][1] + c[1][2]*t + c[1][3]/t + c[1][4]*t*t + c[1][5]/t/t + c[1][6]*pbar + c[1][7]*pbar*t + c[1][8]*pbar/t + c[1][9]*pbar/t/t ;
	molarv[2] = c[2][1] + c[2][2]*t + c[2][3]/t + c[2][4]*t*t + c[2][5]/(630.0 - t) + c[2][6]*pbar + c[2][7]*pbar*log(t) + c[2][8]*pbar/t + c[2][9]*pbar/(630.0 - t) + c[2][10]*pbar*pbar/(630.0 - t)/(630.0 - t) ;
	molarv[3] = c[3][1] + c[3][2]*t + c[3][3]/t + c[3][4]*t*t + c[3][5]/t/t + c[3][6]*pbar + c[3][7]*pbar*t + c[3][8]*pbar/t + c[3][9]*pbar*pbar/t ;
	molarv[4] = c[4][1] + c[4][2]*t + c[4][3]/t + c[4][4]*t*t + c[4][5]/(630.0 - t) + c[4][6]*pbar + c[4][7]*pbar/(630.0 - t) + c[4][8]*pbar*pbar/t ;

}

void
MoskitoMixture2P::GAS_MV_DERIV(const Real & t, const Real & p, Real molarv_deriv[]) const
{

  double c[5][12], molarv_plus[5], molarv_minus[5];
	int i;

	GAS_MV(t+0.001*t, p, molarv_plus);
	GAS_MV(t-0.001*t, p, molarv_minus);

	//cout<<t<<"   "<<molarv_plus[1]<<"   "<<molarv_minus[1]<<endl;

	for(i=1 ; i<=4 ; i++)
	{
		molarv_deriv[i] = (molarv_plus[i] - molarv_minus[i])/2.0/0.001/t;
	}

}

void
MoskitoMixture2P::GAS_LD(const Real & t, Real molarv_vol[]) const
{

  double c[5][5], sum, tc;
	int i,j ;

	c[1][1] = 32.98;
	c[1][2] = 1.648E-1;
	c[1][3] = -1.278E-3;
	c[1][4] = 4.62E-6;

	c[2][1] = 37.51;
	c[2][2] = -9.585E-2;
	c[2][3] = 8.74E-4;
	c[2][4] = -5.044E-7;

	c[3][1] = 0.0;
	c[3][2] = 0.0;
	c[3][3] = 0.0;
	c[3][4] = 0.0;

	c[4][1] = 32.19;
	c[4][2] = 1.13E-1;
	c[4][3] = -6.901E-4;
	c[4][4] = 2.679E-6;


	tc = t - 273.15;

	for(i=1 ; i<=4 ; i++)
	{
		molarv_vol[i] = 0.0;
		for(j=1 ; j<=4 ; j++)
		{
			molarv_vol[i] += c[i][j]*pow(tc, j -1);
		}
		molarv_vol[i] = molarv_vol[i]*0.000001;
	}

}

void
MoskitoMixture2P::enti(const Real & non_aqueous_massf, const Real &  gas_prop, const Real &  liq_prop, Real & ent) const
{

  ent = non_aqueous_massf*gas_prop + (1.0 - non_aqueous_massf)*liq_prop;

}

void
MoskitoMixture2P::LIQ_CONDUCTIVITY(const Real & t, Real & l_con, const Real s_massf[], const Real massf[]) const
{

  double w_lambda, shift[5], sum, frac[5], st;
	int i;

	w_lambda = -0.92247 + 2.8395*(t/273.15) - 1.800*pow(t/273.15, 2) + 0.52577*pow(t/273.15, 3) - 0.07344*pow(t/273.15, 4) ;

	shift[1] = -1.9004e-3;
	shift[2] = -1.9004e-3;
	shift[3] = -2.0856e-3;
	shift[4] = -2.0856e-3;

	sum = s_massf[1] + s_massf[2] + s_massf[3] + s_massf[4] + massf[1];

	st = 0.0;

	for(i=1 ; i<=4 ; i++)
	{
		frac[i] = s_massf[i]/sum;
		st += frac[i]*100.0*shift[i];
	}
	//cout<<st<<endl;

	l_con = w_lambda*(1.0 + st);

}
