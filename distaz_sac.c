/** 
 * @file   distaz_sac.c
 * 
 * @brief   Compute distance and azimuth
 * 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <errno.h>


#include "mach.h"
#include "ucf.h"

/** 
 * Compute the distance and azimuth between locations
 * 
 * @param the 
 *    Event latitude, North positive
 * @param phe 
 *    Event longitude, East positive
 * @param ths 
 *    Array of station latitudes
 * @param phs 
 *    Array of station longitudes
 * @param ns 
 *    Length of \p ths and \p phs
 * @param dist 
 *    Output array of epicentral distances in km
 * @param az 
 *    Output array of azimuths in degrees
 * @param baz 
 *    Output array of back azimuths in degrees
 * @param xdeg 
 *    Output array of great circle arc lengths in degrees
 * @param nerr 
 *    Error Return Flag
 *    - 0 on Success
 *    - Non-Zero on Error
 *    - 904 Internal consistency checks failed
 *
 * @note Calculations are based upon the reference spheroid of 1968 and
 *       are defined by the major radius (RAD) and the flattening (FL).
 * 
 * @note Distance is computed using Rudoe's forumula given in GEODESY,
 *       secontion 2.15(b). 
 *       (There is some numerical problem with the following formulae.
 *        If the station is in the southern hemisphere and the event in
 *        in the northern, these equations give the longer, not the
 *        shorter distance between the two locations.  Since the
 *        equations are fairly messy, the simplist solution is to reverse
 *        the meanings of the two locations for this case.)
 *
 * @note Geodesy, G. Bomford, Carlendon Press, 4th ed, 1980, Section 2.15(b)
 *       page 121
 *
 * @bug Documentation for computing distance and azimumth should be added
 * 
 * @bug Variables here are declared static, which is really not what they
 *      should be declared.  The actual values are difficult to identify
 *      what they are referring to as well.
 *
 * @date   830603:  Fixed bug with negative station latiudes.
 * @date   810000:  Original version.
 *
 */

double powi(double b, long x)
{
  double temp;
  long i;


  if ( b == 0.0 ) return( (double) 0.0 );
  if ( x == 0L ) return( (double) 1.0 ) ;

  if ( x > 0L ) {
    temp = b;
    for ( i = x-1; i > 0L; i-- ) temp *= b;
    return temp;
  }

  if ( x < 0L ) {
    temp = 1.0 / b;
    for ( i = x+1; i < 0L; i++ ) temp *= (1.0/b);
    return temp;
  }
}

int main(int argc, char *argv[])
{
	double the, phe;
	double ths, phs;
	if (argc != 5) {
	printf("Usage: distaz_sac evt_lat evt_lon sta_lat sta_lon\n");
	printf("       Returns:  Delta Baz Az\n");
	exit(1);
	} else {
	the = atof(argv[1]);
	phe = atof(argv[2]);
        ths = atof(argv[3]);
        phs = atof(argv[4]);
	}

        double     dist;
        double     az;
        double     baz;
        double     xdeg;
	int laz, lbaz, ldist, lxdeg;
	double a, a1, a12, a12bot, a12top, al, b, b0, b1, c, c0, c1, c2, 
	 c4, cosa12, costhi, costhk, d, d1, dl, du, e, e1, e1p1, e2, e3, 
	 ec2, el, f, f1, g, g1, h, h1, onemec2, p1, p2, pdist, pherad, 
	 phsrad, sc, sd, sina12, sinthi, sinthk, sqrte1p1, ss, t1, t2, 
	 tanthi, tanthk, temp, therad, thg, thsrad, u1, u1bot, u2, u2bot, 
	 u2top, v1, v2, x2, y2, z1, z2;


	static double rad = 6378.160;  /* Earth Radius */
	static double fl = 0.00335293; /* Earth Flattening */
	static double twopideg = 360.; /* Two Pi in Degrees */
	static double c00 = 1.;
	static double c01 = 0.25;
	static double c02 = -4.6875e-02;
	static double c03 = 1.953125e-02;
	static double c21 = -0.125;
	static double c22 = 3.125e-02;
	static double c23 = -1.46484375e-02;
	static double c42 = -3.90625e-03;
	static double c43 = 2.9296875e-03;
	static double degtokm = 111.3199;  /* Conversion from degrees to km */


	/* - Initialize. */
	ec2 = 2.*fl - fl*fl;
	onemec2 = 1. - ec2;
	/* eps = 1. + ec2/onemec2; */

	/* - Convert event location to radians.
	 *   (Equations are unstable for latidudes of exactly 0 degrees.) */


	temp = the;
	if( temp == 0. )
	    temp = 1.0e-08;
	therad = TORAD*temp;
	pherad = TORAD*phe;
        

	/* - Must convert from geographic to geocentric coordinates in order
	 *   to use the spherical trig equations.  This requires a latitude
	 *   correction given by: 1-EC2=1-2*FL+FL*FL */
	if ( the == 90 || the == -90 ) {/* special attention at the poles */
	    thg = the*TORAD ;		/* ... to avoid division by zero. */
        } else {
          thg = atan( onemec2 * tan(therad) );
        }
	d = sin( pherad );
	e = -cos( pherad );
	f = -cos( thg );
	c = sin( thg );
	a = f*e;
	b = -f*d;
	g = -c*e;
	h = c*d;

	    /* -- Convert to radians. */
	    temp = ths;
	    if( temp == 0. )
		temp = 1.0e-08;
	    thsrad = TORAD*temp;
	    phsrad = TORAD*phs;

	    /* -- Calculate some trig constants. */
	    if ( ths == 90 || ths == -90 )
		thg = ths * TORAD ;
	    else
		thg = atan( onemec2*tan( thsrad ) );
	    d1 = sin( phsrad );
	    e1 = -cos( phsrad );
	    f1 = -cos( thg );
	    c1 = sin( thg );
	    a1 = f1*e1;
	    b1 = -f1*d1;
	    g1 = -c1*e1;
	    h1 = c1*d1;
	    sc = a*a1 + b*b1 + c*c1;

	    /* - Spherical trig relationships used to compute angles. */

	    sd = 0.5*sqrt( (powi(a - a1,2) + powi(b - b1,2) + powi(c - 
		 c1,2))*(powi(a + a1,2) + powi(b + b1,2) + powi(c + c1,2)) );
	    xdeg = atan2( sd, sc )*TODEG;
	    if( xdeg < 0. )
	        xdeg = xdeg + twopideg;

	    ss = powi(a1 - d,2) + powi(b1 - e,2) + powi(c1,2) - 2.;
	    sc = powi(a1 - g,2) + powi(b1 - h,2) + powi(c1 - f,2) - 2.;
	    az = atan2( ss, sc )*TODEG;
	    if( az < 0. )
	        az = az + twopideg;

	    ss = powi(a - d1,2) + powi(b - e1,2) + powi(c,2) - 2.;
	    sc = powi(a - g1,2) + powi(b - h1,2) + powi(c - f1,2) - 2.;
	    baz = atan2( ss, sc )*TODEG;
	    if( baz < 0. )
	        baz = baz + twopideg;

	    /* - Now compute the distance between the two points using Rudoe's
	     *   formula given in GEODESY, section 2.15(b).
	     *   (There is some numerical problem with the following formulae.
	     *   If the station is in the southern hemisphere and the event in
	     *   in the northern, these equations give the longer, not the
	     *   shorter distance between the two locations.  Since the
	     *   equations are fairly messy, the simplist solution is to reverse
	     *   the meanings of the two locations for this case.) */

		if( thsrad > 0. ){
		    t1 = thsrad;
		    p1 = phsrad;
		    t2 = therad;
		    p2 = pherad;

		    /* special attention at the poles to avoid atan2 troubles 
		       and division by zero. */
                    if ( the == 90.0 ) {
			costhk = 0.0 ;
			sinthk = 1.0 ;
			tanthk = FLT_MAX ;
		    }
                    else if ( the == -90.0 ) {
			costhk = 0.0 ;
			sinthk = -1.0 ;
			tanthk = -FLT_MAX ;
		    }
                    else {
			costhk = cos( t2 );
			sinthk = sin( t2 );
			tanthk = sinthk/costhk;
                    }

		    /* special attention at the poles continued. */
		    if ( ths == 90.0 ) {
                        costhi = 0.0 ;
                        sinthi = 1.0 ;
                        tanthi = FLT_MAX ;
                    }
		    else if ( ths == -90.0 ) {
                        costhi = 0.0 ;
                        sinthi = -1.0 ;
                        tanthi = -FLT_MAX ;
		    }
		    else {
			costhi = cos( t1 );
			sinthi = sin( t1 );
			tanthi = sinthi/costhi;
		    } 
		}
		else{
		    t1 = therad;
		    p1 = pherad;
		    t2 = thsrad;
		    p2 = phsrad;

		    /* more special attention at the poles */
                    if ( ths == 90.0 ) {
                        costhk = 0.0 ;
                        sinthk = 1.0 ;
                        tanthk = FLT_MAX ;
                    }
                    else if ( ths == -90.0 ) {
                        costhk = 0.0 ;
                        sinthk = -1.0 ;
                        tanthk = -FLT_MAX ;
                    }
		    else {
                        costhk = cos( t2 );
                        sinthk = sin( t2 );
                        tanthk = sinthk/costhk;
                    } 

		    /* more special attention at the poles continued */
                    if ( the == 90.0 ) {
                        costhi = 0.0 ;
                        sinthi = 1.0 ;
                        tanthi = FLT_MAX ;
                    }
                    else if ( the == -90.0 ) {
                        costhi = 0.0 ;
                        sinthi = -1.0 ;
                        tanthi = -FLT_MAX ;
                    }
		    else {
                        costhi = cos( t1 );
                        sinthi = sin( t1 );
                        tanthi = sinthi/costhi;
                    }

		}

		el = ec2/onemec2;
		e1 = 1. + el;
		al = tanthi/(e1*tanthk) + ec2*sqrt( (e1 + powi(tanthi,2))/
		 (e1 + powi(tanthk,2)) );
		dl = p1 - p2;
		a12top = sin( dl );
		a12bot = (al - cos( dl ))*sinthk;

		/* Rewrote these three lines with help from trig identities.  maf 990415 */
		a12 = atan2( a12top, a12bot );
		cosa12 = cos( a12 );
		sina12 = sin( a12 );

		/*cosa12 = sqrt ( a12bot*a12bot / ( a12bot*a12bot + a12top*a12top ) ) ;
		sina12 = sqrt ( a12top*a12top / ( a12bot*a12bot + a12top*a12top ) ) ; */

		e1 = el*(powi(costhk*cosa12,2) + powi(sinthk,2));
		e2 = e1*e1;
		e3 = e1*e2;
		c0 = c00 + c01*e1 + c02*e2 + c03*e3;
		c2 = c21*e1 + c22*e2 + c23*e3;
		c4 = c42*e2 + c43*e3;
		v1 = rad/sqrt( 1. - ec2*powi(sinthk,2) );
		v2 = rad/sqrt( 1. - ec2*powi(sinthi,2) );
		z1 = v1*(1. - ec2)*sinthk;
		z2 = v2*(1. - ec2)*sinthi;
		x2 = v2*costhi*cos( dl );
		y2 = v2*costhi*sin( dl );
		e1p1 = e1 + 1.;
		sqrte1p1 = sqrt( e1p1 );
		u1bot = sqrte1p1*cosa12;
		u1 = atan2( tanthk, u1bot );
		u2top = v1*sinthk + e1p1*(z2 - z1);
		u2bot = sqrte1p1*(x2*cosa12 - y2*sinthk*sina12);
		u2 = atan2( u2top, u2bot );
		b0 = v1*sqrt( 1. + el*powi(costhk*cosa12,2) )/e1p1;
		du = u2 - u1;
		pdist = b0*(c2*(sin( 2.*u2 ) - sin( 2.*u1 )) + c4*(sin( 4.*
		 u2 ) - sin( 4.*u1 )));
		dist = fabs( b0*c0*du + pdist );

printf("%12.8f, %12.8f, %12.8f\n", dist, baz, az);
}
