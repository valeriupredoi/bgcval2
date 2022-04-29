import numpy as np


def O2sat(pt, ps):
    '''
    Title  : Calculates O2 saturation at 1 atm pressure
    Author : Andrew Yool
    Date   : 14/10/04 (revised 08/07/11)
    
    This subroutine calculates the oxygen saturation concentration
    at 1 atmosphere pressure in mol/m3 given ambient temperature
    and salinity.  This formulation is (ostensibly) taken from
    Garcia & Gordon (1992, L&O, 1307-1312).  The function works
    in the range -1.9 <= T <= 40, 0 <= S <= 42.
    
    Function inputs are :
    nc_T       temperature (degrees C)
    nc_S       salinity (PSU)
    (*) o2_sato  oxygen saturation (mol/m3)
    
    Where (*) is the function output (note its units). 
    
    Check value : T = 10, S = 35, oxy_sato = 0.282015 mol/m3
    '''
    #    ptmin = -1.9
    #    ptmax = 40.
    # Check values in range:
    #    if pt.min() < ptmin or pt.max()>ptmax:
    #   	print "O2 saturation: Temperature outside range. \tmin T.:", pt.min(),'<',ptmin,'\tor max:',pt.max(),'>',ptmax
    #  	pt = np.ma.clip(pt,ptmin,ptmax)

    #    psmin = 0.
    #   psmax = 42.
    # #   # Check values in range:
    #    if ps.min() < psmin or ps.max()>psmax:
    #    	print "O2 saturation: SALINITY outside range. \tmin S.:", ps.min(),'<',psmin,'\tor max:',ps.max(),'>',psmax
    #    	ps = np.ma.clip(ps,psmin,psmax)
    pt = np.ma.array(pt)
    ps = np.ma.array(ps)

    m = pt.mask + ps.mask

    pt = np.ma.masked_where(m, pt)
    ps = np.ma.masked_where(m, ps)

    a0 = 2.00907
    a1 = 3.22014
    a2 = 4.05010
    a3 = 4.94457
    a4 = -2.56847E-1
    a5 = 3.88767
    b0 = -6.24523E-3
    b1 = -7.37614E-3
    b2 = -1.03410E-2
    b3 = -8.17083E-3
    c0 = -4.88682E-7
    #
    tt = 298.15 - pt
    tk = 273.15 + pt
    ts = np.ma.log(tt / tk)
    ts2 = np.ma.power(ts, 2.)
    ts3 = np.ma.power(ts, 3.)
    ts4 = np.ma.power(ts, 4.)
    ts5 = np.ma.power(ts, 5.)
    ans1 = a0 + a1 * ts + a2 * ts2 + a3 * ts3 + a4 * ts4 + a5 * ts5 + ps * (
        b0 + b1 * ts + b2 * ts2 + b3 * ts3) + c0 * (ps * ps)
    ans2 = np.ma.exp(ans1)

    #
    #Convert from ml/l to mmol/m3
    #
    o2_sato = (ans2 / 22391.6) * 1000000.0
    return o2_sato


def testO2sat():
    print("Starting O2 test")
    testvalue = O2sat(10., 35.)
    readlvalue = 282.015
    print(testvalue, "should be", readlvalue)
    if (abs(testvalue - readlvalue) / (readlvalue + testvalue)) > 0.001:
        print("AOU:\ttestO2sat:\t FAIL:\t test failed.", testvalue, '!=',
              readlvalue)
        assert 0


def AOU(Temperature, Salinity, Oxygen):
    '''
    calculate apparentl oxygen usage using temperature, salinity and oxygen.
    '''
    return O2sat(Temperature, Salinity) - Oxygen


if __name__ == "__main__":
    testO2sat()
