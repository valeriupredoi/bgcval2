#
# Copyright 2015, Plymouth Marine Laboratory
#
# This file is part of the bgc-val library.
#
# bgc-val is free software: you can redistribute it and/or modify it
# under the terms of the Revised Berkeley Software Distribution (BSD) 3-clause license.

# bgc-val is distributed in the hope that it will be useful, but
# without any warranty; without even the implied warranty of merchantability
# or fitness for a particular purpose. See the revised BSD license for more details.
# You should have received a copy of the revised BSD license along with bgc-val.
# If not, see <http://opensource.org/licenses/BSD-3-Clause>.
#
# Address:
# Plymouth Marine Laboratory
# Prospect Place, The Hoe
# Plymouth, PL1 3DH, UK
#
# Email:
# ledm@pml.ac.uk
#


def shifttimes(mdata, jobID, year0=False):

    times, datas = [], []
    try:
        t0 = float(sorted(mdata.keys())[0])
    except:
        return times, datas

    if year0 == 'piControl':
        for t in sorted(mdata.keys()):
            if jobID == 'u-ar766':
                t1 = t + (2594 - 1850)  #-1869
            else:
                t1 = t
            times.append(float(t1))
            datas.append(mdata[t])
        return times, datas

    if year0 == 'juggling2':
        for t in sorted(mdata.keys()):
            if jobID == 'u-an869': t1 = t - (2061 - 56)  #-1862
            if jobID == 'u-ao586': t1 = t - (2561 - 556)  #-1869
            if t1 < 540: continue
            if t1 > 700: continue
            times.append(float(t1))
            datas.append(mdata[t])
        return times, datas

    if year0 == 'FullSpinUp':
        for t in sorted(mdata.keys()):
            if jobID == 'u-ak900': t1 = t - 2106
            if jobID == 'u-an869': t1 = t - (2061 - 56) + (3996 - 2106)  #-1862
            if jobID in [
                    'u-ar538',
            ]:
                t1 = t - (2061 - 56) + (3996 - 2106)  #-1862
            if jobID in [
                    'u-ar977',
            ]:
                t1 = t - (2061 - 56) + (3996 - 2106)  #-1862

            if jobID == 'u-ao586':
                t1 = t - (2561 - 556) + (3996 - 2106)  #-1869
            if jobID in [
                    'u-ar783',
            ]: t1 = t + (4965. - 2108.)

            print(jobID, t1, t, t0)
            times.append(float(t1))
            datas.append(mdata[t])
        return times, datas

    if year0 == 'UKESMv1SpinUp':
        for t in sorted(mdata.keys()):
            if jobID == 'u-aj588': t1 = t - 2406
            if jobID == 'u-ak900': t1 = t - 2106
            if jobID == 'u-an869': t1 = t - (2061 - 56) + (3996 - 2106)  #-1862
            if jobID in [
                    'u-ar538',
            ]:
                t1 = t - (2061 - 56) + (3996 - 2106)  #-1862
            #if jobID in ['u-ar977',]: t1 = t - (2061-56) + (3996-2106)      #-1862
            if jobID in [
                    'u-ar783',
            ]: t1 = t + (4965. - 2108.)
            if jobID in ['u-au835', 'u-av450']:
                t1 = t + 4965. - 2108. + 33.  # +233.

            if jobID == 'u-aj588' and t1 > 0.: continue
            if jobID == 'u-ak900' and t1 > 1947.: continue
            if jobID == 'u-an869' and t1 > 4537.: continue
            if jobID == 'u-ar538' and t1 > 4966.: continue
            if jobID == 'u-au835' and t1 > 5269.: continue
            times.append(float(t1))
            datas.append(mdata[t])

        print('shifttimes:\t', year0, jobID, min(times), max(times))
        return times, datas

    if year0 == 'UKESM_CN_control':
        for t in sorted(mdata.keys()):
            # align with u-aw310
            if jobID == 'u-av651':
                if t < 2418: continue
                if t > 2568: continue
                t1 = t - (2418. - 1850.) - 50.

            if jobID == 'u-aw310':  # this is the time basis.
                if t < 1850: continue
                if t > 2070: continue
                t1 = t  #- (2418. - 1850.)

            if jobID == 'u-ay124':
                #                                if t < 1850: continue
                #				if t > 1900: continue
                t1 = t - 29.  # + 20. - 50.

            if jobID == 'u-ay694':
                #                                if t < 1850: continue
                #                                if t > 1930: continue
                t1 = t + 20.

            times.append(float(t1))
            datas.append(mdata[t])
        return times, datas

    if year0 == 'control_2100':
        for t in sorted(mdata.keys()):
            if t < 2090: continue
            if t > 2200: continue
            times.append(float(t))
            datas.append(mdata[t])
        return times, datas

    if year0 in [
            'Historical',
    ]:
        for t in sorted(mdata.keys()):
            if t < 1850: continue
            if t > 2020: continue
            times.append(float(t))
            datas.append(mdata[t])
        return times, datas

    # if year0 in ['UKESM11_Fast_piControl_2',]:
    # Lee: I suggest you adjust the dates so 1991 in ck416 becomes 1850
    # 2743 in by230 also becomes 1850

    for t in sorted(mdata.keys()):
        if jobID == 'u-by230':
            t1 = t - (2743. - 1850.)
            times.append(float(t1))
            datas.append(mdata[t])
        elif jobID == 'u-ck416':
            t1 = t - (1991. - 1850.)
            times.append(float(t1))
            datas.append(mdata[t])
        else:
            assert False
    return times, datas

    if year0 in [
            'UKESM11_historical1',
    ]:
        for t in sorted(mdata.keys()):
            t1 = t
            if jobID == 'u-by230':
                t1 = t - 2811. + 1850.
            if t1 < 1800.: continue
            if t1 > 2050.: continue
            times.append(float(t1))
            datas.append(mdata[t])
        return times, datas

    if year0 in [
            '1849-1901',
    ]:
        for t in sorted(mdata.keys()):
            if t < 1849: continue
            if t > 1901: continue
            times.append(float(t))
            datas.append(mdata[t])
        return times, datas

    if year0 in [
            'AlignToDECK',
            'HistoricalDECK',
            'AlignToDECK1600',
            'AlignToDECK1930',
            'HistoricalDECK1930',
            'AlignToDECK1600-1930',
            'AlignToDECK1950',
            'HistoricalDECK1950',
            'AlignToDECK1600-1950',
            'AlignToDECK2000',
            'HistoricalDECK2000',
            'AlignToDECK1600-2000',
            'AlignToDECK2020',
            'HistoricalDECK2020',
            'AlignToDECK1600-2020',
            'AlignToDECK2050',
            'HistoricalDECK2050',
            'AlignToDECK1600-2050',
            'AlignToDECK2100',
            'HistoricalDECK2100',
            'AlignToDECK1600-2100',
            'AlignToDECK2150',
            'HistoricalDECK2150',
            'AlignToDECK1600-2150',
            'AlignToDECK2200',
            'HistoricalDECK2200',
            'AlignToDECK1600-2200',
            'AlignToDECK2250',
            'HistoricalDECK2250',
            'AlignToDECK1600-2250',
            'AlignToDECK2400',
            'HistoricalDECK2400',
            'AlignToDECK1600-2400',
            'EnsembleAlign',
    ]:
        for t in sorted(mdata.keys()):
            t1 = t
            if jobID in [
                    'u-ar766',
            ]: t1 = t
            elif jobID in [
                    'u-av651',
            ]: t1 = t - 618.
            elif jobID in [
                    'u-av450',
            ]: t1 = t - 527.
            elif jobID in ['u-aq853']: t1 = t - 743.

            elif jobID in ['u-ar783']:
                if t < 2108. - 1.: continue
                if t > 2341. + 1.: continue
                t1 = t - 2341. + 1690.
            elif jobID in ['u-au835']:
                if t < 2308. - 1.: continue
                if t > 2378. + 1.: continue
                t1 = t - 2405. + 1850. - 63.
            elif jobID in ['u-av472']:
                if t < 2378. - 1.: continue
                if t > 2405. + 1.: continue
                t1 = t - 2405. + 1850. - 63.
            elif jobID in ['u-bp705']:
                t1 = t + 50.

            # early cut:
            if year0 in [
                    'AlignToDECK', 'AlignToDECK1930', 'HistoricalDECK1930',
                    'AlignToDECK2020', 'AlignToDECK2050', 'AlignToDECK2150'
            ]:
                if t1 < 1799.: continue

            if year0 in [
                    'AlignToDECK2020',
                    'EnsembleAlign',
            ]:
                if t1 < 1849.: continue

            if year0 in [
                    'AlignToDECK1600-1930', 'AlignToDECK1600',
                    'AlignToDECK1600-2100'
            ]:
                if t1 < 1599.: continue

        # Mid poing cut:
            if jobID in [
                    'u-av651',
            ] and t1 > 1851: continue

            # late cut
            if year0 in [
                    'AlignToDECK1930', 'AlignToDECK1600-1930',
                    'HistoricalDECK1930'
            ]:
                if t1 > 1931.: continue

            if year0 in [
                    'AlignToDECK1950', 'AlignToDECK1600-1950',
                    'HistoricalDECK1950'
            ]:
                if t1 > 1951.: continue

            if year0 in [
                    'AlignToDECK2000', 'AlignToDECK1600-2000',
                    'HistoricalDECK2000'
            ]:
                if t1 > 2001.: continue

            if year0 in [
                    'AlignToDECK2020', 'AlignToDECK1600-2020',
                    'HistoricalDECK2020'
            ]:
                if t1 > 2021.: continue

            if year0 in [
                    'AlignToDECK2050', 'AlignToDECK1600-2050',
                    'HistoricalDECK2050'
            ]:
                if t1 > 2051.: continue

            if year0 in [
                    'AlignToDECK2100', 'AlignToDECK1600-2100',
                    'HistoricalDECK2100'
            ]:
                if t1 > 2101.: continue

            if year0 in [
                    'AlignToDECK2150', 'AlignToDECK1600-2150',
                    'HistoricalDECK2150'
            ]:
                if t1 > 2151.: continue

            if year0 in [
                    'AlignToDECK2200', 'AlignToDECK1600-2200',
                    'HistoricalDECK2200'
            ]:
                if t1 > 2201.: continue

            if year0 in [
                    'AlignToDECK2250', 'AlignToDECK1600-2250',
                    'HistoricalDECK2250'
            ]:
                if t1 > 2251.: continue

            if year0 in [
                    'AlignToDECK2400', 'AlignToDECK1600-2400',
                    'HistoricalDECK2400'
            ]:
                if t1 > 2401.: continue

            times.append(float(t1))
            datas.append(mdata[t])
        #print year0, jobID,'\t',min(mdata.keys()),max(mdata.keys()), '--->',min(times),max(times)
        return times, datas

    if year0 in [
            'ControlAligned',
    ]:
        histruns = {
            'u-aw331': 1850,  #	UKESM1 first historical member (1850)
            'u-ax195': 1880,  #	UKESM1 second Historical member (1880)
            'u-ax589': 1960,  #	UKESM1 third Historical member (1960)
            'u-ax718': 1922,  #	UKESM1 fourth Historical run (1922)
            'u-ay078': 2020,  #	UKESM1 fifth Historical run (2020)
            'u-ay167': 2050,  #      UKESM1 sixth Historical run (2050)
            'u-ay491': 1995,  #	UKESM1 seventh Historical run (1995)
            'u-az942':
            1995,  #     UKESM1 seventh Historical run (1995) restarted
            'u-az021': 1850,
            'u-az417': 1922,
            'u-az418': 1960,
            'u-az513': 2020,  #	UKESM1 fifth Historical run (2020)
            'u-az515': 2050,  #      UKESM1 sixth Historical run (2050)
            'u-az524': 1995,  #	UKESM1 seventh Historical run (1995)
            'u-bb075': 1960,  #
            'u-bb446': 1960,  #
            'u-bb448': 1960,  #
            'u-bb277': 2395,  #	UKESM1 seventh Historical run (1995)
            'u-bc179': 2250,  #      UKESM1 Historical run (2250)
            'u-bc292': 2165,  #      UKESM1 Historical run(2165)
            'u-bc370': 2120,  #      UKESM1 Historical run(2120)
            'u-bc470': 2285,  #      UKESM1 Historical run(2285)
            'u-bd288': 2340,  #	UKESM1 Historical run (2340)
            'u-bd416': 2460,  #	UKESM1 Historical run (2460)
            'u-bd483': 2200,  #	UKESM1 Historical run (2200)
            'u-bf647': 2619,
            'u-bf656': 2716,
            'u-bf703': 2760,
            #'u-bf705': 2815, # dead
            'u-bh162': 2815,
            'u-bf935': 2565,
            'u-bh100': 2685,
            'u-bh101': 2745,
        }

        for t in sorted(mdata.keys()):
            if jobID in [
                    'u-aw310',
            ]: t1 = t
            elif jobID in [
                    'u-av651',
            ]: t1 = t - 618.
            elif jobID in list(histruns.keys()):
                t1 = t + histruns[jobID] - 1850.
            else:
                print("shiftimes.py:\tControlAligned:\t, job not recognised:",
                      jobID, "set to 1850")
                t1 = t
            times.append(float(t1))
            datas.append(mdata[t])
        return times, datas

    if year0 == 'new_emissions':
        for t in sorted(mdata.keys()):
            if jobID in [
                    'u-az508',
            ]: t1 = t - 110
            else: t1 = t
            times.append(float(t1))
            datas.append(mdata[t])
        return times, datas

    if type(year0) in [
            type('str'),
    ]:
        if year0[:10] == 'hist_vs_pi':
            y = float(year0[11:])
            diff = y - 1850.
            for t in sorted(mdata.keys()):
                if jobID in ['u-aw310']: t1 = t - diff
                else: t1 = t

                if t1 > 2020.: continue
                if t1 < 1830.: continue

                times.append(float(t1))
                datas.append(mdata[t])
            return times, datas

    if year0 == 'Staggered':
        for t in sorted(mdata.keys()):
            if jobID in [
                    'u-au362',
            ]: t1 = t + 40.
            elif jobID in [
                    'u-au364',
            ]: t1 = t + 45.
            elif jobID in [
                    'u-au365',
            ]: t1 = t + 50.
            else: t1 = t
            if t1 < 2215.: continue

            times.append(float(t1))
            datas.append(mdata[t])
        return times, datas

    if year0 in [
            'OOvFC',
            'OOvFC1',
            'OOvFC2',
    ]:
        for t in sorted(mdata.keys()):
            if jobID == 'u-ak900': t1 = t - 2106
            if jobID == 'u-an869': t1 = t - (2061 - 56) + (3996 - 2106)  #-1862
            if jobID == 'u-ar538': t1 = t - (2061 - 56) + (3996 - 2106)  #-1862
            if jobID == 'u-ao586':
                t1 = t - (2561 - 556) + (3996 - 2106)  #-1869
            if jobID in [
                    'u-ar783',
            ]: t1 = t + (4965. - 2108.)
            if year0 in [
                    'OOvFC',
                    'OOvFC2',
            ]:
                if t1 < 4940: continue
                if t1 > 5110: continue
            if year0 in [
                    'OOvFC1',
            ]:
                if t1 < 2390: continue
                if t1 > 2860: continue

            print(jobID, t1, t, t0)
            times.append(float(t1))
            datas.append(mdata[t])
        return times, datas

    if year0 == 'Drift':
        for t in sorted(mdata.keys()):
            #                        if jobID in ['u-as462','u-as643',]:	t1 = t  +3312.
            if jobID == 'u-ar977': t1 = t
            else: t1 = t + 3312.
            if t1 < 5400.: continue
            if t1 > 5600.: continue

            times.append(float(t1))
            datas.append(mdata[t])
        return times, datas

    if year0 == 'Drift2':
        for t in sorted(mdata.keys()):
            #                        if jobID in ['u-as462','u-as643',]:    t1 = t  +3312.
            if jobID == 'u-ar977': t1 = t
            else: t1 = t + 3312.
            if t1 < 5400.: continue
            if t1 > 5550.: continue
            times.append(float(t1))
            datas.append(mdata[t])
        return times, datas

    if year0 == 'Drift3':
        for t in sorted(mdata.keys()):
            #                        if jobID in ['u-as462','u-as643',]:    t1 = t  +3312.
            if jobID == 'u-ar977': t1 = t
            else: t1 = t + 3312.
            if t1 < 5405.: continue
            if t1 > 5565.: continue
            if jobID == 'u-as462':
                if t1 > 5438: continue

            if jobID == 'u-as858':
                if t1 > 5485: continue

            times.append(float(t1))
            datas.append(mdata[t])
        return times, datas

    if year0 == 'UKESM0.8':
        for t in sorted(mdata.keys()):
            if jobID == 'u-am004': t1 = t - 1978 - 203
            if jobID == 'u-ao365': t1 = t - 1978 - 33
            if jobID == 'u-ao837': t1 = t - 1978
            if t1 < 0: continue

            print(jobID, t1, t)
            times.append(float(t1))
            datas.append(mdata[t])
        return times, datas

    if year0 == 'UKESM_0.9.1':
        for t in sorted(mdata.keys()):
            if jobID == 'u-ar379': t1 = t + 458.
            else: t1 = t
            times.append(float(t1))
            datas.append(mdata[t])
        return times, datas

    if year0 == 'UKESM_0.9.2':
        for t in sorted(mdata.keys()):
            #                     if jobID == 'u-ar379': t1 = t + 458.
            #                     else: t1=t

            times.append(float(t1))
            datas.append(mdata[t])
        return times, datas

    if year0 == 'u-ai567-minus3':
        t0 = float(sorted(mdata.keys())[0])
        if jobID == 'u-ai567': t0 = t0 + 3
        for t in sorted(mdata.keys()):
            times.append(float(t) - t0)
            datas.append(mdata[t])
        return times, datas

    if year0 == 'ignoreStart':
        t0 = float(sorted(mdata.keys())[0])
        for t in sorted(mdata.keys()):
            if float(t) < 4600.: continue
            times.append(float(t))
            datas.append(mdata[t])
        return times, datas

    if year0 == '2000-2250':
        for t in sorted(mdata.keys()):
            if float(t) < 2000.: continue
            if float(t) > 2250.: continue
            times.append(float(t))
            datas.append(mdata[t])
        return times, datas
    if year0 == '4830-5000':
        for t in sorted(mdata.keys()):
            if float(t) < 4830.: continue
            if float(t) > 5000.: continue
            times.append(float(t))
            datas.append(mdata[t])
        return times, datas
    if year0 in [
            '4945-5110', '4945-5110i', '4945-5110ii', '4945-5110iii',
            '4800-5100'
    ]:
        if jobID in [
                'u-an869', 'u-ar783', 'u-ar538', 'u-ar766', 'u-at629',
                'u-at793', 'u-at628', 'u-at760', 'u-as462', 'u-as858',
                'u-at572', 'u-au027', 'u-au835', 'u-au521', 'u-au756',
                'u-au828', 'u-av450', 'u-av472'
        ]:
            for t in sorted(mdata.keys()):
                if jobID in [
                        'u-at629',
                        'u-at793',
                        'u-at628',
                        'u-at760',
                        'u-as462',
                        'u-as858',
                        'u-at572',
                        'u-au027',
                        'u-au521',
                        'u-au756',
                        'u-au828',
                ]:
                    t1 = t + 2838.
                    #if jobID in ['u-at629','u-at793','u-at628','u-at760',]:#'u-as462','u-as858',]:
                    #				t1 = t + 2794.
                    #elif jobID in ['u-as462','u-as858',]:    t1 = t  +2858.:
                elif jobID in [
                        'u-ar783',
                ]:
                    t1 = t + 4965. - 2108.
                elif jobID in ['u-au835', 'u-av450', 'u-av472']:
                    t1 = t + 4965. - 2108. + 33.  # +233.

                elif jobID in [
                        'u-ar766',
                ]:
                    t1 = t + 3095.
                else:
                    t1 = t

                if year0 in [
                        '4945-5110',
                        '4945-5110i',
                        '4945-5110ii',
                        '4945-5110iii',
                ]:
                    if float(t1) < 4945.: continue
                    if float(t1) > 5410.: continue
                if year0 in [
                        '4800-5100',
                ]:
                    if float(t1) < 4780.: continue
                    #if float(t1) >5410. : continue

                if jobID == 'u-as462':
                    if t > 2126.: continue
                if jobID == 'u-as858':
                    if t > 2175.: continue

#        if t1 > 5485 -3312. + 2794.: continue

                times.append(float(t1))
                datas.append(mdata[t])

        if jobID in [
                'u-am064',
                'u-am927',
                'u-aq853',
                'u-am064i',
                'u-am927i',
                'u-aq853i',
                'u-am064ii',
                'u-am927ii',
                'u-aq853ii',
                'u-am064iii',
                'u-am927iii',
                'u-aq853iii',
        ]:
            for t in sorted(mdata.keys()):
                if jobID in [
                        'u-am064i',
                        'u-am927i',
                        'u-aq853i',
                ]:
                    t1 = t + 2835  # 2130-2230
                if jobID in [
                        'u-am064ii',
                        'u-am927ii',
                        'u-aq853ii',
                ]:
                    t1 = t + 2675  # 2290-2390
                if jobID in [
                        'u-am064iii',
                        'u-am927iii',
                        'u-aq853iii',
                ]:
                    t1 = t + 2515  # 2450-2550

                #if year0 == '4945-5110i':	 t1 = t+2835	# 2130-2230
                #if year0 == '4945-5110ii':      t1 = t+2675	# 2290-2390
                #if year0 == '4945-5110iii':     t1 = t+2515	# 2450-2550
                if jobID == 'u-am927iii':
                    if t1 > 4994: continue
                #if year0 in ['4945-5110','4945-5110i','4945-5110ii','4945-5110iii',]:
                if float(t1) < 4945.: continue
                if float(t1) > 5110.: continue

                times.append(float(t1))
                datas.append(mdata[t])
        return times, datas

    if year0 in [
            'ransom',
    ]:
        if jobID in ['u-am927i', 'u-am927ii', 'u-aq853iii', 'u-aq853']:
            for t in sorted(mdata.keys()):
                if jobID in [
                        'u-am927i',
                ]: t1 = t + 2835  # 2130-2230
                if jobID in [
                        'u-am927ii',
                        'u-aq853iii',
                        'u-aq853',
                ]:
                    t1 = t + 2535  # 2290-2390

                if float(t1) < 4945.: continue
                if float(t1) > 5245.: continue

                times.append(float(t1))
                datas.append(mdata[t])
        else:
            for t in sorted(mdata.keys()):
                if jobID in [
                        'u-at629',
                        'u-at793',
                        'u-at628',
                        'u-at760',
                        'u-as462',
                        'u-as858',
                        'u-at572',
                        'u-au027',
                        'u-au521',
                        'u-au756',
                        'u-au828',
                ]:
                    t1 = t + 2838.
                elif jobID in [
                        'u-ar783',
                ]:
                    t1 = t + 4965. - 2108.
                elif jobID in ['u-au835', 'u-av450', 'u-av472']:
                    t1 = t + 4965. - 2108. + 33.  # +233.
                elif jobID in [
                        'u-ar766',
                ]:
                    t1 = t + 3095.
                else:
                    t1 = t

                if float(t1) < 4945.: continue
                if jobID == 'u-as462' and t > 2126.: continue
                if jobID == 'u-as858' and t > 2175.: continue
                times.append(float(t1))
                datas.append(mdata[t])
        return times, datas
    if year0 in [
            'ransom2',
    ]:
        if jobID in ['u-am927i', 'u-am927ii', 'u-aq853iii', 'u-aq853']:
            for t in sorted(mdata.keys()):
                if jobID in [
                        'u-am927i',
                        'u-am927ii',
                        'u-aq853iii',
                        'u-aq853',
                ]:
                    t1 = t + 2835  # 2130-2230

                if float(t1) < 4945.: continue
                if jobID in [
                        'u-am927i',
                        'u-am927ii',
                ] and t1 > 5315.:
                    continue

                times.append(float(t1))
                datas.append(mdata[t])
        else:
            for t in sorted(mdata.keys()):
                if jobID in [
                        'u-at629',
                        'u-at793',
                        'u-at628',
                        'u-at760',
                        'u-as462',
                        'u-as858',
                        'u-at572',
                        'u-au027',
                        'u-au521',
                        'u-au756',
                        'u-au828',
                ]:
                    t1 = t + 2838.
                elif jobID in [
                        'u-ar783',
                ]:
                    t1 = t + 4965. - 2108.
                elif jobID in ['u-au835', 'u-av450', 'u-av472']:
                    t1 = t + 4965. - 2108. + 33.  # +233.
                elif jobID in [
                        'u-ar766',
                ]:
                    t1 = t + 3095.
                else:
                    t1 = t

                if float(t1) < 4945.: continue
                if jobID == 'u-as462' and t > 2126.: continue
                if jobID == 'u-as858' and t > 2175.: continue
                times.append(float(t1))
                datas.append(mdata[t])
        return times, datas

    if year0 in [
            'Strattrop',
            'Strattrop_fromStart',
    ]:
        for t in sorted(mdata.keys()):
            #if jobID == 'u-as462':
            #        if t > 2126.: continue
            #if jobID == 'u-as858':
            #        if t > 2175.: continue

            if jobID in [
                    'u-at629',
                    'u-at793',
                    'u-at628',
                    'u-at760',
                    'u-as462',
                    'u-as858',
                    'u-at572',
                    'u-au027',
                    'u-au521',
                    'u-au756',
                    'u-au828',
            ]:
                t1 = t + 2838.
            elif jobID in [
                    'u-ar783',
            ]:
                t1 = t + 4965. - 2108.
            elif jobID in ['u-au835', 'u-av450', 'u-av472', 'u-av651']:
                t1 = t + 4965. - 2108. + 33.  # +233.
            elif jobID in [
                    'u-av937',
                    'u-aw072',
                    'u-aw310',
            ]:  #DECK RUNS
                t1 = t + 4965. - 2108. + 33. + 589  # +233.

            elif jobID in [
                    'u-ar766',
            ]:
                t1 = t + 3095.
            else:
                print(jobID)
                assert 0
                t1 = t + 2838.

            if year0 in [
                    'Strattrop',
            ] and float(t1) < 5100.: continue
            #if float(t1) >5410. : continue

            times.append(float(t1))
            datas.append(mdata[t])
        print(year0, jobID, '\t', min(mdata.keys()), max(mdata.keys()), '--->',
              min(times), max(times))
        return times, datas

    if year0 == 'from2228':
        buff = 0
        for t in sorted(mdata.keys()):
            if float(t) < 2228. - buff: continue
            times.append(float(t))
            datas.append(mdata[t])
        return times, datas
    if year0 == 'from2265':
        buff = 0
        for t in sorted(mdata.keys()):
            if float(t) < 2265. - buff: continue
            times.append(float(t))
            datas.append(mdata[t])
        return times, datas
    if year0 == 'from4950':
        buff = 0
        for t in sorted(mdata.keys()):
            if float(t) < 4950. - buff: continue
            times.append(float(t))
            datas.append(mdata[t])
        return times, datas

    if year0 == 'from2350':
        buff = 0
        for t in sorted(mdata.keys()):
            if float(t) < 2350. - buff: continue
            times.append(float(t))
            datas.append(mdata[t])
        return times, datas

    if year0 == 'from2240':
        buff = 0
        for t in sorted(mdata.keys()):
            if float(t) < 2240. - buff: continue
            times.append(float(t))
            datas.append(mdata[t])
        return times, datas

    if year0 == 'from6400':
        buff = 1
        for t in sorted(mdata.keys()):
            if float(t) < 6400. - buff: continue
            times.append(float(t))
            datas.append(mdata[t])
        return times, datas

    if year0 == '2000-2600normu-ak900':
        for t in sorted(mdata.keys()):
            if jobID == 'u-ak900': t = t - 1937.
            if float(t) < 2000.: continue
            if float(t) > 2600.: continue
            times.append(float(t))
            datas.append(mdata[t])
        return times, datas

    if year0 == '2500-3000':
        for t in sorted(mdata.keys()):
            if float(t) < 2500.: continue
            if float(t) > 3000.: continue
            times.append(float(t))
            datas.append(mdata[t])
        return times, datas

    if year0 == '1950-2020':
        for t in sorted(mdata.keys()):
            if float(t) < 1950.: continue
            if float(t) > 2020.: continue
            times.append(float(t))
            datas.append(mdata[t])
        return times, datas

    if year0 == '2300-2400':
        for t in sorted(mdata.keys()):
            if float(t) < 2300.: continue
            if float(t) > 2400.: continue
            times.append(float(t))
            datas.append(mdata[t])
        return times, datas

    if year0 == '2300-2500':
        for t in sorted(mdata.keys()):
            if float(t) < 2300.: continue
            if float(t) > 2500.: continue
            times.append(float(t))
            datas.append(mdata[t])
        return times, datas

    if year0 == '2300-2700':
        for t in sorted(mdata.keys()):
            if float(t) < 2300.: continue
            if float(t) > 2700.: continue
            times.append(float(t))
            datas.append(mdata[t])
        return times, datas

    if year0 == '2300-2900':
        for t in sorted(mdata.keys()):
            if float(t) < 2300.: continue
            if float(t) > 2900.: continue
            times.append(float(t))
            datas.append(mdata[t])
        return times, datas

    if year0 == 'fast4':
        for t in sorted(mdata.keys()):
            if jobID == 'u-bw717':
                t1 = t - 94.
            else:
                t1 = t
            if float(t) > 2400.: continue
            times.append(float(t1))
            datas.append(mdata[t])
        return times, datas

    if year0 == 'N48-ORCA1':
        for t in sorted(mdata.keys()):
            t1 = t
            if jobID == 'u-aw310':
                t1 = t - (1960 - 1850)
            if float(t1) < 1850.: continue
            if float(t1) > 2020.: continue
            times.append(float(t1))
            datas.append(mdata[t])
        return times, datas

#####
# Set all jobs to start at time zero.
    if year0 in [
            'True',
            True,
            'First100Years',
    ]:
        for t in sorted(mdata.keys()):
            if year0 == 'First100Years' and float(t) - t0 > 100.: continue
            times.append(float(t) - t0)
            datas.append(mdata[t])
        return times, datas

    ######
    # No year shift requested, returning sorted arrays.
    times = sorted(mdata.keys())
    datas = [mdata[t] for t in times]
    return times, datas
