
#!/usr/bin/python


# =========================================================================================
# imports
#
#from __future__ import division
#import re
import pandas as pd
import os
#import stat
#import csv
import sys
import numpy as np
from numpy import fft
import math
import random
import scipy
from scipy.integrate import quad
from scipy.integrate import fixed_quad
import argparse
version = "devel"
pi = math.pi
sys.dont_write_bytecode = True
# =========================================================================================
"""
 This is a public beta release (version 0.97) of a 1D and 2D NUS scheduling routine.
 This version is dated and released April 23, 2018, and supercedes v0.96.

 Changes in  v 0.98: Added 1D backfilling, adjusted example schedules, added repeat testing

 Changes in v 0.97: refactored code to allow easier imports and calling programmatically;
                    no changes made to the algorithm (AMM).
 Changes in v 0.96: runs on all platforms; all-in-one file; 1D schedules can have any or no suffix;
                    minor typo fixed in a warning message; removed some commented code; examples below
                    have been tweaked

 Copyright 2017, 2018 D. Levi Craft; David Rovnyak, Virginia G. Rovnyak

 This program is distributed under the terms of the GNU General Public License (http://www.gnu.org/licenses/),
 which is also provided at the bottom of this file.

 We thank Adam Schuyler of University of Connecticut Health Center, Department of Molecular Biology and
 Biophysics Farmington Ave, for developing the parser used within this code.

 This can and probably does have bugs, but has been found to meet expectations under a variety of conditions.

 USAGE

   example 1D:
       python qsched_098.py -t quant-sin -d 1024 -b 256 -o test.txt -x 1.5 -e 3.0 -ot varian
       python qsched_098.py -t quant-poly -d 64 -b 16 -o test.txt -x 1.5 -e 3.14 -ot bruker
       python qsched_098.py -t linear -d 128 -b 16 -o linear.txt -x 1 -e 2 -ot varian -l 0.7
       python qsched_098py -t noweight -d 32 -b 16 -o constant_time.txt -x .3 -e 3.14 -ot jeol

   example 2D:
       python qsched_098.py -t quant-poly quant-poly -d 96 40 -b 24 12 -o test.txt -j 0.8 -x 1.5 1.5 -e 3.14 1.0 -ot varian --inclusion 0
       python qsched_098.py -t quant-exp quant-exp -d 512 512 -b 48 48 -o test.txt -j 0.7 -x 1.5 1.5 -e 2.0 2.0 -ot jeol  --inclusion 0
       python qsched_098.py -t quant-exp quant-sin -d 64 32 -b 24 16 -o test.txt -j 0.7 -x 1. 1. -e 2.5 2.5 -ot bruker --inclusion 1 --backfill 0

  For additional directions, see the provided user guide (qschedXXX_Guide_Distrib.docx)
"""
# =========================================================================================

def qsched(dict_par):
    # dict_par are the parsed arguments as a dictionary

    if dict_par['output_type'] == '0-start':
        output_type = 0
    elif dict_par['output_type'] == '1-start':
        output_type = 1
    elif dict_par['output_type'] == '1-start':
        output_type = 1
    # =========================================================================================

    inclusion = dict_par['inclusion']
    backfill = dict_par['backfill']
    append_corner = dict_par['appendcorner']
    file_name = dict_par['out'] if 'out' in dict_par else None
    if 'guassian' in dict_par['type']:
        lw = dict_par['linewidth']

    #
    # Generate the quantiles
    # dimt1t2 is a list that has t1 quants first, t2 quants second
    #
    dimt1t2 = []
    dim = 0
    tolerance = 10e-10  # required for sine

    while dim < len(dict_par['dims']):
        number_quantile = dict_par['bins'][dim]
        range_quantile = dict_par['dims'][dim]
        type = dict_par['type'][dim]
        if type == 'noweight':
            time = 0.1
        elif type != 'noweight':
            time = dict_par['evolution'][dim]
            if time <= 0.5:
                print("CAUTION: Choosing e < 0.5 for weighted sampling can produce highly non-random schedules. Use a stronger bias or the noweight option instead.")

        bias = dict_par['bias'][dim]
        a = 0
        j1d = dict_par['jitter2d']
        if dim == 0:
            nt1 = number_quantile
            uni1 = range_quantile

        number_areas = number_quantile - 1

        if dim == 1:
            nt2 = number_quantile
            uni2 = range_quantile
            alpha = dict_par['jitter2d']
            if (alpha <0.01) or (alpha>0.99):
                raise ValueError('Jitter for 2D quantiles must be between 0.01 and 0.99; we recommend 0.7')

        if type == 'quant-sin':
            def y(x):
                if bias not in [1, 1.5, 2]:
                    raise ValueError('For quant-sin bias must be 1, 1.5, or 2')
                if bias == 1:
                    return (1 - scipy.sin(x))
                elif bias == 1.5:
                    return ((1 - scipy.sin(x)) ** 1.8)
                elif bias == 2:
                    return ((0.85 - scipy.sin(x)) ** 2)
            time /= 2
        elif type == 'quant-poly':
            if bias not in [1, 1.5, 2]:
                raise ValueError('For quant-poly bias must be 1, 1.5, or 2')
            if bias == 1:
                time /= pi
                def y(x):
                    return ((1 - x) ** 2)
            elif bias == 1.5:
                time = (0.95 * time) / pi
                def y(x):
                    return ((0.95 - x) ** 3)
            elif bias == 2:
                time = (0.92 * time) / pi
                def y(x):
                    return ((0.92 - x) ** 4)
        elif type == 'noweight':
            def y(x):
                return scipy.exp(-x)
        elif type == 'quant-exp':
            def y(x):
                return scipy.exp(-bias * x)
        elif type == 'guassian':
            a_guass = ((math.pi)**2 * (bias * lw)**2)/(4 * scipy.log(2))
            time *= pi
            def y(x):
                return (scipy.exp(-a_guass * x**2))
        elif type == 'linear':
            lin_decay = dict_par['linear'][dim]
            time = lin_decay
            def y(x):
                return (1 - x)
        A_tot, err = quad(y, 0, time)

        # The decision whether or not to remove the following two lines --> leaning towards them being needed
        if type == 'linear':
           time = 1

        # generate quantiles
        while a <= number_areas:
            if a == 0:
                quant = 1
                list = [quant, ]
                dimt1t2.append(quant)
                a += 1
            elif a != 0:
                A_each = a * A_tot / number_areas
                A_calc = 0
                rect_guess = A_each  # find equivalent rectangle first
                b = 1                # then use simple trapezoid to improve first guess
                while b < 100:
                    rect_guess = 2.0 * A_each / (y(rect_guess) + y(0))
                    b +=  1
                t_guess = rect_guess

                while (abs(A_calc - A_each) > tolerance) and (t_guess < time):
                    A_gaus = fixed_quad(y, 0, t_guess, n=16)
                    A_calc = A_gaus[0]
                    del A_gaus
                    delta_t = (A_each - A_calc) / y(t_guess)
                    # delta_t = abs(A_calc - A_each) / y(t_guess)
                    t_guess = t_guess + delta_t

                if t_guess > time:  # problem with some funcs is that very tiny error in integrals
                    t_guess = time  # forces algorithm to search beyond end point; have to catch and stop

                raw_quant = range_quantile * t_guess / (time)
               #if (type == 'noweight') & (len(list) < (number_quantile - 1)):
                if (len(list) < (number_quantile - 1)):
                    shift = j1d * (range_quantile / number_quantile) * (random.randint(-1,1))
                    raw_quant += shift
                quant = int(round(raw_quant))
                if quant == 0:
                    quant += 1
                elif quant > range_quantile:
                    quant = range_quantile
                elif quant != 0:
                    quant = quant
                while True:
                    try:
                        list.index(quant) > -1
                        quant += 1
                    except ValueError:
                        break
                list.append(quant)
                # print(raw_quant)
                dimt1t2.append(quant)
                a += 1
        dim += 1

    # Temporary fix to quantiles being jittered past range_quantile
    # Theoretically this should only occur for noweight
    # However need to create a suitable way to handle when running a 2D sched


    #
    # Used for checking for points exceeding respective range
    #

    if len(dict_par['dims']) == 1:
        for i in dimt1t2:
            if i > range_quantile:
                temporary = i
                while temporary in dimt1t2:
                    temporary += -1
                dimt1t2.remove(i)
                dimt1t2.append(temporary)
    elif len(dict_par['dims']) == 2:
        dimension1 = dimt1t2[:-dict_par['bins'][1]]
        dimension2 = dimt1t2[-dict_par['bins'][1]:]
        for i in dimension1:
            if i > dict_par['dims'][0]:
                temporary = i
                while temporary in dimension1:
                    temporary += -1
                dimension1.remove(i)
                dimension1.append(temporary)
        for i in dimension2:
            if i > dict_par['dims'][1]:
                temporary = i
                while temporary in dimension2:
                    temporary += -1
                dimension2.remove(i)
                dimension2.append(temporary)
        dimension1.sort()
        dimension2.sort()
        dimt1t2 = dimension1 + dimension2
    #
    #
    #

    # If this is a 1-D Schedule...Write it Out
    # Add linearization check here

    if len(dict_par['dims']) == 1:

        # write the schedule
        schedule = []

        dimt1t2.sort()
        if output_type == 0:
                difference = 1
        elif output_type == 1:
                difference = 0
        """
        f = open(file_name,"w") if not file_name is None else sys.stdout
        i = 0
        while i < len(dimt1t2):
                nus1 = dimt1t2[i] - difference
                nus = str(nus1)
                f.write(nus + '\n')
                i+=1
        """

	### Add linearization to dimt1t2 which starts at 1###
        back1d = dict_par['backfill']
        if back1d == 1:
            print('1D Schedule: no linear backfilling')
        elif back1d >1:
            print(dict_par['backfill'])
            print(dimt1t2)
            j=0
            while j < back1d:
                if dimt1t2[j] == j+1:
                    j=j+1
                else:
                    dimt1t2.insert(j,j+1) #don't j++ here because it checks itself on next loop

        schedule = [x-difference for x in dimt1t2]
        print(schedule)
            #print(schedule)
            #schedule.insert(1,-1)
            #print(schedule)

        # for 1D sched, make a psf
        listZF = []
        for i in range(0, range_quantile):
            listZF.append(0)
        for i in range(0, number_quantile):
            temp = dimt1t2[i] - 1
            listZF[temp] = 1
        psf = fft.fft(listZF)
        dwell = 1
        freq =fft.fftfreq(range_quantile, dwell)

        # make psf filename if this is 1D schedule
        """
        if file_name is None:
            file_3 = None
        else:
            if '.' in file_name:
                prefix = file_name.rsplit('.')[0]
                suffix = file_name.rsplit('.')[1]
                #file_3 = file_name.replace('.txt', '_psf.txt')
                file_3 = prefix + '_psf.' + suffix
            else:
                file_3 = file_name +'.psf'
        """

        # optional: this should give peak to sidelobe ratio (PSR)
        # psr_calc = psf.real
        # peak = max(psr_calc)
        # psr_calc_temp = np.setdiff1d(psr_calc, max(psr_calc))
        # lobe = max(psr_calc_temp)
        # psr_ratio = str(peak / lobe)
        # print(psr_ratio)

        """
        f_3 = open(file_3, "a") if not file_3 is None else sys.stdout
        i = 0
        #have to fix sw/2 shift that nu,py does
        while i < len(psf.real):
            psr_val = psf.real[i]
            psr_val = str(psr_val)
            if len(psf.real)/2 -i > 0:
                f_3.write(str(len(psf.real)/2 - i) + '\t')
            else:
                f_3.write(str(3*len(psf.real)/2 - i) + '\t')
            f_3.write(psr_val + '\n')
            i+=1
        raise SystemExit
        """
        def mkpsf(i, val):
            mul = 1 if len(psf.real)/2 -i > 0 else 3
            return (mul * len(psf.real)/2 - i, val)

        psf = [mkpsf(i, val) for i, val in enumerate(psf.real)]

        return schedule, psf

    #
    #
    # Make a shifted dimt1t2
    #
    #
    shiftdimt1t2 = []
    quant_t1t2 = []
    n = 0
    t = 0
    while n < nt1 + nt2:
        t = dimt1t2[n] - 1
        shiftdimt1t2.append(t)
        n += 1

    #
    #
    # Make the midpoints
    #
    #
    i = 0
    midt1t2 = []
    while i < (nt1 - 1):
        j = 0
        while j < (nt2 - 1):
            midt1 = (shiftdimt1t2[i] + shiftdimt1t2[i + 1])/2
            midt2 = (shiftdimt1t2[nt1 + j] + shiftdimt1t2[nt1 + j + 1])/2
            midt1t2.append(midt1)
            midt1t2.append(midt2)
            j += 1
        i += 1


    #
    #
    # jittering: alpha is box size (.1 means 10% jitter inside quantile)
    #
    #
    l = 0
    jitter = []
    tmp = 0
    i = 0
    l = 0
    while i < (nt1 - 1):
        j = 0
        while j < (nt2 - 1):
            boxt1 = (shiftdimt1t2[i + 1] - shiftdimt1t2[i])
            samp = random.random()
            interval = ((samp - 0.5) * alpha * boxt1) + midt1t2[l]
            point1 = round(interval)
            point1 = abs(point1)
            l += 1
            boxt2 = (shiftdimt1t2[nt1 + j + 1] - shiftdimt1t2[nt1 + j])
            samp = random.random()
            interval = ((samp - 0.5) * alpha * boxt2) + midt1t2[l]
            point2 = round(interval)
            point2 = abs(point2)
            jitter.append(point1)
            jitter.append(point2)
            j += 1
            l += 1
        i += 1



    #
    # Now combine the quantiles and the interior points; stay shifted
    #
    if inclusion == 0:
        jitter.append(shiftdimt1t2[0])                 # bottom  left corner
        jitter.append(shiftdimt1t2[nt1])
        jitter.append(shiftdimt1t2[nt1 - 1])           # bottom right corner
        jitter.append(0)
        jitter.append(0)                               # top left corner
        jitter.append(shiftdimt1t2[nt1 + nt2 - 1])
        #jitter.append(shiftdimt1t2[nt1-1])             # top right corner
        #jitter.append(shiftdimt1t2[nt1+nt2-1])
    else:
        j = 0
        k = 0
        while j == 0:
            while k < nt1:
                jitter.append(shiftdimt1t2[k])
                jitter.append(0)
                k += 1
            j = 1
        j = 0
        k = 0
        while j == 0:
            while k < nt2:
                jitter.append(0)
                jitter.append(shiftdimt1t2[nt1 + k])
                k += 1
            j = 1

        jitter.append(shiftdimt1t2[0])              # bottom  left corner
        jitter.append(shiftdimt1t2[nt1])
        #jitter.append(shiftdimt1t2[nt1-1])             # top right corner
        #jitter.append(shiftdimt1t2[nt1+nt2-1])

    #
    #   Have to check for multiples (not just duplicates)
    #   recall jitter is now edge and interiors points
    #   Handle multiples by popping each match off of jitter
    #
    dup = 0
    duplicate = []
    tick1 = 0
    tick2 = 0
    while tick1 < (len(jitter) - 4):      # allows jitter to change size
        id1 = jitter[tick1]               # get first x value
        id2 = jitter[tick1 + 1]           # get first y value
        tick2 = tick1 + 2                 # get ready for next x value
        while tick2 < (len(jitter) - 1):  # loop over all x values
            if (id1 == jitter[tick2]) & (id2 == jitter[tick2 + 1]):  # brute force
                dup += 1
                duplicate.append(id1)
                duplicate.append(id2)
                jitter.pop(tick2)  # multiples now get nuked
                jitter.pop(tick2)  # remember tick2 now points to next
                tick2 -= 2         # reset tick2 if we popped jitter
            tick2 += 2             # loop over all x values
        tick1 += 2                 # first pair done; advance to next pair


    #
    #
    # Make matrix of all cartesian coordinates
    # Note this will still be shifted
    #
    #
    cartGrid = []
    nusGrid = []
    x1 = 0
    y1 = 0
    while x1 < uni1:
        while y1 < uni2:
            cartGrid.append(x1)
            cartGrid.append(y1)
            y1 += 1
        y1 = 0
        x1 += 1

    for x in cartGrid:
        nusGrid.append(x)

    #
    # Jitter has all unique points from quantiles, edges and interior
    # We know that we need 2*nt1*nt2 final points so we need to backfill
    # the difference between length of jitter and 2*nt1*nt2
    #
    # Flag the non-picked points in cartGrid - stay shifted
    #        (0 to uni1-1,  0 to uni2-1)
    # Do this by comparing jitter to cartGrid
    #   Turn each one in to -1,-1 to be popped later
    #
    #
    tick1 = 0
    while tick1 < (len(jitter) - 2):   # search every pair in jitter
        id1 = jitter[tick1]            # get x value
        id2 = jitter[tick1 + 1]        # get y value
        pt = 2 * id1 * uni2 + 2 * id2  # the value of x,y in jitter is used to find it in cartGrid
        pt = int(pt)                   # make integer to use as list locale
        cartGrid[pt] = -1              # set it to (-1,-1) to be removed next
        cartGrid[pt + 1] = -1          #
        tick1 += 2                     # first pair done; advance to next pair


    #
    # And remove flagged points from cartGrid
    #
    tick1 = 0
    while tick1 < (len(cartGrid) - 2):
        if (cartGrid[tick1] == -1) & (cartGrid[tick1 + 1] == -1):
            cartGrid.pop(tick1)
            cartGrid.pop(tick1)
            tick1 -= 2
        tick1 += 2


    #
    #
    # Find norm of each of normalized unpicked points,
    #
    #
    norm = []
    idx = 0
    while idx < (len(cartGrid) - 1):
        temp = ((cartGrid[idx] / uni1)**2 + (cartGrid[idx + 1] / uni2)**2)**(0.5)
        norm.append(temp)
        idx += 2


    #
    #
    # Algorithm  to backfill
    #       get min value of norm
    #       go through cartGrid until you find it
    #             append the point to jitter
    #             pop it from norm (if the min occurs twice the next instant will be found)
    #

    if (backfill ==1):
        fill = (2 * nt2 * nt1 - len(jitter)) / 2

        while fill > 0:
            idx = 0
            nearest = min(norm)
            # print "check nearest" , nearest
            while idx < (len(cartGrid) - 2):
                temp = ((cartGrid[idx] / uni1) ** 2 + (cartGrid[idx + 1] / uni2) ** 2) ** (0.5)
                if temp == nearest:
                    jitter.append(cartGrid[idx])
                    jitter.append(cartGrid[idx + 1])
                    idxtmp = int(idx / 2)
                    norm.pop(idxtmp)  # careful, now norm and cartGrid don't align
                    cartGrid.pop(idx)  # have to pop cartGrid too to realign
                    cartGrid.pop(idx)
                    idx = len(cartGrid)  # exit loop
                else:
                    idx += 2
            fill -= 1


    #
    #
    # Schedule is done, but the points are in wrong order and are still zero-origin (shifted)
    #
    #    Remember we made copy of cartGrid above to use now
    #
    #    Similar to before, flag points in nusGrid with (-1,-1)
    #
    tick1 = 0
    while tick1 < (len(jitter) - 1):   # search every pair in jitter
        id1 = jitter[tick1]            # get x value
        id2 = jitter[tick1 + 1]        # get y value
        pt = 2 * id1 * uni2 + 2 * id2  # the value of x,y in jitter is used to find it in cartGrid
        pt = int(pt)                   # make integer to use as list locale
        nusGrid[pt] = -1               # set it to (-1,-1) to be removed next
        nusGrid[pt + 1] = -1           #
        tick1 += 2                     # first pair done; advance to next pair


    #
    #
    # Now go through nusGrid to pick out points in order
    #
    #
    idx = 0
    nusFinal = []
    while idx < len(nusGrid) - 1:
        getx = nusGrid[idx]
        gety = nusGrid[idx + 1]
        if (getx == -1) & (gety == -1):
            xpt = int(idx / (2 * uni2))
            ypt = int((idx - 2 * xpt * uni2) / 2)
            nusFinal.append(xpt + output_type)  # the +1 shifts back to 1-origin for spectrometer
            nusFinal.append(ypt + output_type)
        idx += 2

    """
    f = open(file_name,"a")
    i = 0

    while i < len(nusFinal):
        nus_1 = nusFinal[i]
        nus_2 = nusFinal[i+1]
        nus1 = str(nus_1)
        nus2 = str(nus_2)
        nus = nus1 + "   " + nus2
        f.write(nus + '\n')
        i+=2

    length = len(nusFinal)
    if (nusFinal[length - 1] != (shiftdimt1t2[nt1+nt2-1] + output_type)) & (append_corner == 1):
        nus_1 = shiftdimt1t2[nt1-1] + output_type
        nus_2 = shiftdimt1t2[nt1+nt2-1] + output_type
        nus1 = str(nus_1)
        nus2 = str(nus_2)
        nus = nus1 + "   " + nus2
        f.write(nus + '\n')
    f.close()
    #print datetime.now() - startTime
    """
    schedule = [(nusFinal[i], nusFinal[i+1]) for i in range(0, len(nusFinal), 2)]
    length = len(nusFinal)
    if (nusFinal[length - 1] != (shiftdimt1t2[nt1+nt2-1] + output_type)) & \
       (append_corner == 1):
        schedule.append( (shiftdimt1t2[nt1-1] + output_type,
                          shiftdimt1t2[nt1+nt2-1] + output_type))

    return schedule

###
### Analyze repeated sequences in 1D Sampling Schedules 
###    D. Rovnyak, L. Cullen, Bucknell University
###    June 2022
###


#make our files for later
results_hist = open("results_hist.txt","w")
results      = open("results.txt","w")


def find_longest(schedlist, filter_repeat, appendseq):
#expand schedule to fill in zeros

    sched = np.array(schedlist)

    x = sched.size
    y0=int(sched[0])
    yfinal = int(sched[x-1])
    print('Sanity check: First = ' + str(y0) + ' Last = ' + str(yfinal))


# cute way to set sampled points to 1
    sched_expand = np.zeros(yfinal+1)
    x=0
    for x in range(sched.size):
        index = sched[x]
        sched_expand[index]=1

    check = 0
    last_array = filter_repeat
    last_length = filter_repeat.size
    last_count = 0
    loop_counter = 0

    arrayX = np.array([])
    arrayY = np.array([])

    a=0
    while loop_counter >= 0:
        iterations = sched_expand.size - filter_repeat.size + 1
        hits = np.zeros(filter_repeat.size + 1)

#run through the iterations (note j and i start at 0)
        for j in range(iterations):
            temp = 0
            for i in range(filter_repeat.size):
                temp += bool(sched_expand[j+i] == filter_repeat[i])
            hits[temp] += 1.0
			# to do: capture where sequences are found, only find position of longest

		#count = hits[len(filter_repeat)] + hits[0]
        count = hits[len(filter_repeat)]

		#make the histogram file
        i=0
        while i<len(hits):
            tmp = hits[i]
            tmpstring = str(tmp)
            ival = str(i)
            results_hist.write(ival + '  ' + tmpstring + '\n')
            i=i+1

	#write out the count
        results.write(str(a) + '\t' + str(count) +'\n')
        arrayY = np.append(arrayY,count)
        a = a +1

	#finding the longest here
        if count == 0:
            loop_counter = -1
        else:
	#save copy of the sequence that had at least one match 
            last_array = filter_repeat
            last_length = filter_repeat.size
            last_count = count
	#now prepare the next sequence
            check += 1
            lenArray = len(filter_repeat)

            if str(appendseq) == "special":
                if filter_repeat[lenArray-1] == 0 :
                    filter_repeat = np.append(filter_repeat,1)
                else:
                    filter_repeat = np.append(filter_repeat,0)
            else:
                filter_repeat=np.append(filter_repeat,appendseq)

	##
	## now we know the longest sequence, but it could occur once or
	## many times; so we could write a function that finds all of the places
	## where that sequence was; let' sit on that; repeats generally occur early,
	## but not always.
	##


    if check == 0:
        print('\n No occurrences of minimal filter found:' + str(last_array))
    else:
        print('\nFound ' + str(last_count) + ' occurrences of length: ' + str(last_length) + '\n' + str(last_array))
    print('/-----------------------------------------------------/')

	#print(arrayX)
	#print(arrayY)

    return arrayY


# =========================================================================================


if __name__=="__main__":
    # for CLI
    def argin2dict(argin, parser):
        if type(argin) is dict:
            return (argin)

        if type(argin) is str:
            args = parser.parse_args(list(argin.split()))
        else:
            args = argin

        dict_par = vars(args)
        dict_par = {key: dict_par[key] for key in dict_par if dict_par[key] != None}
        return (dict_par)

    def ensure_dir(dir_path):
        # DESC: recursive directory creation
        if not os.path.isdir(dir_path):
            if os.path.exists(dir_path):
                print('ERROR: can not create directory, name in use')
            else:
                os.makedirs(dir_path)

    # Note that linear does not require bias, so cannot make bias required=True
    # Also, linear does not require an evolution time, but still ensuring that it is required
    # Catch command line parameters - do sanity checks
    #
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--type", required=True, type = str, nargs = '*',
                        help='type of schedule generator(quant-exp, quant-poly, quant-sin, guass, linear, noweight)')
    parser.add_argument("-d", "--dims", required=True, type=int, nargs='*',
                        help="indel dimensions")
    parser.add_argument("-o", "--out", type=str, default=None,
                        help="name of output file (if ommited print to stdout)")
    parser.add_argument('-j', "--jitter2d", type=float,
                        help="percent jitter for 2D quantiles (0.7 recommended; between .01-.99)")
    parser.add_argument('-b', "--bins", required=True, type=int, nargs='*',
                        help="number of points to be sampled along each dimension")
    parser.add_argument('-x', "--bias", type=float, nargs='*',
                        help="bias strength on exponential")
    parser.add_argument("-e", "--evolution", type=float, nargs='*',
                        help="T2 evolution time for sampling in each indirect dimension")
    parser.add_argument("-ot", "--output_type", required=True, type=str,
                        help='select output type: 0-start or 1-start')
    parser.add_argument("--inclusion", type=int, default=1,
                        help='0 removes this feature, 1 ensures edge forcing')
    parser.add_argument("--backfill", type=int, default=1,
                        help='0 removes this feature, 1 backfills 2D schedules, >1 backfills 1D schedules to value')
    parser.add_argument('-lw',"--linewidth", type=float, nargs = '*',
                        help ="required for guassian option, please enter your linewidth")
    parser.add_argument("-l", "--linear", type=float, nargs ='*',
                        help = 'percent of linear sampling')
    parser.add_argument("--appendcorner", type=int, default=1,
                        help='0 removes this feature, 1 appends top right corner')
    # =========================================================================================
    # This parser is given by Adam Schuyler

    def proc_par(dict_par):
        if 'out' in dict_par:
            dict_par['out'] = os.path.abspath(dict_par['out'])
            ensure_dir(os.path.dirname(dict_par['out']))
        return (dict_par)

    dict_par = proc_par(argin2dict(parser.parse_args(), parser))

    #print(dict_par)
    #If 1D, run qsched here
    if len(dict_par['dims']) == 1:
        schedule, psf = qsched(dict_par)
        print('REPEATS OF TWO')
        filter = np.array([1,0,1])
        addon = np.array([0,1])
        col1 = find_longest(schedule, filter, addon)
        print('REPEATS OF THREE')
        filter = np.array([1,0,0,1])
        addon = np.array([0,0,1])
        col2 = find_longest(schedule, filter, addon)
        print('REPEATS OF FOUR')
        filter = np.array([1,0,0,0,1])
        addon = np.array([0,0,0,1])
        col3 = find_longest(schedule, filter, addon)
        print('REPEATS OF FIVE')
        filter = np.array([1,0,0,0,0,1])
        addon = np.array([0,0,0,0,1])
        col4 = find_longest(schedule, filter, addon)

        print(col1,col2,col3,col4)
       ##
       ## Format Output to Excel with Pandas
       ##

        writer =  pd.ExcelWriter("RLC.xlsx", engine='xlsxwriter')

        b=0
        c=0
        allColumns = ["101","1001","10001", "100001"]
        allColumnsData = [col1, col2, col3, col4]
        for  h in allColumns:
            d =  {h: allColumnsData[c]}
            df = pd.DataFrame(data=d)
            df.to_excel(writer, sheet_name= 'RLC', startcol=b)
            b=b+3
            c=c+1

        workbook = writer.book
        worksheet = writer.sheets['RLC']
        writer.save()



        file_name = dict_par['out'] if 'out' in dict_par else None

        # write the schedule
        with open(file_name,"a") if not file_name is None else sys.stdout as f:
            for val in schedule:
                f.write("{}\n".format(val))

        # write the psf
        if file_name is None:
            file_3 = None
        else:
            if '.' in file_name:
                prefix = file_name.rsplit('.')[0]
                suffix = file_name.rsplit('.')[1]
                #file_3 = file_name.replace('.txt', '_psf.txt')
                file_3 = prefix + '_psf.' + suffix
            else:
                file_3 = file_name +'.psf'

        with open(file_3, "w") if not file_3 is None else sys.stdout as f_3:
            for data in psf:
                f_3.write("{}\t{}\n".format(*data))
        results_hist.close()
        results.close()

        exit()

    elif len(dict_par['dims']) == 2:
        schedule = qsched(dict_par)

        file_name = dict_par['out'] if 'out' in dict_par else None
        # write the schedule
        with open(file_name,"w") if not file_name is None else sys.stdout as f:
            for val in schedule:
                f.write("{}   {}\n".format(*val))
    else:
        raise NotImplementedError("Only 1 and 2D is supported!")



if False:
    pass
    """
    The GPL License itself


                       GNU GENERAL PUBLIC LICENSE
                          Version 3, 29 June 2007

    Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
    Everyone is permitted to copy and distribute verbatim copies
    of this license document, but changing it is not allowed.

                               Preamble

     The GNU General Public License is a free, copyleft license for
     software and other kinds of works.

     The licenses for most software and other practical works are designed
    to take away your freedom to share and change the works.  By contrast,
    the GNU General Public License is intended to guarantee your freedom to
    share and change all versions of a program--to make sure it remains free
    software for all its users.  We, the Free Software Foundation, use the
    GNU General Public License for most of our software; it applies also to
    any other work released this way by its authors.  You can apply it to
    your programs, too.

     When we speak of free software, we are referring to freedom, not
    price.  Our General Public Licenses are designed to make sure that you
    have the freedom to distribute copies of free software (and charge for
    them if you wish), that you receive source code or can get it if you
    want it, that you can change the software or use pieces of it in new
    free programs, and that you know you can do these things.

     To protect your rights, we need to prevent others from denying you
    these rights or asking you to surrender the rights.  Therefore, you have
    certain responsibilities if you distribute copies of the software, or if
    you modify it: responsibilities to respect the freedom of others.

     For example, if you distribute copies of such a program, whether
    gratis or for a fee, you must pass on to the recipients the same
    freedoms that you received.  You must make sure that they, too, receive
    or can get the source code.  And you must show them these terms so they
    know their rights.

     Developers that use the GNU GPL protect your rights with two steps:
    (1) assert copyright on the software, and (2) offer you this License
    giving you legal permission to copy, distribute and/or modify it.

     For the developers' and authors' protection, the GPL clearly explains
    that there is no warranty for this free software.  For both users' and
    authors' sake, the GPL requires that modified versions be marked as
    changed, so that their problems will not be attributed erroneously to
    authors of previous versions.

     Some devices are designed to deny users access to install or run
    modified versions of the software inside them, although the manufacturer
    can do so.  This is fundamentally incompatible with the aim of
    protecting users' freedom to change the software.  The systematic
    pattern of such abuse occurs in the area of products for individuals to
    use, which is precisely where it is most unacceptable.  Therefore, we
    have designed this version of the GPL to prohibit the practice for those
    products.  If such problems arise substantially in other domains, we
    stand ready to extend this provision to those domains in future versions
    of the GPL, as needed to protect the freedom of users.

     Finally, every program is threatened constantly by software patents.
    States should not allow patents to restrict development and use of
    software on general-purpose computers, but in those that do, we wish to
    avoid the special danger that patents applied to a free program could
    make it effectively proprietary.  To prevent this, the GPL assures that
    patents cannot be used to render the program non-free.

     The precise terms and conditions for copying, distribution and
    modification follow.

                          TERMS AND CONDITIONS

     0. Definitions.

     "This License" refers to version 3 of the GNU General Public License.

     "Copyright" also means copyright-like laws that apply to other kinds of
    works, such as semiconductor masks.

     "The Program" refers to any copyrightable work licensed under this
    License.  Each licensee is addressed as "you".  "Licensees" and
    "recipients" may be individuals or organizations.

     To "modify" a work means to copy from or adapt all or part of the work
    in a fashion requiring copyright permission, other than the making of an
    exact copy.  The resulting work is called a "modified version" of the
    earlier work or a work "based on" the earlier work.

     A "covered work" means either the unmodified Program or a work based
    on the Program.

     To "propagate" a work means to do anything with it that, without
    permission, would make you directly or secondarily liable for
    infringement under applicable copyright law, except executing it on a
    computer or modifying a private copy.  Propagation includes copying,
    distribution (with or without modification), making available to the
    public, and in some countries other activities as well.

     To "convey" a work means any kind of propagation that enables other
    parties to make or receive copies.  Mere interaction with a user through
    a computer network, with no transfer of a copy, is not conveying.

     An interactive user interface displays "Appropriate Legal Notices"
    to the extent that it includes a convenient and prominently visible
    feature that (1) displays an appropriate copyright notice, and (2)
    tells the user that there is no warranty for the work (except to the
    extent that warranties are provided), that licensees may convey the
    work under this License, and how to view a copy of this License.  If
    the interface presents a list of user commands or options, such as a
    menu, a prominent item in the list meets this criterion.

     1. Source Code.

     The "source code" for a work means the preferred form of the work
    for making modifications to it.  "Object code" means any non-source
    form of a work.

     A "Standard Interface" means an interface that either is an official
    standard defined by a recognized standards body, or, in the case of
    interfaces specified for a particular programming language, one that
    is widely used among developers working in that language.

     The "System Libraries" of an executable work include anything, other
    than the work as a whole, that (a) is included in the normal form of
    packaging a Major Component, but which is not part of that Major
    Component, and (b) serves only to enable use of the work with that
    Major Component, or to implement a Standard Interface for which an
    implementation is available to the public in source code form.  A
    "Major Component", in this context, means a major essential component
    (kernel, window system, and so on) of the specific operating system
    (if any) on which the executable work runs, or a compiler used to
    produce the work, or an object code interpreter used to run it.

     The "Corresponding Source" for a work in object code form means all
    the source code needed to generate, install, and (for an executable
    work) run the object code and to modify the work, including scripts to
    control those activities.  However, it does not include the work's
    System Libraries, or general-purpose tools or generally available free
    programs which are used unmodified in performing those activities but
    which are not part of the work.  For example, Corresponding Source
    includes interface definition files associated with source files for
    the work, and the source code for shared libraries and dynamically
    linked subprograms that the work is specifically designed to require,
    such as by intimate data communication or control flow between those
    subprograms and other parts of the work.

     The Corresponding Source need not include anything that users
    can regenerate automatically from other parts of the Corresponding
    Source.

     The Corresponding Source for a work in source code form is that
    same work.

    2. Basic Permissions.

     All rights granted under this License are granted for the term of
    copyright on the Program, and are irrevocable provided the stated
    conditions are met.  This License explicitly affirms your unlimited
    permission to run the unmodified Program.  The output from running a
    covered work is covered by this License only if the output, given its
    content, constitutes a covered work.  This License acknowledges your
    rights of fair use or other equivalent, as provided by copyright law.

     You may make, run and propagate covered works that you do not
    convey, without conditions so long as your license otherwise remains
    in force.  You may convey covered works to others for the sole purpose
    of having them make modifications exclusively for you, or provide you
    with facilities for running those works, provided that you comply with
    the terms of this License in conveying all material for which you do
    not control copyright.  Those thus making or running the covered works
    for you must do so exclusively on your behalf, under your direction
    and control, on terms that prohibit them from making any copies of
    your copyrighted material outside their relationship with you.

     Conveying under any other circumstances is permitted solely under
    the conditions stated below.  Sublicensing is not allowed; section 10
    makes it unnecessary.

     3. Protecting Users' Legal Rights From Anti-Circumvention Law.

     No covered work shall be deemed part of an effective technological
    measure under any applicable law fulfilling obligations under article
    11 of the WIPO copyright treaty adopted on 20 December 1996, or
    similar laws prohibiting or restricting circumvention of such
    measures.

     When you convey a covered work, you waive any legal power to forbid
    circumvention of technological measures to the extent such circumvention
    is effected by exercising rights under this License with respect to
    the covered work, and you disclaim any intention to limit operation or
    modification of the work as a means of enforcing, against the work's
    users, your or third parties' legal rights to forbid circumvention of
    technological measures.

     4. Conveying Verbatim Copies.

     You may convey verbatim copies of the Program's source code as you
    receive it, in any medium, provided that you conspicuously and
    appropriately publish on each copy an appropriate copyright notice;
    keep intact all notices stating that this License and any
    non-permissive terms added in accord with section 7 apply to the code;
    keep intact all notices of the absence of any warranty; and give all
    recipients a copy of this License along with the Program.

     You may charge any price or no price for each copy that you convey,
    and you may offer support or warranty protection for a fee.

    5. Conveying Modified Source Versions.

     You may convey a work based on the Program, or the modifications to
    produce it from the Program, in the form of source code under the
    terms of section 4, provided that you also meet all of these conditions:

       a) The work must carry prominent notices stating that you modified
       it, and giving a relevant date.

       b) The work must carry prominent notices stating that it is
       released under this License and any conditions added under section
       7.  This requirement modifies the requirement in section 4 to
       "keep intact all notices".

      c) You must license the entire work, as a whole, under this
       License to anyone who comes into possession of a copy.  This
       License will therefore apply, along with any applicable section 7
       additional terms, to the whole of the work, and all its parts,
       regardless of how they are packaged.  This License gives no
       permission to license the work in any other way, but it does not
       invalidate such permission if you have separately received it.

       d) If the work has interactive user interfaces, each must display
       Appropriate Legal Notices; however, if the Program has interactive
       interfaces that do not display Appropriate Legal Notices, your
       work need not make them do so.

     A compilation of a covered work with other separate and independent
    works, which are not by their nature extensions of the covered work,
    and which are not combined with it such as to form a larger program,
    in or on a volume of a storage or distribution medium, is called an
    "aggregate" if the compilation and its resulting copyright are not
    used to limit the access or legal rights of the compilation's users
    beyond what the individual works permit.  Inclusion of a covered work
    in an aggregate does not cause this License to apply to the other
    parts of the aggregate.

     6. Conveying Non-Source Forms.

     You may convey a covered work in object code form under the terms
    of sections 4 and 5, provided that you also convey the
    machine-readable Corresponding Source under the terms of this License,
    in one of these ways:

       a) Convey the object code in, or embodied in, a physical product
       (including a physical distribution medium), accompanied by the
       Corresponding Source fixed on a durable physical medium
       customarily used for software interchange.

       b) Convey the object code in, or embodied in, a physical product
       (including a physical distribution medium), accompanied by a
       written offer, valid for at least three years and valid for as
       long as you offer spare parts or customer support for that product
       model, to give anyone who possesses the object code either (1) a
       copy of the Corresponding Source for all the software in the
       product that is covered by this License, on a durable physical
       medium customarily used for software interchange, for a price no
       more than your reasonable cost of physically performing this
       conveying of source, or (2) access to copy the
       Corresponding Source from a network server at no charge.

       c) Convey individual copies of the object code with a copy of the
       written offer to provide the Corresponding Source.  This
       alternative is allowed only occasionally and noncommercially, and
       only if you received the object code with such an offer, in accord
       with subsection 6b.

       d) Convey the object code by offering access from a designated
       place (gratis or for a charge), and offer equivalent access to the
       Corresponding Source in the same way through the same place at no
       further charge.  You need not require recipients to copy the
       Corresponding Source along with the object code.  If the place to
       copy the object code is a network server, the Corresponding Source
       may be on a different server (operated by you or a third party)
       that supports equivalent copying facilities, provided you maintain
       clear directions next to the object code saying where to find the
       Corresponding Source.  Regardless of what server hosts the
       Corresponding Source,  you remain obligated to ensure that it is
       available for as long as needed to satisfy these requirements.

       e) Convey the object code using peer-to-peer transmission, provided
       you inform other peers where the object code and Corresponding
       Source of the work are being offered to the general public at no
       charge under subsection 6d.

     A separable portion of the object code, whose source code is excluded
    from the Corresponding Source as a System Library, need not be
    included in conveying the object code work.

     A "User Product" is either (1) a "consumer product", which means any
    tangible personal property which is normally used for personal, family,
    or household purposes, or (2) anything designed or sold for incorporation
    into a dwelling.  In determining whether a product is a consumer product,
    doubtful cases shall be resolved in favor of coverage.  For a particular
    product received by a particular user, "normally used" refers to a
    typical or common use of that class of product, regardless of the status
    of the particular user or of the way in which the particular user
    actually uses, or expects or is expected to use, the product.  A product
    is a consumer product regardless of whether the product has substantial
    commercial, industrial or non-consumer uses, unless such uses represent
    the only significant mode of use of the product.

     "Installation Information" for a User Product means any methods,
    procedures, authorization keys, or other information required to install
    and execute modified versions of a covered work in that User Product from
    a modified version of its Corresponding Source.  The information must
    suffice to ensure that the continued functioning of the modified object
    code is in no case prevented or interfered with solely because
    modification has been made.

     If you convey an object code work under this section in, or with, or
    specifically for use in, a User Product, and the conveying occurs as
    part of a transaction in which the right of possession and use of the
    User Product is transferred to the recipient in perpetuity or for a
    fixed term (regardless of how the transaction is characterized), the
    Corresponding Source conveyed under this section must be accompanied
    by the Installation Information.  But this requirement does not apply
    if neither you nor any third party retains the ability to install
    modified object code on the User Product (for example, the work has
    been installed in ROM).

     The requirement to provide Installation Information does not include a
    requirement to continue to provide support service, warranty, or updates
    for a work that has been modified or installed by the recipient, or for
    the User Product in which it has been modified or installed.  Access to a
    network may be denied when the modification itself materially and
    adversely affects the operation of the network or violates the rules and
    protocols for communication across the network.

     Corresponding Source conveyed, and Installation Information provided,
    in accord with this section must be in a format that is publicly
    documented (and with an implementation available to the public in
    source code form), and must require no special password or key for
    unpacking, reading or copying.

     7. Additional Terms.

     "Additional permissions" are terms that supplement the terms of this
    License by making exceptions from one or more of its conditions.
    Additional permissions that are applicable to the entire Program shall
    be treated as though they were included in this License, to the extent
    that they are valid under applicable law.  If additional permissions
    apply only to part of the Program, that part may be used separately
    under those permissions, but the entire Program remains governed by
    this License without regard to the additional permissions.

     When you convey a copy of a covered work, you may at your option
    remove any additional permissions from that copy, or from any part of
    it.  (Additional permissions may be written to require their own
    removal in certain cases when you modify the work.)  You may place
    additional permissions on material, added by you to a covered work,
    for which you have or can give appropriate copyright permission.
       Notwithstanding any other provision of this License, for material you
    add to a covered work, you may (if authorized by the copyright holders of
    that material) supplement the terms of this License with terms:

       a) Disclaiming warranty or limiting liability differently from the
       terms of sections 15 and 16 of this License; or

       b) Requiring preservation of specified reasonable legal notices or
       author attributions in that material or in the Appropriate Legal
       Notices displayed by works containing it; or

       c) Prohibiting misrepresentation of the origin of that material, or
       requiring that modified versions of such material be marked in
       reasonable ways as different from the original version; or

       d) Limiting the use for publicity purposes of names of licensors or
       authors of the material; or

       e) Declining to grant rights under trademark law for use of some
       trade names, trademarks, or service marks; or

       f) Requiring indemnification of licensors and authors of that
       material by anyone who conveys the material (or modified versions of
       it) with contractual assumptions of liability to the recipient, for
       any liability that these contractual assumptions directly impose on
       those licensors and authors.

     All other non-permissive additional terms are considered "further
    restrictions" within the meaning of section 10.  If the Program as you
    received it, or any part of it, contains a notice stating that it is
    governed by this License along with a term that is a further
    restriction, you may remove that term.  If a license document contains
    a further restriction but permits relicensing or conveying under this
    License, you may add to a covered work material governed by the terms
    of that license document, provided that the further restriction does
    not survive such relicensing or conveying.

     If you add terms to a covered work in accord with this section, you
    must place, in the relevant source files, a statement of the
    additional terms that apply to those files, or a notice indicating
    where to find the applicable terms.

     Additional terms, permissive or non-permissive, may be stated in the
    form of a separately written license, or stated as exceptions;
    the above requirements apply either way.

     8. Termination.

     You may not propagate or modify a covered work except as expressly
    provided under this License.  Any attempt otherwise to propagate or
    modify it is void, and will automatically terminate your rights under
    this License (including any patent licenses granted under the third
    paragraph of section 11).

     However, if you cease all violation of this License, then your
    license from a particular copyright holder is reinstated (a)
    provisionally, unless and until the copyright holder explicitly and
    finally terminates your license, and (b) permanently, if the copyright
    holder fails to notify you of the violation by some reasonable means
    prior to 60 days after the cessation.

     Moreover, your license from a particular copyright holder is
    reinstated permanently if the copyright holder notifies you of the
    violation by some reasonable means, this is the first time you have
    received notice of violation of this License (for any work) from that
    copyright holder, and you cure the violation prior to 30 days after
    your receipt of the notice.

     Termination of your rights under this section does not terminate the
    licenses of parties who have received copies or rights from you under
    this License.  If your rights have been terminated and not permanently
    reinstated, you do not qualify to receive new licenses for the same
    material under section 10.

     9. Acceptance Not Required for Having Copies.

     You are not required to accept this License in order to receive or
    run a copy of the Program.  Ancillary propagation of a covered work
    occurring solely as a consequence of using peer-to-peer transmission
    to receive a copy likewise does not require acceptance.  However,
    nothing other than this License grants you permission to propagate or
    modify any covered work.  These actions infringe copyright if you do
    not accept this License.  Therefore, by modifying or propagating a
    covered work, you indicate your acceptance of this License to do so.

     10. Automatic Licensing of Downstream Recipients.

     Each time you convey a covered work, the recipient automatically
    receives a license from the original licensors, to run, modify and
    propagate that work, subject to this License.  You are not responsible
    for enforcing compliance by third parties with this License.

     An "entity transaction" is a transaction transferring control of an
    organization, or substantially all assets of one, or subdividing an
    organization, or merging organizations.  If propagation of a covered
    work results from an entity transaction, each party to that
    transaction who receives a copy of the work also receives whatever
    licenses to the work the party's predecessor in interest had or could
    give under the previous paragraph, plus a right to possession of the
    Corresponding Source of the work from the predecessor in interest, if
    the predecessor has it or can get it with reasonable efforts.

     You may not impose any further restrictions on the exercise of the
    rights granted or affirmed under this License.  For example, you may
    not impose a license fee, royalty, or other charge for exercise of
    rights granted under this License, and you may not initiate litigation
    (including a cross-claim or counterclaim in a lawsuit) alleging that
    any patent claim is infringed by making, using, selling, offering for
    sale, or importing the Program or any portion of it.

     11. Patents.

     A "contributor" is a copyright holder who authorizes use under this
    License of the Program or a work on which the Program is based.  The
    work thus licensed is called the contributor's "contributor version".

     A contributor's "essential patent claims" are all patent claims
    owned or controlled by the contributor, whether already acquired or
    hereafter acquired, that would be infringed by some manner, permitted
    by this License, of making, using, or selling its contributor version,
    but do not include claims that would be infringed only as a
    consequence of further modification of the contributor version.  For
    purposes of this definition, "control" includes the right to grant
    patent sublicenses in a manner consistent with the requirements of
    this License.

     Each contributor grants you a non-exclusive, worldwide, royalty-free
    patent license under the contributor's essential patent claims, to
    make, use, sell, offer for sale, import and otherwise run, modify and
    propagate the contents of its contributor version.

     In the following three paragraphs, a "patent license" is any express
    agreement or commitment, however denominated, not to enforce a patent
    (such as an express permission to practice a patent or covenant not to
    sue for patent infringement).  To "grant" such a patent license to a
    party means to make such an agreement or commitment not to enforce a
    patent against the party.

     If you convey a covered work, knowingly relying on a patent license,
    and the Corresponding Source of the work is not available for anyone
    to copy, free of charge and under the terms of this License, through a
    publicly available network server or other readily accessible means,
    then you must either (1) cause the Corresponding Source to be so
    available, or (2) arrange to deprive yourself of the benefit of the
    patent license for this particular work, or (3) arrange, in a manner
    consistent with the requirements of this License, to extend the patent
    license to downstream recipients.  "Knowingly relying" means you have
    actual knowledge that, but for the patent license, your conveying the
    covered work in a country, or your recipient's use of the covered work
    in a country, would infringe one or more identifiable patents in that
    country that you have reason to believe are valid.

     If, pursuant to or in connection with a single transaction or
    arrangement, you convey, or propagate by procuring conveyance of, a
    covered work, and grant a patent license to some of the parties
    receiving the covered work authorizing them to use, propagate, modify
    or convey a specific copy of the covered work, then the patent license
    you grant is automatically extended to all recipients of the covered
    work and works based on it.

     A patent license is "discriminatory" if it does not include within
    the scope of its coverage, prohibits the exercise of, or is
    conditioned on the non-exercise of one or more of the rights that are
    specifically granted under this License.  You may not convey a covered
    work if you are a party to an arrangement with a third party that is
    in the business of distributing software, under which you make payment
    to the third party based on the extent of your activity of conveying
    the work, and under which the third party grants, to any of the
    parties who would receive the covered work from you, a discriminatory
    patent license (a) in connection with copies of the covered work
    conveyed by you (or copies made from those copies), or (b) primarily
    for and in connection with specific products or compilations that
    contain the covered work, unless you entered into that arrangement,
    or that patent license was granted, prior to 28 March 2007.

     Nothing in this License shall be construed as excluding or limiting
    any implied license or other defenses to infringement that may
    otherwise be available to you under applicable patent law.

     12. No Surrender of Others' Freedom.

     If conditions are imposed on you (whether by court order, agreement or
    otherwise) that contradict the conditions of this License, they do not
    excuse you from the conditions of this License.  If you cannot convey a
    covered work so as to satisfy simultaneously your obligations under this
    License and any other pertinent obligations, then as a consequence you may
    not convey it at all.  For example, if you agree to terms that obligate you
    to collect a royalty for further conveying from those to whom you convey
    the Program, the only way you could satisfy both those terms and this
    License would be to refrain entirely from conveying the Program.

     13. Use with the GNU Affero General Public License.

     Notwithstanding any other provision of this License, you have
    permission to link or combine any covered work with a work licensed
    under version 3 of the GNU Affero General Public License into a single
    combined work, and to convey the resulting work.  The terms of this
    License will continue to apply to the part which is the covered work,
    but the special requirements of the GNU Affero General Public License,
    section 13, concerning interaction through a network will apply to the
    combination as such.

     14. Revised Versions of this License.

     The Free Software Foundation may publish revised and/or new versions of
    the GNU General Public License from time to time.  Such new versions will
    be similar in spirit to the present version, but may differ in detail to
    address new problems or concerns.

     Each version is given a distinguishing version number.  If the
    Program specifies that a certain numbered version of the GNU General
    Public License "or any later version" applies to it, you have the
    option of following the terms and conditions either of that numbered
    version or of any later version published by the Free Software
    Foundation.  If the Program does not specify a version number of the
    GNU General Public License, you may choose any version ever published
    by the Free Software Foundation.

     If the Program specifies that a proxy can decide which future
    versions of the GNU General Public License can be used, that proxy's
    public statement of acceptance of a version permanently authorizes you
    to choose that version for the Program.

     Later license versions may give you additional or different
    permissions.  However, no additional obligations are imposed on any
    author or copyright holder as a result of your choosing to follow a
    later version.

     15. Disclaimer of Warranty.

     THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
    APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT
    HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY
    OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,
    THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
    PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM
    IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF
    ALL NECESSARY SERVICING, REPAIR OR CORRECTION.

     16. Limitation of Liability.

     IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
    WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS
    THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY
    GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE
    USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF
    DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD
    PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS),
    EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
    SUCH DAMAGES.

     17. Interpretation of Sections 15 and 16.

     If the disclaimer of warranty and limitation of liability provided
    above cannot be given local legal effect according to their terms,
    reviewing courts shall apply local law that most closely approximates
    an absolute waiver of all civil liability in connection with the
    Program, unless a warranty or assumption of liability accompanies a
    copy of the Program in return for a fee.

                        END OF TERMS AND CONDITIONS

               How to Apply These Terms to Your New Programs

     If you develop a new program, and you want it to be of the greatest
    possible use to the public, the best way to achieve this is to make it
    free software which everyone can redistribute and change under these terms.

     To do so, attach the following notices to the program.  It is safest
    to attach them to the start of each source file to most effectively
    state the exclusion of warranty; and each file should have at least
    the "copyright" line and a pointer to where the full notice is found.

       <one line to give the program's name and a brief idea of what it does.>
       Copyright (C) <year>  <name of author>

       This program is free software: you can redistribute it and/or modify
       it under the terms of the GNU General Public License as published by
       the Free Software Foundation, either version 3 of the License, or
       (at your option) any later version.

       This program is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Also add information on how to contact you by electronic and paper mail.

     If the program does terminal interaction, make it output a short
    notice like this when it starts in an interactive mode:

       <program>  Copyright (C) <year>  <name of author>
       This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
       This is free software, and you are welcome to redistribute it
       under certain conditions; type `show c' for details.

    The hypothetical commands `show w' and `show c' should show the appropriate
    parts of the General Public License.  Of course, your program's commands
    might be different; for a GUI interface, you would use an "about box".

     You should also get your employer (if you work as a programmer) or school,
    if any, to sign a "copyright disclaimer" for the program, if necessary.
    For more information on this, and how to apply and follow the GNU GPL, see
    <http://www.gnu.org/licenses/>.

     The GNU General Public License does not permit incorporating your program
    into proprietary programs.  If your program is a subroutine library, you
    may consider it more useful to permit linking proprietary applications with
    the library.  If this is what you want to do, use the GNU Lesser General
    Public License instead of this License.  But first, please read
    <http://www.gnu.org/philosophy/why-not-lgpl.html>.
    """
