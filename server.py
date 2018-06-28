from __future__ import division
from flask import Flask, request, jsonify
from flask_cors import CORS
import json
import collections


#from qsched_097 import qsched as qs
import re
import os
import stat
import csv
import sys
#import numpy as np
from numpy import fft
import math
import random
import scipy
from scipy.integrate import quad
from scipy.integrate import fixed_quad

qapp = Flask(__name__)

cors = CORS(qapp,
            resources={
                r"/*":{"origins": "*"}
                    #{"origins": "http://localhost:3000"}
                }
            )
# =========================================================================================

def qs(dict_par):
    #from 0.97
    # dict_par are the parsed arguments as a dictionary

    if dict_par['output_type'] == 'bruker':
        output_type = 0
    elif dict_par['output_type'] == 'jeol':
        output_type = 1
    elif dict_par['output_type'] == 'varian':
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
                if (type == 'noweight') & (len(list) < (number_quantile - 1)):
                    shift = bias * (range_quantile / number_quantile) * (random.randint(-1,1))
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
    # Think further to handling a 3D sched...


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

    if len(dict_par['dims']) == 1:

        # write the schedule
        schedule = []

        dimt1t2.sort()
        if output_type == 0:
                difference = 1
        elif output_type == 1:
                difference = 0
        """
        f = open(file_name,"a") if not file_name is None else sys.stdout
        i = 0
        while i < len(dimt1t2):
                nus1 = dimt1t2[i] - difference
                nus = str(nus1)
                f.write(nus + '\n')
                i+=1
        """
        schedule = [x-difference for x in dimt1t2]


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

# =========================================================================================
def numify(typestr, value):
    if typestr == 'int':
        if isinstance(value, collections.Sequence) and not isinstance(value, str):
            return [numify(typestr, x) for x in value]
        else:
            return int(value)
    elif typestr == 'float':
        if isinstance(value, collections.Sequence) and not isinstance(value, str):
            return [numify(typestr, x) for x in value]
        else:
            return float(value)
    else:
        raise NotImplementedError('type ({}) isn\'t supported'.format(typestr))
@qapp.route('/', methods=["GET"])
def root():
    return jsonify({"status": "OK"})

@qapp.route('/qsched', methods=["POST"])
def qsched():
    qapp.logger.info("form: " + str(request.data))

    try:

        #args = {k:request.form[k] for k in request.form.keys()}
        args = json.loads(request.data)

        numtypes = {
            'dims': 'int',
            'jitter2d': 'float',
            'bins': 'int',
            'bias': 'float',
            'evolution': 'float',
            'linewidth': 'float',
            'linear': 'float'
        }

        twodlist = ['type', 'dims', 'bins', 'bias', 'evolution',
                  'linewidth', 'linear']
        for n in twodlist:
            if n in args:
                args[n] = args[n].replace(","," ").split(" ")

        bools = ['inclusion', 'backfill', 'appendcorner']
        for n in bools:
            if n in args:
                args[n] = 1 if args[n] else 0

        for nums in numtypes.keys():
            if nums in args:
                args[nums] = numify(numtypes[nums], args[nums])

        qapp.logger.info(args)

        return jsonify(qs(args))
    except Exception as x:
        qapp.logger.warning(x)
        return jsonify({"error":"{}".format(x),
                        "args": args})
