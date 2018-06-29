from flask import Flask, request, jsonify
from flask_cors import CORS
import json
import collections

# load the qsched_097 module.
from qsched_097 import qsched as qs

qapp = Flask(__name__)

# opne up to the world.
cors = CORS(qapp,
            resources={
                r"/*":{"origins": "*"}
                    #{"origins": "http://localhost:3000"}
                }
            )
# =========================================================================================

def numify(typestr, value):
    # given a typestring convert the value appropriately.
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
    args = None
    try:
        #args = {k:request.form[k] for k in request.form.keys()}
        body_data = request.data
        qapp.logger.info('REQUEST: ' + str(body_data))
        args = json.loads(body_data)

        # defines how each argument needs to be treated
        numtypes = {
            'dims': 'int',
            'jitter2d': 'float',
            'bins': 'int',
            'bias': 'float',
            'evolution': 'float',
            'linewidth': 'float',
            'linear': 'float'
        }

        # defines the array-like arguments
        twodlist = ['type', 'dims', 'bins', 'bias', 'evolution',
                  'linewidth', 'linear']
        for n in twodlist:
            if n in args:
                args[n] = args[n].replace(","," ").split(" ")

        # and the boolean 0/1 arguments
        bools = ['inclusion', 'backfill', 'appendcorner']
        for n in bools:
            if n in args:
                args[n] = 1 if args[n] else 0

        # numify does the hard work converting strings to their correct types
        # while preserving arras for array-like args.
        for nums in numtypes.keys():
            if nums in args:
                args[nums] = numify(numtypes[nums], args[nums])

        #qapp.logger.info(args)

        # run the scheduler and pass the results back to the client as JSON
        return jsonify(qs(args))
    except Exception as x:
        # if something went wrong, log it and return it to the client as well.
        # it's usually a type error.
        qapp.logger.warning(x)
        return jsonify({"error":"{}".format(x),
                        "args": args})
