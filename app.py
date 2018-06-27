from flask import Flask, request, jsonify
from flask_cors import CORS
import json
import collections

# rename import to avoid name conflict
from qsched_097 import qsched as qs

app = Flask(__name__)

cors = CORS(app,
            resources={
                r"/*":
                    {"origins": "http://localhost:3000"}
                }
            )

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

@app.route('/qsched', methods=["POST"])
def qsched():
    app.logger.info("form: " + str(request.data))

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

        app.logger.info(args)

        return jsonify(qs(args))
    except Exception as x:
        app.logger.warning(x)
        return jsonify({"error":"{}".format(x),
                        "args": args})
