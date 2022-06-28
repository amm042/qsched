## about this project

This project provides a serverless web application for the qsched algorithm.

### architecture

The scheduling algorithm is wrapped with Python Flask and runs on the server side.
A web application presents a form for collecting arguments and then makes an
AJAX call to the server to produce output.

#### web frontend

The web frontend is implemented using [ReactJS](https://reactjs.org/). The
project was created with `create-react-app` in the folder `./qsched`. The
form is generated by the Qsched in the file `./qsched/src/qsched/qsched.js`.
This file also stores the URL to the hosted server application. Currently
this is hosted on AWS Lambda. The web frontend itself is deployed and hosted
on AWS S3.

If you submit a pull request to this repo, I will push changes to the production
environment.

#### server
The scheduling algorithm is from qsched_097.py. Since it runs in python with
SciPy/NumPy, it is run on the backend (server). To expose this we wrapped
qsched_097 with a Python Flask wrapper with entry point server.qsched().

To run the server locally for development, create a python virtual environment
and install the modules from requirements.txt.

```
$ python3 -m venv new_venv
$ source new_venv/bin/activate
$ pip install -r requirements.txt
```

Now you are ready to run the server locally.
```
$ export FLASK_APP=server.py
$ export FLASK_ENV=development
$ flask run
```

If you open a browser at [http://localhost:5000](http://localhost:5000) you should see:
```
{
  "status": "OK"
}
```

The scheduling algorithm takes input as a JSON body to a POST request to the `/qsched` endpoint. For example,
```
$ curl -H "Content-Type: application/json" \
  --request POST \
  --data '{"showModal":false, "type":"quant-sin", "dims":"1024", "jitter2d":"0.7",  "bins":"128", "bias":"1.0", "evolution":"2.5", "output_type":"varian",  "inclusion":true, "backfill":true, "linewidth":"1.0", "linear":"0.7",  "appendcorner":true}' \
  http://localhost:5000/qsched
```
will return a JSON object containing the schedule at index 0 and psf data at
index 1.

The Flask server passes the arguments directly to the qsched algorithm as if
they were passed on the command line, so generating 2D schedules can be done by
passing arguments for additional dimensions on the command line as is done
but the original python script.

#### AWS Lambda server

The Flask application is deployed to AWS Lambda via [Zappa](https://github.com/Miserlou/Zappa). The file `zappa_settings.json`
configures Zappa (you need to first setup your AWS credentials).

NOTE: I could not get this working on AWS when `"slim_handler": true`. There
are dependency problems with SciPy/NumPy in this case.

#### Local uWSGI server

In 2022 we transitioned from AWS lambda to a locally hosted solution using uWSGI and flask directly.

see [uWSG](https://flask.palletsprojects.com/en/2.1.x/deploying/uwsgi/).

We have uwsgi installed so, it just has to be run with the qsched flask app:

```
$ uwsgi --http 0.0.0.0:4002 --master -p 2 --enable-threads -w wsgi_qsched:app
```

## about qsched

This is a public beta release of a 1D and 2D NUS scheduling routine.

Copyright 2017, 2018 D. Levi Craft; David Rovnyak, Virginia G. Rovnyak

This program is distributed under the terms of the GNU General Public License (http://www.gnu.org/licenses/),
which is also provided at the bottom of qsched_xxx.py.

We thank Adam Schuyler of University of Connecticut Health Center, Department of Molecular Biology and
Biophysics Farmington Ave, for developing the parser used within this code.

This can and probably does have bugs, but has been found to meet expectations under a variety of conditions.
