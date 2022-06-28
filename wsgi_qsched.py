# https://flask.palletsprojects.com/en/2.1.x/deploying/proxy_fix/
from werkzeug.middleware.proxy_fix import ProxyFix

from server import qapp as app

app.wsgi_app = ProxyFix(
    app.wsgi_app, x_for=1, x_proto=1, x_host=1, x_prefix=1
)

