from flask import Flask, render_template, request
from flask_cors import CORS
import os
import datetime

# Import Blueprints
from modules.GUI import GUI_bp
from modules.Extend import Extend_bp
from modules.Sidechain_Swap import Sidechain_Swap_bp  # Importing the same blueprint as in Extend.py
#from modules.image_upload import image_upload_bp




def create_app():
    """Factory function to create and configure the Flask app."""
    app = Flask(__name__, template_folder='templates', static_folder='static')
    CORS(app, resources={r"/*": {"origins": "*"}})


    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    UPLOAD_FOLDER = os.path.join(BASE_DIR, "upload")

    os.makedirs(UPLOAD_FOLDER, exist_ok=True)
    app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER


    # Register Blueprints
    app.register_blueprint(GUI_bp, url_prefix='/GUI')
    app.register_blueprint(Extend_bp, url_prefix='/Extend')
    app.register_blueprint(Sidechain_Swap_bp, url_prefix='/Sidechain_Swap')

    @app.before_request
    def log_request():
        # Log the access
        ip_address = request.remote_addr  # Get the user's IP address
        requested_url = request.url  # The URL being accessed
        timestamp = datetime.datetime.now()  # The time of the request

        #print(f"Page accessed: {requested_url}")
        #print(f"IP Address: {ip_address}")
        #print(f"Timestamp: {timestamp}")

    # Home route
    @app.route('/')
    def home():
        return render_template('index.html')

    return app

if __name__ == '__main__':
    app = create_app()
    app.run(debug=True, port=5000, use_reloader=False)
else:
    app = create_app()