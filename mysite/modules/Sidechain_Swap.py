from flask import Blueprint, render_template, request, jsonify, send_from_directory
from werkzeug.utils import secure_filename
import os
from scripts.itp_parser import Itp_parser
from scripts.itpX2 import POLY_Extend
import json

#Blueprint
Sidechain_Swap_bp = Blueprint('Sidechain_Swap', __name__, url_prefix='/Sidechain_Swap')

UPLOAD_FOLDER = os.path.join(os.getcwd(), 'upload')
ALLOWED_EXTENSIONS = {'itp', 'gro'}



os.makedirs(UPLOAD_FOLDER, exist_ok=True)  # Create folder if it doesn't exist


# Define a route for this page
@Sidechain_Swap_bp.route('/')
def new_page():
    return render_template('Sidechain_Swap.html')



# Helper function to check allowed file extensions
def allowed_file(filename):
    """Check if the file has an allowed extension."""
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@Sidechain_Swap_bp.route('/upload', methods=['POST'])
def upload_file():
    """Handle file uploads."""
    if 'file' not in request.files:
        return jsonify({'error': 'No file part'}), 400
    file = request.files['file']
    if file.filename == '':
        return 'No selected file', 400
    if file and allowed_file(file.filename):
        filename = secure_filename(file.filename)
        file.save(os.path.join(UPLOAD_FOLDER, filename))
        result = load_files()  # Process files
        return result, 200
    else:
        return jsonify({'error': 'No selected file'}), 400


@Sidechain_Swap_bp.route('/upload/<filename>', methods=['GET'])
def serve_file(filename):
    json_dir = UPLOAD_FOLDER
    file_path = os.path.join(json_dir, filename)

    if os.path.exists(file_path):
        return send_from_directory(json_dir, filename)
    else:
        print(f"File not found: {file_path}")  # Debugging info
        return jsonify({'error': 'File not found'}), 404



def load_files():
    """Process uploaded ITP and GRO files."""
    directory = UPLOAD_FOLDER
    gro_file = []
    itp_file = []

    for filename in os.listdir(directory):
        path = os.path.join(directory, filename)
        if filename.endswith('.gro'):
            gro_file.append(path)
        elif filename.endswith('.itp'):
            itp_file.append(path)

    if len(gro_file) == 0 or len(itp_file) == 0:
        return 'Please Upload Other File'
    elif len(gro_file) > 1 or len(itp_file) > 1:
        print('Warning: multiple files found')

    try:
        itp = Itp_parser(itp_file[0])
        itp.load_gro(gro_file[0])
        #for file in itp_file + gro_file:
        #    os.remove(file)

        json_file_name = itp_file[0].replace(".itp", ".json")
        json_string = dict(itp)
        json_string['coordinates'] = itp.coordinates
        json_string = json.dumps(json_string)

        with open(json_file_name, "w") as file:
            json.dump(json.loads(json_string), file, indent=0)


        return f'Files uploaded successfully'
    except Exception as e:
        print(e)
        return f'Error loading ITP and GRO to JSON: {str(e)}'
