from flask import Blueprint, render_template, request, jsonify, send_from_directory
from werkzeug.utils import secure_filename
import os
from scripts.itp_parser import Itp_parser
import json
import uuid
import zipfile
from flask import current_app

# Blueprint
GUI_bp = Blueprint('GUI', __name__, url_prefix='/GUI')

ALLOWED_EXTENSIONS = {'itp', 'gro'}

def get_session_folder():
    session_id = str(uuid.uuid4())
    base_upload = current_app.config['UPLOAD_FOLDER']
    path = os.path.join(base_upload, 'GUI', session_id)
    os.makedirs(path, exist_ok=True)
    return path, session_id

# Define a route for this page
@GUI_bp.route('/')
def new_page():
    return render_template('GUI.html')

# Helper function to check allowed file extensions
def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@GUI_bp.route('/upload', methods=['POST'])
def upload_file():
    """Handle file uploads with optional session_id from frontend."""
    if 'file' not in request.files:
        return jsonify({'error': 'No file part'}), 400

    file = request.files['file']
    if file.filename == '':
        return jsonify({'error': 'No selected file'}), 400

    if file and allowed_file(file.filename):
        # Accept or generate session ID
        session_id = request.form.get('session_id')
        if session_id:
            upload_path = os.path.join(current_app.config['UPLOAD_FOLDER'], 'GUI', session_id)
        else:
            upload_path, session_id = get_session_folder()

        os.makedirs(upload_path, exist_ok=True)

        filename = secure_filename(file.filename)
        file_path = os.path.join(upload_path, filename)
        file.save(file_path)

        result = load_files(upload_path)

        return jsonify({
            'message': result,
            'session_id': session_id
        }), 200
    else:
        return jsonify({'error': 'File not allowed'}), 400
    
@GUI_bp.route('/upload/<session_id>/<filename>', methods=['GET'])
def serve_file(session_id, filename):
    json_dir = os.path.join(current_app.config['UPLOAD_FOLDER'], 'GUI', session_id)
    file_path = os.path.join(json_dir, filename)

    if os.path.exists(file_path):
        return send_from_directory(json_dir, filename, as_attachment=True)
    else:
        return jsonify({'error': 'File not found'}), 404

def load_files(directory):
    gro_file = []
    itp_file = []

    for filename in os.listdir(directory):
        path = os.path.join(directory, filename)
        if filename.endswith('.gro'):
            gro_file.append(path)
        elif filename.endswith('.itp'):
            itp_file.append(path)

    if len(gro_file) == 0 or len(itp_file) == 0:
        return 'Please upload both .gro and .itp files'
    elif len(gro_file) > 1 or len(itp_file) > 1:
        print('Warning: multiple files found')

    try:
        itp = Itp_parser(itp_file[0])
        itp.load_gro(gro_file[0])

        json_file_name = os.path.join(directory, "processed.json")
        json_data = dict(itp)
        json_data['coordinates'] = itp.coordinates
        with open(json_file_name, "w") as f:
            json.dump(json_data, f, indent=0)

        return 'Files uploaded and processed successfully'
    except Exception as e:
        print(e)
        return f'Error loading ITP and GRO to JSON: {str(e)}'