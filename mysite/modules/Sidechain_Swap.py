from flask import Blueprint, render_template, request, jsonify, send_from_directory
from werkzeug.utils import secure_filename
import os
from scripts.itp_parser import Itp_parser
from scripts.SwapX import POLY_Swap
import json
import uuid
import zipfile
from flask import current_app
#Blueprint
Sidechain_Swap_bp = Blueprint('Sidechain_Swap', __name__, url_prefix='/Sidechain_Swap')

ALLOWED_EXTENSIONS = {'itp', 'gro'}

def get_session_folder():
    session_id = str(uuid.uuid4())
    base_upload = current_app.config['UPLOAD_FOLDER']
    path = os.path.join(base_upload, 'Swap', session_id)
    os.makedirs(path, exist_ok=True)
    return path, session_id


# Define a route for this page
@Sidechain_Swap_bp.route('/')
def new_page():
    return render_template('Sidechain_Swap.html')

# Helper function to check allowed file extensions
def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@Sidechain_Swap_bp.route('/upload', methods=['POST'])
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
            upload_path = os.path.join(current_app.config['UPLOAD_FOLDER'], 'Swap', session_id)
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
    
@Sidechain_Swap_bp.route('/upload/<session_id>/<filename>', methods=['GET'])
def serve_file(session_id, filename):
    json_dir = os.path.join(current_app.config['UPLOAD_FOLDER'], 'Swap', session_id)
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

@Sidechain_Swap_bp.route('/check_swap', methods=['GET'])
def check_extend():
    result = {'status': 'success', 'message': 'Swap function is ready to proceed'}
    return jsonify(result), 200

@Sidechain_Swap_bp.route('/swap_action', methods=['POST'])
def swap_action():
    try:
        data = request.get_json()
        print(f"data received: {data}")

        sc_1  = list(data.get('sc_1'))
        sc_rep = list(data.get('sc_rep'))

        sc_1 = to_pairs(list(sc_1)) #converts to format POLY_SWAP likes
        sc_rep = to_pairs(list(sc_rep)) #converts to format POLY_SWAP likes

        hs_on_monomer = list(data.get('hs_on_monomer'))
        hs_ontrimer = list(data.get('hs_ontrimer'))
        

        itp_file_path_tri = data.get('itp_file_path_tri')
        gro_file_path_tri = data.get('gro_file_path_tri')
        itp_file_path_1mer = data.get('itp_file_path_1mer')
        gro_file_path_1mer = data.get('gro_file_path_1mer')

        session_id = data.get('session_id')
        
        session_dir =os.path.join(current_app.config['UPLOAD_FOLDER'], 'Swap', session_id)

        
        Swap(
            itp_file_path_tri,
            gro_file_path_tri,
            itp_file_path_1mer,
            gro_file_path_1mer,
            sc_1,
            sc_rep,
            hs_on_monomer,
            hs_ontrimer,
            session_dir=session_dir
        )


        zip_and_cleanup(session_dir, output_zip_name='Results.zip')

        return jsonify({
            'message': 'successfully Swapped!',
            'files': ['Results.zip'],
            'session_id': session_id
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500

def Swap(itp_file_path_tri,gro_file_path_tri,itp_file_path_1mer,gro_file_path_1mer,sc_1,sc_rep,hs_on_monomer,hs_ontrimer,session_dir):
    print(sc_1,sc_rep,hs_on_monomer,hs_ontrimer,itp_file_path_1mer,gro_file_path_1mer,itp_file_path_tri,gro_file_path_tri)
    try:
        
        POLY_Swap(
            itp_file_path_tri=os.path.join(session_dir,itp_file_path_tri),
            gro_file_path_tri=os.path.join(session_dir,gro_file_path_tri),
            itp_file_path_1mer= os.path.join(session_dir,itp_file_path_1mer),
            gro_file_path_1mer= os.path.join(session_dir,gro_file_path_1mer),
            sc_1=sc_1,
            sc_rep=sc_rep,
            hs_on_monomer=hs_on_monomer,
            hs_ontrimer=hs_ontrimer,
            Out_Name=os.path.join(session_dir, "Swapped_POLYMERS")
        )
    except Exception as e:
        print(e)
        return f'Error during extension: {str(e)}'



def zip_and_cleanup(directory_path, output_zip_name='all_files.zip'):
    output_zip_path = os.path.join(directory_path, output_zip_name)
    
    # Create the zip file and add all files
    with zipfile.ZipFile(output_zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for filename in os.listdir(directory_path):
            file_path = os.path.join(directory_path, filename)
            if os.path.isfile(file_path):
                zipf.write(file_path, arcname=filename)
    
    # Delete all files except the zip file itself
    for filename in os.listdir(directory_path):
        file_path = os.path.join(directory_path, filename)
        if os.path.isfile(file_path) and filename != output_zip_name:
            os.remove(file_path)




# Convert flat list → list of 2-tuples
def to_pairs(lst):
    if len(lst) % 2 != 0:
        raise ValueError("Sidechain values must contain an even number of entries.")
    return [(lst[i], lst[i+1]) for i in range(0, len(lst), 2)]