from flask import Blueprint, render_template, request, jsonify, send_from_directory
from werkzeug.utils import secure_filename
import os
from scripts.itp_parser import Itp_parser
from scripts.itpX2 import POLY_Extend
import json

#Blueprint
Extend_bp = Blueprint('Extend', __name__, url_prefix='/Extend')

UPLOAD_FOLDER = os.path.join(os.getcwd(), 'upload')
ALLOWED_EXTENSIONS = {'itp', 'gro'}



os.makedirs(UPLOAD_FOLDER, exist_ok=True)  # Create folder if it doesn't exist


# Define a route for this page
@Extend_bp.route('/')
def new_page():
    return render_template('Extend.html')



# Helper function to check allowed file extensions
def allowed_file(filename):
    """Check if the file has an allowed extension."""
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@Extend_bp.route('/upload', methods=['POST'])
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


@Extend_bp.route('/upload/<filename>', methods=['GET'])
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

@Extend_bp.route('/check_extend', methods=['GET'])
def check_extend():
    # Perform the necessary logic or return information about the extension status
    # For example, if you're checking whether some task has completed:
    result = {'status': 'success', 'message': 'Extend function is ready to proceed'}
    return jsonify(result), 200

@Extend_bp.route('/extend_action', methods=['POST'])
def extend_action():
    try:
        # Get data from the request (the data sent from JavaScript)
        data = request.get_json()
        print(F"data received: {data}")
        # Access the variables from the data dictionary
        #NOTE this was made before switching dihedrals to bonds, things are mislabled
        dihedral1 = list(data.get('dihedral1'))  #really bonds
        dihedral2 =  list(data.get('dihedral2')) #Really Bonds
        propagation  = list(data.get('propagation')) #Really starting and ending hydorgen
        Mon_num  = int(data.get('Mon_num'))

        print(dihedral1,dihedral2,propagation,Mon_num)


        # Trigger your function or any logic you want
        Extend(Bridge_1=dihedral1, Bridge_2=dihedral2, H_Start=propagation[0], H_End=propagation[1], Repeat=Mon_num)

        return jsonify({
            'message': 'successfully extended!',
            'files': ['EXTENDED_POLYMER.itp', 'EXTENDED_POLYMER.gro']
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# Define your backend function here
def Extend(Bridge_1, Bridge_2, H_Start, H_End, Repeat):
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
        print(itp_file[0],gro_file[0], Repeat,H_Start,H_End, Bridge_1,Bridge_2,"EXTENDED_POLYMER")
        POLY_Extend(ITP_FILE=itp_file[0],GRO_FILE=gro_file[0], Repeat = Repeat, H_Start=H_Start,H_End=H_End, Bridge_1=Bridge_1,Bridge_2=Bridge_2,Out_Name=os.path.join(UPLOAD_FOLDER, "EXTENDED_POLYMER"))



    except Exception as e:
        print(e)
        return f'Error loading ITP and GRO to JSON: {str(e)}'

