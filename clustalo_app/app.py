from flask import Flask, render_template, request, send_file
import requests
import subprocess
import os
import uuid
import re

app = Flask(__name__)

# Ensure the temp directory exists
TEMP_DIR = os.path.join(os.path.dirname(os.path.abspath(__name__)), 'temp')
os.makedirs(TEMP_DIR, exist_ok=True)

def detect_and_fetch_input(user_text):

    user_text = user_text.strip()
    
    # 1. Is it FASTA? (Starts with >)
    if user_text.startswith(">"):
        return user_text, None
        
    # 2. It is a list of IDs. Split strictly by newlines.
    lines = [line.strip() for line in user_text.splitlines() if line.strip()]
    fasta_result = ""

    for line in lines:
        # Enforce the "One per line / No spaces" rule
        if ' ' in line or ',' in line:
            return None, f"Format Error: Multiple items detected on one line ('{line}'). Please enter exactly one ID per line with no spaces or commas."
            
        query = line.upper()
        
        # Check PDB ID 
        if len(query) == 4 and query.isalnum():
            response = requests.get(f"https://www.rcsb.org/fasta/entry/{query}")
            if response.status_code == 200 and ">" in response.text:
                fasta_result += response.text + "\n"
            else:
                return None, f"PDB Error: Non-existent or invalid PDB ID '{query}'."
                
        # Check UniProt ID 
        elif len(query) >= 6 and query.isalnum():
            response = requests.get(f"https://www.uniprot.org/uniprot/{query}.fasta")
            if response.status_code == 200 and ">" in response.text:
                fasta_result += response.text + "\n"
            else:
                return None, f"UniProt Error: Non-existent or invalid UniProt ID '{query}'."
                
        # If it matches neither length/format
        else:
            return None, f"Unrecognized ID Error: '{query}' does not match standard PDB (4 chars) or UniProt (6+ chars) formats."
                
    return fasta_result, None

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/align', methods=['POST'])
def align():
    fasta_data = ""
    error_msg = None
    
    job_title = request.form.get('job_title', '').strip()
    if not job_title:
        job_title = "Untitled Alignment"
        
    # 1. Parse File Upload
    if 'file_upload' in request.files and request.files['file_upload'].filename != '':
        fasta_data = request.files['file_upload'].read().decode('utf-8')
        if not fasta_data.strip().startswith(">"):
            error_msg = "File Error: The uploaded file does not appear to be a valid FASTA format."
            
    # 2. Parse Text Input using our strict function
    elif 'sequence_input' in request.form and request.form['sequence_input'].strip() != '':
        raw_text = request.form['sequence_input']
        fasta_data, fetch_error = detect_and_fetch_input(raw_text)
        if fetch_error:
            error_msg = fetch_error
    else:
        error_msg = "Input Error: Please paste sequences/IDs or upload a file."

    # If any specific error was caught above, stop and show it
    if error_msg:
        return render_template('index.html', error=error_msg)
    
    # Final check: Does it have at least 2 sequences to align?
    if not fasta_data or fasta_data.count(">") < 2:
         return render_template('index.html', error="Execution Error: Clustal Omega requires at least two valid sequences to perform an alignment.")

    # Fasta specific errors check
    has_sequence = False
    for line in fasta_data.strip().splitlines():
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        clean_line = line.replace(" ", "")
        if not re.match(r'^[A-Za-z\-\*]+$', clean_line):
            return render_template('index.html', error=f"FASTA Format Error: Invalid character detected in sequence line: '{line[:30]}...'. Sequences must contain only letters.")
        has_sequence = True
        
    if not has_sequence:
        return render_template('index.html', error="FASTA Format Error: Sequence headers ('>') were found, but the actual sequence data is missing.")
    # ------------------------------------------

    # 3. Prepare files
    job_id = str(uuid.uuid4())
    input_file = os.path.join(TEMP_DIR, f"{job_id}_in.fasta")
    output_file = os.path.join(TEMP_DIR, f"{job_id}_out.aln")
    
    with open(input_file, "w") as f:
        f.write(fasta_data)
        
    outfmt = request.form.get('outfmt', 'clustal')
    
    # Advanced Options
    iterations = request.form.get('iterations', '0')
    dealign = request.form.get('dealign')
    full_matrix = request.form.get('full_matrix')
    
    # 4. Execute ClustalO
    cmd = ["clustalo", "-i", input_file, "-o", output_file, "--outfmt", outfmt, "--force"]
    
    if iterations != '0':
        cmd.extend(["--iter", iterations])
    if dealign == 'yes':
        cmd.append("--dealign")
    if full_matrix == 'yes':
        cmd.append("--full")

    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        return render_template('index.html', error=f"ClustalO System Error: {e.stderr}")

    # 5. Read Results
    with open(output_file, "r") as f:
        alignment_result = f.read()

    return render_template('result.html', alignment=alignment_result, job_id=job_id, job_title=job_title, outfmt=outfmt)

@app.route('/download/<job_id>')
def download(job_id):
    output_file = os.path.join(TEMP_DIR, f"{job_id}_out.aln")
    title = request.args.get('title', 'alignment').replace(" ", "_")
    
    if os.path.exists(output_file):
        return send_file(output_file, as_attachment=True, download_name=f"{title}.aln")
    return "File not found.", 404

if __name__ == '__main__':
    app.run(debug=True)