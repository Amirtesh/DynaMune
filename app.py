from flask import Flask, render_template, request, jsonify, send_file, session
import os
import subprocess
import pandas as pd
from pathlib import Path
import shutil
import uuid
import zipfile
import io
from datetime import datetime
import json
from Bio.SeqUtils.ProtParam import ProteinAnalysis

app = Flask(__name__)
app.secret_key = 'vaccine_builder_secret_key_2025'
app.config['SESSION_COOKIE_SAMESITE'] = 'Lax'
app.config['SESSION_COOKIE_HTTPONLY'] = True

# Store active session ID temporarily
current_prediction_session = None

# Base directories
BASE_DIR = Path(__file__).parent
NETBCE_DIR = BASE_DIR / 'NetBCE'
NETBCE_PREDICTION_DIR = NETBCE_DIR / 'prediction'
NETBCE_RESULT_DIR = NETBCE_DIR / 'result'
RESULTS_DIR = BASE_DIR / 'results'

# Adjuvant sequences dictionary
ADJUVANTS = {
    'β-defensin': 'GIINTLQKYYCRVRGGRCAVLSCLPKEEQIGKCSTRGRKCCRRKK',
    'PADRE': 'AKFVAAWTLKAAA',
    'CTB': 'MTPQNITDLCAEYHNTQIHTLTKSEKQLYVENGKTLVQRIKDFLRNQLSVEINCTFPNKQAINMKDTLRIAYLTEAKVEKLCVWNNKTPHAIAAISMAN',
    '50S_Ribosomal': 'MAKETAAAKFERQHMDSHLKSGEALKKELQAAIDYINGKEGAEVIHLVNNKDVAKLSEKAKDAQAEALQKKLEQALQQKVESLKEAFNASLDTAAAIATVKKDFEGVAMTTVEESLKVAAE',
    'RS09': 'APPHALS',
    'HBHA': 'MTEQQWNFAGIEAAASAIQGNVTSIHSLLDEGKQSLTKLAAAWGGSGSEAYQGVQQKWDATATELNNALQNLARTISEAGQAMASTEGNVTGMFA',
    'Flagellin': 'STNVKSSQDSKSDDLSNTLTQAAIDNKMLST',
    'Hsp70': 'SGVPLPPGMPGGAPPPGPAA',
    'HABA': 'MQMKKLLEAALAAAAP'
}

# Ensure directories exist
RESULTS_DIR.mkdir(exist_ok=True)
NETBCE_PREDICTION_DIR.mkdir(exist_ok=True)
NETBCE_RESULT_DIR.mkdir(exist_ok=True)


def get_session_dir():
    """Get or create a unique session directory for this user's job"""
    if 'session_id' not in session:
        # Generate unique session ID
        session['session_id'] = f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
    
    session_dir = RESULTS_DIR / session['session_id']
    session_dir.mkdir(parents=True, exist_ok=True)
    return session_dir


def cleanup_old_files():
    """Clean up old result files before starting a new prediction"""
    try:
        # Clean up session-specific directory if starting fresh
        if 'session_id' in session:
            session_dir = RESULTS_DIR / session['session_id']
            if session_dir.exists():
                shutil.rmtree(session_dir)
        
        # DON'T remove session_id - keep it for status tracking
        # session.pop('session_id', None)  # COMMENTED OUT
        
        # Also clean up old temp files in base directory
        files_to_remove = [
            RESULTS_DIR / 'B_Cell_Epitopes.csv',
            RESULTS_DIR / 'mhc1_prediction.csv',
            RESULTS_DIR / 'mhc2_prediction.csv',
            BASE_DIR / 'mhc1_prediction.csv',
            BASE_DIR / 'mhc2_prediction.csv',
            BASE_DIR / 'test.fasta',
            BASE_DIR / 'alleles_mhc1.txt',
            BASE_DIR / 'alleles_mhc2.txt',
            NETBCE_PREDICTION_DIR / 'test.fasta',
            NETBCE_RESULT_DIR / 'NetBCE_predictions.tsv'
        ]
        
        for file_path in files_to_remove:
            if file_path.exists():
                file_path.unlink()
    except Exception as e:
        print(f"Cleanup warning: {e}")


@app.route('/')
def index():
    """Main page for input"""
    return render_template('index.html', adjuvants=ADJUVANTS.keys())


@app.route('/session_info')
def session_info():
    """Get current session information - creates session if doesn't exist"""
    # This is called before prediction starts, so we create session here
    session_dir = get_session_dir()  # This creates session_id if it doesn't exist
    session.modified = True
    
    return jsonify({
        'session_id': session['session_id'],
        'session_path': str(session_dir.relative_to(RESULTS_DIR))
    })


@app.route('/check_prediction_status')
def check_prediction_status():
    """Check which prediction files exist"""
    try:
        # First try global variable (for active prediction)
        global current_prediction_session
        session_id = current_prediction_session
        
        # Fallback to Flask session
        if not session_id and 'session_id' in session:
            session_id = session['session_id']
        
        if not session_id:
            return jsonify({
                'bcell_done': False,
                'mhc1_done': False,
                'mhc2_done': False
            })
        
        session_dir = RESULTS_DIR / session_id
        
        if not session_dir.exists():
            return jsonify({
                'bcell_done': False,
                'mhc1_done': False,
                'mhc2_done': False
            })
        
        bcell_file = session_dir / 'B_Cell_Epitopes.csv'
        mhc1_file = session_dir / 'mhc1_prediction.csv'
        mhc2_file = session_dir / 'mhc2_prediction.csv'
        
        bcell_exists = bcell_file.exists()
        bcell_size = bcell_file.stat().st_size if bcell_exists else 0
        mhc1_exists = mhc1_file.exists()
        mhc1_size = mhc1_file.stat().st_size if mhc1_exists else 0
        mhc2_exists = mhc2_file.exists()
        mhc2_size = mhc2_file.stat().st_size if mhc2_exists else 0
        
        status = {
            'bcell_done': bcell_exists and bcell_size > 100,
            'mhc1_done': mhc1_exists and mhc1_size > 100,
            'mhc2_done': mhc2_exists and mhc2_size > 100
        }
        
        print(f"Status: bcell={status['bcell_done']}({bcell_size}), mhc1={status['mhc1_done']}({mhc1_size}), mhc2={status['mhc2_done']}({mhc2_size})")
        
        return jsonify(status)
    except Exception as e:
        print(f"Status check error: {e}")
        return jsonify({
            'bcell_done': False,
            'mhc1_done': False,
            'mhc2_done': False
        })
    except Exception as e:
        print(f"Error checking status: {e}")
        import traceback
        traceback.print_exc()
        return jsonify({
            'bcell_done': False,
            'mhc1_done': False,
            'mhc2_done': False
        })


@app.route('/predict', methods=['POST'])
def predict():
    """Run all predictions pipeline"""
    try:
        # Cleanup old files
        cleanup_old_files()
        
        # Get input data
        fasta_content = request.form.get('fasta_sequence', '').strip()
        mhc1_alleles = request.form.get('mhc1_alleles', '').strip()
        mhc2_alleles = request.form.get('mhc2_alleles', '').strip()
        
        # Validate inputs
        if not fasta_content:
            return jsonify({'success': False, 'error': 'Please provide a protein sequence in FASTA format'})
        
        if not mhc1_alleles:
            return jsonify({'success': False, 'error': 'Please provide MHC-I alleles'})
        
        if not mhc2_alleles:
            return jsonify({'success': False, 'error': 'Please provide MHC-II alleles'})
        
        # Get session directory for this job
        session_dir = get_session_dir()
        session_id = session['session_id']
        
        # Store session ID globally for status checks
        global current_prediction_session
        current_prediction_session = session_id
        
        session.permanent = True
        session.modified = True
        print(f"Using session directory: {session_dir}")
        print(f"Session ID: {session_id}")
        
        # Step 1: B-cell Epitope Prediction (NetBCE)
        print("Step 1: Running NetBCE prediction...")
        netbce_fasta = NETBCE_PREDICTION_DIR / 'test.fasta'
        with open(netbce_fasta, 'w') as f:
            f.write(fasta_content)
        
        # Run NetBCE prediction
        try:
            result = subprocess.run(
                ['python3', 'NetBCE_prediction.py', '-f', 'test.fasta', '-o', '../result/'],
                cwd=NETBCE_PREDICTION_DIR,
                capture_output=True,
                text=True,
                timeout=2000
            )
            
            if result.returncode != 0:
                return jsonify({'success': False, 'error': f'NetBCE prediction failed: {result.stderr}'})
            
            # Read NetBCE predictions and save to results
            netbce_output = NETBCE_RESULT_DIR / 'NetBCE_predictions.tsv'
            if netbce_output.exists():
                # NetBCE output has no header, columns are: sequence_id, peptide, position, score
                df_bcell = pd.read_csv(netbce_output, sep='\t', header=None, 
                                      names=['Sequence_ID', 'Peptide', 'Position', 'Score'])
                
                # Sort by score (descending - higher is better)
                df_bcell = df_bcell.sort_values('Score', ascending=False)
                
                # Save to session directory as CSV
                df_bcell.to_csv(session_dir / 'B_Cell_Epitopes.csv', index=False)
                
                # Clean up NetBCE files
                if netbce_fasta.exists():
                    netbce_fasta.unlink()
                if netbce_output.exists():
                    netbce_output.unlink()
                # Also clean up the HTML report if it exists
                report_file = NETBCE_RESULT_DIR / 'report.html'
                if report_file.exists():
                    report_file.unlink()
            else:
                return jsonify({'success': False, 'error': 'NetBCE did not generate predictions'})
                
        except subprocess.TimeoutExpired:
            return jsonify({'success': False, 'error': 'NetBCE prediction timed out'})
        except Exception as e:
            return jsonify({'success': False, 'error': f'NetBCE prediction error: {str(e)}'})
        
        # Step 2: MHC-I Prediction
        print("Step 2: Running MHC-I prediction...")
        test_fasta = BASE_DIR / 'test.fasta'
        with open(test_fasta, 'w') as f:
            f.write(fasta_content)
        
        mhc1_alleles_file = BASE_DIR / 'alleles_mhc1.txt'
        with open(mhc1_alleles_file, 'w') as f:
            f.write(mhc1_alleles)
        
        try:
            result = subprocess.run(
                ['python3', 'mhc1.py', 'test.fasta', 'alleles_mhc1.txt'],
                cwd=BASE_DIR,
                capture_output=True,
                text=True,
                timeout=2000
            )
            
            if result.returncode != 0:
                return jsonify({'success': False, 'error': f'MHC-I prediction failed: {result.stderr}'})
            
            # Move and sort MHC-I predictions
            mhc1_output = BASE_DIR / 'mhc1_prediction.csv'
            if mhc1_output.exists():
                df_mhc1 = pd.read_csv(mhc1_output)
                # Sort by netmhcpan_el_score (descending)
                if 'netmhcpan_el_score' in df_mhc1.columns:
                    df_mhc1 = df_mhc1.sort_values('netmhcpan_el_score', ascending=False)
                
                df_mhc1.to_csv(session_dir / 'mhc1_prediction.csv', index=False)
            else:
                return jsonify({'success': False, 'error': 'MHC-I prediction did not generate output'})
                
        except subprocess.TimeoutExpired:
            return jsonify({'success': False, 'error': 'MHC-I prediction timed out'})
        except Exception as e:
            return jsonify({'success': False, 'error': f'MHC-I prediction error: {str(e)}'})
        
        # Step 3: MHC-II Prediction
        print("Step 3: Running MHC-II prediction...")
        mhc2_alleles_file = BASE_DIR / 'alleles_mhc2.txt'
        with open(mhc2_alleles_file, 'w') as f:
            f.write(mhc2_alleles)
        
        try:
            result = subprocess.run(
                ['python3', 'mhc2.py', 'test.fasta', 'alleles_mhc2.txt'],
                cwd=BASE_DIR,
                capture_output=True,
                text=True,
                timeout=2000
            )
            
            if result.returncode != 0:
                return jsonify({'success': False, 'error': f'MHC-II prediction failed: {result.stderr}'})
            
            # Move and sort MHC-II predictions
            mhc2_output = BASE_DIR / 'mhc2_prediction.csv'
            if mhc2_output.exists():
                df_mhc2 = pd.read_csv(mhc2_output)
                # Sort by netmhciipan_el_score (descending)
                if 'netmhciipan_el_score' in df_mhc2.columns:
                    df_mhc2 = df_mhc2.sort_values('netmhciipan_el_score', ascending=False)
                
                df_mhc2.to_csv(session_dir / 'mhc2_prediction.csv', index=False)
            else:
                return jsonify({'success': False, 'error': 'MHC-II prediction did not generate output'})
                
        except subprocess.TimeoutExpired:
            return jsonify({'success': False, 'error': 'MHC-II prediction timed out'})
        except Exception as e:
            return jsonify({'success': False, 'error': f'MHC-II prediction error: {str(e)}'})
        
        return jsonify({'success': True, 'message': 'All predictions completed successfully!'})
        
    except Exception as e:
        return jsonify({'success': False, 'error': f'Unexpected error: {str(e)}'})


@app.route('/save_selected_epitopes', methods=['POST'])
def save_selected_epitopes():
    """Save selected epitopes to session for screening"""
    try:
        data = request.get_json()
        
        # Save to session
        session['selected_bcell'] = data.get('bcell', [])
        session['selected_mhc1'] = data.get('mhc1', [])
        session['selected_mhc2'] = data.get('mhc2', [])
        
        return jsonify({'success': True})
        
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@app.route('/results')
def results():
    """Display prediction results"""
    try:
        # Get session directory
        session_dir = get_session_dir()
        
        # Read all prediction files
        bcell_file = session_dir / 'B_Cell_Epitopes.csv'
        mhc1_file = session_dir / 'mhc1_prediction.csv'
        mhc2_file = session_dir / 'mhc2_prediction.csv'
        
        # Check if files exist
        if not bcell_file.exists() or not mhc1_file.exists() or not mhc2_file.exists():
            return render_template('results.html', error='Prediction files not found. Please run predictions first.')
        
        # Read data
        df_bcell = pd.read_csv(bcell_file).head(50)
        df_mhc1 = pd.read_csv(mhc1_file).head(50)
        df_mhc2 = pd.read_csv(mhc2_file).head(50)
        
        # Convert to JSON for frontend
        bcell_data = df_bcell.to_dict('records')
        mhc1_data = df_mhc1.to_dict('records')
        mhc2_data = df_mhc2.to_dict('records')
        
        return render_template('results.html', 
                             bcell_data=bcell_data,
                             mhc1_data=mhc1_data,
                             mhc2_data=mhc2_data,
                             bcell_columns=df_bcell.columns.tolist(),
                             mhc1_columns=df_mhc1.columns.tolist(),
                             mhc2_columns=df_mhc2.columns.tolist())
        
    except Exception as e:
        return render_template('results.html', error=f'Error loading results: {str(e)}')


@app.route('/download/<filename>')
def download_file(filename):
    """Download prediction files"""
    try:
        session_dir = get_session_dir()
        file_path = session_dir / filename
        if file_path.exists():
            return send_file(file_path, as_attachment=True, download_name=filename)
        else:
            return jsonify({'error': 'File not found'}), 404
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/screening')
def screening():
    """Render screening page for allergenicity and toxicity prediction"""
    return render_template('screening.html')


@app.route('/run_screening', methods=['POST'])
def run_screening():
    """Run AlgPred2.0 and ToxinPred3 on selected epitopes"""
    try:
        session_dir = get_session_dir()
        
        # Get epitopes from session
        bcell_epitopes = session.get('selected_bcell', [])
        mhc1_epitopes = session.get('selected_mhc1', [])
        mhc2_epitopes = session.get('selected_mhc2', [])
        
        if not bcell_epitopes and not mhc1_epitopes and not mhc2_epitopes:
            return jsonify({'success': False, 'error': 'No epitopes found in session'})
        
        results = {
            'success': True,
            'bcell': [],
            'mhc1': [],
            'mhc2': []
        }
        
        # Process B-cell epitopes
        if bcell_epitopes:
            bcell_file = session_dir / 'bcell.seq'
            with open(bcell_file, 'w') as f:
                for epitope in bcell_epitopes:
                    f.write(f"{epitope}\n")
            
            # Run AlgPred2.0
            algpred_output = session_dir / 'bcell_algpred.csv'
            result_alg = subprocess.run(
                ['algpred2', '-i', str(bcell_file), '-d', '2', '-m', '2', '-o', str(algpred_output)],
                cwd=str(session_dir),
                capture_output=True,
                text=True,
                timeout=300
            )
            print(f"AlgPred2 B-cell - Return code: {result_alg.returncode}")
            print(f"AlgPred2 B-cell - stdout: {result_alg.stdout}")
            print(f"AlgPred2 B-cell - stderr: {result_alg.stderr}")
            
            # Run ToxinPred3
            toxinpred_output = session_dir / 'bcell_toxinpred.csv'
            result_tox = subprocess.run(
                ['toxinpred3', '-i', str(bcell_file), '-d', '2', '-m', '2', '-o', str(toxinpred_output)],
                cwd=str(session_dir),
                capture_output=True,
                text=True,
                timeout=300
            )
            print(f"ToxinPred3 B-cell - Return code: {result_tox.returncode}")
            print(f"ToxinPred3 B-cell - stdout: {result_tox.stdout}")
            print(f"ToxinPred3 B-cell - stderr: {result_tox.stderr}")
            
            # Check if output files exist
            print(f"AlgPred output exists: {algpred_output.exists()}")
            print(f"ToxinPred output exists: {toxinpred_output.exists()}")
            
            # Parse results - pass epitope list for index mapping
            algpred_results = parse_screening_results(algpred_output, bcell_epitopes)
            toxinpred_results = parse_screening_results(toxinpred_output, bcell_epitopes)
            
            print(f"AlgPred results: {algpred_results}")
            print(f"ToxinPred results: {toxinpred_results}")
            
            for i, epitope in enumerate(bcell_epitopes):
                results['bcell'].append({
                    'sequence': epitope,
                    'algpred': algpred_results.get(epitope, 'N/A'),
                    'toxinpred': toxinpred_results.get(epitope, 'N/A')
                })
        
        # Process MHC-I epitopes
        if mhc1_epitopes:
            mhc1_file = session_dir / 'mhc1.seq'
            with open(mhc1_file, 'w') as f:
                for epitope in mhc1_epitopes:
                    f.write(f"{epitope}\n")
            
            # Run AlgPred2.0
            algpred_output = session_dir / 'mhc1_algpred.csv'
            result_alg = subprocess.run(
                ['algpred2', '-i', str(mhc1_file), '-d', '2', '-m', '2', '-o', str(algpred_output)],
                cwd=str(session_dir),
                capture_output=True,
                text=True,
                timeout=300
            )
            
            # Run ToxinPred3
            toxinpred_output = session_dir / 'mhc1_toxinpred.csv'
            result_tox = subprocess.run(
                ['toxinpred3', '-i', str(mhc1_file), '-d', '2', '-m', '2', '-o', str(toxinpred_output)],
                cwd=str(session_dir),
                capture_output=True,
                text=True,
                timeout=300
            )
            
            # Parse results - pass epitope list for index mapping
            algpred_results = parse_screening_results(algpred_output, mhc1_epitopes)
            toxinpred_results = parse_screening_results(toxinpred_output, mhc1_epitopes)
            
            for i, epitope in enumerate(mhc1_epitopes):
                results['mhc1'].append({
                    'sequence': epitope,
                    'algpred': algpred_results.get(epitope, 'N/A'),
                    'toxinpred': toxinpred_results.get(epitope, 'N/A')
                })
        
        # Process MHC-II epitopes
        if mhc2_epitopes:
            mhc2_file = session_dir / 'mhc2.seq'
            with open(mhc2_file, 'w') as f:
                for epitope in mhc2_epitopes:
                    f.write(f"{epitope}\n")
            
            # Run AlgPred2.0
            algpred_output = session_dir / 'mhc2_algpred.csv'
            result_alg = subprocess.run(
                ['algpred2', '-i', str(mhc2_file), '-d', '2', '-m', '2', '-o', str(algpred_output)],
                cwd=str(session_dir),
                capture_output=True,
                text=True,
                timeout=300
            )
            
            # Run ToxinPred3
            toxinpred_output = session_dir / 'mhc2_toxinpred.csv'
            result_tox = subprocess.run(
                ['toxinpred3', '-i', str(mhc2_file), '-d', '2', '-m', '2', '-o', str(toxinpred_output)],
                cwd=str(session_dir),
                capture_output=True,
                text=True,
                timeout=300
            )
            
            # Parse results - pass epitope list for index mapping
            algpred_results = parse_screening_results(algpred_output, mhc2_epitopes)
            toxinpred_results = parse_screening_results(toxinpred_output, mhc2_epitopes)
            
            for i, epitope in enumerate(mhc2_epitopes):
                results['mhc2'].append({
                    'sequence': epitope,
                    'algpred': algpred_results.get(epitope, 'N/A'),
                    'toxinpred': toxinpred_results.get(epitope, 'N/A')
                })
        
        return jsonify(results)
        
    except subprocess.TimeoutExpired:
        return jsonify({'success': False, 'error': 'Screening timeout. Please try with fewer epitopes.'})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


def parse_screening_results(csv_file, epitope_list):
    """Parse AlgPred2.0 or ToxinPred3 CSV output and map to actual epitope sequences"""
    results = {}
    try:
        if not csv_file.exists():
            print(f"CSV file does not exist: {csv_file}")
            return results
        
        import csv
        with open(csv_file, 'r') as f:
            reader = csv.DictReader(f)
            for i, row in enumerate(reader):
                # Get prediction from the CSV
                prediction = row.get('Prediction', row.get('prediction', row.get('result', row.get('Result', 'N/A'))))
                
                # Map by index - CSV has Seq_1, Seq_2, etc. which correspond to index 0, 1, 2
                if i < len(epitope_list) and prediction:
                    actual_epitope = epitope_list[i]
                    results[actual_epitope] = prediction.strip()
                    
    except Exception as e:
        print(f"Error parsing {csv_file}: {e}")
        import traceback
        traceback.print_exc()
    
    return results


@app.route('/save_screened_epitopes', methods=['POST'])
def save_screened_epitopes():
    """Save screened epitopes to session for construct building"""
    try:
        data = request.get_json()
        
        # Update session with screened epitopes
        session['selected_bcell'] = data.get('bcell', [])
        session['selected_mhc1'] = data.get('mhc1', [])
        session['selected_mhc2'] = data.get('mhc2', [])
        
        return jsonify({'success': True})
        
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@app.route('/construct', methods=['GET', 'POST'])
def construct():
    """Build vaccine construct"""
    if request.method == 'GET':
        return render_template('construct.html', adjuvants=ADJUVANTS)
    
    try:
        # Get selected epitopes, adjuvant, and design
        data = request.get_json()
        bcell_epitopes = data.get('bcell_epitopes', [])
        mhc1_epitopes = data.get('mhc1_epitopes', [])  # CTL epitopes
        mhc2_epitopes = data.get('mhc2_epitopes', [])  # HTL epitopes
        adjuvant_name = data.get('adjuvant', '')
        design = data.get('design', 'design1')
        
        if not adjuvant_name or adjuvant_name not in ADJUVANTS:
            return jsonify({'success': False, 'error': 'Please select a valid adjuvant'})
        
        if not bcell_epitopes and not mhc1_epitopes and not mhc2_epitopes:
            return jsonify({'success': False, 'error': 'Please select at least one epitope'})
        
        # Build construct based on selected design
        construct_parts = []
        
        # Add adjuvant
        construct_parts.append(ADJUVANTS[adjuvant_name])
        construct_parts.append('EAAAK')
        
        # Design 1: CTL → HTL → B-cell (Sequential)
        if design == 'design1':
            if mhc1_epitopes:
                for i, epitope in enumerate(mhc1_epitopes):
                    construct_parts.append(epitope)
                    if i < len(mhc1_epitopes) - 1:
                        construct_parts.append('AAY')
            if mhc2_epitopes:
                if mhc1_epitopes:
                    construct_parts.append('GPGPG')
                for i, epitope in enumerate(mhc2_epitopes):
                    construct_parts.append(epitope)
                    if i < len(mhc2_epitopes) - 1:
                        construct_parts.append('GPGPG')
            if bcell_epitopes:
                if mhc1_epitopes or mhc2_epitopes:
                    construct_parts.append('KK')
                for i, epitope in enumerate(bcell_epitopes):
                    construct_parts.append(epitope)
                    if i < len(bcell_epitopes) - 1:
                        construct_parts.append('KK')
        
        # Design 2: HTL → CTL → B-cell (Sequential)
        elif design == 'design2':
            if mhc2_epitopes:
                for i, epitope in enumerate(mhc2_epitopes):
                    construct_parts.append(epitope)
                    if i < len(mhc2_epitopes) - 1:
                        construct_parts.append('GPGPG')
            if mhc1_epitopes:
                if mhc2_epitopes:
                    construct_parts.append('AAY')
                for i, epitope in enumerate(mhc1_epitopes):
                    construct_parts.append(epitope)
                    if i < len(mhc1_epitopes) - 1:
                        construct_parts.append('AAY')
            if bcell_epitopes:
                if mhc1_epitopes or mhc2_epitopes:
                    construct_parts.append('KK')
                for i, epitope in enumerate(bcell_epitopes):
                    construct_parts.append(epitope)
                    if i < len(bcell_epitopes) - 1:
                        construct_parts.append('KK')
        
        # Design 3: B-cell → HTL → CTL (Sequential)
        elif design == 'design3':
            if bcell_epitopes:
                for i, epitope in enumerate(bcell_epitopes):
                    construct_parts.append(epitope)
                    if i < len(bcell_epitopes) - 1:
                        construct_parts.append('KK')
            if mhc2_epitopes:
                if bcell_epitopes:
                    construct_parts.append('AAY')
                for i, epitope in enumerate(mhc2_epitopes):
                    construct_parts.append(epitope)
                    if i < len(mhc2_epitopes) - 1:
                        construct_parts.append('GPGPG')
            if mhc1_epitopes:
                if bcell_epitopes or mhc2_epitopes:
                    construct_parts.append('AAY')
                for i, epitope in enumerate(mhc1_epitopes):
                    construct_parts.append(epitope)
                    if i < len(mhc1_epitopes) - 1:
                        construct_parts.append('AAY')
        
        # Design 4: CTL-HTL-B alternating pattern
        elif design == 'design4':
            max_len = max(len(mhc1_epitopes) if mhc1_epitopes else 0,
                         len(mhc2_epitopes) if mhc2_epitopes else 0,
                         len(bcell_epitopes) if bcell_epitopes else 0)
            for i in range(max_len):
                if mhc1_epitopes and i < len(mhc1_epitopes):
                    construct_parts.append(mhc1_epitopes[i])
                    if (mhc2_epitopes and i < len(mhc2_epitopes)) or (bcell_epitopes and i < len(bcell_epitopes)):
                        construct_parts.append('AAY')
                if mhc2_epitopes and i < len(mhc2_epitopes):
                    construct_parts.append(mhc2_epitopes[i])
                    if bcell_epitopes and i < len(bcell_epitopes):
                        construct_parts.append('GPGPG')
                if bcell_epitopes and i < len(bcell_epitopes):
                    construct_parts.append(bcell_epitopes[i])
                    if i < max_len - 1:
                        construct_parts.append('KK')
        
        # Design 5: HTL-B-CTL alternating pattern
        elif design == 'design5':
            max_len = max(len(mhc1_epitopes) if mhc1_epitopes else 0,
                         len(mhc2_epitopes) if mhc2_epitopes else 0,
                         len(bcell_epitopes) if bcell_epitopes else 0)
            for i in range(max_len):
                if mhc2_epitopes and i < len(mhc2_epitopes):
                    construct_parts.append(mhc2_epitopes[i])
                    if (bcell_epitopes and i < len(bcell_epitopes)) or (mhc1_epitopes and i < len(mhc1_epitopes)):
                        construct_parts.append('GPGPG')
                if bcell_epitopes and i < len(bcell_epitopes):
                    construct_parts.append(bcell_epitopes[i])
                    if mhc1_epitopes and i < len(mhc1_epitopes):
                        construct_parts.append('KK')
                if mhc1_epitopes and i < len(mhc1_epitopes):
                    construct_parts.append(mhc1_epitopes[i])
                    if i < max_len - 1:
                        construct_parts.append('AAY')
        
        # Design 6: CTL-B-HTL alternating pattern
        elif design == 'design6':
            max_len = max(len(mhc1_epitopes) if mhc1_epitopes else 0,
                         len(mhc2_epitopes) if mhc2_epitopes else 0,
                         len(bcell_epitopes) if bcell_epitopes else 0)
            for i in range(max_len):
                if mhc1_epitopes and i < len(mhc1_epitopes):
                    construct_parts.append(mhc1_epitopes[i])
                    if (bcell_epitopes and i < len(bcell_epitopes)) or (mhc2_epitopes and i < len(mhc2_epitopes)):
                        construct_parts.append('AAY')
                if bcell_epitopes and i < len(bcell_epitopes):
                    construct_parts.append(bcell_epitopes[i])
                    if mhc2_epitopes and i < len(mhc2_epitopes):
                        construct_parts.append('KK')
                if mhc2_epitopes and i < len(mhc2_epitopes):
                    construct_parts.append(mhc2_epitopes[i])
                    if i < max_len - 1:
                        construct_parts.append('GPGPG')
        
        # Design 7: B-CTL-HTL alternating pattern
        elif design == 'design7':
            max_len = max(len(mhc1_epitopes) if mhc1_epitopes else 0,
                         len(mhc2_epitopes) if mhc2_epitopes else 0,
                         len(bcell_epitopes) if bcell_epitopes else 0)
            for i in range(max_len):
                if bcell_epitopes and i < len(bcell_epitopes):
                    construct_parts.append(bcell_epitopes[i])
                    if (mhc1_epitopes and i < len(mhc1_epitopes)) or (mhc2_epitopes and i < len(mhc2_epitopes)):
                        construct_parts.append('KK')
                if mhc1_epitopes and i < len(mhc1_epitopes):
                    construct_parts.append(mhc1_epitopes[i])
                    if mhc2_epitopes and i < len(mhc2_epitopes):
                        construct_parts.append('AAY')
                if mhc2_epitopes and i < len(mhc2_epitopes):
                    construct_parts.append(mhc2_epitopes[i])
                    if i < max_len - 1:
                        construct_parts.append('GPGPG')
        
        # Join all parts with hyphens for display
        final_construct = '-'.join(construct_parts)
        
        # Also create clean sequence (no hyphens) for FASTA
        clean_sequence = ''.join(construct_parts)
        
        # Save construct to session
        session['vaccine_construct'] = clean_sequence
        session['adjuvant_used'] = adjuvant_name
        
        # Save FASTA file immediately to results directory
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        fasta_content = f">Vaccine_Construct_Adjuvant_{adjuvant_name}_{timestamp}\n"
        
        # Split sequence into 80-character lines (FASTA standard)
        for i in range(0, len(clean_sequence), 80):
            fasta_content += clean_sequence[i:i+80] + '\n'
        
        # Save to session directory
        session_dir = get_session_dir()
        output_file = session_dir / 'vaccine_construct.fasta'
        with open(output_file, 'w') as f:
            f.write(fasta_content)
        
        return jsonify({
            'success': True,
            'construct': final_construct,
            'clean_sequence': clean_sequence,
            'length': len(clean_sequence)
        })
        
    except Exception as e:
        return jsonify({'success': False, 'error': f'Construction error: {str(e)}'})


@app.route('/download_construct')
def download_construct():
    """Download vaccine construct as FASTA"""
    try:
        vaccine_sequence = session.get('vaccine_construct', '')
        adjuvant_used = session.get('adjuvant_used', 'unknown')
        
        if not vaccine_sequence:
            return jsonify({'error': 'No vaccine construct found'}), 404
        
        # Create FASTA content
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        fasta_content = f">Vaccine_Construct_Adjuvant_{adjuvant_used}_{timestamp}\n"
        
        # Split sequence into 80-character lines (FASTA standard)
        for i in range(0, len(vaccine_sequence), 80):
            fasta_content += vaccine_sequence[i:i+80] + '\n'
        
        # Save to session directory
        session_dir = get_session_dir()
        output_file = session_dir / 'vaccine_construct.fasta'
        with open(output_file, 'w') as f:
            f.write(fasta_content)
        
        return send_file(output_file, as_attachment=True, download_name='vaccine_construct.fasta')
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/protparam', methods=['POST'])
def protparam():
    """Analyze vaccine construct using ProtParam"""
    try:
        data = request.get_json()
        sequence = data.get('sequence', '').replace('-', '')  # Remove any linker separators
        
        if not sequence:
            return jsonify({'success': False, 'error': 'No sequence provided'})
        
        # Initialize ProteinAnalysis
        protein_analysis = ProteinAnalysis(sequence)
        
        # Get amino acid counts
        aa_counts = protein_analysis.count_amino_acids()
        aa_percent = protein_analysis.get_amino_acids_percent()
        
        # Calculate properties
        molecular_weight = protein_analysis.molecular_weight()
        aromaticity = protein_analysis.aromaticity()
        instability_index = protein_analysis.instability_index()
        isoelectric_point = protein_analysis.isoelectric_point()
        gravy = protein_analysis.gravy()
        
        # Secondary structure fractions
        sec_struc = protein_analysis.secondary_structure_fraction()
        
        # Molar extinction coefficient
        epsilon = protein_analysis.molar_extinction_coefficient()
        
        # Prepare response
        result = {
            'success': True,
            'molecular_weight': round(molecular_weight, 2),
            'aromaticity': round(aromaticity, 4),
            'instability_index': round(instability_index, 2),
            'isoelectric_point': round(isoelectric_point, 2),
            'gravy': round(gravy, 4),
            'helix_fraction': round(sec_struc[0], 4),
            'turn_fraction': round(sec_struc[1], 4),
            'sheet_fraction': round(sec_struc[2], 4),
            'extinction_coefficient_reduced': epsilon[0],
            'extinction_coefficient_oxidized': epsilon[1],
            'amino_acid_composition': {aa: round(percent, 4) for aa, percent in aa_percent.items()},
            'amino_acid_counts': aa_counts,
            'sequence_length': len(sequence),
            'stability': 'Unstable' if instability_index > 40 else 'Stable'
        }
        
        return jsonify(result)
        
    except Exception as e:
        return jsonify({'success': False, 'error': f'ProtParam analysis error: {str(e)}'})


@app.route('/predict_structure', methods=['POST'])
def predict_structure():
    """Predict 3D structure using ESMFold"""
    try:
        data = request.get_json()
        sequence = data.get('sequence', '')
        
        if not sequence:
            return jsonify({'success': False, 'error': 'No sequence provided'})
        
        # Get session directory
        session_dir = get_session_dir()
        
        # Save sequence to temporary file
        sequence_file = session_dir / 'sequence.txt'
        with open(sequence_file, 'w') as f:
            f.write(sequence)
        
        # Run model.py script
        import subprocess
        result = subprocess.run(
            ['python3', 'model.py', str(sequence_file)],
            cwd=BASE_DIR,
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )
        
        # Check if PDB file was created
        pdb_file = BASE_DIR / 'predicted_structure.pdb'
        if pdb_file.exists():
            # Read PDB content
            with open(pdb_file, 'r') as f:
                pdb_content = f.read()
            
            # Save to session directory for later download
            output_pdb = session_dir / 'predicted_structure.pdb'
            shutil.copy(pdb_file, output_pdb)
            # Clean up temp file
            pdb_file.unlink()
            
            return jsonify({
                'success': True,
                'pdb_content': pdb_content,
                'message': 'Structure prediction completed successfully'
            })
        else:
            return jsonify({
                'success': False,
                'error': 'Structure prediction failed. No PDB file generated.',
                'output': result.stdout,
                'error_output': result.stderr
            })
    
    except subprocess.TimeoutExpired:
        return jsonify({'success': False, 'error': 'Structure prediction timed out (5 minutes)'})
    except Exception as e:
        return jsonify({'success': False, 'error': f'Structure prediction error: {str(e)}'})


@app.route('/download_pdb')
def download_pdb():
    """Download predicted structure as PDB file"""
    try:
        session_dir = get_session_dir()
        pdb_file = session_dir / 'predicted_structure.pdb'
        if pdb_file.exists():
            return send_file(pdb_file, as_attachment=True, download_name='vaccine_structure.pdb')
        else:
            return jsonify({'error': 'No PDB file found'}), 404
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/prody_nma')
def prody_nma_page():
    """Render ProDy NMA analysis page"""
    return render_template('prody_nma.html')


@app.route('/run_prody_nma', methods=['POST'])
def run_prody_nma():
    """Run ProDy Normal Mode Analysis"""
    try:
        data = request.get_json()
        method = data.get('method', 'ANM')
        cutoff = data.get('cutoff', 10.0)  # Default changed to 10.0
        n_modes = data.get('n_modes', 20)
        conformers = data.get('conformers', 10)
        interface_cutoff = data.get('interface_cutoff', 10.0)  # New parameter
        
        # Validate inputs
        valid_methods = ['ANM', 'GNM']
        if method not in valid_methods:
            return jsonify({'success': False, 'error': f'Invalid method. Choose from {valid_methods}'})
        
        # Validate cutoff
        try:
            cutoff = float(cutoff)
            if cutoff < 5.0 or cutoff > 30.0:
                return jsonify({'success': False, 'error': 'Cutoff distance must be between 5.0 and 30.0 Å'})
        except ValueError:
            return jsonify({'success': False, 'error': 'Invalid cutoff value'})
        
        # Validate interface_cutoff
        try:
            interface_cutoff = float(interface_cutoff)
            if interface_cutoff < 5.0 or interface_cutoff > 30.0:
                return jsonify({'success': False, 'error': 'Interface cutoff must be between 5.0 and 30.0 Å'})
        except ValueError:
            return jsonify({'success': False, 'error': 'Invalid interface cutoff value'})
        
        try:
            n_modes = int(n_modes)
            if n_modes < 1 or n_modes > 100:
                return jsonify({'success': False, 'error': 'Number of modes must be between 1 and 100'})
        except ValueError:
            return jsonify({'success': False, 'error': 'Invalid number of modes'})
        
        try:
            conformers = int(conformers)
            if conformers < 1 or conformers > 50:
                return jsonify({'success': False, 'error': 'Number of conformers must be between 1 and 50'})
        except ValueError:
            return jsonify({'success': False, 'error': 'Invalid number of conformers'})
        
        # Get session directory
        session_dir = get_session_dir()
        
        # Check if PDB file exists
        pdb_file = session_dir / 'predicted_structure.pdb'
        if not pdb_file.exists():
            return jsonify({'success': False, 'error': 'No PDB file found. Please predict structure first.'})
        
        # Build command with required flags
        cmd = ['python3', str(BASE_DIR / 'prody_nma.py'), str(pdb_file), 
               '--method', method, '--nmodes', str(n_modes), '--conformers', str(conformers),
               '--cutoff', str(cutoff), '--interface-cutoff', str(interface_cutoff),
               '--ensemble-pca', '--allosteric', '--clustenmd']  # Always include these 3 flags
        
        # Run ProDy NMA analysis with output in session directory
        # Change working directory to session dir to keep all output there
        import subprocess
        result = subprocess.run(
            cmd,
            cwd=str(session_dir),
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout
        )
        
        # The output directory is created by prody_nma.py as {pdb_basename}_NMA_PCA_results
        pdb_basename = pdb_file.stem  # 'predicted_structure'
        nma_outdir = session_dir / f'{pdb_basename}_NMA_PCA_results'
        
        # Check for output files (including new ones)
        comprehensive_plot = nma_outdir / 'comprehensive_analysis.png'
        mode_shapes_plot = nma_outdir / 'mode_shapes_combined.png'
        pca_comprehensive = nma_outdir / 'pca_comprehensive_analysis.png'
        anm_pca_comparison = nma_outdir / 'anm_pca_comparison.png'
        clustenmd_analysis = nma_outdir / 'clustenmd_analysis.png'
        
        # Debug: List all files in output directory
        debug_info = {
            'nma_outdir': str(nma_outdir),
            'nma_outdir_exists': nma_outdir.exists(),
            'files_in_dir': []
        }
        if nma_outdir.exists():
            debug_info['files_in_dir'] = [f.name for f in nma_outdir.iterdir()]
        
        # Check if analysis completed successfully
        if comprehensive_plot.exists() and mode_shapes_plot.exists():
            return jsonify({
                'success': True,
                'message': 'ProDy NMA analysis completed successfully',
                'has_plots': True,
                'output': result.stdout,
                'result_dir': f'{pdb_basename}_NMA_PCA_results',
                'has_pca_comprehensive': pca_comprehensive.exists(),
                'has_anm_pca_comparison': anm_pca_comparison.exists(),
                'has_clustenmd': clustenmd_analysis.exists(),
                'debug': debug_info
            })
        else:
            return jsonify({
                'success': False,
                'error': f'ProDy NMA analysis completed but output plots not found.',
                'output': result.stdout,
                'error_output': result.stderr,
                'debug': debug_info
            })
    
    except subprocess.TimeoutExpired:
        return jsonify({'success': False, 'error': 'ProDy NMA analysis timed out (10 minutes)'})
    except Exception as e:
        return jsonify({'success': False, 'error': f'ProDy NMA error: {str(e)}'})


@app.route('/download_nma/<filename>')
def download_nma(filename):
    """Download ProDy NMA plot files"""
    try:
        # Allow specific filenames for security
        allowed_files = ['comprehensive_analysis.png', 'mode_shapes_combined.png',
                        'pca_comprehensive_analysis.png', 'anm_pca_comparison.png', 
                        'clustenmd_analysis.png']
        
        if filename not in allowed_files:
            return jsonify({'error': 'Invalid filename'}), 400
        
        session_dir = get_session_dir()
        pdb_basename = 'predicted_structure'
        nma_dir = session_dir / f'{pdb_basename}_NMA_PCA_results'
        file_path = nma_dir / filename
        
        if file_path.exists():
            return send_file(file_path, mimetype='image/png')
        else:
            return jsonify({'error': f'File not found: {filename}'}), 404
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/download_nma_zip')
def download_nma_zip():
    """Download entire NMA results folder as ZIP"""
    try:
        session_dir = get_session_dir()
        pdb_basename = 'predicted_structure'
        nma_dir = session_dir / f'{pdb_basename}_NMA_PCA_results'
        
        if not nma_dir.exists():
            return jsonify({'error': 'NMA results not found'}), 404
        
        # Create zip file in session directory
        zip_path = session_dir / 'NMA_results.zip'
        
        import zipfile
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, dirs, files in os.walk(nma_dir):
                for file in files:
                    file_path = Path(root) / file
                    arcname = file_path.relative_to(session_dir)
                    zipf.write(file_path, arcname)
        
        return send_file(zip_path, as_attachment=True, download_name='ProDy_NMA_Results.zip')
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/standalone_nma')
def standalone_nma_page():
    """Render standalone NMA analysis page"""
    return render_template('standalone_nma.html')


@app.route('/upload_standalone_pdb', methods=['POST'])
def upload_standalone_pdb():
    """Upload PDB file for standalone NMA analysis"""
    try:
        if 'pdb_file' not in request.files:
            return jsonify({'success': False, 'error': 'No file uploaded'})
        
        file = request.files['pdb_file']
        if file.filename == '':
            return jsonify({'success': False, 'error': 'No file selected'})
        
        if not file.filename.endswith('.pdb'):
            return jsonify({'success': False, 'error': 'Invalid file type. Please upload a PDB file'})
        
        # Get session directory
        session_dir = get_session_dir()
        
        # Save uploaded PDB file
        pdb_filename = 'uploaded_structure.pdb'
        pdb_path = session_dir / pdb_filename
        file.save(str(pdb_path))
        
        return jsonify({
            'success': True,
            'filename': file.filename,
            'pdb_path': str(pdb_path)
        })
        
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@app.route('/run_standalone_nma', methods=['POST'])
def run_standalone_nma():
    """Run ProDy NMA on uploaded PDB file"""
    try:
        data = request.get_json()
        pdb_path = data.get('pdb_path', '')
        method = data.get('method', 'ANM')
        cutoff = data.get('cutoff')
        n_modes = data.get('n_modes', 20)
        conformers = data.get('conformers', 10)
        interface_cutoff = data.get('interface_cutoff', 10.0)
        
        # Validate inputs
        if not pdb_path or not Path(pdb_path).exists():
            return jsonify({'success': False, 'error': 'PDB file not found'})
        
        valid_methods = ['ANM', 'GNM']
        if method not in valid_methods:
            return jsonify({'success': False, 'error': f'Invalid method. Choose from {valid_methods}'})
        
        # Validate cutoff if provided, default to 10.0
        if cutoff is None:
            cutoff = 10.0
        try:
            cutoff = float(cutoff)
            if cutoff < 5.0 or cutoff > 30.0:
                return jsonify({'success': False, 'error': 'Cutoff distance must be between 5.0 and 30.0 Å'})
        except ValueError:
            return jsonify({'success': False, 'error': 'Invalid cutoff value'})
        
        # Validate interface_cutoff
        try:
            interface_cutoff = float(interface_cutoff)
            if interface_cutoff < 5.0 or interface_cutoff > 30.0:
                return jsonify({'success': False, 'error': 'Interface cutoff must be between 5.0 and 30.0 Å'})
        except ValueError:
            return jsonify({'success': False, 'error': 'Invalid interface cutoff value'})
        
        try:
            n_modes = int(n_modes)
            if n_modes < 1 or n_modes > 100:
                return jsonify({'success': False, 'error': 'Number of modes must be between 1 and 100'})
        except ValueError:
            return jsonify({'success': False, 'error': 'Invalid number of modes'})
        
        try:
            conformers = int(conformers)
            if conformers < 1 or conformers > 50:
                return jsonify({'success': False, 'error': 'Number of conformers must be between 1 and 50'})
        except ValueError:
            return jsonify({'success': False, 'error': 'Invalid number of conformers'})
        
        # Get session directory
        session_dir = Path(pdb_path).parent
        
        # Build command - always include ensemble-pca, allosteric, and clustenmd flags
        cmd = ['python3', str(BASE_DIR / 'prody_nma.py'), pdb_path, 
               '--method', method, '--nmodes', str(n_modes), '--conformers', str(conformers),
               '--cutoff', str(cutoff), '--interface-cutoff', str(interface_cutoff),
               '--ensemble-pca', '--allosteric', '--clustenmd']
        
        # Run ProDy NMA analysis
        import subprocess
        result = subprocess.run(
            cmd,
            cwd=str(session_dir),
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout
        )
        
        # The output directory is created by prody_nma.py as {pdb_basename}_NMA_PCA_results
        pdb_basename = Path(pdb_path).stem  # 'uploaded_structure'
        nma_outdir = session_dir / f'{pdb_basename}_NMA_PCA_results'
        
        # Check for output files
        comprehensive_plot = nma_outdir / 'comprehensive_analysis.png'
        mode_shapes_plot = nma_outdir / 'mode_shapes_combined.png'
        pca_comprehensive = nma_outdir / 'pca_comprehensive_analysis.png'
        anm_pca_comparison = nma_outdir / 'anm_pca_comparison.png'
        clustenmd_analysis = nma_outdir / 'clustenmd_analysis.png'
        
        # Debug: List all files in output directory
        debug_info = {
            'nma_outdir': str(nma_outdir),
            'nma_outdir_exists': nma_outdir.exists(),
            'files_in_dir': []
        }
        if nma_outdir.exists():
            debug_info['files_in_dir'] = [f.name for f in nma_outdir.iterdir()]
        
        # Check if analysis completed successfully
        if comprehensive_plot.exists() and mode_shapes_plot.exists():
            # Store NMA output directory in session for downloads
            session['standalone_nma_dir'] = str(nma_outdir)
            
            return jsonify({
                'success': True,
                'message': 'ProDy CG-NMA analysis completed successfully',
                'has_plots': True,
                'has_pca_comprehensive': pca_comprehensive.exists(),
                'has_anm_pca_comparison': anm_pca_comparison.exists(),
                'has_clustenmd': clustenmd_analysis.exists(),
                'output': result.stdout,
                'result_dir': f'{pdb_basename}_NMA_PCA_results',
                'debug': debug_info
            })
        else:
            return jsonify({
                'success': False,
                'error': f'ProDy CG-NMA analysis completed but output plots not found.',
                'output': result.stdout,
                'error_output': result.stderr,
                'debug': debug_info
            })
    
    except subprocess.TimeoutExpired:
        return jsonify({'success': False, 'error': 'ProDy CG-NMA analysis timed out (10 minutes)'})
    except Exception as e:
        return jsonify({'success': False, 'error': f'ProDy CG-NMA error: {str(e)}'})


@app.route('/download_standalone_nma/<filename>')
def download_standalone_nma(filename):
    """Download standalone NMA plot files"""
    try:
        # Only allow specific filenames for security
        allowed_files = ['comprehensive_analysis.png', 'mode_shapes_combined.png',
                        'pca_comprehensive_analysis.png', 'anm_pca_comparison.png', 
                        'clustenmd_analysis.png']
        if filename not in allowed_files:
            return jsonify({'error': 'Invalid filename'}), 400
        
        nma_dir_str = session.get('standalone_nma_dir', None)
        if not nma_dir_str:
            return jsonify({'error': 'No NMA results found'}), 404
        
        nma_dir = Path(nma_dir_str)
        file_path = nma_dir / filename
        
        if file_path.exists():
            return send_file(file_path, mimetype='image/png')
        else:
            return jsonify({'error': f'File not found: {filename}'}), 404
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/download_standalone_nma_zip')
def download_standalone_nma_zip():
    """Download entire standalone NMA results folder as ZIP"""
    try:
        nma_dir_str = session.get('standalone_nma_dir', None)
        if not nma_dir_str:
            return jsonify({'error': 'No NMA results found'}), 404
        
        nma_dir = Path(nma_dir_str)
        if not nma_dir.exists():
            return jsonify({'error': 'NMA results not found'}), 404
        
        session_dir = nma_dir.parent
        
        # Create zip file in session directory
        zip_path = session_dir / 'Standalone_NMA_Results.zip'
        
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, dirs, files in os.walk(nma_dir):
                for file in files:
                    file_path = Path(root) / file
                    arcname = file_path.relative_to(session_dir)
                    zipf.write(file_path, arcname)
        
        return send_file(zip_path, as_attachment=True, download_name='Standalone_NMA_Results.zip')
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/anm_comparison')
def anm_comparison_page():
    """Render ANM comparison analysis page"""
    # Force create a new session for each visit
    session.pop('session_id', None)
    session.pop('anm_comparison_dir', None)
    get_session_dir()
    return render_template('anm_comparison.html')


@app.route('/upload_anm_comparison_pdbs', methods=['POST'])
def upload_anm_comparison_pdbs():
    """Upload apo and complex PDB files for ANM comparison"""
    try:
        if 'apo_pdb' not in request.files or 'complex_pdb' not in request.files:
            return jsonify({'success': False, 'error': 'Both apo and complex PDB files are required'})
        
        apo_file = request.files['apo_pdb']
        complex_file = request.files['complex_pdb']
        
        if apo_file.filename == '' or complex_file.filename == '':
            return jsonify({'success': False, 'error': 'Please select both files'})
        
        if not apo_file.filename.endswith('.pdb') or not complex_file.filename.endswith('.pdb'):
            return jsonify({'success': False, 'error': 'Invalid file type. Please upload PDB files'})
        
        # Get session directory
        session_dir = get_session_dir()
        
        # Save uploaded PDB files
        apo_path = session_dir / 'apo.pdb'
        complex_path = session_dir / 'complex.pdb'
        apo_file.save(str(apo_path))
        complex_file.save(str(complex_path))
        
        return jsonify({
            'success': True,
            'apo_filename': apo_file.filename,
            'complex_filename': complex_file.filename,
            'apo_path': str(apo_path),
            'complex_path': str(complex_path)
        })
        
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@app.route('/run_anm_comparison', methods=['POST'])
def run_anm_comparison():
    """Run ANM comparison analysis between apo and complex structures"""
    try:
        data = request.get_json()
        apo_path = data.get('apo_path', '')
        complex_path = data.get('complex_path', '')
        apo_chain = data.get('apo_chain', 'A')
        complex_chain = data.get('complex_chain', 'B')
        n_modes = data.get('n_modes', 20)
        
        # Validate inputs
        if not apo_path or not Path(apo_path).exists():
            return jsonify({'success': False, 'error': 'Apo PDB file not found'})
        
        if not complex_path or not Path(complex_path).exists():
            return jsonify({'success': False, 'error': 'Complex PDB file not found'})
        
        # Validate chain IDs
        if not apo_chain or len(apo_chain) != 1:
            return jsonify({'success': False, 'error': 'Apo chain ID must be a single character'})
        
        if not complex_chain or len(complex_chain) != 1:
            return jsonify({'success': False, 'error': 'Complex chain ID must be a single character'})
        
        # Validate n_modes
        try:
            n_modes = int(n_modes)
            if n_modes < 1 or n_modes > 100:
                return jsonify({'success': False, 'error': 'Number of modes must be between 1 and 100'})
        except ValueError:
            return jsonify({'success': False, 'error': 'Invalid number of modes'})
        
        # Get session directory
        session_dir = Path(apo_path).parent
        
        # Output directory
        output_dir = session_dir / 'anm_comparison'
        
        # Check if anm_comparison.py exists
        anm_script = BASE_DIR / 'anm_comparison.py'
        if not anm_script.exists():
            return jsonify({
                'success': False,
                'error': f'ANM comparison script not found at {anm_script}'
            })
        
        # Build command
        cmd = ['python3', str(anm_script),
               '--apo', str(apo_path),
               '--complex', str(complex_path),
               '--c_apo', apo_chain,
               '--c_complex', complex_chain,
               '--modes', str(n_modes),
               '--outdir', str(output_dir)]
        
        # Debug: log the command
        print(f"Running ANM comparison command: {' '.join(cmd)}")
        print(f"Working directory: {session_dir}")
        
        # Run ANM comparison analysis
        import subprocess
        result = subprocess.run(
            cmd,
            cwd=str(session_dir),
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout
        )
        
        # Debug: log output
        print(f"Return code: {result.returncode}")
        print(f"STDOUT: {result.stdout}")
        print(f"STDERR: {result.stderr}")
        
        if result.returncode != 0:
            # Try to get available chains from the PDB files
            try:
                from prody import parsePDB
                apo_struct = parsePDB(str(apo_path))
                complex_struct = parsePDB(str(complex_path))
                
                apo_chains = set([chain.getChid() for chain in apo_struct.iterChains()])
                complex_chains = set([chain.getChid() for chain in complex_struct.iterChains()])
                
                chain_info = f"\nAvailable chains:\n  Apo ({Path(apo_path).name}): {sorted(apo_chains)}\n  Complex ({Path(complex_path).name}): {sorted(complex_chains)}"
                error_msg = f'ANM comparison analysis failed. {chain_info}\n\nYou specified: apo_chain={apo_chain}, complex_chain={complex_chain}\n\nError output:\n{result.stderr}'
            except:
                chain_info = ""
                error_msg = f'ANM comparison analysis failed.\n\nError output:\n{result.stderr}'
            
            return jsonify({
                'success': False,
                'error': error_msg,
                'output': result.stdout,
                'error_output': result.stderr,
                'command': ' '.join(cmd)
            })
        
        # Check if output directory exists
        print(f"\n=== ANM COMPARISON DEBUG ===")
        print(f"Session dir: {session_dir}")
        print(f"Output dir: {output_dir}")
        print(f"Output dir exists: {output_dir.exists()}")
        
        if output_dir.exists():
            plot_files = list(output_dir.glob('*.png'))
            print(f"Plot files found: {plot_files}")
            print(f"Plot file names: {[f.name for f in plot_files]}")
            print(f"=== END DEBUG ===\n")
            
            # Store the output directory name for downloads
            session['anm_comparison_dir'] = output_dir.name
            
            return jsonify({
                'success': True,
                'message': 'ANM comparison analysis completed successfully',
                'output': result.stdout,
                'plot_count': len(plot_files),
                'plots': [f.name for f in plot_files]
            })
        else:
            # List what's actually in session_dir
            session_contents = list(session_dir.iterdir()) if session_dir.exists() else []
            return jsonify({
                'success': False,
                'error': 'Analysis completed but output directory not found',
                'output': result.stdout,
                'session_dir': str(session_dir),
                'session_contents': [str(f) for f in session_contents]
            })
            
    except subprocess.TimeoutExpired:
        return jsonify({
            'success': False,
            'error': 'Analysis timed out. Try reducing the number of modes.'
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@app.route('/download_anm_comparison/<filename>')
def download_anm_comparison(filename):
    """Download a specific ANM comparison plot file"""
    try:
        session_id = session.get('session_id')
        if not session_id:
            return "Session expired", 404
        
        session_dir = RESULTS_DIR / session_id
        dir_name = session.get('anm_comparison_dir', 'anm_comparison')
        output_dir = session_dir / dir_name
        
        file_path = output_dir / filename
        
        if not file_path.exists():
            return f"File not found: {filename}", 404
        
        return send_file(
            str(file_path),
            as_attachment=True,
            download_name=filename
        )
    except Exception as e:
        return str(e), 500


@app.route('/download_anm_comparison_zip')
def download_anm_comparison_zip():
    """Download all ANM comparison results as a ZIP file"""
    try:
        session_id = session.get('session_id')
        if not session_id:
            return "Session expired", 404
        
        session_dir = RESULTS_DIR / session_id
        dir_name = session.get('anm_comparison_dir', 'anm_comparison')
        output_dir = session_dir / dir_name
        
        if not output_dir.exists():
            return "Results directory not found", 404
        
        # Create a ZIP file in memory
        memory_file = io.BytesIO()
        with zipfile.ZipFile(memory_file, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for file_path in output_dir.rglob('*'):
                if file_path.is_file():
                    arcname = file_path.relative_to(output_dir)
                    zipf.write(file_path, arcname)
        
        memory_file.seek(0)
        
        return send_file(
            memory_file,
            mimetype='application/zip',
            as_attachment=True,
            download_name='anm_comparison_results.zip'
        )
    except Exception as e:
        return str(e), 500


@app.route('/deformability')
def deformability_page():
    """Render deformability analysis page"""
    return render_template('deformability.html')


@app.route('/upload_deformability_pdbs', methods=['POST'])
def upload_deformability_pdbs():
    """Upload reference and target PDB files for deformability analysis"""
    try:
        if 'reference_pdb' not in request.files or 'target_pdb' not in request.files:
            return jsonify({'success': False, 'error': 'Both reference and target PDB files are required'})
        
        ref_file = request.files['reference_pdb']
        target_file = request.files['target_pdb']
        
        if ref_file.filename == '' or target_file.filename == '':
            return jsonify({'success': False, 'error': 'Please select both files'})
        
        if not ref_file.filename.endswith('.pdb') or not target_file.filename.endswith('.pdb'):
            return jsonify({'success': False, 'error': 'Invalid file type. Please upload PDB files'})
        
        # Get session directory
        session_dir = get_session_dir()
        
        # Save uploaded PDB files
        ref_path = session_dir / 'reference.pdb'
        target_path = session_dir / 'target.pdb'
        ref_file.save(str(ref_path))
        target_file.save(str(target_path))
        
        return jsonify({
            'success': True,
            'reference_filename': ref_file.filename,
            'target_filename': target_file.filename,
            'ref_path': str(ref_path),
            'target_path': str(target_path)
        })
        
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@app.route('/run_deformability_analysis', methods=['POST'])
def run_deformability_analysis():
    """Run deformability analysis comparing two PDB structures"""
    try:
        data = request.get_json()
        ref_path = data.get('ref_path', '')
        target_path = data.get('target_path', '')
        ref_chain = data.get('ref_chain', 'A')
        target_chain = data.get('target_chain', 'B')
        n_modes = data.get('n_modes', 20)
        
        # Validate inputs
        if not ref_path or not Path(ref_path).exists():
            return jsonify({'success': False, 'error': 'Reference PDB file not found'})
        
        if not target_path or not Path(target_path).exists():
            return jsonify({'success': False, 'error': 'Target PDB file not found'})
        
        # Validate chain IDs
        if not ref_chain or len(ref_chain) != 1:
            return jsonify({'success': False, 'error': 'Reference chain ID must be a single character'})
        
        if not target_chain or len(target_chain) != 1:
            return jsonify({'success': False, 'error': 'Target chain ID must be a single character'})
        
        # Validate n_modes
        try:
            n_modes = int(n_modes)
            if n_modes < 1 or n_modes > 100:
                return jsonify({'success': False, 'error': 'Number of modes must be between 1 and 100'})
        except ValueError:
            return jsonify({'success': False, 'error': 'Invalid number of modes'})
        
        # Get session directory
        session_dir = Path(ref_path).parent
        
        # Check if deformation_analysis.py exists
        deformation_script = BASE_DIR / 'deformation_analysis.py'
        if not deformation_script.exists():
            return jsonify({
                'success': False,
                'error': f'Deformation analysis script not found at {deformation_script}'
            })
        
        # Build command
        cmd = ['python3', str(deformation_script),
               '--reference', str(ref_path),
               '--target', str(target_path),
               '--ref_chain', ref_chain,
               '--target_chain', target_chain,
               '--nmodes', str(n_modes)]
        
        # Debug: log the command
        print(f"Running deformability analysis command: {' '.join(cmd)}")
        print(f"Working directory: {session_dir}")
        
        # Run deformability analysis
        import subprocess
        result = subprocess.run(
            cmd,
            cwd=str(session_dir),
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout
        )
        
        # Debug: log output
        print(f"Return code: {result.returncode}")
        print(f"STDOUT: {result.stdout}")
        print(f"STDERR: {result.stderr}")
        
        if result.returncode != 0:
            # Try to get available chains from the PDB files
            try:
                from prody import parsePDB
                ref_struct = parsePDB(str(ref_path))
                target_struct = parsePDB(str(target_path))
                
                ref_chains = set([chain.getChid() for chain in ref_struct.iterChains()])
                target_chains = set([chain.getChid() for chain in target_struct.iterChains()])
                
                chain_info = f"\nAvailable chains:\n  Reference ({ref_path.name}): {sorted(ref_chains)}\n  Target ({target_path.name}): {sorted(target_chains)}"
                error_msg = f'Deformability analysis failed. {chain_info}\n\nYou specified: ref_chain={ref_chain}, target_chain={target_chain}\n\nError output:\n{result.stderr}'
            except:
                chain_info = ""
                error_msg = f'Deformability analysis failed.\n\nError output:\n{result.stderr}'
            
            return jsonify({
                'success': False,
                'error': error_msg,
                'output': result.stdout,
                'error_output': result.stderr,
                'command': ' '.join(cmd)
            })
        
        # The script creates output in {reference_basename}_Deformation_Analysis/
        # Since we save as 'reference.pdb', the output should be 'reference_Deformation_Analysis'
        output_dir = session_dir / 'reference_Deformation_Analysis'
        
        # Debug: check what exists
        print(f"\n=== DEFORMABILITY DEBUG ===")
        print(f"Session dir: {session_dir}")
        print(f"Session dir exists: {session_dir.exists()}")
        print(f"Expected output dir: {output_dir}")
        print(f"Output dir exists: {output_dir.exists()}")
        
        if session_dir.exists():
            print(f"Session dir contents: {list(session_dir.iterdir())}")
            # Also check for any _Deformation_Analysis directories
            deform_dirs = [d for d in session_dir.iterdir() if d.is_dir() and '_Deformation_Analysis' in d.name]
            print(f"Directories with '_Deformation_Analysis': {deform_dirs}")
        
        if output_dir.exists():
            plot_files = list(output_dir.glob('*.png'))
            print(f"Plot files found: {plot_files}")
            print(f"Plot file names: {[f.name for f in plot_files]}")
            print(f"=== END DEBUG ===\n")
            
            # Store the output directory name for downloads
            session['deformability_dir'] = output_dir.name
            
            return jsonify({
                'success': True,
                'message': 'Deformability analysis completed successfully',
                'output': result.stdout,
                'plot_count': len(plot_files),
                'plots': [f.name for f in plot_files]
            })
        else:
            # List what's actually in session_dir
            session_contents = list(session_dir.iterdir()) if session_dir.exists() else []
            return jsonify({
                'success': False,
                'error': 'Analysis completed but output directory not found',
                'output': result.stdout,
                'error_output': result.stderr,
                'session_dir': str(session_dir),
                'session_contents': [f.name for f in session_contents]
            })
    
    except subprocess.TimeoutExpired:
        return jsonify({'success': False, 'error': 'Deformability analysis timed out (10 minutes)'})
    except Exception as e:
        return jsonify({'success': False, 'error': f'Deformability analysis error: {str(e)}'})


@app.route('/download_deformability/<filename>')
def download_deformability(filename):
    """Download deformability analysis plot files"""
    try:
        session_dir = get_session_dir()
        # Get the stored output directory name from session
        dir_name = session.get('deformability_dir', 'reference_Deformation_Analysis')
        output_dir = session_dir / dir_name
        file_path = output_dir / filename
        
        if file_path.exists() and file_path.suffix == '.png':
            return send_file(file_path, mimetype='image/png')
        else:
            return jsonify({'error': f'File not found: {filename}'}), 404
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/download_deformability_zip')
def download_deformability_zip():
    """Download entire deformability analysis results as ZIP"""
    try:
        session_dir = get_session_dir()
        # Get the stored output directory name from session
        dir_name = session.get('deformability_dir', 'reference_Deformation_Analysis')
        output_dir = session_dir / dir_name
        
        if not output_dir.exists():
            return jsonify({'error': 'Deformability analysis results not found'}), 404
        
        # Create zip file
        zip_path = session_dir / 'Deformability_Analysis.zip'
        
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, dirs, files in os.walk(output_dir):
                for file in files:
                    file_path = Path(root) / file
                    arcname = file_path.relative_to(session_dir)
                    zipf.write(file_path, arcname)
        
        return send_file(zip_path, as_attachment=True, download_name='Deformability_Analysis.zip')
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)

