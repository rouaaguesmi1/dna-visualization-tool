"""
Professional DNA Visualization and Analysis Backend Server
Provides a REST API for advanced DNA sequence analysis using Biopython.
"""
import io
import re
import tempfile
import json
import os
from typing import Dict, List, Any

# --- Imports for the corrected 3D download feature ---
from flask import Flask, request, jsonify, send_file, after_this_request
import shutil # Used for reliably deleting directories
# ----------------------------------------------------

from flask_cors import CORS
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
from Bio.PDB import PDBList

# --- Hide Biopython's normal "partial codon" warning ---
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
# --------------------------------------------------------

# --- Flask App Initialization ---
app = Flask(__name__)
CORS(app)

# --- Core DNA Analysis Logic ---
class DNAAnalyzer:
    def __init__(self):
        self.sequence = ""
        self.original_sequence = ""
        self.sequence_id = "unknown"
        self.sequence_description = "Unknown sequence"

    def get_stats(self):
        """Helper function to get stats for the current sequence."""
        if not self.sequence:
            return {'length': 0, 'gc_content': 0.0}
        return {
            'length': len(self.sequence),
            'gc_content': round(gc_fraction(self.sequence) * 100, 2)
        }

    def load_sequence_from_string(self, content: str, file_format: str) -> Dict[str, Any]:
        try:
            handle = io.StringIO(content)
            records = list(SeqIO.parse(handle, file_format))
            if not records:
                raise ValueError(f"No valid {file_format} records found.")
            
            first_record = records[0]
            self.sequence = re.sub(r'[^ATCG]', '', str(first_record.seq).upper())
            if not self.sequence:
                raise ValueError("Sequence contains no valid DNA bases.")
            
            self.original_sequence = self.sequence
            self.sequence_id = first_record.id
            self.sequence_description = first_record.description
            
            return {
                "id": self.sequence_id,
                "description": self.sequence_description,
                "sequence": self.sequence,
                "length": len(self.sequence),
                "gc_content": round(gc_fraction(self.sequence) * 100, 2),
                "stats": self.get_stats()
            }
        except Exception as e:
            raise ValueError(f"Error parsing sequence: {str(e)}")

    def load_sequence_direct(self, sequence: str) -> Dict[str, Any]:
        """Load sequence directly from string input"""
        clean_sequence = re.sub(r'[^ATCG]', '', sequence.upper())
        if not clean_sequence:
            raise ValueError("Sequence contains no valid DNA bases.")
        
        self.sequence = clean_sequence
        self.original_sequence = clean_sequence
        self.sequence_id = "direct_input"
        self.sequence_description = "Directly input sequence"
        
        return {
            "id": self.sequence_id,
            "description": self.sequence_description,
            "sequence": self.sequence,
            "length": len(self.sequence),
            "gc_content": round(gc_fraction(self.sequence) * 100, 2),
            "stats": self.get_stats()
        }

    def find_orfs(self, min_length: int = 100) -> List[Dict[str, Any]]:
        """Find Open Reading Frames in the sequence"""
        orfs = []
        seq_obj = Seq(self.sequence)
        
        for strand, nuc_seq in [(1, seq_obj), (-1, seq_obj.reverse_complement())]:
            for frame in range(3):
                translated = nuc_seq[frame:].translate(to_stop=True)
                
                for match in re.finditer(r'M[^*]*', str(translated)):
                    orf_start_aa = match.start()
                    orf_end_aa = match.end()
                    orf_length_aa = orf_end_aa - orf_start_aa
                    orf_length_nt = orf_length_aa * 3
                    
                    if orf_length_nt >= min_length:
                        if strand == 1:
                            nt_start = frame + (orf_start_aa * 3)
                            nt_end = nt_start + orf_length_nt
                        else:
                            nt_end = len(self.sequence) - (frame + (orf_start_aa * 3))
                            nt_start = nt_end - orf_length_nt
                        
                        orfs.append({
                            'start': nt_start,
                            'end': nt_end,
                            'strand': strand,
                            'frame': frame + 1,
                            'length': orf_length_nt,
                            'protein_length': orf_length_aa,
                            'sequence': str(nuc_seq[nt_start:nt_end])
                        })
        
        return sorted(orfs, key=lambda x: x['length'], reverse=True)

    def find_cpg_islands(self, window_size: int = 200, min_gc: float = 0.5, min_obs_exp: float = 0.6) -> List[Dict[str, Any]]:
        """Find CpG islands in the sequence"""
        cpg_islands = []
        
        for i in range(len(self.sequence) - window_size + 1):
            window = self.sequence[i:i+window_size]
            gc_content = gc_fraction(window)
            
            if gc_content > min_gc:
                c_count = window.count('C')
                g_count = window.count('G')
                cg_count = window.count('CG')
                
                if c_count > 0 and g_count > 0:
                    observed_expected = (cg_count * window_size) / (c_count * g_count)
                    
                    if observed_expected > min_obs_exp:
                        cpg_islands.append({
                            'start': i,
                            'end': i + window_size,
                            'length': window_size,
                            'gc_content': round(gc_content * 100, 2),
                            'cpg_count': cg_count,
                            'observed_expected': round(observed_expected, 3)
                        })
        
        merged_islands = []
        for island in cpg_islands:
            if not merged_islands or island['start'] > merged_islands[-1]['end']:
                merged_islands.append(island)
            else:
                last_island = merged_islands[-1]
                last_island['end'] = max(last_island['end'], island['end'])
                last_island['length'] = last_island['end'] - last_island['start']
        
        return merged_islands

    def find_promoters(self) -> List[Dict[str, Any]]:
        """Find potential promoter regions (TATA boxes and other motifs)"""
        promoters = []
        expanded_tata = r'TATA[AT][AT][AG]?'
        
        for match in re.finditer(expanded_tata, self.sequence, re.IGNORECASE):
            start = match.start()
            end = match.end()
            promoter_start = max(0, start - 200)
            promoter_end = min(len(self.sequence), end + 50)
            
            promoters.append({
                'start': promoter_start,
                'end': promoter_end,
                'tata_start': start,
                'tata_end': end,
                'tata_sequence': match.group(),
                'type': 'TATA_box'
            })
        
        return promoters

    def apply_mutation(self, mut_type: str, pos: int, sequence: str = "", new_base: str = "", insert_seq: str = ""):
        """Apply mutation to the sequence"""
        if not (0 <= pos < len(self.sequence)):
            raise ValueError(f"Position {pos} is out of bounds (sequence length: {len(self.sequence)}).")
        
        if mut_type == 'SNP':
            if not new_base or new_base not in 'ATCG':
                raise ValueError("SNP requires a valid new base (A, T, C, G).")
            self.sequence = self.sequence[:pos] + new_base + self.sequence[pos + 1:]
            
        elif mut_type == 'insertion':
            insert_sequence = re.sub(r'[^ATCG]', '', (insert_seq or new_base).upper())
            if not insert_sequence:
                raise ValueError("Insert sequence must contain valid DNA bases (A, T, C, G).")
            self.sequence = self.sequence[:pos] + insert_sequence + self.sequence[pos:]
            
        elif mut_type == 'deletion':
            if pos >= len(self.sequence):
                raise ValueError("Cannot delete beyond sequence length.")
            self.sequence = self.sequence[:pos] + self.sequence[pos + 1:]
            
        else:
            raise ValueError(f"Unknown mutation type '{mut_type}'. Use 'SNP', 'insertion', or 'deletion'.")

    def reset_sequence(self):
        """Reset sequence to original"""
        self.sequence = self.original_sequence

    def get_sample_sequence(self) -> str:
        """Generate a sample DNA sequence for testing"""
        return ("ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCT"
                "GACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAA"
                "AGGTAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGA"
                "ACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGT"
                "CCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAAAACCATTAGCGGC"
                "CAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTTGACGGGACTCG"
                "CCGCCGCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTTGCCC"
                "AAATAAAACATGTCCTGCATGGCATTAGTTTGTTGGTTCAATGGTGGCAAATACACCACCGGAG"
                "GTGGAGATGATTTACGGTCGTCGCATGTATGCAGCGGGTGGCGATTTTGCCGATTTGGCCGCC"
                "CGTGATGATGTTGATGTTGCCGAACGCCATCGTGATGTTGAAGCCGATAAACGTCACATAACC"
                "TGGGATAATGTTAATGGGCAACGTGGTGGCGGTGGCCGTGATGGCGGTCACGGTGATGGCGGT"
                "TGGGACGCTCGTGATGTGGACGGCGGTGATGTTGATGATGTCGGCGGTGATGCCGATGACGGC"
                "CGTGATGGTGCCGGTCGTGATGGCGGCCGTGATGGCGGTGATGCCGGCGATGGTGGCGGCCGT"
                "GATGGCGCCGGTGATGGCGGTGATGCCGGTGATGTTGATGATGGCGGTGATGCCGGCTAG")

# Initialize analyzer
analyzer = DNAAnalyzer()

def api_response(success, data=None, error=None):
    """Standardized API response format"""
    return jsonify({"success": success, "data": data, "error": error})

# --- API Endpoints ---

@app.route('/api/upload', methods=['POST'])
def upload_sequence_endpoint():
    """Upload and load DNA sequence from file or direct input."""
    try:
        if 'file' in request.files:
            file = request.files['file']
            if file.filename == '':
                return api_response(False, error="No file selected.")
            
            content = file.read().decode('utf-8')
            fmt = 'genbank' if file.filename.lower().endswith(('.gb', '.gbk')) else 'fasta'
            result = analyzer.load_sequence_from_string(content, fmt)
            
        elif request.json and 'sequence' in request.json:
            sequence = request.json['sequence'].strip()
            if not sequence:
                return api_response(False, error="Empty sequence provided.")
            
            if sequence.startswith('>'):
                result = analyzer.load_sequence_from_string(sequence, 'fasta')
            else:
                result = analyzer.load_sequence_direct(sequence)
        else:
            return api_response(False, error="No file or sequence data provided.")
        
        return api_response(True, data=result)
        
    except Exception as e:
        return api_response(False, error=str(e))


@app.route('/api/download/pdb', methods=['POST'])
def download_pdb_endpoint():
    """
    Downloads a 3D structure file (.pdb) from the RCSB PDB database.
    This version correctly handles temporary file cleanup on all OS, including Windows.
    """
    data = request.json
    if not data or 'pdb_id' not in data:
        return api_response(False, error="Request must include a 'pdb_id'."), 400

    pdb_id = data['pdb_id'].strip().lower()
    if len(pdb_id) != 4 or not re.match(r'^[a-z0-9]+$', pdb_id):
        return api_response(False, error="Invalid PDB ID format. Must be 4 alphanumeric characters."), 400

    # Create a temporary directory that we will manage manually to avoid file lock issues.
    temp_dir = tempfile.mkdtemp()
    
    try:
        pdbl = PDBList()
        # Download the PDB file into our temporary directory.
        file_path = pdbl.retrieve_pdb_file(pdb_id, pdir=temp_dir, file_format="pdb")

        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Could not retrieve file for PDB ID '{pdb_id}'. The ID may be invalid or obsolete.")

        # Define a cleanup function to run AFTER the file has been sent.
        @after_this_request
        def cleanup(response):
            try:
                shutil.rmtree(temp_dir)
            except Exception as error:
                app.logger.error(f"Error cleaning up temp directory {temp_dir}: {error}")
            return response

        # Send the file. The 'cleanup' function will run automatically afterwards.
        return send_file(
            file_path,
            as_attachment=True,
            download_name=f"{pdb_id}.pdb",
            mimetype='chemical/x-pdb'
        )
    except Exception as e:
        # If an error happens before sending the file, clean up the directory immediately.
        shutil.rmtree(temp_dir)
        app.logger.error(f"Failed to download PDB for {pdb_id}: {e}")
        return api_response(False, error=f"An error occurred during download: {str(e)}"), 500


@app.route('/api/sample', methods=['GET'])
def load_sample_endpoint():
    """Load a sample DNA sequence for testing"""
    try:
        sample_seq = analyzer.get_sample_sequence()
        result = analyzer.load_sequence_direct(sample_seq)
        result['description'] = "Sample E. coli gene sequence"
        result['id'] = "sample_ecoli"
        return api_response(True, data=result)
    except Exception as e:
        return api_response(False, error=str(e))

@app.route('/api/analyze', methods=['POST'])
def analyze_sequence_endpoint():
    """Analyze the loaded DNA sequence"""
    if not analyzer.sequence:
        return api_response(False, error="No sequence loaded.")
    try:
        options = request.json.get('options', {}) if request.json else {}
        min_orf_length = options.get('min_orf_length', 100)
        orfs = analyzer.find_orfs(min_orf_length)
        cpg_islands = analyzer.find_cpg_islands()
        promoters = analyzer.find_promoters()
        results = {
            "orfs": orfs,
            "cpg_islands": cpg_islands,
            "promoters": promoters,
            "sequence_stats": analyzer.get_stats()
        }
        return api_response(True, data=results)
    except Exception as e:
        return api_response(False, error=f"Analysis failed: {e}")

@app.route('/api/mutate', methods=['POST'])
def mutate_sequence_endpoint():
    """Apply mutation to the loaded sequence"""
    if not analyzer.sequence:
        return api_response(False, error="No sequence loaded.")
    try:
        data = request.json
        if not data:
            return api_response(False, error="No mutation data provided.")
        
        mut_type = data.get('type')
        position = data.get('position')
        new_base = data.get('new_base', '')
        sequence = data.get('sequence', '')
        insert_seq = data.get('insert_seq', '')
        
        if position is None:
            return api_response(False, error="Position is required.")
        position = int(position)
        
        analyzer.apply_mutation(mut_type, position, sequence, new_base, insert_seq)
        
        result = {
            "new_sequence": analyzer.sequence,
            "stats": analyzer.get_stats()
        }
        return api_response(True, data=result)
    except Exception as e:
        return api_response(False, error=str(e))

@app.route('/api/reset', methods=['POST'])
def reset_mutations_endpoint():
    """Reset sequence to original state"""
    if not analyzer.original_sequence:
        return api_response(False, error="No original sequence to reset to.")
    try:
        analyzer.reset_sequence()
        result = {
            "sequence": analyzer.sequence,
            "stats": analyzer.get_stats()
        }
        return api_response(True, data=result)
    except Exception as e:
        return api_response(False, error=str(e))

@app.route('/api/export/fasta', methods=['POST'])
def export_fasta_endpoint():
    """Export current sequence as FASTA file"""
    if not analyzer.sequence:
        return api_response(False, error="No sequence to export.")
    try:
        data = request.json if request.json else {}
        seq_id = data.get('id', analyzer.sequence_id)
        description = data.get('description', analyzer.sequence_description)
        
        record = SeqRecord(Seq(analyzer.sequence), id=seq_id, description=description)
        
        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fasta") as f:
            SeqIO.write(record, f, 'fasta')
            temp_filename = f.name
        
        @after_this_request
        def cleanup(response):
            try:
                os.remove(temp_filename)
            except Exception as error:
                app.logger.error(f"Error cleaning up FASTA file: {error}")
            return response

        return send_file(
            temp_filename, 
            as_attachment=True, 
            download_name=f"{seq_id}.fasta",
            mimetype='text/plain'
        )
    except Exception as e:
        return api_response(False, error=f"Export failed: {str(e)}")

@app.route('/api/export/data', methods=['POST'])
def export_analysis_data_endpoint():
    """Export analysis data as JSON"""
    if not analyzer.sequence:
        return api_response(False, error="No sequence data to export.")
    try:
        analysis_data = {
            "orfs": analyzer.find_orfs(),
            "cpg_islands": analyzer.find_cpg_islands(),
            "promoters": analyzer.find_promoters()
        }
        
        export_data = {
            "sequence": {
                "id": analyzer.sequence_id,
                "description": analyzer.sequence_description,
                "sequence": analyzer.sequence,
                "original_sequence": analyzer.original_sequence,
                "stats": analyzer.get_stats()
            },
            "analysis": analysis_data,
            "metadata": {"version": "1.0"}
        }
        return api_response(True, data=export_data)
    except Exception as e:
        return api_response(False, error=f"Data export failed: {str(e)}")

@app.route('/api/sequence/stats', methods=['GET'])
def get_sequence_stats_endpoint():
    """Get current sequence statistics"""
    if not analyzer.sequence:
        return api_response(False, error="No sequence loaded.")
    try:
        stats = analyzer.get_stats()
        stats.update({
            "id": analyzer.sequence_id,
            "description": analyzer.sequence_description,
            "has_mutations": analyzer.sequence != analyzer.original_sequence
        })
        return api_response(True, data=stats)
    except Exception as e:
        return api_response(False, error=str(e))

@app.route('/api/health', methods=['GET'])
def health_check():
    """Health check endpoint"""
    return api_response(True, data={"status": "healthy", "version": "1.o"})

# Error handlers
@app.errorhandler(404)
def not_found(error):
    return api_response(False, error="Endpoint not found"), 404

@app.errorhandler(500)
def internal_error(error):
    return api_response(False, error="Internal server error"), 500

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)