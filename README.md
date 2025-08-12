# DNA Visualization & Analysis System

A comprehensive bioinformatics tool for 3D DNA visualization and sequence analysis with interactive features for mutation simulation and annotation.
<<img width="1865" height="887" alt="image" src="https://github.com/user-attachments/assets/2df3e398-5265-4334-8310-9cf99f038fef" />



## 🧬 Features

### 3D Visualization
- **Accurate double helix rendering** with molecular geometry
- **Interactive camera controls** (zoom, rotate, pan)
- **Real-time mutation visualization**
- **Color-coded base pairs** (A=Red, T=Teal, C=Blue, G=Yellow)
- **Annotation highlighting** for ORFs, CpG islands, and promoters

### Sequence Analysis
- **GC content calculation**
- **Open Reading Frame (ORF) detection**
- **CpG island identification**
- **Promoter region prediction**
- **Molecular weight calculation**
- **Complement and reverse complement generation**

### Mutation Simulation
- **Single Nucleotide Polymorphisms (SNPs)**
- **Insertions and deletions**
- **Instant visual updates**
- **Mutation tracking and reset functionality**

### File Format Support
- **FASTA** (.fasta, .fa)
- **GenBank** (.gbk, .gb)
- **Plain text** sequences
- **Export capabilities** (PNG, FASTA, JSON)

## 🚀 Installation

### Prerequisites
- Python 3.8 or higher
- Modern web browser with WebGL support

### Step 1: Clone or Download Files
Create a new directory and save the following files:
- `app.py` (Backend server)
- `index.html` (Frontend interface)
- `requirements.txt` (Python dependencies)

### Step 2: Install Python Dependencies
```bash
# Create virtual environment (recommended)
python -m venv dna_analyzer
source dna_analyzer/bin/activate  # On Windows: dna_analyzer\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Step 3: Start the Backend Server
```bash
python app.py
```
The server will start on `http://localhost:5000`

### Step 4: Open the Frontend
Open `index.html` in your web browser, or serve it with a simple HTTP server:
```bash
# Option 1: Direct file opening
# Simply double-click index.html

# Option 2: Local server (recommended)
python -m http.server 8080
# Then visit http://localhost:8080
```

## 📊 Usage Guide

### Loading DNA Sequences

#### Method 1: File Upload
1. Click **"Upload FASTA/GenBank"** button
2. Select your sequence file (.fasta, .fa, .gbk, .gb)
3. The sequence will be automatically loaded and visualized

#### Method 2: Direct Input
1. Paste your DNA sequence in the text area
2. Click **"Load Sequence"**
3. Only A, T, C, G characters are accepted

#### Method 3: Sample Data
1. Click **"Load Sample Data"** for a demonstration
2. Uses E. coli lacZ promoter region with 1000+ base pairs

### 3D Visualization Controls

#### Mouse Controls
- **Left click + drag**: Rotate the DNA model
- **Mouse wheel**: Zoom in/out
- **Right panel sliders**: 
  - Zoom Level: Fine-tune camera distance
  - Helix Pitch: Adjust spacing between base pairs
  - Base Size: Scale nucleotide spheres

#### View Controls
- **Reset View**: Return to default camera position
- **Toggle Rotation**: Enable/disable automatic rotation

### Sequence Analysis

#### Running Analysis
1. Load a DNA sequence
2. Click **"Run Analysis"** button
3. Results appear in the statistics panel and annotations list

#### Analysis Features
- **Basic Statistics**: Length, GC content, AT content, molecular weight
- **ORF Detection**: Finds open reading frames in all 6 reading frames
- **CpG Islands**: Identifies regions with high CG dinucleotide frequency
- **Promoter Prediction**: Locates potential TATA box sequences

### Mutation Simulation

#### Applying Mutations
1. Select mutation type:
   - **SNP**: Single nucleotide change
   - **Insertion**: Add nucleotides
   - **Deletion**: Remove nucleotides
2. Enter position (0-indexed)
3. Choose new base (for SNP/insertion)
4. Click **"Apply Mutation"**

#### Visual Updates
- Mutations are instantly reflected in the 3D model
- Statistics update automatically
- Re-run analysis to see effect on features

#### Reset Functionality
- Click **"Reset"** to restore original sequence
- All mutations are cleared

### Color Coding System

#### Base Pairs
- 🔴 **Red**: Adenine (A)
- 🟢 **Teal**: Thymine (T)  
- 🔵 **Blue**: Cytosine (C)
- 🟡 **Yellow**: Guanine (G)

#### Annotations (when analysis is run)
- 🟣 **Purple**: Open Reading Frames (ORFs)
- 🌸 **Pink**: CpG Islands
- 🟢 **Green**: Promoter Regions

### Export Options

#### PNG Export
- Captures current 3D view as high-resolution image
- Useful for presentations and publications

#### FASTA Export
- Downloads sequence in standard FASTA format
- Includes any applied mutations

#### Analysis Export
- Saves complete analysis results as JSON
- Contains sequence, statistics, and all annotations

## 🧪 Sample Datasets

### Included Sample
The built-in sample contains:
- **E. coli lacZ promoter region** (1000+ bp)
- Multiple ORFs in different reading frames
- CpG islands and regulatory elements
- TATA box sequences for promoter prediction

### Creating Custom Datasets
For testing, you can use sequences from:
- **NCBI GenBank**: https://www.ncbi.nlm.nih.gov/genbank/
- **EMBL-EBI**: https://www.ebi.ac.uk/ena/
- **UniProt**: https://www.uniprot.org/

## 🔧 Technical Architecture

### Backend (Python/Flask)
- **Flask**: Web framework for REST API
- **Biopython**: Sequence analysis and file parsing
- **NumPy**: Numerical computations for statistics

### Frontend (JavaScript/Three.js)
- **Three.js**: 3D visualization and rendering
- **WebGL**: Hardware-accelerated graphics
- **Responsive CSS**: Modern UI with gradient backgrounds

### API Endpoints
- `POST /api/upload`: Upload sequence files or text
- `POST /api/analyze`: Run comprehensive analysis
- `POST /api/mutate`: Apply mutations to sequence
- `GET /api/sequence`: Retrieve current sequence
- `GET /api/sample`: Load demonstration data
- `POST /api/export/fasta`: Export sequence file

## 🚀 Performance Optimization

### Visualization Limits
- **Maximum bases rendered**: 1000 (for smooth performance)
- **Automatic downsampling**: For longer sequences
- **LOD (Level of Detail)**: Reduces complexity at distance

### Memory Management
- **Efficient geometry reuse**: Minimizes WebGL buffer usage
- **Progressive loading**: Large files processed in chunks
- **Garbage collection**: Proper cleanup of Three.js objects

## 🔬 Scientific Accuracy

### Molecular Geometry
- **Helix parameters**: 10.5 base pairs per turn
- **Rise per base**: 3.4 Å (scaled for visualization)
- **Major/minor grooves**: Represented through backbone positioning
- **Hydrogen bonds**: Visualized as connecting cylinders

### Analysis Algorithms
- **ORF detection**: Standard start/stop codon recognition
- **CpG islands**: Sliding window with O/E ratio calculation
- **GC content**: Precise nucleotide counting
- **Promoter prediction**: TATA box pattern matching

## 🐛 Troubleshooting

### Common Issues

#### Server Won't Start
```bash
# Check Python version
python --version  # Should be 3.8+

# Install missing dependencies
pip install -r requirements.txt

# Check port availability
netstat -an | grep 5000
```

#### CORS Errors
- Ensure backend server is running on port 5000
- Try serving frontend with local HTTP server instead of file://

#### Performance Issues
- Reduce sequence length (< 1000 bp for best performance)
- Lower quality settings in browser if needed
- Close other GPU-intensive applications

#### WebGL Errors
- Update graphics drivers
- Enable hardware acceleration in browser
- Try different browser (Chrome/Firefox recommended)

### Browser Compatibility
- ✅ **Chrome 80+**: Full support
- ✅ **Firefox 75+**: Full support  
- ✅ **Safari 13+**: Full support
- ✅ **Edge 80+**: Full support

## 📚 API Documentation

### Upload Sequence
```javascript
POST /api/upload
Content-Type: application/json

{
    "sequence": "ATCGATCGATCG"
}
```

### Run Analysis
```javascript
POST /api/analyze
Content-Type: application/json

{
    "options": {
        "min_orf_length": 100
    }
}
```

### Apply Mutation
```javascript
POST /api/mutate
Content-Type: application/json

{
    "type": "SNP",
    "position": 10,
    "new_base": "G"
}
```

## 📖 Educational Use

### Learning Objectives
This tool helps students and researchers understand:
- **DNA structure**: Double helix geometry and base pairing
- **Sequence analysis**: ORF detection and feature annotation
- **Mutations**: Effects of genetic variations
- **Bioinformatics**: Computational sequence analysis

### Classroom Integration
- **Molecular biology courses**: Visual DNA structure
- **Bioinformatics training**: Hands-on sequence analysis
- **Research projects**: Mutation effect visualization
- **Presentations**: High-quality 3D renderings

## 🤝 Contributing

### Code Structure
- `app.py`: Flask backend with analysis algorithms
- `index.html`: Complete frontend with Three.js integration
- Modular design allows easy feature additions

### Future Enhancements
- **Protein translation**: Visualize amino acid sequences
- **Multiple sequence alignment**: Compare sequences
- **Real-time collaboration**: Share sessions
- **Advanced export**: GLTF 3D model export
- **Mobile support**: Touch controls optimization

## 📄 License

This project is intended for educational and research purposes. Please ensure proper attribution when using in academic work.

## 🆘 Support

For issues or questions:
1. Check the troubleshooting section above
2. Verify all dependencies are installed correctly
3. Test with the included sample data first
4. Check browser console for error messages

---

**Enjoy exploring the fascinating world of DNA with interactive 3D visualization!** 🧬✨
