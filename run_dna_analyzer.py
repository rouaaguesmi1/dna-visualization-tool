# run_dna_analyzer.py
"""
Quick Start Script for the Professional DNA Visualization System
Automatically sets up the environment and launches the application.
"""

import subprocess
import sys
import os
import webbrowser
import time
from pathlib import Path

def check_python_version():
    """Checks for a compatible Python version (3.8+)."""
    print("--- Checking Prerequisites ---")
    if sys.version_info < (3, 8):
        print(f"❌ ERROR: Python 3.8 or higher is required. You are using {sys.version.split()[0]}.")
        return False
    print(f"✅ Python version {sys.version.split()[0]} is compatible.")
    return True

def install_dependencies():
    """Installs required Python packages from requirements.txt."""
    requirements_file = Path("requirements.txt")
    if not requirements_file.exists():
        print(f"❌ ERROR: '{requirements_file}' not found. Cannot install dependencies.")
        return False
        
    print("\n--- Installing Dependencies ---")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", requirements_file])
        print("✅ Dependencies installed successfully.")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ ERROR: Failed to install dependencies. See error below:\n{e}")
        print("Please try running this command manually in your terminal: pip install -r requirements.txt")
        return False

def check_files():
    """Checks if all required application files exist."""
    required_files = ['app.py', 'index.html']
    missing_files = [file for file in required_files if not Path(file).exists()]
    
    if missing_files:
        print(f"❌ ERROR: Missing required files: {', '.join(missing_files)}")
        return False
    
    print("✅ All required application files found.")
    return True

def start_backend():
    """Starts the Flask backend server as a background process."""
    print("\n--- Starting Backend Server ---")
    try:
        # Use Popen to run the server in a non-blocking background process
        process = subprocess.Popen([sys.executable, "app.py"])
        print("✅ Backend server is starting up on http://localhost:5000")
        return process
    except Exception as e:
        print(f"❌ ERROR: Failed to start the backend server: {e}")
        return None

def open_frontend():
    """Opens the frontend HTML file in the default web browser."""
    print("\n--- Opening Frontend Interface ---")
    # Wait a moment for the server to initialize
    time.sleep(3)
    
    frontend_path = Path("index.html").resolve()
    frontend_url = frontend_path.as_uri()
    
    try:
        webbrowser.open(frontend_url)
        print(f"✅ Frontend opened in your browser. URL: {frontend_url}")
    except Exception as e:
        print(f"❌ WARNING: Could not automatically open the browser: {e}")
        print(f"Please manually open this file in your browser: {frontend_path}")

def main():
    """Main function to orchestrate the application startup."""
    print("===============================================")
    print("  DNA Visualization & Analysis System - Startup")
    print("===============================================")
    
    if not (check_python_version() and check_files() and install_dependencies()):
        sys.exit(1) # Exit if prerequisites or installation fail

    backend_process = start_backend()
    if not backend_process:
        sys.exit(1) # Exit if backend fails to start

    open_frontend()
    
    print("\n===============================================")
    print("   🚀 Application is now running! 🚀")
    print("===============================================")
    print("Backend API is active at: http://localhost:5000")
    print("To stop the application, close this terminal window or press Ctrl+C.")
    
    try:
        # Wait for the backend process to complete (it won't, until terminated)
        backend_process.wait()
    except KeyboardInterrupt:
        print("\n\n--- Shutting Down ---")
        backend_process.terminate()
        backend_process.wait() # Ensure the process is fully terminated
        print("✅ Backend server has been stopped.")
        print("===============================================")

if __name__ == "__main__":
    main()
