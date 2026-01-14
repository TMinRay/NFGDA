# NFGDA (Neuro-Fuzzy Gust Front Detection Algorithm)

This repository provides tools and scripts for the **Neuro-Fuzzy Gust Front Detection Algorithm** (NFGDA). Follow the instructions below to set up and run the system.

## Setup Instructions

### 1. Clone the repository

```bash
git clone git@github.com:TMinRay/NFGDA.git
cd NFGDA
```

### 2. Set up venv

Create and activate a virtual environment, then install dependencies:

```bash
# (Optional) deactivate any existing virtual environment
deactivate

# Create a new virtual environment
python3.12 -m venv ~/nfgda

# Activate the virtual environment
source ~/nfgda/bin/activate

# Install the package in editable mode
python -m pip install -e .
```

### 3. (Optional) Configure for an event or a different radar site

Edit the `scripts/NFGDA.ini` to select the radar site and time range.

#### Real-time forecasting

```ini
radar_id = KABX
custom_start_time = None
custom_end_time =   None
```

#### Historic event analysis
```ini
radar_id = KABX
custom_start_time = 2023,07,03,01,15,23
custom_end_time =   2023,07,03,03,18,58
```
#### Configure output and runtime directories
```ini
export_preds_dir     = ./runtime/nfgda_detection/
export_forecast_dir  = ./runtime/forecast/
V06_dir              = ./runtime/V06/
```

### 4. Run NFGDA Server

```bash
cd scripts
python NFGDA_Host.py
```
