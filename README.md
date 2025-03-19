# NFGDA (Neuro-Fuzzy Gust Front Detection Algorithm)

This repository provides tools and scripts for the **Neuro-Fuzzy Gust Front Detection Algorithm** (NFGDA). Follow the instructions below to set up and run the system.

## Setup Instructions

### 1. Clone the repository

```bash
git clone git@github.com:TMinRay/NFGDA.git
```

### 2. Set up paths

Run the following command to configure necessary paths:

```bash
cd NFGDA/main
bash set_path.sh
```

### 3. Download or copy data

You have two options to obtain the required data:

#### Option 1: Copy data to the `V06` directory

If you already have the data, copy it to the `V06` directory.

#### Option 2: Download the data

Run the script `download_nexrad_l2_data.py` to automatically download the data:

```bash
python download_nexrad_l2_data.py
```

### 4. Edit the INI configuration file

Open the `main/NFGDA.ini` file and set the `case_name` to your desired case. Example:

```ini
case_name = KABX20200705_21
```

### 5. Convert the data

Navigate to the `python_accessories` directory and run the following script to convert the data into the required format:

```bash
cd ../python_accessories
python NF01_convert_V06_to_mat.py
```

#### Optional: Labeling

If needed, you can copy the labels to the following directory:

```bash
cp [case_id]_labels/ V06/[case_id]/
```

Alternatively, you can run the label script with GUI manual labeling:

```bash
python NF07_label_evalbox.py
```

### 6. Running the analysis

You can choose between running the analysis with **MATLAB** or **Python**.

#### Option 1: MATLAB

Run the following MATLAB script in `main` for analysis:

```matlab
NF_single_run.m
```

#### Option 2: Python

Alternatively, run the Python script in `main`:

```bash
python NFGDA.py
```

#### Additional MATLAB Visualization

To visualize the results, run the following MATLAB script:

```matlab
NFFIG_v2.m
```

---
For any issues or further information, feel free to open an issue on this repository.
