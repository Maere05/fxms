# Force Validation Notebook Setup

This folder contains the Python helper code and notebook for organizing FluidX3D force-validation runs.

## What to Install

The easiest setup is:

1. Install Python 3.11 or newer from https://www.python.org/downloads/.
2. In VS Code, install these extensions:
   - Python
   - Jupyter
3. Open this repo folder in VS Code.
4. Open a terminal in the repo root and run:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r py\requirements.txt
```

If PowerShell refuses to run `Activate.ps1`, allow locally-created scripts for your user account once:

```powershell
Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy RemoteSigned
```

Then open `py/force_validation_lab.ipynb` and select the `.venv` kernel when VS Code asks.

## What the Notebook Does

- Runs local benchmark commands such as `bin/FluidX3D.exe skijumper 512 0`.
- Keeps human labels for runs.
- Loads CSV files from `bin/export/force_validation/`.
- Plots force convergence and coefficient trends.
- Includes a dry-run remote execution template for a future Linux server workflow.

Remote execution is intentionally a scaffold for now. Once server hostname, user, repo path, and copy method are known, fill in the remote config cell in the notebook.
