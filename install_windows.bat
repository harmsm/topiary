@echo off
echo Setting up the topiary environment...

:: 1. Create the virtual environment folder named .topiary
python -m venv .topiary

:: 2. Check if the environment was created successfully
if not exist ".topiary\Scripts\activate" (
    echo.
    echo Error: Python was not found or failed to create the environment.
    exit /b
)

echo Activating environment and installing topiary...

:: 3. Run the activation and installation in the current shell
:: We use 'call' so the batch script continues after activation
call .topiary\Scripts\activate

:: 4. Upgrade pip and install the current directory
python -m pip install --upgrade pip
pip install .

echo.
echo Installation complete! 
echo To use topiary in the future, run: .topiary\Scripts\activate
echo 
echo Please make sure that muscle and ncbi tools are installed in the $PATH. 
echo See https://topiary-asr.readthedocs.io/en/latest/installation.html for
echo details.
