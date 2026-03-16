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

:: 4. Upgrade pip and install build dependencies
python -m pip install --upgrade pip
pip install setuptools Cython^<3.0 wheel

:: 5. Install ete4 dependencies manually since we'll use --no-build-isolation
:: These are base dependencies for building ete4
pip install numpy

:: 6. Install ete4 without build isolation to avoid Cython 3.0 path issues
pip install "ete4<4.4.0" --no-build-isolation

:: 7. Install the current directory
pip install .

:: 8. Install muscle 5.1
echo Downloading and installing muscle 5.1...
powershell -Command "Invoke-WebRequest -Uri https://github.com/rcedgar/muscle/releases/download/v5.1/muscle5.1.win64.exe -OutFile .topiary\Scripts\muscle.exe"

:: 9. Install NCBI BLAST+ 2.16.0
echo Downloading and installing NCBI BLAST+ 2.16.0...
powershell -Command "Invoke-WebRequest -Uri https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/ncbi-blast-2.16.0+-win64.zip -OutFile blast.zip"
powershell -Command "Expand-Archive -Path blast.zip -DestinationPath .topiary\_blast_temp"
xcopy /Y /E ".topiary\_blast_temp\ncbi-blast-2.16.0+\bin\*" ".topiary\Scripts\"
rmdir /S /Q ".topiary\_blast_temp"
del blast.zip

echo.
echo Installation complete! 
echo To use topiary in the future, run: .topiary\Scripts\activate
echo.
echo Muscle 5 and NCBI BLAST+ tools have been installed into your virtual environment.
echo.
