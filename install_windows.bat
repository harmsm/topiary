@echo off
echo Setting up the topiary conda environment...

:: 1. Create the conda environment from environment_windows.yml
call conda env create -f environment_windows.yml -p .topiary -y
if %errorlevel% neq 0 (
    echo.
    echo Error: Failed to create conda environment.
    exit /b %errorlevel%
)

echo Activating environment and installing topiary...

:: 2. Activate the environment
:: Use 'call' so the batch script continues after activation
call conda activate .\.topiary
if %errorlevel% neq 0 (
    echo.
    echo Error: Failed to activate conda environment.
    exit /b %errorlevel%
)

:: 3. Install the current directory
pip install .
if %errorlevel% neq 0 exit /b %errorlevel%

:: 4. Install muscle 5.1
echo Downloading and installing muscle 5.1...
powershell -Command "Invoke-WebRequest -Uri https://github.com/rcedgar/muscle/releases/download/v5.1/muscle5.1.win64.exe -OutFile .topiary\Scripts\muscle.exe"
if %errorlevel% neq 0 exit /b %errorlevel%

:: 5. Install NCBI BLAST+ 2.16.0
echo Downloading and installing NCBI BLAST+ 2.16.0...
powershell -Command "Invoke-WebRequest -Uri https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/ncbi-blast-2.16.0+-win64.zip -OutFile blast.zip"
if %errorlevel% neq 0 exit /b %errorlevel%
powershell -Command "Expand-Archive -Path blast.zip -DestinationPath .topiary\_blast_temp"
if %errorlevel% neq 0 exit /b %errorlevel%
xcopy /Y /E ".topiary\_blast_temp\ncbi-blast-2.16.0+\bin\*" ".topiary\Scripts\"
if %errorlevel% neq 0 exit /b %errorlevel%
rmdir /S /Q ".topiary\_blast_temp"
if %errorlevel% neq 0 exit /b %errorlevel%
del blast.zip
if %errorlevel% neq 0 exit /b %errorlevel%

echo.
echo Installation complete! 
echo To use topiary in the future, run: conda activate .\.topiary
echo.
echo Muscle 5 and NCBI BLAST+ tools have been installed into your virtual environment.
echo.
