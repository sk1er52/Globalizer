@echo off
setlocal

set START_DIR=%cd%
set ROOT_DIR=%~dp0\..


cd /d "%ROOT_DIR%"

if not exist "build_64" mkdir build_64
cd build_64

call conda init

echo [1/5] Creating a Conda Environment...
call conda create -p "%ROOT_DIR%build_64\Globalizer_env" python=3.11 -y

echo [2/5] activate a Conda Environment...
call conda activate "%ROOT_DIR%build_64\Globalizer_env"

echo [3/5] Installing the library...
call conda install -p "%ROOT_DIR%build_64\Globalizer_env" numpy -y
call conda install -p "%ROOT_DIR%build_64\Globalizer_env" jsonschema -y

echo [4/5] CMake Configuration...
call cmake -G "Visual Studio 17 2022" -DGLOBALIZER_BUILD_PROBLEMS=ON -DBUILD_ALL_TASK=ON -DGLOBALIZER_MAX_DIMENSION=130 -DGLOBALIZER_MAX_Number_Of_Function=70 ..

if %errorlevel% neq 0 goto error

echo [5/5] Opening Visual Studio...
if exist "globalizer.sln" (
    start "" "globalizer.sln"
) else (
    echo Error: globalizer.sln not found!
)

cd /d "%START_DIR%"
echo.

exit /b 0

