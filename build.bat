@echo off
setlocal EnableDelayedExpansion
cd /d "%~dp0"

set "VCVARS=C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvarsall.bat"
set "NINJA=C:\Program Files\Microsoft Visual Studio\18\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja\ninja.exe"
set "BUILD_DIR=build_msvc"

if not exist "%VCVARS%" (
    echo ERROR: vcvarsall.bat not found at %VCVARS%
    pause & exit /b 1
)
if not exist "%NINJA%" (
    echo ERROR: ninja.exe not found at %NINJA%
    pause & exit /b 1
)

call "%VCVARS%" x64 >nul 2>&1
if errorlevel 1 ( echo ERROR: vcvarsall.bat failed. & pause & exit /b 1 )

if not exist "%BUILD_DIR%" mkdir "%BUILD_DIR%"
cd "%BUILD_DIR%"

set VCPKG_ROOT=
cmake .. -G "Ninja" -DCMAKE_MAKE_PROGRAM="%NINJA%" -DCMAKE_BUILD_TYPE=Release
if errorlevel 1 ( echo ERROR: CMake configure failed. & cd .. & pause & exit /b 1 )

cmake --build . --config Release
if errorlevel 1 ( echo ERROR: Build failed. & cd .. & pause & exit /b 1 )

cd ..
echo.
echo BUILD SUCCESSFUL
echo Run: build_msvc\langevin_bench.exe
pause
