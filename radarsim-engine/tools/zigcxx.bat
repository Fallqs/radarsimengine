@echo off
rem CMake compiler wrapper for the conda zig toolchain (self-contained
rem mingw-w64 target; no Visual Studio needed). Used on machines without
rem MSVC. Requires the conda env with the zig package:
rem   conda install -n <env> -c conda-forge zig
rem
rem   cmake -G Ninja -DCMAKE_CXX_COMPILER=tools/zigcxx.bat -DCMAKE_RC_COMPILER=tools/zigwindres.bat
rem
rem ZIGCXX_ENV points at the conda env; defaults to %CONDA_PREFIX%\..\envs\rsimdev.
if "%ZIGCXX_ENV%"=="" set ZIGCXX_ENV=%CONDA_PREFIX%\..\envs\rsimdev
set ZIGBIN=%ZIGCXX_ENV%\Library\bin
set PREFIX=%ZIGCXX_ENV%
"%ZIGBIN%\x86_64-w64-mingw32-zig-cxx.exe" -target x86_64-windows-gnu %*
