@echo off
rem RC compiler counterpart to zigcxx.bat (see its header).
if "%ZIGCXX_ENV%"=="" set ZIGCXX_ENV=%CONDA_PREFIX%\..\envs\rsimdev
set ZIGBIN=%ZIGCXX_ENV%\Library\bin
set PREFIX=%ZIGCXX_ENV%
"%ZIGBIN%\x86_64-w64-mingw32-zig-windres.exe" %*
