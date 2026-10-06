@echo off
setlocal
cd /d "%~dp0"
where py >nul 2>nul
if not errorlevel 1 (
    py -3 simulador_celeste.py
) else (
    python simulador_celeste.py
)
if errorlevel 1 (
    echo.
    echo Confira se o Python 3.9 ou superior com Tcl/Tk esta disponivel.
    echo Leia o arquivo LEIA-ME.md para as instrucoes.
    pause
)
endlocal
