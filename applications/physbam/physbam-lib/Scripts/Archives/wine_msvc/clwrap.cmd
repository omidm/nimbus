@echo off
set VCHOME=z:\n\curvature\data\msvc8
set PATH=%PATH%;%VCHOME%\VC\redist\x86\Microsoft.VC80.CRT
set PATH=%PATH%;%VCHOME%\Common7\IDE
set PATH=%PATH%;%VCHOME%\VC\bin

@set INCLUDE=%VCHOME%\VC\ATLMFC\INCLUDE;%VCHOME%\VC\INCLUDE;%VCHOME%\VC\PlatformSDK\include;%VCHOME%\SDK\v2.0\include;%INCLUDE%
@set LIB=%VCHOME%\VC\ATLMFC\LIB;%VCHOME%\VC\LIB;%VCHOME%\VC\PlatformSDK\lib;%VCHOME%\SDK\v2.0\lib;%LIB%
@set LIBPATH=C:\WINDOWS\Microsoft.NET\Framework\v2.0.50727;%VCHOME%\VC\ATLMFC\LIB
%*
