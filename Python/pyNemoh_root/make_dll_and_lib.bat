cd dll_build
::mingw32-make clean
::mingw32-make
:: Alternative to lib.exe /verbose /DEF:libnemoh.def /MACHINE:X64 /OUT:libnemoh.lib as lib.exe is vs
dlltool -llibnemoh.lib -Dlibnemoh.dll -dlibnemoh.def
copy libnemoh.lib ..
copy libnemoh.dll ..\..\pyNemoh
cd..
