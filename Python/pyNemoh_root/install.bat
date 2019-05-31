 copy ..\..\Fortran\Build\fNemoh_Build\libnemoh.dll .\pyNemoh\libnemoh.dll
del libnemoh.lib
dlltool -d..\..\Fortran\Build\fNemoh_Build\libnemoh.def -D.\pyNemoh\libnemoh.dll -llibnemoh.lib
python.exe setup.py cleanall
python.exe setup.py build_ext --inplace
python.exe setup.py install
pause

