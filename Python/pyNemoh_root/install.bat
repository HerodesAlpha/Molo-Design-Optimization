echo off
::set "PATH=%PATH%;C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build"
@CALL C:\Users\eison\AppData\Local\Continuum\anaconda3\Scripts\activate.bat C:\Users\eison\AppData\Local\Continuum\anaconda3
::@CALL  conda activate py27

copy .\dll_build\libnemoh.def .
copy  .\dll_build\libnemoh.dll .\pyNemoh
python.exe setup.py cleanall
python.exe setup.py build_ext --inplace
python.exe setup.py install


