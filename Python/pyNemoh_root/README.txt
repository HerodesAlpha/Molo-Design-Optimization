cp ../../Fortran/Build/fNemoh_Build/libnemoh.dll ./pyNemoh/libnemoh.dll
dlltool -d../../Fortran/Build/fNemoh_Build/libnemoh.def -D./pyNemoh/libnemoh.dll -llibnemoh.lib