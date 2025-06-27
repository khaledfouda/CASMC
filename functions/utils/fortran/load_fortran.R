# gfortran -shared -fPIC -o suvC.so suvC.f
#gfortran -shared -fPIC -o suvC.so suvC.f90

dyn.load(paste0(path_to_code,"/utils/fortran/suvC.so"))
#dyn.load("./functions/SoftImpute_Incomplete/src/sparse_prod.so")
#dyn.unload("./functions/SoftImpute_Incomplete/src/sparse_prod.so")

