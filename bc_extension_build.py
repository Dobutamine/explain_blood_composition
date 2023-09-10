from cffi import FFI
ffibuilder = FFI()

ffibuilder.cdef("""
                typedef struct {
                    int valid_ab;
                    int valid_o2;
                    double po2;
                    double so2;
                    double ph;
                    double pco2;
                    double hco3;
                    double be;
                    double sid_app;
                    double steps_ab;
                    double steps_o2;
                } bloodResult;

                bloodResult GetBloodComposition(
                    double _to2,
                    double _tco2,
                    double _sid,
                    double _albumin,
                    double _phosphates,
                    double _uma,
                    double _hemoglobin,
                    double _dpg,
                    double _temp
                );     
                """)

ffibuilder.set_source("_blood_composition",  # name of the output C extension
"""
    #include "blood_composition.h"
""",
    sources=['blood_composition.c'],   # includes pi.c as additional sources
    libraries=['m'])    # on Unix, link with the math library

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)