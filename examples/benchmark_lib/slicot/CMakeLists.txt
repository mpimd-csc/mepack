SET (SRC
	MB01RD.f
	MB01RW.f
	SG03AD.f
	SG03AY.f
	SG03AX.f
	SB04PY.f
	SB04PX.f
	SB03MD.f
	SB03MW.f
	SB03MV.f
	SB03MX.f
	SB03MY.f
    MB02UU.f
    MB02UV.f
    select.f
    delctg.f
    )

IF( CMAKE_Fortran_COMPILER_ID STREQUAL "XL")
	SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qfixed -qnosave")
ENDIF()
add_library (slicot OBJECT ${SRC})

