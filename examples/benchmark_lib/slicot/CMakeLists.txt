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
	SB04MD.f
	SB04ND.f
	SB04QD.f
	SB04RD.f
	SB04MY.f
	SB04MU.f
	SB04MR.f
	SB04MW.f
	SB04NV.f
	SB04NW.f
	SB04NX.f
	SB04NY.f
	SB04QR.f
	SB04QU.f
	SB04QY.f
	SB04RV.f
	SB04RW.f
	SB04RX.f
	SB04RY.f
    MB02UU.f
    MB02UV.f
    select.f
    delctg.f
    dgegs.f
    )

IF( CMAKE_Fortran_COMPILER_ID STREQUAL "XL")
	SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qfixed -qnosave")
ENDIF()
add_library (slicot OBJECT ${SRC})

