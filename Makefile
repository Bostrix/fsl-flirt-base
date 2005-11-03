# Makefile for FLIRT

include ${FSLCONFDIR}/default.mk

PROJNAME = flirt

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB}

LIBS = -lnewimage -lmiscmaths -lprob -lfslio -lniftiio -lznz -lnewmat -lutils -lm -lz

FL_OBJS = globaloptions.o flirt.o 
C_OBJS = convert_xfm.o
A_OBJS = avscale.o
R_OBJS = rmsdiff.o
TAL_OBJS = tal2imgcoord.o
EPI_OBJS = img2talcoord.o
IMG_OBJS = img2imgcoord.o
X_OBJS = applyxfm4D.o
P_OBJS = pointflirt.o
M_OBJS = makerot.o
MID_OBJS = midtrans.o
NL_OBJS = nonlin.o

RUNTCLS = Flirt InvertXFM ApplyXFM InvertMEDxXFM ConcatXFM Nudge
XFILES = flirt convert_xfm avscale rmsdiff tal2imgcoord img2talcoord \
	img2imgcoord applyxfm4D pointflirt makerot midtrans
TESTXFILES = nonlin
HFILES =
SCRIPTS = extracttxt pairreg fixxfm standard_space_roi

all:	${XFILES} schedule

schedule:
	@if [ ! -d ${DESTDIR}/etc ] ; then ${MKDIR} ${DESTDIR}/etc ; ${CHMOD} g+w ${DESTDIR}/etc ; fi
	${CP} -rf flirtsch ${DESTDIR}/etc

flirt:    	${FL_OBJS}
	        $(CXX)  ${CXXFLAGS} ${LDFLAGS} -o $@ ${FL_OBJS} ${LIBS}


convert_xfm:    ${C_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${C_OBJS} ${LIBS}

avscale:        ${A_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${A_OBJS} ${LIBS}

rmsdiff:        ${R_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${R_OBJS} ${LIBS}

tal2imgcoord:	${TAL_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${TAL_OBJS} ${LIBS}

img2talcoord:	${EPI_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${EPI_OBJS} ${LIBS}

img2imgcoord:	${IMG_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${IMG_OBJS} ${LIBS}

applyxfm4D:	${X_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${X_OBJS} ${LIBS}

pointflirt:	${P_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${P_OBJS} ${LIBS}

makerot:	${M_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${M_OBJS} ${LIBS}

midtrans:	${MID_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${MID_OBJS} ${LIBS}

nonlin:		${NL_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${NL_OBJS} ${LIBS}
