# Makefile for the MEDx version of FLIRT: by Mark Jenkinson, 21/9/99

include ${FSLDIR}/etc/default.mk

DBGFLAGS = -g

PROJNAME = flirt

USRINCFLAGS =  -I${INCDIR}/newmat -I${INCDIR}/avwio \
		-I${INCDIR}/miscmaths -I${INCDIR}/newimage
LIBS = -lnewimage -lmiscmaths -lavwio -lnewmat -lm 
INHERITVARS = "USRLDFLAGS=${USRLDFLAGS}" "USRINCFLAGS=${USRINCFLAGS}"

FL_OBJS = globaloptions.o flirt.o 
C_OBJS = convert_xfm.o
A_OBJS = avscale.o
R_OBJS = rmsdiff.o
TAL_OBJS = tal2imgcoord.o
EPI_OBJS = img2talcoord.o
IMG_OBJS = img2imgcoord.o
X_OBJS = applyxfm4D.o
P_OBJS = pointflirt.o

RUNTCLS = Flirt InvertXFM ApplyXFM
XFILES = flirt convert_xfm avscale rmsdiff tal2imgcoord img2talcoord \
	img2imgcoord applyxfm4D
TESTXFILES = 
HFILES =
SCRIPTS = extracttxt pairreg fixxfm

all:	${XFILES} schedule

schedule:
	/bin/cp -rf flirtsch ${FSLDIR}/etc

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




