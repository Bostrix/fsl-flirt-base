# Makefile for the MEDx version of FLIRT: by Mark Jenkinson, 21/9/99

include ${FSLDIR}/etc/default.mk

DBGFLAGS = -g

PROJNAME = flirt

USRINCFLAGS =  -I${INCDIR}/newmat -I${INCDIR}/avwio -I${INCDIR}/mjimage_dev \
		-I${INCDIR}/miscmaths
LIBS = -lmjimage_dev -lmiscmaths -lavwio -lnewmat -lm 
INHERITVARS = "USRLDFLAGS=${USRLDFLAGS}" "USRINCFLAGS=${USRINCFLAGS}"

FL_OBJS = globaloptions.o costfns.o flirt.o 
T_OBJS = globaloptions.o costfns.o testcode.o
C_OBJS = convert_xfm.o
A_OBJS = avscale.o
R_OBJS = rmsdiff.o
TAL_OBJS = tal2epicoord.o
EPI_OBJS = epi2talcoord.o
IMG_OBJS = img2imgcoord.o

RUNTCLS = Flirt
XFILES = flirt convert_xfm avscale rmsdiff tal2epicoord epi2talcoord img2imgcoord
TESTXFILES = testcode
HFILES =
SCRIPTS = extract pairreg pairregrc1 pairregrc2 pairregrc3 fixxfm

all:	${XFILES}


flirt:    	${FL_OBJS}
	        $(CXX)  ${CXXFLAGS} ${LDFLAGS} -o $@ ${FL_OBJS} ${LIBS}


testcode:	${T_OBJS}
		$(CXX)  ${CXXFLAGS} ${LDFLAGS} -o $@ ${T_OBJS} ${LIBS}

convert_xfm:    ${C_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${C_OBJS} ${LIBS}

avscale:        ${A_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${A_OBJS} ${LIBS}

rmsdiff:        ${R_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${R_OBJS} ${LIBS}

tal2epicoord:	${TAL_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${TAL_OBJS} ${LIBS}

epi2talcoord:	${EPI_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${EPI_OBJS} ${LIBS}

img2imgcoord:	${IMG_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${IMG_OBJS} ${LIBS}




