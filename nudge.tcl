#{{{ copyright and setup

#   nudge - interactive (if slow) alignment correction
#
#   Stephen Smith, FMRIB Image Analysis Group
#
#   Copyright (C) 2005-2006 University of Oxford
#
#   TCLCOPYRIGHT

source [ file dirname [ info script ] ]/fslstart.tcl

set VARS(history) {}

#}}}
#{{{ nudge

proc nudge { w } {
    global FSLDIR VARS PWD nvars TMP

    toplevel $w
    wm title      $w "Nudge"
    wm iconname   $w "Nudge"
    wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm

    set TMP [ exec sh -c "${FSLDIR}/bin/tmpnam /tmp/nudge" ]

    #{{{ select images and initial transform

frame $w.inputs -relief raised -borderwidth 1

FSLFileEntry $w.inputs.input \
	-variable nvars(input) \
	-pattern "IMAGE" \
	-directory $PWD \
	-label "Input image" \
	-title "Select the input image" \
	-width 40 \
	-filterhist VARS(history) \
        -command "nudge_inputrange"

FSLFileEntry $w.inputs.reference \
	-variable nvars(reference) \
	-pattern "IMAGE" \
	-directory $PWD \
	-label "Reference image" \
	-title "Select the reference image" \
	-width 40 \
	-filterhist VARS(history)

FSLFileEntry $w.inputs.initial_xfm \
        -variable nvars(initial_xfm) \
        -pattern "*.mat" \
        -directory $PWD \
        -label "Initial transformation (optional)" \
	-title "Select the initial transformation matrix" \
	-width 40 \
	-filterhist VARS(history)

pack $w.inputs.input $w.inputs.reference $w.inputs.initial_xfm -in $w.inputs -padx 5 -pady 5 -side top -anchor w -expand yes

#}}}
    #{{{ nudge

frame $w.nudge -relief raised -borderwidth 1

#{{{ rot

frame $w.nudge.rot

label $w.nudge.rot.label -text "Rotation (degrees):"

set nvars(xrot) 0
set nvars(yrot) 0
set nvars(zrot) 0
set nvars(rotinc) 5

tixControl $w.nudge.rot.x -label "X " -variable nvars(xrot) -step $nvars(rotinc) -selectmode immediate -options { entry.width 5 }
tixControl $w.nudge.rot.y -label "Y " -variable nvars(yrot) -step $nvars(rotinc) -selectmode immediate -options { entry.width 5 }
tixControl $w.nudge.rot.z -label "Z " -variable nvars(zrot) -step $nvars(rotinc) -selectmode immediate -options { entry.width 5 }
tixControl $w.nudge.rot.inc -label "increment " -variable nvars(rotinc) -step 1 -selectmode immediate -options { entry.width 5 } -command "nudge_update $w"

pack $w.nudge.rot.label $w.nudge.rot.x $w.nudge.rot.y $w.nudge.rot.z $w.nudge.rot.inc -in $w.nudge.rot -padx 5 -side left -expand yes

#}}}
#{{{ trans

frame $w.nudge.trans

label $w.nudge.trans.label -text "Translation (mm):"

set nvars(xtrans) 0
set nvars(ytrans) 0
set nvars(ztrans) 0
set nvars(transinc) 5

tixControl $w.nudge.trans.x -label "X " -variable nvars(xtrans) -step $nvars(transinc) -selectmode immediate -options { entry.width 5 }
tixControl $w.nudge.trans.y -label "Y " -variable nvars(ytrans) -step $nvars(transinc) -selectmode immediate -options { entry.width 5 }
tixControl $w.nudge.trans.z -label "Z " -variable nvars(ztrans) -step $nvars(transinc) -selectmode immediate -options { entry.width 5 }
tixControl $w.nudge.trans.inc -label "increment " -variable nvars(transinc) -step 1 -selectmode immediate -options { entry.width 5 } -command "nudge_update $w"

pack $w.nudge.trans.label $w.nudge.trans.x $w.nudge.trans.y $w.nudge.trans.z $w.nudge.trans.inc -in $w.nudge.trans -padx 5 -side left -expand yes

#}}}
#{{{ scale

frame $w.nudge.scale

label $w.nudge.scale.label -text "Scaling (X):"

set nvars(xscale) 1
set nvars(yscale) 1
set nvars(zscale) 1
set nvars(scaleinc) .01

tixControl $w.nudge.scale.x -label "X " -variable nvars(xscale) -step $nvars(scaleinc) -selectmode immediate -options { entry.width 5 }
tixControl $w.nudge.scale.y -label "Y " -variable nvars(yscale) -step $nvars(scaleinc) -selectmode immediate -options { entry.width 5 }
tixControl $w.nudge.scale.z -label "Z " -variable nvars(zscale) -step $nvars(scaleinc) -selectmode immediate -options { entry.width 5 }
tixControl $w.nudge.scale.inc -label "increment " -variable nvars(scaleinc) -step .01 -selectmode immediate -options { entry.width 5 } -command "nudge_update $w"

pack $w.nudge.scale.label $w.nudge.scale.x $w.nudge.scale.y $w.nudge.scale.z $w.nudge.scale.inc -in $w.nudge.scale -padx 5 -side left -expand yes

#}}}

pack $w.nudge.rot $w.nudge.trans $w.nudge.scale -in $w.nudge -padx 5 -pady 5 -side top -anchor w

#}}}
    #{{{ view options

frame $w.fslview -relief raised -borderwidth 1

label $w.fslview.label -text "FSLView options for reference and input"

entry $w.fslview.reference -textvariable nvars(fslview_reference) -width 50

set nvars(fslview_input) "-l Red-Yellow -b 3,10"
entry $w.fslview.input -textvariable nvars(fslview_input) -width 50

pack $w.fslview.label $w.fslview.reference $w.fslview.input -in $w.fslview -padx 5 -pady 5 -side top -anchor w -expand yes

#}}}
    #{{{ bottom buttons

frame $w.btns -relief raised -borderwidth 1

button $w.btns.apply  -text "View" -command "nudge_run $w 0"

pack $w.btns.apply -in $w.btns -padx 5 -pady 5 -side left -expand yes

#}}}

    pack $w.inputs $w.nudge $w.fslview $w.btns -in $w -padx 5 -pady 5 -fill x 
}

#}}}
#{{{ nudge_update

proc nudge_update { w dummy } {
    global FSLDIR nvars

    $w.nudge.rot.x configure -step $nvars(rotinc)
    $w.nudge.rot.y configure -step $nvars(rotinc)
    $w.nudge.rot.z configure -step $nvars(rotinc)

    $w.nudge.trans.x configure -step $nvars(transinc)
    $w.nudge.trans.y configure -step $nvars(transinc)
    $w.nudge.trans.z configure -step $nvars(transinc)

    $w.nudge.scale.x configure -step $nvars(scaleinc)
    $w.nudge.scale.y configure -step $nvars(scaleinc)
    $w.nudge.scale.z configure -step $nvars(scaleinc)
}

#}}}
#{{{ nudge_inputrange

proc nudge_inputrange { dummy } {
    global FSLDIR nvars

    set minmax [ fsl:exec "${FSLDIR}/bin/avwstats++ $nvars(input) -l 0.1 -r" ]
    set nvars(fslview_input) "-l Red-Yellow -b [ lindex $minmax 0 ],[ lindex $minmax 1 ]"
}

#}}}
#{{{ nudge_run

proc nudge_run { w dummy } {
    global FSLDIR nvars fslviewpid TMP

    fsl:exec "${FSLDIR}/bin/makerot -o ${TMP}_tmp.xfm --cov=$nvars(reference) -a 1,0,0 -t $nvars(xrot)"
    if { [ file exists $nvars(initial_xfm) ] } {
	fsl:exec "${FSLDIR}/bin/convert_xfm -omat ${TMP}.xfm -concat ${TMP}_tmp.xfm $nvars(initial_xfm)"
    } else {
	fsl:exec "cp ${TMP}_tmp.xfm ${TMP}.xfm"
    }

    fsl:exec "${FSLDIR}/bin/makerot -o ${TMP}_tmp.xfm --cov=$nvars(reference) -a 0,1,0 -t $nvars(yrot)"
    fsl:exec "${FSLDIR}/bin/convert_xfm -omat ${TMP}.xfm -concat ${TMP}_tmp.xfm ${TMP}.xfm"

    fsl:exec "${FSLDIR}/bin/makerot -o ${TMP}_tmp.xfm --cov=$nvars(reference) -a 0,0,1 -t $nvars(zrot)"
    fsl:exec "${FSLDIR}/bin/convert_xfm -omat ${TMP}.xfm -concat ${TMP}_tmp.xfm ${TMP}.xfm"

    set nudgexfm [ open ${TMP}_tmp.xfm "w" ]
    puts $nudgexfm "$nvars(xscale) 0 0 $nvars(xtrans)
0 $nvars(yscale) 0 $nvars(ytrans)
0 0 $nvars(zscale) $nvars(ztrans)
0 0 0 1
"
    close $nudgexfm
    fsl:exec "${FSLDIR}/bin/convert_xfm -omat ${TMP}.xfm -concat ${TMP}_tmp.xfm ${TMP}.xfm"

    fsl:exec "${FSLDIR}/bin/flirt -ref $nvars(reference) -in $nvars(input) -out $TMP -applyxfm -init ${TMP}.xfm"

    if { [ info exists fslviewpid ] } {
	set fslviewpid [ expr $fslviewpid + 1 ]
	catch { exec sh -c "kill -9 $fslviewpid" } errmsg
    }

    set FSLVIEW /usr/local/bin/fslview
    if { [ file exists ${FSLDIR}/bin/fslview ] } {
	set FSLVIEW ${FSLDIR}/bin/fslview
    }

    set fslviewpid [ exec sh -c "$FSLVIEW $nvars(reference) $nvars(fslview_reference) $TMP $nvars(fslview_input)" & ]

    puts "final transform is in ${TMP}.xfm"
}

#}}}
#{{{ call GUI and wait

if { ! [ info exists INGUI ] } {
    wm withdraw .
    nudge .r
    tkwait window .r
}

#}}}

