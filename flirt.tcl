#

# FLIRT - FMRIB's Linear Image Registration Tool
#
# Mark Jenkinson, Stephen Smith and Matthew Webster, FMRIB Image Analysis Group
#
# Copyright (C) 1999-2001 University of Oxford
#
# TCLCOPYRIGHT

#{{{ setup

source [ file dirname [ info script ] ]/fslstart.tcl

catch { exec sh -c "${FSLDIR}/bin/flirt 2>&1 | grep -i 'flirt version' | awk '{ print \$3 }' -" } VERSION

#}}}
#{{{ flirt

proc flirt { w } {

    #{{{ setup main window etc

    global reg entries USER FSLDIR argc argv PWD PADY IGS VERSION gui_ext

    set reg($w,maxnstats) 20

    # ---- Set up Frames ----
    toplevel $w
    wm title $w "FLIRT - FMRIB's Linear Image Registration Tool - v${VERSION}"
    wm iconname $w "FLIRT"
    wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm
    frame $w.f

    TitleFrame $w.f.basic -relief groove 
    set lfbasic [ $w.f.basic getframe ]

	set PADY 3
	set IG  "image"
	set IGS "images"

#}}}
    #{{{ number of secondary images
    TitleFrame $w.f.stats -relief groove 
    set lfstats [ $w.f.stats getframe ]

    set reg($w,nstats) 0
    LabelSpinBox $w.f.nstats -label "Number of secondary $IGS to apply transform to " -textvariable reg($w,nstats) -range " 0 10000 $reg($w,maxnstats) "  -command " $w.f.nstats.spin.e validate; flirt:updatestats $w" -modifycmd  " flirt:updatestats $w"
    pack $w.f.nstats -in $lfstats -side top -anchor w -pady 3 -padx 5

#}}}
    #{{{ standard image

	set entries($w,1) ${FSLDIR}/etc/standard/avg152T1_brain.hdr

    FileEntry  $w.f.ref -textvariable entries($w,1) -label "Reference image   " -title "Select" -width 50 -filedialog directory  -filetypes IMAGE

#}}}

    #{{{ DOF
frame $w.f.dof
label $w.f.dof.label -text "  Model/DOF (input to ref)"
optionMenu2 $w.f.dof.menu reg($w,dof) 2Dmenu  "   2D to 2D registration" 2D "Rigid Body (3 parameter model)" 3Dmenu "   3D to 3D registration" TRANS "Translation Only (3 parameter model)" 6 "Rigid Body (6 parameter model)" 7 "Global Rescale (7 parameter model)" 9 "Traditional (9 parameter model)" 12 "Affine (12 parameter model)"
pack $w.f.dof.label $w.f.dof.menu -in $w.f.dof -side top -side left

  $w.f.dof.menu entryconfigure 0 -state disabled


#tixOptionMenu $w.f.dof -label "  Model/DOF (input to ref)" -variable reg($w,dof) -options {label.anchor w;menubutton.width 30}
    #$w.f.dof add command 2Dmenu -label "   2D to 2D registration" -state disabled -background "#555555"
    #$w.f.dof add command 2D -label "Rigid Body (3 parameter model)"
    #$w.f.dof add command 3Dmenu -label "   3D to 3D registration" -state disabled -background "#555555"
    #$w.f.dof add command TRANS -label "Translation Only (3 parameter model)"
    #$w.f.dof add command 6 -label "Rigid Body (6 parameter model)"
    #$w.f.dof add command 7 -label "Global Rescale (7 parameter model)"
    #$w.f.dof add command 9 -label "Traditional (9 parameter model)"
    #$w.f.dof add command 12 -label "Affine (12 parameter model)"
    set reg($w,dof) 12

#}}}

    #{{{ input/high res image
    FileEntry  $w.f.test -textvariable entries($w,2) -label "Input image   " -title "Select" -width 50 -filedialog directory  -filetypes IMAGE

#}}}
    #{{{ low res image
    FileEntry  $w.f.test2 -textvariable entries($w,3) -label "Low res image   " -title "Select" -width 50 -filedialog directory  -filetypes IMAGE

# DOF 2
frame $w.f.doftwo
label $w.f.doftwo.label -text "  Model/DOF (lowres to highres)"
optionMenu2 $w.f.doftwo.menu reg($w,doftwo) 2Dmenu  "   2D to 2D registration" 2D "Rigid Body (3 parameter model)" 3Dmenu "   3D to 3D registration" TRANS "Translation Only (3 parameter model)" 6 "Rigid Body (6 parameter model)" 7 "Global Rescale (7 parameter model)" 9 "Traditional (9 parameter model)" 12 "Affine (12 parameter model)"
pack $w.f.doftwo.label $w.f.doftwo.menu -in $w.f.doftwo -side top -side left




    #tixOptionMenu $w.f.doftwo -label "  Model/DOF (lowres to highres)" -variable reg($w,doftwo) -options { label.anchor w menubutton.width 30 }
    #$w.f.doftwo add command 2Dmenu -label "   2D to 2D registration" -state disabled -background "#555555"
    #$w.f.doftwo add command 2D -label "Rigid Body (3 parameter model)"
    #$w.f.doftwo add command 3Dmenu -label "   3D to 3D registration" -state disabled -background "#555555"
    #$w.f.doftwo add command TRANS -label "Translation Only (3 parameter model)"
    #$w.f.doftwo add command 6 -label "Rigid Body (6 parameter model)"
    #$w.f.doftwo add command 7 -label "Global Rescale (7 parameter model)"
    #$w.f.doftwo add command 9 -label "Traditional (9 parameter model)"
    #$w.f.doftwo add command 12 -label "Affine (12 parameter model)"
    set reg($w,doftwo) 12

#}}}
    #{{{ mode

    set reg($w,mode) 1

frame $w.f.mode
label $w.f.mode.label -text "Mode " 
optionMenu2 $w.f.mode.menu  reg($w,mode)  -command "flirt:updatemode $w" 1 "Input $IG -> Reference image" 2 "Low res $IG -> High res $IG -> Reference image"
pack $w.f.mode.label $w.f.mode.menu -in $w.f.mode -side top -side left
    #tixOptionMenu $w.f.mode -label "Mode " -variable reg($w,mode) -options { label.anchor w }
    #$w.f.mode add command 1 -label "Input $IG -> Reference image"
    #$w.f.mode add command 2 -label "Low res $IG -> High res $IG -> Reference image"

#}}}
    #{{{ output image

    FileEntry  $w.f.output -textvariable entries($w,4) -label "Output image   " -title "Select" -width 50 -filedialog directory  -filetypes IMAGE


#}}}

    pack $w.f.mode $w.f.ref $w.f.dof $w.f.test $w.f.output -in $lfbasic -side top -anchor w -pady $PADY -padx 5

    #{{{ secondary image(s)

set i 1
while { $i <= $reg($w,maxnstats) } {


    FileEntry  $w.f.second$i -textvariable entries($w,[ expr $i + 4 ]) -label "Secondary image $i" -title "Select" -width 50 -filedialog directory  -filetypes IMAGE

    incr i 1
}

#}}}
    pack $w.f.basic $w.f.stats -in $w.f -side top -anchor w -pady 0 -padx 5
    #{{{ advanced options

    # ---- Optional stuff ----

    collapsible frame $w.f.opts -title "Advanced Options"    


NoteBook $w.nb -side top -bd 2 -tabpady {5 10} -arcradius 3
$w.nb insert 0 search  -text "Search"
$w.nb insert 1 cost    -text "Cost Function"
$w.nb insert 2 interp  -text "Interpolation"
$w.nb insert 3 weights -text "Weighting Volumes"
$w.nb raise search
    #{{{ Search

    set lf [$w.nb getframe search]
    set reg($w,searchrxmin) -90
    set reg($w,searchrxmax) 90
    set reg($w,searchrymin) -90
    set reg($w,searchrymax) 90
    set reg($w,searchrzmin) -90
    set reg($w,searchrzmax) 90
    frame $w.search
    frame $w.searchf
    label $w.searchf.angleslabel -text  "Search Angles"
    LabelSpinBox  $w.searchf.rxmin -label "X-axis (degrees): min " -textvariable reg($w,searchrxmin) -range {-180.0 180 1 } 
    LabelSpinBox  $w.searchf.rxmax -label  "  max" -textvariable reg($w,searchrxmax) -range {-180.0 180 1 }  
    frame $w.searchf.rx 
    pack $w.searchf.rxmin $w.searchf.rxmax -in $w.searchf.rx  -side left
    LabelSpinBox  $w.searchf.rymin -label "Y-axis (degrees): min " -textvariable reg($w,searchrymin) -range {-180.0 180 1 } 
    LabelSpinBox  $w.searchf.rymax -label  "  max" -textvariable reg($w,searchrymax) -range {-180.0 180 1 }   
    frame $w.searchf.ry
    pack $w.searchf.rymin $w.searchf.rymax -in $w.searchf.ry  -side left
    LabelSpinBox  $w.searchf.rzmin -label "Z-axis (degrees): min " -textvariable reg($w,searchrzmin) -range {-180.0 180 1 } 
    LabelSpinBox  $w.searchf.rzmax -label  "  max" -textvariable reg($w,searchrzmax) -range {-180.0 180 1 }  
    frame $w.searchf.rz
    pack $w.searchf.rzmin $w.searchf.rzmax -in $w.searchf.rz  -side left

    frame  $w.search.range 
    label  $w.search.range.label -text "Mode " 
    optionMenu2 $w.search.range.menu reg($w,search) -command "flirt:updatesearch $w" 0 "Already virtually aligned (no search)" 1 "Not aligned, but same orientation" 2 "Incorrectly oriented"
    pack $w.search.range.label $w.search.range.menu -in $w.search.range -side top -side left

    set reg($w,search) 1
    set reg($w,disablesearch_yn) 0

    pack $w.searchf.angleslabel $w.searchf.rx $w.searchf.ry $w.searchf.rz -in $w.searchf -side top \
	    -anchor w -padx 3 -pady 3

    pack $w.search.range $w.searchf -in $w.search -side top -anchor w -padx 3 -pady 3

    pack $w.search -in $lf -side top -anchor w



#}}}
    #{{{ Cost Function

    set costlf [$w.nb getframe cost]

    radiobutton $w.corratio -text "Correlation Ratio" \
	    -variable reg($w,cost) -value corratio -anchor w -command "flirt:updatecost $w $costlf"
    radiobutton $w.mutualinfo -text "Mutual Information" \
	    -variable reg($w,cost) -value mutualinfo -anchor w -command "flirt:updatecost $w $costlf"
    radiobutton $w.nmi -text "Normalised Mutual Information" \
	    -variable reg($w,cost) -value normmi -anchor w -command "flirt:updatecost $w $costlf"
    radiobutton $w.normcorr -text "Normalised Correlation (intra-modal)" \
	    -variable reg($w,cost) -value normcorr -anchor w -command "flirt:updatecost $w $costlf"
    radiobutton $w.leastsq -text "Least Squares (intra-modal)" \
	    -variable reg($w,cost) -value leastsq -anchor w -command "flirt:updatecost $w $costlf"
    set reg($w,bins) 256
    LabelSpinBox  $w.bins -label "Number of Histogram Bins " -textvariable reg($w,bins) -range {1 5000 1 } 
    # ---- pack ----
    pack $w.corratio $w.mutualinfo $w.nmi $w.normcorr $w.leastsq $w.bins -in $costlf -side top -anchor w -padx 3
    set reg($w,cost) corratio

#}}}
    #{{{ Interpolation

    set interplf [$w.nb getframe interp]

    label $w.interpbanner -text "Final Interpolation Method (Reslice Only)"
    radiobutton $w.trilinear -text "Tri-Linear" \
	    -variable reg($w,interp) -value trilinear -anchor w -command "flirt:updateinterp $w $interplf"
    radiobutton $w.nearestneighbour -text "Nearest Neighbour" \
	    -variable reg($w,interp) -value nearestneighbour -anchor w -command "flirt:updateinterp $w $interplf"
    radiobutton $w.sinc -text "Sinc" \
	    -variable reg($w,interp) -value sinc -anchor w -command "flirt:updateinterp $w $interplf"
    set reg($w,sincwidth) 7
    LabelSpinBox $w.sincwidth -label " Width of Sinc Window (full width - voxels)" -textvariable reg($w,sincwidth) -range {1 5000 1 } 
    frame $w.swinopt
    label $w.swinbanner -text "Sinc Window Options"
    radiobutton $w.rectangular -text "Rectangular" \
	    -variable reg($w,sincwindow) -value rectangular -anchor w
    radiobutton $w.hanning -text "Hanning" \
	    -variable reg($w,sincwindow) -value hanning -anchor w
    radiobutton $w.blackman -text "Blackman" \
	    -variable reg($w,sincwindow) -value blackman -anchor w
    set reg($w,sincwindow) hanning
    
    # ---- pack ----
    pack $w.interpbanner $w.trilinear -in $interplf -side top -anchor w -padx 3
    pack $w.nearestneighbour $w.sinc -in $interplf -side top -anchor w -padx 3
    set reg($w,interp) trilinear

    pack $w.swinbanner -in $w.swinopt -side top -anchor w -padx 3
    pack $w.rectangular $w.hanning $w.blackman -in $w.swinopt -side left -anchor w -padx 3

#}}}
    #{{{ Weightings

    set weightlf [$w.nb getframe weights]

    set entries($w,35) ""

    FileEntry  $w.wgt -textvariable entries($w,35) -label "Reference weighting   " -title "Select" -width 50 -filedialog directory  -filetypes IMAGE

    set entries($w,36) ""

    FileEntry  $w.iwgt -textvariable entries($w,36) -label "Input weighting   " -title "Select" -width 50 -filedialog directory  -filetypes IMAGE

    set entries($w,37) ""

    FileEntry  $w.iwgt2 -textvariable entries($w,37) -label  "Low res weighting   " -title "Select" -width 50 -filedialog directory  -filetypes IMAGE

    # ---- pack ----
    pack $w.wgt $w.iwgt -in $weightlf -side top -anchor w -padx 3 -pady $PADY

#}}}

    frame $w.f.advopts
    pack $w.nb -in $w.f.advopts -side top
    pack $w.f.advopts -in $w.f.opts.b -side left -padx 8 -pady 6 -expand yes -fill both
    pack $w.f.opts -in $w.f -side left -padx 5 -pady 5

#}}}
    #{{{ buttons

    # ---- Button Frame ----

    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
    
    button $w.apply     -command "flirt:apply $w keep" \
	    -text "Go" -width 5
    bind $w.apply <Return> {
	[winfo toplevel %W].apply invoke
    }

    button $w.cancel    -command "flirt:destroy $w" \
	    -text "Exit" -width 5
    bind $w.cancel <Return> {
	[winfo toplevel %W].cancel invoke
    }

    button $w.help -command "FmribWebHelp file: ${FSLDIR}/doc/flirt/index.html" \
	    -text "Help" -width 5

    #{{{ Utils

menubutton $w.utils -text "Utils" -menu $w.utils.menu -relief raised -bd 2

menu $w.utils.menu

$w.utils.menu add command -label "Apply FLIRT transform" -command { exec sh -c "${FSLDIR}/bin/ApplyXFM$gui_ext" & }
$w.utils.menu add command -label "Concat FLIRT transforms" -command { exec sh -c "${FSLDIR}/bin/ConcatXFM$gui_ext" & }
$w.utils.menu add command -label "Invert FLIRT transform" -command { exec sh -c "${FSLDIR}/bin/InvertXFM$gui_ext" & }

#}}}

    pack $w.btns.b -side bottom -fill x
    pack $w.apply $w.cancel $w.help $w.utils -in $w.btns.b \
	    -side left -expand yes -padx 3 -pady 10 -fill y
    
    pack $w.f $w.btns -expand yes -fill both

#}}}
}

#}}}
#{{{ flirt:apply

proc flirt:apply { w dialog } {
    global reg entries

    if { $reg($w,nstats) == 0 } {
	set statslist "grot"
    } else {
	set statslist ""
    }
    set i 1
    while { $i <= $reg($w,nstats) } {
	lappend statslist $entries($w,[ expr $i + 4 ])
	incr i 1
    }

    set status [ flirt:proc $reg($w,mode) $entries($w,1) $entries($w,2) $entries($w,3) $reg($w,nstats) $statslist $entries($w,4) $reg($w,dof) $reg($w,doftwo) $reg($w,bins) $reg($w,searchrxmin) $reg($w,searchrxmax) $reg($w,searchrymin) $reg($w,searchrymax) $reg($w,searchrzmin) $reg($w,searchrzmax) $reg($w,disablesearch_yn) $reg($w,cost) $reg($w,interp) $reg($w,sincwidth) $reg($w,sincwindow) $entries($w,35) $entries($w,36) $entries($w,37) 1 ]

    update idletasks
    
    # destroy if the OK button was used AND a normal exit occurred
    if { $dialog == "destroy" && $status == 0 } {
	flirt:destroy $w
    }
}

#}}}
#{{{ flirt:destroy

proc flirt:destroy { w } {
    destroy $w
}

#}}}
#{{{ flirt:updatemode

proc flirt:updatemode { w } {
    global reg PADY IGS

    if { $reg($w,mode) == 1 } {
	$w.f.dof.label configure -text "  Model/DOF (input to ref)"
	$w.f.test.frame.label configure -text "Input image"
	$w.iwgt.frame.label configure -text "Input weighting"
	$w.f.nstats configure -label "Number of secondary $IGS to apply transform to "
	pack forget $w.f.test2
	pack forget $w.f.doftwo
	pack forget $w.iwgt2
#	pack $w.wgt -in [$w.nb subwidget weights] -side top -anchor w -padx 3
#	pack $w.iwgt -in [$w.nb subwidget weights] -side top -anchor w -padx 3 -pady $PADY
    } else {
	$w.f.dof.label configure -text "  Model/DOF (highres to ref)"
	$w.f.test.frame.label configure -text "High res image"
	$w.iwgt.frame.label configure -text "High res weighting"
	$w.f.nstats configure -label "Number of secondary $IGS to apply combined transform to "
	pack $w.f.doftwo $w.f.test2 -in $w.f -side top -anchor w -pady $PADY -padx 5 -after $w.f.test
	pack $w.iwgt2 -in [$w.nb getframe weights] -side top -anchor w -padx 3 -pady $PADY
#	pack forget $w.wgt
#	pack forget $w.iwgt
    }
}

#}}}
#{{{ flirt:updatestats

proc flirt:updatestats { w } {
    global reg PADY

    set i 1
    while { $i <= $reg($w,maxnstats) } {
	pack forget $w.f.second$i
	incr i 1
    }

    if { $reg($w,nstats) > 0 } {
	pack $w.f.second1 -in $w.f -side top -anchor w -pady $PADY -padx 5 -after $w.f.nstats
	$w.f.second1 configure -label "Secondary image"
    }
    set i 2
    while { $i <= $reg($w,nstats) } {
	$w.f.second1 configure -label "Secondary image 1"
	pack $w.f.second$i -in $w.f -side top -anchor w -pady $PADY -padx 5 -after $w.f.second[ expr $i - 1 ]
	incr i 1
    }
}

#}}}
#{{{ flirt:updatesearch

proc flirt:updatesearch { w } {
    global reg

    if { $reg($w,search) == 0 } {
	set reg($w,disablesearch_yn) 1
    } else {
	set reg($w,disablesearch_yn) 0

	if { $reg($w,search) == 1 } {
	    set reg($w,searchrxmin) -90
	    set reg($w,searchrxmax) 90
	    set reg($w,searchrymin) -90
	    set reg($w,searchrymax) 90
	    set reg($w,searchrzmin) -90
	    set reg($w,searchrzmax) 90
	} else {
	    set reg($w,searchrxmin) -180
	    set reg($w,searchrxmax) 180
	    set reg($w,searchrymin) -180
	    set reg($w,searchrymax) 180
	    set reg($w,searchrzmin) -180
	    set reg($w,searchrzmax) 180
	}
    }

    if { $reg($w,disablesearch_yn) } {
	pack forget $w.searchf
    } else {
	pack $w.searchf -in $w.search -side top -anchor w -padx 3 -pady 3
    }
}

#}}}
#{{{ flirt:updatecost

proc flirt:updatecost { w costlf } {
    global reg

    if { [ string match $reg($w,cost) "normcorr" ] == 1  || 
         [ string match $reg($w,cost) "leastsq" ] == 1 } {
	pack forget $w.bins
    } else {
	pack $w.bins -in $costlf -side top -anchor w -padx 3
    }
}

#}}}
#{{{ flirt:updateinterp

proc flirt:updateinterp { w interplf } {
    global reg

    if { [ string match $reg($w,interp) "sinc" ] == 1 } {
	pack $w.swinopt -in $interplf -side top -anchor w -padx 40
	pack $w.sincwidth -in $interplf -side top -anchor w -padx 40
    } else {
	pack forget $w.swinopt
	pack forget $w.sincwidth
    }
}

#}}}

#{{{ start TK window

wm withdraw .
flirt .rename
tkwait window .rename

#}}}

