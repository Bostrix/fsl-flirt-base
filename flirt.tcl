#

# FLIRT - FMRIB's Linear Image Registration Tool
#
# Mark Jenkinson and Stephen Smith, FMRIB Image Analysis Group
#
# Copyright (C) 1999-2000 University of Oxford
#
# TCLCOPYRIGHT

#{{{ setup

source [ file dirname [ info script ] ]/fslstart.tcl

set VARS(history) {}

#}}}

proc flirt { w } {

    #{{{ setup main window etc

    global reg entries USER FSLDIR INMEDX argc argv PWD PADY IGS

    set reg($w,maxnstats) 20

    # ---- Set up Frames ----
    toplevel $w
    wm title $w "FLIRT - FMRIB's Linear Image Registration Tool - v2.2"
    wm iconname $w "FLIRT"
    wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm
    tixBalloon    $w.bhelp
    frame $w.f

    tixLabelFrame $w.f.basic
    set lfbasic [ $w.f.basic subwidget frame ]

    if { $INMEDX } {
	set PADY 0
	set IG  "image/group"
	set IGS "images/groups"
    } else {
	set PADY 3
	set IG  "image"
	set IGS "images"
    }

#}}}
    #{{{ number of secondary images

    tixLabelFrame $w.f.stats
    set lfstats [ $w.f.stats subwidget frame ]
    tixControl $w.f.nstats -label "Number of secondary $IGS to apply transform to " \
	-variable reg($w,nstats) -step 1 -min 0 -max $reg($w,maxnstats) -command "flirt:updatestats $w"
    set reg($w,nstats) 0
    pack $w.f.nstats -in $lfstats -side top -anchor w -pady 3 -padx 5

#}}}
    #{{{ standard image

set entries($w,1) ${FSLDIR}/etc/standard/avg152T1_brain.hdr

    if { $INMEDX } {
        frame $w.f.ref
        label $w.f.reftxt -width 23 -text "Reference image"
        entry $w.f.refvar -textvariable entries($w,1) -width 30
        button $w.f.refsel -text "Select" -command "SelectPage:Dialog $w 1 0 30 entries"
        pack $w.f.reftxt $w.f.refvar $w.f.refsel -in $w.f.ref -padx 3 -pady 0 -side left
    } else {
	FSLFileEntry $w.f.ref \
		-variable entries($w,1) \
		-pattern "*.hdr" \
		-directory $PWD \
		-label "Reference image   " \
		-labelwidth 18 \
		-title "Select" \
		-width 50 \
		-filterhist VARS(history)
    }

#}}}
    #{{{ input/high res image

    if { $INMEDX } {
        frame $w.f.test
        label $w.f.testtxt -width 23 -text "Input image/group"
        entry $w.f.testvar -textvariable entries($w,2) -width 30
        button $w.f.testsel -text "Select" -command "SelectPage:Dialog $w 2 0 30 entries"
        pack $w.f.testtxt $w.f.testvar $w.f.testsel -in $w.f.test -padx 3 -pady 0 -side left
    } else {
	FSLFileEntry $w.f.test \
		-variable entries($w,2) \
		-pattern "*.hdr" \
		-directory $PWD \
		-label "Input image   " \
		-labelwidth 18 \
		-title "Select" \
		-width 50 \
		-filterhist VARS(history)
    }

#}}}
    #{{{ low res image

    if { $INMEDX } {
        frame $w.f.test2
        label $w.f.testtxt2 -width 23 -text "Low res image/group"
        entry $w.f.testvar2 -textvariable entries($w,3) -width 30
        button $w.f.testsel2 -text "Select" -command "SelectPage:Dialog $w 3 0 30 entries"
        pack $w.f.testtxt2 $w.f.testvar2 $w.f.testsel2 -in $w.f.test2 -padx 3 -pady 0 -side left
    } else {
	FSLFileEntry $w.f.test2 \
		-variable entries($w,3) \
		-pattern "*.hdr" \
		-directory $PWD \
		-label "Low res image   " \
		-labelwidth 18 \
		-title "Select" \
		-width 50 \
		-filterhist VARS(history)
    }

#}}}
    #{{{ mode

    set reg($w,mode) 1
    tixOptionMenu $w.f.mode -label "Mode " \
	    -variable reg($w,mode) \
	    -command "flirt:updatemode $w" \
	    -options {
	label.anchor w
    }
    $w.f.mode add command 1 -label "Input $IG -> Reference image"
    $w.f.mode add command 2 -label "Low res $IG -> High res $IG -> Reference image"

#}}}
    pack $w.f.mode $w.f.ref $w.f.test -in $lfbasic -side top -anchor w -pady $PADY -padx 5
    #{{{ output image

    if { ! $INMEDX } {
	FSLFileEntry $w.f.output \
		-variable entries($w,4) \
		-pattern "*.hdr" \
		-directory $PWD \
		-label "Output image   " \
		-labelwidth 18 \
		-title "Select" \
		-width 50 \
		-filterhist VARS(history)
	pack $w.f.output -in $lfbasic -side top -anchor w -pady $PADY -padx 5
    } else {
	set entries($w,4) /tmp/grot
    }

#}}}
    #{{{ secondary image(s)

set i 1
while { $i <= $reg($w,maxnstats) } {

    if { $INMEDX } {
	frame $w.f.second$i
	label $w.f.secondtxt($i) -width 23 -text "Secondary image/group $i"
	entry $w.f.secondvar($i) -textvariable entries($w,[ expr $i + 4 ]) -width 30
        button $w.f.secondsel($i) -text "Select" -command "SelectPage:Dialog $w [ expr $i + 4 ] 0 2 entries"
	pack $w.f.secondtxt($i) $w.f.secondvar($i) $w.f.secondsel($i) -in $w.f.second$i -padx 0 -pady 0 -side left
    } else {
	FSLFileEntry $w.f.second$i \
		-variable entries($w,[ expr $i + 4 ]) \
		-pattern "*.hdr" \
		-directory $PWD \
		-label "Secondary image $i" \
		-labelwidth 18 \
		-title "Select" \
		-width 50 \
		-filterhist VARS(history)
    }

    incr i 1
}

#}}}
    #{{{ DOF

    tixLabelFrame $w.f.doff
    set lfdoff [ $w.f.doff subwidget frame ]

    tixOptionMenu $w.f.dof -label "Model/DOF " \
	    -variable reg($w,dof) \
	    -options {
	label.anchor w
	menubutton.width 30
    }
    $w.f.dof add command 3D -label "         3D Models" -state disabled
    $w.f.dof add command 6 -label "Rigid Body (6 parameter model)"
    $w.f.dof add command 7 -label "Global Rescale (7 parameter model)"
    $w.f.dof add command 9 -label "Traditional (9 parameter model)"
    $w.f.dof add command 12 -label "Affine (12 parameter model)"
    set reg($w,dof) 12

#}}}
    pack $w.f.dof -in $lfdoff -side top -anchor w -padx 5 -pady 1
    pack $w.f.basic $w.f.stats $w.f.doff -in $w.f -side top -anchor w -pady 0 -padx 5
    #{{{ advanced options

    # ---- Optional stuff ----

    collapsible frame $w.f.opts -title "Advanced Options"    

    tixNoteBook $w.nb -ipadx 5 -ipady 5

    $w.nb add search -label "Search"
    $w.nb add anglerep -label "Angle Rep."
    $w.nb add cost -label "Cost Function"
    
    # ---- Search ----
    set lf [$w.nb subwidget search]

    frame $w.search
    frame $w.searchf
    label $w.searchf.angleslabel -text  "Search Angles"

    frame $w.searchf.rx
    tixControl $w.searchf.rxmin -label "X-axis (degrees): min " \
	    -variable reg($w,searchrxmin) -step 1 -min -180 -max 180 \
	    -selectmode immediate
    
    tixControl $w.searchf.rxmax -label "  max" -variable reg($w,searchrxmax) \
	    -step 1 -min -180 -max 180 -selectmode immediate
    pack $w.searchf.rxmin $w.searchf.rxmax -in $w.searchf.rx  -side left

    frame $w.searchf.ry
    tixControl $w.searchf.rymin -label "Y-axis (degrees): min " \
	    -variable reg($w,searchrymin) -step 1 -min -180 -max 180 \
	    -selectmode immediate
    tixControl $w.searchf.rymax -label "  max" -variable reg($w,searchrymax) \
	    -step 1 -min -180 -max 180 -selectmode immediate
    pack $w.searchf.rymin $w.searchf.rymax -in $w.searchf.ry  -side left


    frame $w.searchf.rz
    tixControl $w.searchf.rzmin -label "Z-axis (degrees): min " \
	    -variable reg($w,searchrzmin) -step 1 -min -180 -max 180 \
	    -selectmode immediate
    tixControl $w.searchf.rzmax -label "  max" -variable reg($w,searchrzmax) \
	    -step 1 -min -180 -max 180 -selectmode immediate
    pack $w.searchf.rzmin $w.searchf.rzmax -in $w.searchf.rz  -side left

    global reg($w,search)
    tixOptionMenu $w.search.range -label "Images: " -variable reg($w,search) -command "flirt:updatesearch $w"
    $w.search.range add command 0 -label "Already virtually aligned (no search)"
    $w.search.range add command 1 -label "Not aligned, but same orientation"
    $w.search.range add command 2 -label "Incorrectly oriented"
    set reg($w,search) 1
    set reg($w,disablesearch_yn) 0

    pack $w.searchf.angleslabel $w.searchf.rx $w.searchf.ry $w.searchf.rz -in $w.searchf -side top \
	    -anchor w -padx 3 -pady 3

    pack $w.search.range $w.searchf -in $w.search -side top -anchor w -padx 3 -pady 3

    pack $w.search -in $lf -side top -anchor w

    set reg($w,searchrxmin) -90
    set reg($w,searchrxmax) 90
    set reg($w,searchrymin) -90
    set reg($w,searchrymax) 90
    set reg($w,searchrzmin) -90
    set reg($w,searchrzmax) 90


    # ---- Angle Representation ----
    set lf [$w.nb subwidget anglerep]

    label $w.anglebanner -text "Internal Angle Representation"
    radiobutton $w.euler -text "Euler Angles" \
	    -variable reg($w,anglerep) -value euler

    radiobutton $w.quaternion -text "Quaternions" \
	    -variable reg($w,anglerep) -value quaternion


    # ---- pack ----
    pack $w.anglebanner $w.euler $w.quaternion -in $lf -side top -padx 3 -pady 3
    set reg($w,anglerep) euler



    # ---- Cost Function ----
    set costlf [$w.nb subwidget cost]

    radiobutton $w.corratio -text "Correlation Ratio" \
	    -variable reg($w,cost) -value corratio -anchor w -command "flirt:updatecost $w $costlf"
    radiobutton $w.mutualinfo -text "Mutual Information" \
	    -variable reg($w,cost) -value mutualinfo -anchor w -command "flirt:updatecost $w $costlf"
    radiobutton $w.woods -text "Woods Function" \
	    -variable reg($w,cost) -value woods -anchor w -command "flirt:updatecost $w $costlf"
    radiobutton $w.normcorr -text "Normalised Correlation" \
	    -variable reg($w,cost) -value normcorr -anchor w -command "flirt:updatecost $w $costlf"

    tixControl $w.bins -label "Number of Histogram Bins " \
	    -variable reg($w,bins) -step 1 -min 1 -max 5000 -selectmode immediate
    set reg($w,bins) 256
    
    # ---- pack ----
    pack $w.corratio $w.mutualinfo $w.woods $w.normcorr $w.bins -in $costlf -side top -anchor w -padx 3
    set reg($w,cost) corratio

    frame $w.f.advopts
    pack $w.nb -in $w.f.advopts -side top
    pack $w.f.advopts -in $w.f.opts.b -side left -padx 8 -pady 6 -expand yes -fill both


    pack $w.f.opts -in $w.f -side left -padx 5 -pady 5

#}}}
    #{{{ buttons

    # ---- Button Frame ----

    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
    
    button $w.ok \
	    -text "OK" -width 5 \
	    -command "flirt:apply $w destroy"
    bind $w.ok <Return> {
        [winfo toplevel %W].ok invoke
    }
    
    button $w.apply     -command "flirt:apply $w keep" \
	    -text "Apply" -width 5
    bind $w.apply <Return> {
	[winfo toplevel %W].apply invoke
    }

    button $w.cancel    -command "flirt:destroy $w" \
	    -text "Cancel" -width 5
    bind $w.cancel <Return> {
	[winfo toplevel %W].cancel invoke
    }

    button $w.help -command "FmribWebHelp file: ${FSLDIR}/doc/flirt/index.html" \
	    -text "Help" -width 5

    pack $w.btns.b -side bottom -fill x
    pack $w.ok $w.apply $w.cancel $w.help -in $w.btns.b \
	    -side left -expand yes -padx 3 -pady 10 -fill y
    
    pack $w.f $w.btns -expand yes -fill both

#}}}
}

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

    set status [ flirt_proc $reg($w,mode) $entries($w,1) $entries($w,2) $entries($w,3) $reg($w,nstats) $statslist $entries($w,4) $reg($w,dof) $reg($w,bins) $reg($w,searchrxmin) $reg($w,searchrxmax) $reg($w,searchrymin) $reg($w,searchrymax) $reg($w,searchrzmin) $reg($w,searchrzmax) $reg($w,disablesearch_yn) $reg($w,anglerep) $reg($w,cost) 1 ]

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

proc flirt:updatemode { w junk } {
    global reg INMEDX PADY IGS

    if { $reg($w,mode) == 1 } {
	if { $INMEDX } {
	    $w.f.testtxt configure -text "Input image/group"
	} else {
	    $w.f.test.frame.label configure -text "Input image"
	}
	$w.f.nstats configure -label "Number of secondary $IGS to apply transform to "
	pack forget $w.f.test2
    } else {
	if { $INMEDX } {
	    $w.f.testtxt configure -text "High res image/group"
	} else {
	    $w.f.test.frame.label configure -text "High res image"
	}
	$w.f.nstats configure -label "Number of secondary $IGS to apply combined transform to "
	pack $w.f.test2 -in $w.f -side top -anchor w -pady $PADY -padx 5 -after $w.f.test
    }
}

#}}}
#{{{ flirt:updatestats

proc flirt:updatestats { w junk } {
    global reg INMEDX PADY

    set i 1
    while { $i <= $reg($w,maxnstats) } {
	pack forget $w.f.second$i
	incr i 1
    }

    if { $reg($w,nstats) > 0 } {
	pack $w.f.second1 -in $w.f -side top -anchor w -pady $PADY -padx 5 -after $w.f.nstats
	if { $INMEDX } {
	    $w.f.secondtxt(1) configure -text "Secondary image/group"
	} else {
	    $w.f.second1.frame.label configure -text "Secondary image"
	}
    }
    set i 2
    while { $i <= $reg($w,nstats) } {
	if { $INMEDX } {
	    $w.f.secondtxt(1) configure -text "Secondary image/group 1"
	} else {
	    $w.f.second1.frame.label configure -text "Secondary image 1"
	}
	pack $w.f.second$i -in $w.f -side top -anchor w -pady $PADY -padx 5 -after $w.f.second[ expr $i - 1 ]
	incr i 1
    }
}

#}}}
#{{{ flirt:updatesearch

proc flirt:updatesearch { w dummy } {
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

    if { [ string match $reg($w,cost) "normcorr" ] == 1 } {
	pack forget $w.bins
    } else {
	pack $w.bins -in $costlf -side top -anchor w -padx 3
    }
}

#}}}

wm withdraw .
flirt .rename
tkwait window .rename
