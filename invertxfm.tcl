#

# InvertXFM - the GUI for convert_xfm
#
# Mark Jenkinson and Stephen Smith, FMRIB Image Analysis Group
#
# Copyright (C) 2001 University of Oxford
#
# TCLCOPYRIGHT

#{{{ setup

source [ file dirname [ info script ] ]/fslstart.tcl

set VARS(history) {}

#}}}

proc invertxfm { w } {

    #{{{ setup main window etc

    global reg entries USER FSLDIR INMEDX argc argv PWD PADY IGS

    # ---- Set up Frames ----
    toplevel $w
    wm title $w "InvertXFM - front end for convert_xfm"
    wm iconname $w "InvertXFM"
    wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm
    tixBalloon    $w.bhelp
    frame $w.f

    tixLabelFrame $w.f.basic
    set lfbasic [ $w.f.basic subwidget frame ]

    if { $INMEDX } {
	set PADY 0
    } else {
	set PADY 3
    }

#}}}

    #{{{ input transform

    set entries($w,1) ""

    if { $INMEDX } {
        frame $w.f.xfm
        label $w.f.xfmtxt -width 23 -text "Transformation File for A to B"
        entry $w.f.xfmvar -textvariable entries($w,1) -width 30
        button $w.f.xfmsel -text "Select" -command "SelectPage:Dialog $w 1 0 30 entries"
        pack $w.f.xfmtxt $w.f.xfmvar $w.f.xfmsel -in $w.f.xfm -padx 3 -pady 0 -side left
    } else {
	FSLFileEntry $w.f.xfm \
		-variable entries($w,1) \
		-pattern "*.hdr" \
		-directory $PWD \
		-label "Transformation File for A to B   " \
		-labelwidth 27 \
		-title "Select" \
		-width 41 \
		-filterhist VARS(history)
    }

#}}}

    #{{{ standard image

    set entries($w,2) ""

    if { $INMEDX } {
        frame $w.f.ref
        label $w.f.reftxt -width 23 -text "Volume A"
        entry $w.f.refvar -textvariable entries($w,2) -width 30
        button $w.f.refsel -text "Select" -command "SelectPage:Dialog $w 1 0 30 entries"
        pack $w.f.reftxt $w.f.refvar $w.f.refsel -in $w.f.ref -padx 3 -pady 0 -side left
    } else {
	FSLFileEntry $w.f.ref \
		-variable entries($w,2) \
		-pattern "*.hdr" \
		-directory $PWD \
		-label "Volume A   " \
		-labelwidth 18 \
		-title "Select" \
		-width 50 \
		-filterhist VARS(history)
    }

#}}}
    #{{{ input/high res image

    if { $INMEDX } {
        frame $w.f.test
        label $w.f.testtxt -width 23 -text "Volume B"
        entry $w.f.testvar -textvariable entries($w,3) -width 30
        button $w.f.testsel -text "Select" -command "SelectPage:Dialog $w 2 0 30 entries"
        pack $w.f.testtxt $w.f.testvar $w.f.testsel -in $w.f.test -padx 3 -pady 0 -side left
    } else {
	FSLFileEntry $w.f.test \
		-variable entries($w,3) \
		-pattern "*.hdr" \
		-directory $PWD \
		-label "Volume B   " \
		-labelwidth 18 \
		-title "Select" \
		-width 50 \
		-filterhist VARS(history)
    }

#}}}

    #{{{ input transform

    set entries($w,4) ""

    tixLabelFrame $w.f.outxfm
    set lfoxfm [ $w.f.outxfm subwidget frame ]

    label $w.oxfmbanner -text "Save Inverse Transform (B to A)"

    if { $INMEDX } {
        frame $w.f.oxfm
        label $w.f.oxfmtxt -width 23 -text "Filename"
        entry $w.f.oxfmvar -textvariable entries($w,4) -width 30
        button $w.f.oxfmsel -text "Select" -command "SelectPage:Dialog $w 1 0 30 entries"
        pack $w.f.oxfmtxt $w.f.oxfmvar $w.f.oxfmsel -in $w.f.oxfm -padx 3 -pady 0 -side left
    } else {
	FSLFileEntry $w.f.oxfm \
		-variable entries($w,4) \
		-pattern "*.hdr" \
		-directory $PWD \
		-label "Filename   " \
		-labelwidth 27 \
		-title "Select" \
		-width 41 \
		-filterhist VARS(history)
    }

    pack $w.oxfmbanner -in $lfoxfm -padx 3 -pady $PADY -side top -anchor w
    pack $w.f.oxfm -in $lfoxfm -padx 3 -pady $PADY -side left -anchor w
#}}}



    # ---- packing ---- #

    pack $w.f.xfm $w.f.ref $w.f.test $w.f.outxfm -in $lfbasic -side top -anchor w -pady $PADY -padx 5
    pack $w.f.basic -in $w.f -side top -anchor w -pady 0 -padx 5

    #{{{ buttons

    # ---- Button Frame ----

    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
    
    button $w.ok \
	    -text "OK" -width 5 \
	    -command "InvertXFM:apply $w destroy"
    bind $w.ok <Return> {
        [winfo toplevel %W].ok invoke
    }
    
    button $w.apply     -command "InvertXFM:apply $w keep" \
	    -text "Apply" -width 5
    bind $w.apply <Return> {
	[winfo toplevel %W].apply invoke
    }

    button $w.cancel    -command "InvertXFM:destroy $w" \
	    -text "Cancel" -width 5
    bind $w.cancel <Return> {
	[winfo toplevel %W].cancel invoke
    }

    button $w.help -command "FmribWebHelp file: ${FSLDIR}/doc/flirt/overview.html" \
	    -text "Help" -width 5

    pack $w.btns.b -side bottom -fill x
    pack $w.ok $w.apply $w.cancel $w.help -in $w.btns.b \
	    -side left -expand yes -padx 3 -pady 10 -fill y
    
    pack $w.f $w.btns -expand yes -fill both

#}}}
}

#{{{ InvertXFM:apply

proc InvertXFM:apply { w dialog } {
    global reg entries

    set status [ invertxfm_proc $entries($w,1) $entries($w,2) $entries($w,3) $entries($w,4) 1 ]

    update idletasks
    
    # destroy if the OK button was used AND a normal exit occurred
    if { $dialog == "destroy" && $status == 0 } {
	InvertXFM:destroy $w
    }
}

#}}}
#{{{ InvertXFM:destroy

proc InvertXFM:destroy { w } {
    destroy $w
}

#}}}




wm withdraw .
invertxfm .rename
tkwait window .rename

