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
#{{{ invertxfm

proc invertxfm { w } {

    #{{{ setup main window etc

    global entries FSLDIR PWD

    # ---- Set up Frames ----
    toplevel $w
    wm title $w "InvertXFM"
    wm iconname $w "InvertXFM"
    wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm
    tixBalloon    $w.bhelp
    frame $w.f

#}}}

    tixLabelFrame $w.f.input -label "Input"
    set lfinput [ $w.f.input subwidget frame ]
    #{{{ input transform

set entries($w,1) ""

FSLFileEntry $w.f.xfm \
	-variable entries($w,1) \
	-pattern "*.mat" \
	-directory $PWD \
	-label "Transformation File for A to B   " \
		-labelwidth 29 \
		-title "Select" \
		-width 40 \
		-filterhist VARS(history)

#}}}
    #{{{ A image

set entries($w,2) ""

FSLFileEntry $w.f.ref \
	-variable entries($w,2) \
	-pattern "*.hdr" \
	-directory $PWD \
	-label "Volume A   " \
	-labelwidth 29 \
	-title "Select" \
	-width 40 \
	-filterhist VARS(history)

#}}}
    #{{{ B image

FSLFileEntry $w.f.test \
	-variable entries($w,3) \
	-pattern "*.hdr" \
	-directory $PWD \
	-label "Volume B   " \
	-labelwidth 29 \
	-title "Select" \
	-width 40 \
	-filterhist VARS(history)

#}}}
    pack $w.f.xfm $w.f.ref $w.f.test -in $lfinput -side top -anchor w -pady 3 -padx 5

    tixLabelFrame $w.f.output -label "Output"
    set lfoutput [ $w.f.output subwidget frame ]
    #{{{ output filename

set entries($w,4) ""

FSLFileEntry $w.f.oxfm \
	-variable entries($w,4) \
	-pattern "*.mat" \
	-directory $PWD \
	-label "Save Inverse Transform (B to A)" \
	-labelwidth 29 \
		-title "Select" \
		-width 40 \
		-filterhist VARS(history)

#}}}
    #{{{ Type of transformation to save

set entries($w,5) n
    
frame $w.f.type
label $w.f.typebanner -text "Type of transform to save"

radiobutton $w.f.a -text "MEDx AlignLinearReslice" -variable entries($w,5) -value a -anchor w
radiobutton $w.f.u -text "MEDx UserTransformation" -variable entries($w,5) -value u -anchor w
radiobutton $w.f.t -text "MEDx IntoTalairachSpace" -variable entries($w,5) -value t -anchor w
radiobutton $w.f.n -text "FLIRT matrix"            -variable entries($w,5) -value n -anchor w

pack $w.f.typebanner -in $w.f.type -side top -anchor w

frame $w.f.type.a
frame $w.f.type.b

pack $w.f.u $w.f.a -in $w.f.type.a -side left
pack $w.f.t $w.f.n -in $w.f.type.b -side left

pack $w.f.type.a $w.f.type.b -in $w.f.type -side top -anchor w

#}}}
    pack $w.f.oxfm $w.f.type -in $lfoutput -side top -anchor w -pady 3 -padx 5

    pack $w.f.input $w.f.output -in $w.f -side top -anchor w -pady 0 -padx 5

    #{{{ buttons

    # ---- Button Frame ----

    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
    
    button $w.go     -command "InvertXFM:go $w" \
	    -text "Go" -width 5

    button $w.apply     -command "InvertXFM:apply $w" \
	    -text "Apply" -width 5

    button $w.cancel    -command "destroy $w" \
	    -text "Exit" -width 5

    button $w.help -command "FmribWebHelp file: ${FSLDIR}/doc/flirt/overview.html" \
	    -text "Help" -width 5

    pack $w.btns.b -side bottom -fill x
    pack $w.go $w.apply $w.cancel $w.help -in $w.btns.b \
	    -side left -expand yes -padx 3 -pady 10 -fill y
    
    pack $w.f $w.btns -expand yes -fill both

#}}}
}

#}}}
#{{{ InvertXFM:apply

proc InvertXFM:go { w } {
    global entries

    catch { InvertXFM:apply $w }
    destroy $w
}

proc InvertXFM:apply { w } {
    global entries

    catch { invertxfm:proc $entries($w,1) $entries($w,2) $entries($w,3) $entries($w,4) $entries($w,5) }

    update idletasks
    puts "Done"
}

#}}}
#{{{ invertxfm:proc

proc invertxfm:proc { transAB volA volB invxfmfilename xfmtype } {

    global FSLDIR

    if { $xfmtype == "n" } {
	set savemedx "-omat $invxfmfilename"
    } else {
	set savemedx "-omedx $invxfmfilename -xfmtype $xfmtype"
    }
    
    set thecommand "${FSLDIR}/bin/convert_xfm -ref $volB -in $volA $savemedx -inverse $transAB"
    puts $thecommand
    catch { exec sh -c $thecommand } errmsg
    puts $errmsg

    if { [ file readable $invxfmfilename ] == 0 } {
	puts "No transformation saved!"
	return 4
    }

    return 0
}

#}}}
#{{{ call GUI

wm withdraw .
invertxfm .rename
tkwait window .rename

#}}}
