#

# Concatxfm - the GUI for convert_xfm
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
#{{{ concatxfm

proc concatxfm { w } {

    #{{{ setup main window etc

    global entries FSLDIR PWD

    # ---- Set up Frames ----
    toplevel $w
    wm title $w "Concatxfm"
    wm iconname $w "Concatxfm"
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
    pack $w.f.xfm -in $lfinput -side top -anchor w -pady 3 -padx 5

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
    pack $w.f.oxfm -in $lfoutput -side top -anchor w -pady 3 -padx 5

    pack $w.f.input $w.f.output -in $w.f -side top -anchor w -pady 0 -padx 5

    #{{{ buttons

    # ---- Button Frame ----

    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
    
    button $w.go     -command "Concatxfm:go $w" \
	    -text "OK" -width 5

    button $w.apply     -command "Concatxfm:apply $w" \
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
#{{{ Concatxfm:apply

proc Concatxfm:go { w } {
    global entries

    catch { Concatxfm:apply $w }
    destroy $w
}

proc Concatxfm:apply { w } {
    global entries

    catch { concatxfm:proc $entries($w,1) $entries($w,4) }

    update idletasks
    puts "Done"
}

#}}}
#{{{ concatxfm:proc

proc concatxfm:proc { transAB invxfmfilename } {

    global FSLDIR

    set thecommand "${FSLDIR}/bin/convert_xfm -matonly -omat $invxfmfilename -inverse $transAB"
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
concatxfm .rename
tkwait window .rename

#}}}
