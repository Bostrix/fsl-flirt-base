#

# InvertXFM - the GUI for convert_xfm

# Mark Jenkinson and Stephen Smith, FMRIB Image Analysis Group

# Copyright (C) 2001 University of Oxford

# TCLCOPYRIGHT


proc invertxfm_proc { transAB volA volB transBA popups } {

    global PXHOME FSLDIR USER MEDXV HOME INMEDX CLUSTERRSH

    invertxfm_run $transAB $volA $volB $transBA u

    return 0

}



proc invertxfm_run { transAB volA volB invxfmfilename xfmtype } {

    # setup vars
    global FSLDIR MEDXV CLUSTERRSH INMEDX
    
    if { "X${xfmtype}X" == "XX" } {
	set xfmtype u
    }
    
    # set up the command to be executed
    set thecommand "$CLUSTERRSH ${FSLDIR}/bin/convert_xfm -ref $volB -in $volA -omedx $invxfmfilename -xfmtype $xfmtype -inverse $transAB"
    puts $thecommand
    catch { exec sh -c $thecommand } errmsg
    puts $errmsg

    set success [ file readable "$invxfmfilename" ]
    if { $success == 0 } {
	puts "No transformation saved!"
	return 4
    }

    return 0
}

