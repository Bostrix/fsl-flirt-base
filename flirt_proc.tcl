#

# FLIRT - FMRIB's Linear Image Registration Tool
#
# Mark Jenkinson and Stephen Smith, FMRIB Image Analysis Group
#
# Copyright (C) 1999-2001 University of Oxford
#
# TCLCOPYRIGHT

#{{{ FixMedxTransform 

proc FixMedxTransform { xfmfname OriginX OriginY OriginZ XPixelSeparation YPixelSeparation ZPixelSeparation } {


    set notfound 0
    set optr [ open ${xfmfname}TMP w ]
    set fptr [ open $xfmfname r ]

    while { ( $notfound == 0 ) && ( [ gets $fptr line ] >= 0 ) } {
	if { [ regexp outputusermatrix $line ] } {
	    set notfound 1
	}
	puts $optr $line
    }

    puts $optr "            $XPixelSeparation"
    puts $optr "            0"
    puts $optr "            0"
    puts $optr "            $OriginX"
    puts $optr "            0"
    puts $optr "            $YPixelSeparation"
    puts $optr "            0"
    puts $optr "            $OriginY"
    puts $optr "            0"
    puts $optr "            0"
    puts $optr "            $ZPixelSeparation"
    puts $optr "            $OriginZ"
    puts $optr "            0"
    puts $optr "            0"
    puts $optr "            0"
    puts $optr "            1"

    for { set i 0 } { $i < 16 } { incr i 1 } {
	gets $fptr line
    }

    while { [ gets $fptr line ] >= 0 } {
	puts $optr $line
    }
    
    close $fptr
    close $optr

    set errrr [ catch { exec sh -c "/bin/mv -f ${xfmfname}TMP ${xfmfname}" } junk ]    
}

#}}}
#{{{ TalairachMedxTransform

proc TalairachMedxTransform { xfmfname } {

    set optr [ open ${xfmfname}TMP w ]
    set fptr [ open $xfmfname r ]

    while { [ gets $fptr line ] >= 0 } {
	if { [ regexp talairachcalibrated $line ] } {
	    puts $optr "        /talairachcalibrated?	true"
	} else {
	    puts $optr $line
	}
    }
    
    close $fptr
    close $optr

    set errrr [ catch { exec sh -c "/bin/mv -f ${xfmfname}TMP ${xfmfname}" } junk ]    
}

#}}}
#{{{ ConcatedTalairachMedxTransform

proc ConcatedTalairachMedxTransform { xfmfname } {

    set optr [ open ${xfmfname}TMP w ]
    set fptr [ open $xfmfname r ]
    set notfixed 1

    while { [ gets $fptr line ] >= 0 } {
	if { $notfixed && [ regexp ">>" $line ] } {
	    puts $optr "        /talairachcalibrated?	true"
	    set notfixed 0
	}
	puts $optr $line
    }
    
    close $fptr
    close $optr

    set errrr [ catch { exec sh -c "/bin/mv -f ${xfmfname}TMP ${xfmfname}" } junk ]    
}

#}}}
#{{{ flirt:proc

proc flirt:proc { regmode refname testname testname2 nstats statslist output dof doftwo bins searchrxmin searchrxmax searchrymin searchrymax searchrzmin searchrzmax disablesearch_yn cost interp sincwidth sincwindow refweight inweight inweight2 popups } {

    global PXHOME FSLDIR USER MEDXV HOME INMEDX CLUSTERRSH

    #{{{ setup options

#{{{ setup MEDx stuff

if { $INMEDX } {

    set TempFileBase [ exec sh -c "${FSLDIR}/bin/tmpnam ${HOME}/.fl" ]

    MxGetCurrentFolder Folder
 
    if { ! [ info exists Folder ] } {
	if { $popups } {
	    LocalErrorBox "A folder must be opened!"
	}
	puts "A folder must be opened!"
	return 1
    }
}

#}}}

set flirtoptions "-bins $bins -cost $cost"
set dofoptions ""
set doftwooptions ""

if { $dof == "2D" } {
    set dofoptions "-2D -dof 12"
} else {
  if { $dof == "TRANS" } {
      set dofoptions "-dof 6 -schedule ${FSLDIR}/etc/flirtsch/xyztrans.sch"
  } else {
    set dofoptions "-dof $dof"
  }
}

if { $doftwo == "2D" } {
    set doftwooptions "-2D -dof 12"
} else {
  if { $doftwo == "TRANS" } {
      set doftwooptions "-dof 6 -schedule ${FSLDIR}/etc/flirtsch/xyztrans.sch"
  } else {
    set doftwooptions "-dof $doftwo"
  }
}


if { $disablesearch_yn } {
    set flirtoptions "$flirtoptions -searchrx 0 0 -searchry 0 0 -searchrz 0 0"
} else {
    set flirtoptions "$flirtoptions -searchrx $searchrxmin $searchrxmax -searchry $searchrymin $searchrymax -searchrz $searchrzmin $searchrzmax"
}

set flirtweights1 ""
set flirtweights2 ""
if { $refweight != "" } {
    if { $INMEDX } {
	MxGetPageByName $Folder $refweight tmpname
	FSLSaveAs $tmpname AVW ${TempFileBase}refweight.hdr true
	set refweight ${TempFileBase}refweight
    }
    set flirtweights1 "$flirtweights1 -refweight $refweight"
}
if { $inweight != "" } {
    if { $INMEDX } {
	MxGetPageByName $Folder $inweight tmpname
	FSLSaveAs $tmpname AVW ${TempFileBase}inweight.hdr true
	set inweight ${TempFileBase}inweight
    }
    set flirtweights1 "$flirtweights1 -inweight $inweight"
    set flirtweights2 "$flirtweights2 -refweight $inweight"
}
if { $inweight2 != "" } {
    if { $INMEDX } {
	MxGetPageByName $Folder $inweight2 tmpname
	FSLSaveAs $tmpname AVW ${TempFileBase}inweight2.hdr true
	set inweight2 ${TempFileBase}inweight2
    }
    set flirtweights2 "$flirtweights2 -inweight $inweight2"
}

set flirtinterp "-interp $interp"
if { $interp == "sinc" } {
    set flirtinterp "$flirtinterp -sincwidth $sincwidth -sincwindow $sincwindow"
}

#}}}

    if { $INMEDX} {
	#{{{ check standard image

    set problem [ MxGetPageByName $Folder $refname refptr ]

    if { $problem == 1 } {
	if { $popups } {
	    LocalErrorBox "Standard Image is not a valid page!"
	}
	puts "Standard Image is not a valid page!"
	return 1
    }

    MxGetPageProperties $refptr Props

    if { [ string match [ keylget Props Type ] Volume ] != 1 } {
	if { $popups } {
	    LocalErrorBox "Standard Image isn't a volume!"
	}
	puts "Standard Image isn't a volume!"
	return 1
    }

#}}}
	#{{{ check image/group 1

    set problem [ MxGetPageByName $Folder $testname testptr ]

    if { $problem == 1 } {
	if { $popups } {
	    LocalErrorBox "Reslice Image/Group is not a valid page!"
	}
	puts "Reslice Image/Group is not a valid page!"
	return 1
    }

    MxSetCurrentPage $testptr
    set testtype [ what_type 0 ]
    if { $testtype != 3 && $testtype != 4 } {
	if { $popups } {
	    LocalErrorBox "Reslice Image/Group isn't a volume or a group!"
	}
	puts "Reslice Image/Group isn't a volume or a group!"
	return 1
    }

    if { $testtype == 3 } {
	set testpages $testptr
	set testvolumes 1
    } else {
	MxGetPagesFromGroup $testptr testpages
	set testvolumes [ llength $testpages ]
    }

#}}}
	#{{{ check secondary (stats) group

    if { $nstats > 0 } {
	set i 1
	while { $i <= $nstats } {
	    set problem [ MxGetPageByName $Folder [ lindex $statslist [ expr $i - 1 ] ] secondptr($i) ]
	    if { $problem == 1 } {
		if { $popups } {
		    LocalErrorBox "Secondary Image/Group(s) not valid page(s)!"
		}
		puts "Secondary Image/Group(s) not valid page(s)!"
		return 1
	    } else {
		MxSetCurrentPage $secondptr($i)
		set secondtype [ what_type 0 ]
		if { $testtype != $secondtype } {
		    if { $popups } {
			LocalErrorBox "All Secondary and High Res entries must either be single volumes or be groups!"
		    }
		    puts "All Secondary and High Res entries must either be single volumes or be groups!"
		    return 1
		}

		if { $testtype == 3 } {
		    set secondpages($i) $secondptr($i)
		    set secondvolumes 1
		} else {
		    MxGetPagesFromGroup $secondptr($i) secondpages($i)
		    set secondvolumes [ llength $secondpages($i) ]
		}
		
		if { $testvolumes != $secondvolumes } {
		    if { $popups } {
			LocalErrorBox "Secondary Image/Group doesn't contain the right number of images!"
		    }
		    puts "Secondary Image/Group doesn't contain the right number of images!"
		    return 1
		}
	    }
	    incr i 1
	}
    }

#}}}
	#{{{ check image/group 2

    if { $regmode == 2 } {
	set problem [ MxGetPageByName $Folder $testname2 testptr2 ]
	if { $problem == 1 } {
	    if { $popups } {
		LocalErrorBox "Low Res Image/Group is not a valid page!"
	    }
	    puts "Low Res Image/Group is not a valid page!"
	    return 1
	} else {
	    MxSetCurrentPage $testptr2
	    set testtype2 [ what_type 0 ]
	    if { $testtype != $testtype2 } {
		if { $popups } {
		    LocalErrorBox "Low Res and High Res entries must either both be single volumes or both be groups!"
		}
		puts "Low Res and High Res entries must either both be single volumes or both be groups!"
		return 1
	    }

	    if { $testtype == 3 } {
		set testpages2 $testptr2
		set testvolumes2 1
	    } else {
		MxGetPagesFromGroup $testptr2 testpages2
		set testvolumes2 [ llength $testpages2 ]
	    }

	    if { $testvolumes != $testvolumes2 } {
		if { $popups } {
		    LocalErrorBox "Low Res Image/Group doesn't contain the right number of images!"
		}
		puts "Low Res Image/Group doesn't contain the right number of images!"
		return 1
	    }
	}
    }

#}}}
	#{{{ register image/group 1 to standard

    set i 0
    set NewList ""
    foreach testpage $testpages {
	set returnval [ flirt:medxrun $refptr $testpage $flirtoptions $dofoptions $flirtweights1 $flirtinterp ${TempFileBase}_${i} ]
	MxGetCurrentPage output
	lappend NewList $output
	incr i 1
    }
    if { $testtype == 4 } {
	MxGroupVolumes $NewList "Registered Group" OutputGroup
    }

#}}}
	#{{{ register image/group 2 to image/group 1 and then apply concated transforms to group 2

    if { $regmode == 2 } {

	set i 0
	set NewList ""
	set NewList2 ""
	foreach testpage2 $testpages2 {

	    set refptr2 [ lindex $testpages $i ]
	    set returnval [ flirt:medxrun $refptr2 $testpage2 $flirtoptions $doftwooptions $flirtweights2 $flirtinterp ${TempFileBase}_level2_${i} ]
	    MxGetCurrentPage output
	    lappend NewList $output
	    MxConcatTransforms ${TempFileBase}_level2_${i}_result.xfm AIR ${TempFileBase}_${i}_result.xfm AIR ${TempFileBase}_combined_${i}_result.xfm
	    MxGetImageProperties $refptr ImageProps1
	    if { [ regexp true [ keylget ImageProps1 TalairachCalibrated ] ] } {
		ConcatedTalairachMedxTransform ${TempFileBase}_combined_${i}_result.xfm
	    }

	    MxSetCurrentPage $testpage2
	    if { $MEDXV == 3.2 } {
		MxApplySavedTransform $testpage2 ${TempFileBase}_combined_${i}_result.xfm AIR resultptr   
	    } else {
		MxApplySavedTransform $testpage2 ${TempFileBase}_combined_${i}_result.xfm Reslice F false resultptr   
	    }
	    MxGetImageProperties $testpage2 ImageProps2
	    name_page $resultptr "Shadow Registration of [keylget ImageProps2 Name] to [keylget ImageProps1 Name]"
	    MxSetCurrentPage $resultptr
	    lappend NewList2 $resultptr

	    incr i 1
	}
	if { $testtype == 4 } {
	    MxGroupVolumes $NewList "Registered Low Res to High Res Group" OutputGroup2
	    MxGroupVolumes $NewList2 "Shadow Registered Low Res Group" ShadowOutputGroup2
	}
    }

#}}}
	#{{{ apply to stats images

    if { $nstats > 0 } {
	set ii 1
	while { $ii <= $nstats } {
	    set i 0
	    set NewList ""
	    foreach secondpage $secondpages($ii) {
		MxSetCurrentPage $secondpage
		if { $regmode == 1 } {
		    set xfmname ${TempFileBase}_${i}_result.xfm
		} else {
		    set xfmname ${TempFileBase}_combined_${i}_result.xfm
		}
		if { $MEDXV == 3.2 } {
		    MxApplySavedTransform $secondpage $xfmname AIR resultptr   
		} else {
		    if { $regmode == 1 } {
			MxApplySavedTransform $secondpage $xfmname AIR F false resultptr
		    } else {
			MxApplySavedTransform $secondpage $xfmname Reslice F false resultptr
		    }
		}
		MxGetImageProperties $refptr ImageProps1
		MxGetImageProperties $secondpage ImageProps2
		name_page $resultptr "Shadow Registration of [keylget ImageProps2 Name] to [keylget ImageProps1 Name]"
		MxSetCurrentPage $resultptr
		lappend NewList $resultptr
		incr i 1
	    }
	    if { $testtype == 4 } {
		MxGroupVolumes $NewList "Shadow Registered Group $ii" ShadowOutputGroup($ii)
	    }
	    incr ii 1
	}
    }

#}}}
	catch { exec sh -c "rm -f ${TempFileBase}*" } errmsg
    }  else {
	#{{{ non-MEDx FLIRTing

set outroot [ file rootname $output ]

if { $regmode == 1 } {
    set thecommand "$CLUSTERRSH ${FSLDIR}/bin/flirt -in $testname -ref $refname -out $output -omat ${outroot}.mat $flirtoptions $dofoptions $flirtweights1 $flirtinterp"
    puts $thecommand
    catch { exec sh -c $thecommand } errmsg
    puts $errmsg
} else {
    set thecommand "$CLUSTERRSH ${FSLDIR}/bin/flirt -in $testname -ref $refname -omat ${outroot}1.mat $flirtoptions $dofoptions $flirtweights1"
    puts $thecommand
    catch { exec sh -c $thecommand } errmsg
    puts $errmsg

    set thecommand "$CLUSTERRSH ${FSLDIR}/bin/flirt -in $testname2 -ref $testname -omat ${outroot}2.mat $flirtoptions $doftwooptions $flirtweights2"
    puts $thecommand
    catch { exec sh -c $thecommand } errmsg
    puts $errmsg

    set thecommand "${FSLDIR}/bin/convert_xfm -matonly -concat ${outroot}1.mat -omat ${outroot}.mat ${outroot}2.mat"
    puts $thecommand
    catch { exec sh -c $thecommand } errmsg
    puts $errmsg

    set thecommand "${FSLDIR}/bin/flirt -in $testname2 -ref $refname -out $output -applyxfm -init ${outroot}.mat $flirtinterp"
    puts $thecommand
    catch { exec sh -c $thecommand } errmsg
    puts $errmsg

}

for { set i 1 } { $i <= $nstats } { incr i 1 } {
    set statsname [ lindex $statslist [ expr $i - 1 ] ]

    set thecommand "${FSLDIR}/bin/flirt -in $statsname -ref $refname -out [ file rootname ${output} ]_shadowreg_[ file tail [ file rootname $statsname ] ] -applyxfm -init ${outroot}.mat $flirtinterp"
    puts $thecommand
    catch { exec sh -c $thecommand } errmsg
    puts $errmsg
}

puts Finished

set returnval 0

#}}}
    }

    return $returnval
}

#}}}
#{{{ flirt:medxrun

proc flirt:medxrun { reffile infile flirtoptions dofoptions flirtweights flirtinterp TempFileBase } {

    # setup vars
    global FSLDIR MEDXV CLUSTERRSH INMEDX
    
    # save images and fix permissions
    set refptr $reffile
    set inptr  $infile
    MxGetCoordinateSystem $refptr OriginX OriginY OriginZ XPixelSeparation YPixelSeparation ZPixelSeparation
    MxGetImageProperties $refptr ImageProps1
    set reffile ${TempFileBase}_ref
    set infile  ${TempFileBase}_in
    set xfmfilename ${TempFileBase}_result.xfm
    FSLSaveAs $refptr AVW ${reffile}.hdr true
    FSLSaveAs $inptr  AVW ${infile}.hdr  true
    catch { exec sh -c "/bin/chmod 755 ${reffile}* ${infile}*" } junk
    
    # set up the command to be executed
    set fullcommand "$CLUSTERRSH ${FSLDIR}/bin/flirt -ref $reffile -in $infile -omedx $xfmfilename $flirtoptions $dofoptions $flirtweights $flirtinterp"

    # run command
    ScriptUpdate "$fullcommand"
    puts $fullcommand
    set result [ catch { exec sh -c $fullcommand } ErrMsg ]
    if { $result != 0 } {
	puts "$ErrMsg"
    }
    CancelScriptUpdate

    set success [ file readable "$xfmfilename" ]
    if { $success == 0 } {
	puts "No transformation saved!"
	return 4
    }

    FixMedxTransform $xfmfilename $OriginX $OriginY $OriginZ $XPixelSeparation $YPixelSeparation $ZPixelSeparation
	
    if { [ regexp true [ keylget ImageProps1 TalairachCalibrated ] ] } {
	TalairachMedxTransform $xfmfilename
    }

    MxSetCurrentPage $inptr

    if { $MEDXV == 3.2 } {
	MxApplySavedTransform $inptr "$xfmfilename" AIR outptr   
    } else {
	MxApplySavedTransform $inptr "$xfmfilename" AIR F false outptr
    }

    MxGetImageProperties $inptr ImageProps2
    name_page $outptr "Registration of [keylget ImageProps2 Name] to [keylget ImageProps1 Name]"
    MxSetCurrentPage $outptr

    return 0
}

#}}}
