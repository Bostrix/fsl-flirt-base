#!/bin/sh

if [ $# -lt 3 ] ; then
  echo "Usage: $0 <inputimage> <refimage> <outputmat> [qform|sform] [qform|sform]"
  exit 0;
fi

im1=`$FSLDIR/bin/remove_ext $1`
im2=`$FSLDIR/bin/remove_ext $2`

# extract the qform and sform and calculate the flirt vox2mm (samp) matrix
nam=""
oldnam=""
for fn in $im1 $im2 ; do
  oldnam=$nam
  nam=`$FSLDIR/bin/tmpnam`;
  $FSLDIR/bin/avwhd $fn | grep qto_xyz | sed 's/qto_xyz:. *//' > $nam.qform.mtx
  $FSLDIR/bin/avwhd $fn | grep sto_xyz | sed 's/sto_xyz:. *//' > $nam.sform.mtx

  nx=`$FSLDIR/bin/avwval $fn dim1`;
  dx=`$FSLDIR/bin/avwval $fn pixdim1`;
  dy=`$FSLDIR/bin/avwval $fn pixdim2`;
  dz=`$FSLDIR/bin/avwval $fn pixdim3`;

  echo "$dx 0 0 0" > $nam.samp.mtx
  echo "0 $dy 0 0" >> $nam.samp.mtx
  echo "0 0 $dz 0" >> $nam.samp.mtx
  echo "0 0 0 1" >> $nam.samp.mtx

  if [ `$FSLDIR/bin/avworient -getorient $fn` = NEUROLOGICAL ] ; then
    nx1=`echo \( $nx - 1 \) | bc -l`;
    echo "-1 0 0 $nx1" > $nam.swapx.mtx
    echo "0 1 0 0" >> $nam.swapx.mtx
    echo "0 0 1 0" >> $nam.swapx.mtx
    echo "0 0 0 1" >> $nam.swapx.mtx
    $FSLDIR/bin/convert_xfm -omat $nam.samp.mtx -concat $nam.samp.mtx $nam.swapx.mtx
  fi
done

nam1=$oldnam
nam2=$nam

# decide which of sform or qform will be used
codetype1=sform
code1=`$FSLDIR/bin/avwval $im1 sform_code`;
if [ X$code1 = X0 ] ; then
  codetype1=qform
  code1=`$FSLDIR/bin/avwval $im1 qform_code`;
fi
codetype2=sform
code2=`$FSLDIR/bin/avwval $im2 sform_code`;
if [ X$code2 = X0 ] ; then
  codetype2=qform
  code2=`$FSLDIR/bin/avwval $im2 qform_code`;
fi

if [ $# -ge 4 ] ; then
  codetype1=$4;
  code1=`$FSLDIR/bin/avwval $im1 ${codetype1}_code`;
fi

if [ $# -ge 5 ] ; then
  codetype2=$5;
  code2=`$FSLDIR/bin/avwval $im2 ${codetype2}_code`;
fi


# default init is the identity
cp $FSLDIR/etc/flirtsch/ident.mat $3

# if the sform or qform codes are non-zero then calculate the init matrix
# formula is: init = (samp2) (qsform2)^{-1} (qsform1) (samp1)^{-1}
if [ X$code1 != X0 ] ; then
  # (it might be good to have manual control to force some of this logic)
  # if [ X$code2 = X$code1 ] ; then  # alternative to only do for same space
  if [ X$code2 != X0 ] ; then
     $FSLDIR/bin/convert_xfm -omat ${nam1}.samp.inv.mtx -inverse ${nam1}.samp.mtx 
     $FSLDIR/bin/convert_xfm -omat $3 -concat ${nam1}.samp.inv.mtx $3
     $FSLDIR/bin/convert_xfm -omat $3 -concat ${nam1}.${codetype1}.mtx $3
     $FSLDIR/bin/convert_xfm -omat ${nam2}.${codetype2}.inv.mtx -inverse ${nam2}.${codetype2}.mtx
     $FSLDIR/bin/convert_xfm -omat $3 -concat ${nam2}.${codetype2}.inv.mtx $3
     $FSLDIR/bin/convert_xfm -omat $3 -concat ${nam2}.samp.mtx $3
  fi
fi

# clean up
if [ X${nam1} != X ] ; then rm -f ${nam1}* ; fi
if [ X${nam2} != X ] ; then rm -f ${nam2}* ; fi
