#!/usr/local/bin/bash

if [ $# -lt 1 ] ; then
  echo "Usage: $0 rcfile"
  exit -1
fi

cat $1 | sed 's/^/  comms.push_back("/' | sed 's/$/");/'

