#!/usr/local/bin/bash

cat defaultschedule.h | grep comms.push_back | sed 's/.*(\"//' | sed 's/\").*$//'

