#!/bin/bash

PROJ_DIR="$HOME/landrace10"
PGGB_DIR="$PROJ_DIR/04_pggb"
POST_PGGB_DIR="$PROJ_DIR/05_post_pggb"
LOG_DIR="$PROJ_DIR/00_logs"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
THREADS=10

{

find  

} > "$LOG_DIR/sv_analysis_$TIMESTAMP.log" 2>&1

