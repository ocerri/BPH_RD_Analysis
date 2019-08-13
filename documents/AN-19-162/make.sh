#!/bin/bash

# config
compile=${1:-0} # compile main document (0) or auxiliary file (1)
tmpdir="tmp"

if (( ${compile} == 0 ))
then
    infile="AN-19-162.tex"
    outfile="${tmpdir}/AN-19-162_temp.pdf"
elif (( ${compile} == 1 ))
then
    infile="AN-19-162_aux.tex"
    outfile="${tmpdir}/AN-19-162_aux_temp.pdf"
else
    echo "Only options are 0 (compile main) or 1 (compile auxiliary)! Exiting..."
    exit
fi

# build PAPER
eval `utils/tdr runtime -sh`
tdr --temp_dir="${tmpdir}" --style=paper b "${infile}"

# detect OS to open PDF file
if [[ "${OSTYPE}" == "linux-gnu" ]] || [[ "${OSTYPE}" == "freebsd"* ]]
then
    xdg-open "${outfile}"
elif [[ "${OSTYPE}" == "darwin"* ]]
then
    open "${outfile}"
elif [[ "${OSTYPE}" == "cygwin" ]]
then
    cygstart "${outfile}"
else
    echo "Not sure which application to use to open PDF files for OS: ${OSTYPE}..."
    echo "Please add your own line in this script!"
fi
