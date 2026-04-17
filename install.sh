#!/usr/bin/env bash
# Check if Makefile, parallel.h, or parameters.F90 exist
# if they do exist, check with user whether to overwrite.
# if don't exist, copy from Scripts

# files to act on
files=("./Makefile" "./program/parameters.F90" "./parallel.h" "./Makefile-tests")
sourcefiles=("./Scripts/Makefile_template" "./Scripts/parameters.F90_template" "./Scripts/parallel.h_template" "./Scripts/Makefile-tests_template")

# Check if source exists before copying
mycopy () {
    if ! [ -f ${sourcefiles[$1]} ]; then
        echo "can't copy ${sourcefiles[$1]}"
        echo "Error, default file not found."
        exit 1
    else
        cp ${sourcefiles[$1]} ${files[$1]} 
    fi
}

# Loop through all files
for i in ${!files[@]}; do
    if [ -f ./${files[$i]} ]; then
        # if file already exists, check with user
        echo "${files[$i]} already present,"
        read -p "    Do you wish to copy ${sourcefiles[$i]} to ${files[$i]}? (yes/no):" yesno
        case ${yesno} in 
            # copy if user wants to
            yes ) echo "    Ok, copying ${sourcefiles[$i]} to ${files[$i]}...";
                mycopy ${i};
                echo "    ...done";;
            # don't copy if user doesn't want to
            no ) echo "    Ok, leaving ${files[$i]} as is.";;
            # don't copy if invalid response
            * ) echo "    invalid response, not copying.";;
        esac
    else
        # if file doesn't exist, copy
        echo "${files[$i]} not present in root of repository, copying..."
        mycopy ${i}
        echo "   ...done"
    fi
done
