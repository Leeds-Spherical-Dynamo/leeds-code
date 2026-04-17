#! /bin/bash


make -f Makefile-tests clean

for CASE in 0 1 2 0_symmetry 1_symmetry 2_symmetry
do
  echo "======================
  Making Christensen case ${CASE}
  ========================="

  make -f Makefile-tests CASE=${CASE}
  if [ $? -eq 0 ]
  then
    echo "Created case ${CASE}"
  else
    echo "Failed to create case ${CASE}" >&2
    exit 1
  fi
  make -f Makefile-tests CASE=${CASE} test
  if [ $? -eq 0 ]
  then
    echo "Installed case ${CASE}"
  else
    echo "Failed to install case ${CASE}" >&2
    exit 1
  fi
  make -f Makefile-tests clean
done

cp ./tests/check_output.* ./test_run/
if [ $? -eq 0 ]
then
  echo "Copied test scripts"
  exit 0
else
  echo "Failed to copy test scripts" >&2
  exit 1
fi
