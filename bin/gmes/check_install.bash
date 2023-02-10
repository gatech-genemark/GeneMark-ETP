#!/bin/bash
# -----------------------
# Alex Lomsadze
# 2020
# GeneMark-ES installation check
# -----------------------

echo "Checking GeneMark-ES installation"

echo "Checking Perl setup"
for value in YAML Hash::Merge Parallel::ForkManager MCE::Mutex threads Thread::Queue Math::Utils
do
  if ( ! perl -M$value -e 1 ) ;
  then
    echo "Error, Perl CPAN library $value not found"
    exit 1
  fi
done
echo "All required Perl modules were found"

echo "Checking GeneMark.hmm setup"
dir=$(dirname $0)

if [ ! -x "$dir/gmhmme3" ] ;
then 
  echo "File gmhmme3 is not an executable or does not exist"
  exit 1
fi
echo "GeneMark.hmm was found"

message=$($dir/gmhmme3 -v 2>&1 >/dev/null)
if [[ $message == *"License key"* ]] ;
then
  echo "Error, installation key is missing or expired"
  exit 1
fi

if [[ $message == *"cannot execute"* ]] ;
then
  echo "Error, mismatch between the binary and OS was detected"
  exit 1
fi
echo "GeneMark.hmm is set"

message=$($dir/gmhmme3 -v)
if [[ $message != *"gmhmme3 : error, sequence file's name is not provided"* ]]
then
  echo "Error, unexpected condition was detected"
  exit 1
fi
echo "GeneMark.hmm is executable"

echo "Performing GeneMark.hmm test run"
if [ -d "GeneMark-E-tests/GeneMark.hmm/" ]
then
  mkdir -p "GeneMark-E-tests/GeneMark.hmm/test"
  $dir/gmhmme3  -o GeneMark-E-tests/GeneMark.hmm/test/genemark.gff3  -m GeneMark-E-tests/GeneMark.hmm/input/athaliana.mod  -f gff3  GeneMark-E-tests/GeneMark.hmm/input/sequence.fasta 

  if [ ! -e "GeneMark-E-tests/GeneMark.hmm/test/genemark.gff3" ]
  then
    echo "Error, GeneMark.hmm run test failed"
    exit 1
  fi
fi


echo "All required components for GeneMark-ES were found"

