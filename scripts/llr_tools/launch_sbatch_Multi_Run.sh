#!/bin/bash
ARGV=`basename -a $1`
set -eu
script_file_name=$(basename "$0")
tput bold;
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
echo "!                                                                       !"
echo "!     Code to plot grid output from the run files                       !"
echo "!     $script_file_name                                        !"
echo "!     [Author]: Frederic Bonnet November 2025                           !"
echo "!     [usage]: bash -s < ./launch_sbatch_Multi_Run.sh   {Input list}    !"
echo "!     [example]: bash -s < ./launch_sbatch_Multi_Run.sh search_string   !"
echo "!                                                                       !"
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
tput sgr0;
#colors
red="tput setaf 1"  ;green="tput setaf 2"  ;yellow="tput setaf 3"
blue="tput setaf 4" ;magenta="tput setaf 5";cyan="tput setaf 6"
white="tput setaf 7";bold=""               ;reset_colors="tput sgr0"
#-------------------------------------------------------------------------------
# Method to check if file exists
#-------------------------------------------------------------------------------
file_exists (){
_input_filename=$1
#-------------------------------------------------------------------------------
file_exists="no"
if [ -f "${_input_filename}" ]
  then
    $white; printf "File                   : "; $bold;
    $yellow; printf '%s'"${_input_filename}"; $green; printf " --->: exist.\n";
    $white; $reset_colors;
    file_exists="yes"
    printf "      file_exists ---->: "; $bold;
    $green; printf "$file_exists\n"; $reset_colors;
  else
    $white; printf "File                   : "; $bold;
    $yellow; printf '%s'"${_input_filename}"; $red;printf " --->: does not exist.\n";
    $white; $reset_colors;
    file_exists="no"
    printf "      file_exists ---->: "; $bold;
    $red; printf "$file_exists\n"; $reset_colors;
fi
}
finish_up_message (){
  #-------------------------------------------------------------------------------
  #End of the script
  echo
  echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
  $cyan; echo `date`; $blue;
  echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
  echo "-                  launch_sbatch_Multi_Run.sh Done.                     -"
  echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
  $reset_colors
}
#-------------------------------------------------------------------------------
# Check if executable exists
#-------------------------------------------------------------------------------
# Application paths
grid_dwf_telos_dir=/home/dp208/dp208/dc-bonn2/SwanSea/SourceCodes/Grid-Telos-Devel/Grid
grid_dwf_telos_build_dir=$grid_dwf_telos_dir/build-Sp2n-withLLR
path_to_test_llr=$grid_dwf_telos_build_dir/tests/llr/

file_exists "${path_to_test_llr}/Test_llr_contrained_action"
echo
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo
if [ "$file_exists" = 'yes' ]
then
  printf "                       : "; $bold;
  $white; printf "YES ---> we can proceed ....\n"; $reset_colors;
  # Submitting the batch script to the slurm queue.
elif [ "$file_exists" = 'no' ]
then
  printf "                       : "; $bold;
  $white; printf "NO  ---> ${path_to_test_llr}/Test_llr_contrained_action does exits, we will exit\n"; $reset_colors;
  finish_up_message;
  exit;
fi;
#-------------------------------------------------------------------------------
# Getting the common code setup and variables, #setting up the environment properly.
#-------------------------------------------------------------------------------
file_in=$1
#-------------------------------------------------------------------------------
# Stripping and creating the files
#-------------------------------------------------------------------------------
filename="${file_in}"
basename="${filename%.sh}"

declare -a mdstep_testcase=(18 22 26 28 30 32 34 36 38 40 42 44 46 48 50)
#-------------------------------------------------------------------------------
# Stripping and creating the files
#-------------------------------------------------------------------------------
echo
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo
H=1
for i in $(seq 0 `expr ${#mdstep_testcase[@]} - 1`)
do
  MD=${mdstep_testcase[$i]}

  #XSecho "Generating file with MDsteps=$MD"

  # Output filename
  outfile="${basename}_MDsteps-${MD}.sh"

  $cyan;printf "File --->: ";$red;printf "${file_in}";$reset_colors;
  $cyan;printf " <---> "; $green;printf "${i}";$reset_colors;
  $cyan;printf " <---> MDsteps --->: "; $magenta;printf "${mdstep_testcase[$i]}";$reset_colors;
  $cyan;printf " <---> Output file --->: ";$yellow;printf "${outfile}\n";
  $white; $reset_colors;

  # Replace MDsteps value and write new file
  sed \
  -e "s/^MDsteps=.*/MDsteps=${MD}/" \
  -e "s/^#SBATCH --job-name=.*/#SBATCH --job-name=${basename}_MDsteps-${MD}/" \
  -e "s/^job_name=.*/job_name=${basename}_MDsteps-${MD}/" \
  -e "s/^#SBATCH --qos=.*/#SBATCH --qos=standard/" \
  "$filename" > "$outfile"

  #sbatch "$outfile";

  index_i=$(printf "%04d" "$i")
  index_H=$(printf "%04d" "$H")
  $cyan; printf "#------>: i #------>: "; $green;   printf "${index_i} ";
  $cyan; printf "#------>: H #------>: "; $magenta; printf "${index_H} ";
  $cyan; printf "#------>: File #--->: "; $red;     printf "${outfile}";
  $cyan; printf "\n"; $white; $reset_colors;

  file_exists "${outfile}"

  if [ "$file_exists" = 'yes' ]
  then
    printf "                       : "; $bold;
    $white; printf "YES ---> sbatch submitting to the queue....\n"; $reset_colors;
    # Submitting the batch script to the slurm queue.

    sbatch "${outfile}" >> ./"Batch_submission.log" &

  elif [ "$file_exists" = 'no' ]
  then
    printf "                       : "; $bold;
    $white; printf "NO  ---> sbatch NO GO.\n"; $reset_colors;
  fi;

  H=$(expr $H + 1)

done
#-------------------------------------------------------------------------------
# Finishing up
#-------------------------------------------------------------------------------
finish_up_message;
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------


