-------make-multiple-jobs:
location: /work/hallc/alphaE/ruonanli/simc_michael_copy/makejobs
scripts: kinXX_vcs/pions_makejobs.pl to generate jobcards
		 swif.pl to submit jobs to the farm

-------generate SIMC input files with random seeds:
location:/work/hallc/alphaE/ruonanli/simc_michael_copy/infiles/gen_pion/vcs_infiles
scripts: random_replace.sh to replace random seeds

First copy initial inp files here, generate files and check. Then move generated correct files to the upper directory.

-------conver root
location:/work/hallc/alphaE/ruonanli/simc_michael_copy/worksim/TBC
script:convert_root.sh & merge_pion_simc.sh

copy .rzdat files here and run script to convert .rzdat to .root
Then merge root files
In order to save spaces, delete individual root files after converting. Do keep .rzdat files.
