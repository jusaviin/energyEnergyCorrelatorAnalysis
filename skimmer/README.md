# Running the skimmer

1) Test the code locally. Compile it with the Makefile by just typing

    make
 
   After compilation, you can run the skimmer:

   ./skimmer

   Running without arguments gives usage instructions for the main program. For the file name list provide a text file with one file name in each row. Check that everything works and that the output file contains the branches you desire.

2) Tar the files for running in lxplus:

    tar -cvzf skimmer.tar.gz skimmer.cxx Makefile

    Copy this file to your running area in lxplus.

3) Set up crab configuration in lxplus. Example of configuration can be found from the crab folder. Send the jobs to crab.

4) Enjoy your newly skimmed files!
