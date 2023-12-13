## GeneMark-EP+ example

Create a test folder and run GeneMark-EP+ with the following command

    mkdir test; cd test
    ../../../gmes_petap.pl --seq ../input/genome.fasta --EP --dbep ../input/proteins.fasta --verbose --cores=8 --max_intergenic 10000 --mask_penalty 0

If everything is configured correctly, the resulting `genemark.gtf` file  in the
`test` folder should match the contents of the `genemark.gtf` in the `output` folder.

Note that the order of predicted genes in the GTD file can differ between the runs.
To compare outputs:

    ../../../compare_intervals_exact.pl  --f1 genemark.gtf  --f2 ../output/genemark.gtf  --v

The expected runtime for this example on a system with 8 CPUs (each 2.40GHz)
and 8GB of RAM is about **7 minutes**.

The example input genome comprises of the last 1Mbp of Arabidopsis thaliana's 5 chromosome.

