## GeneMark.hmm run example

Create a test folder and run GeneMark.hmm with the following command

    mkdir test; cd test
    ../../../gmhmme3  -o genemark.gff3  -m ../input/athaliana.mod  -f gff3  ../input/sequence.fasta

If everything is configured correctly, the resulting `genemark.gff3` file  in the
`test` folder should match the contents of the `genemark.gff3` in the `output` folder.

To compare outputs:

    ../../../compare_intervals_exact.pl  --f1 genemark.gff3  --f2 ../output/genemark.gff3  --v

The expected runtime for this example is under one minute.


