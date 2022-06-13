less -S feature-table.tsv |perl -ne '\@a=split;\$n++;if(/#/){print \$_}else{print "ASV\$n","\t",join("\t",\@a[1..\$\#a]),"\n";}' > new-table.tsv
biom convert -i new-table.tsv -o new.feature-table.biom --to-hdf5 --table-type="ASV table"
less -S dna-sequences.fasta|perl -ne 'chomp;if(/>/){\$n++;print ">ASV\$n\n"}else{print \$_,"\n";}' > new-sequence.fasta