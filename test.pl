#use File::Temp qw(tempdir);

#system("setup cernlib/2005");
system("cd /w/hallc-scifs17exp/alphaE/ruonanli/simc_michael/"); #path for tar rc_externals 
my $infile = "coin_pi0_kin2a_scincut_ytarcorr3.inp"; 
my @in = ($infile);
#my $dir = tempdir( CLEANUP => 1);
open(my $fh, ">in.txt");
	foreach(@in){
		print $fh $_."\n";
}

close ($fh);

#my @output = qx{simc < "@in\n"};
my @output = qx{simc < in.txt};

print "@output\n"
