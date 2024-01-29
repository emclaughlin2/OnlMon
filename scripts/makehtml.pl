#!/usr/bin/perl

use strict;
use warnings;

sub findruns;

my $histodir = sprintf("/sphenix/lustre01/sphnxpro/commissioning/online_monitoring/histograms");
my @subsystems = ("BBCMON", "CEMCMON", "IHCALMON", "INTTMON", "LL1MON", "OHCALMON", "TPCMON", "TPOTMON");

for my $subsys (@subsystems)
{
    my %donehash = ();
    my $rundone = sprintf("%s.done",$subsys);
    if (-f $rundone)
    {
	open(F,"$rundone");
	while (my $run = <F>)
	{
            chomp $run;
	    $donehash{$run} = 1;
	}
	close(F); 
    }
    my $subsysdir = sprintf("%s/%s",$histodir,$subsys);
    my %todoruns = findruns($subsysdir,$subsys, \%donehash );
# store the last run and exclude from run done list so
# we process it again in the next round, just in case we
# transferred the histo files before the last one arrived
    my $lastrun = 0;
    for my $run (sort { $a <=> $b } keys %todoruns)
    {
	if ($run < 80000) # tpot/intt creates large bad runnumbers
	{
	    my $listfile = sprintf("%s.list",$subsys);
	    my $cmd = sprintf("ls -1 /sphenix/lustre01/sphnxpro/commissioning/online_monitoring/histograms/%s/Run_%d-* > %s.list",$subsys,$run,$subsys);
	    print "$cmd\n";
	    system($cmd);
	    my $rootcmd = sprintf("root.exe -q makehtml.C\\\(\\\"%s.list\\\",\\\"%s\\\"\\\)",$subsys,$subsys);
	    print "$rootcmd\n";
	    system($rootcmd);
	    print "handled $run for $subsys\n";
	    $lastrun = $run; #process the last run again next ti
	}
    }
    for my $run (sort { $a <=> $b } keys %todoruns)
    {
	if ($lastrun == $run)
	{
	    next;
	}
	my $updatedone = sprintf("echo %d >> %s",$run,$rundone);
	system($updatedone);
    }

}

sub findruns
{
    my $subsysdir = shift;
    my $subsys = shift;
    my %donehash = %{$_[0]};
    opendir my $dir, $subsysdir  or die "Cannot open directory: $!";
    my @files = readdir $dir;
    closedir $dir;
    my %todoruns = ();
    for my $file (@files)
    {
	if ($file =~ /$subsys/)
	{
#	    print "found $file\n";
	    if ($file =~ /(\S+)_(\d+)-(\S+).*\..*/ )
	    {
		my $run = int($2);
#		print "found run $run\n";
		if (exists $donehash{$run})
		{
		    next;
		}
		if ($run > 80000)
		{
		    next;
		}
		$todoruns{$run} = 1;
	    }
	}
    }
    return %todoruns;
}

